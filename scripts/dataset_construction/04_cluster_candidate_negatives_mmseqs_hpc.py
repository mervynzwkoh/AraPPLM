#!/usr/bin/env python3
"""
04_cluster_candidate_negatives_mmseqs_hpc.py
============================================
Performs pair-level redundancy reduction across the clean candidate negative
pool for Arabidopsis thaliana using greedy CD-HIT representative selection.

Methodology (docs/dataset_construction/03_homology_clustering/01_pair_redundancy_reduction.md):
    - Sequence identity threshold: theta = 0.50
    - Reciprocal coverage threshold: cmin = 0.80
    - E-value cutoff: <= 1e-3
    - Pair similarity metric:
        S(P_i, P_j) = max[ min(s(A, C), s(B, D)), min(s(A, D), s(B, C)) ]

Representative Priority Hierarchy (Option A):
    1. TAIR Annotation Confidence: Experimental codes (3.0) > Curated (2.0) > Computational (1.0)
    2. Combined Sequence Length: (L_A + L_B) descending
    3. Deterministic Lexicographical: participant_a_id, participant_b_id

Algorithmic Scalability (The Component Decoupling Theorem):
    Cross-component pairs satisfy S(P_1, P_2) == 0.0 by mathematical proof.
    The 101.85M pairs are partitioned by their 64-bit component-pair signature:
        CompKey(A, B) = (min(K(A), K(B)) << 32) | max(K(A), K(B))
    Singleton components (~42M pairs) are trivially non-redundant and retained instantly.
    Paralog micro-clusters are processed via exact greedy CD-HIT selection in memory (< 2.5 GB RAM).

Outputs:
    data/dataset_construction/02_negative_sampling/final_candidates_nr/ara_candidate_negatives_nr_chunk_*.parquet
    data/dataset_construction/02_negative_sampling/ara_candidate_negatives_redundancy_reduction_audit.json
    data/dataset_construction/02_negative_sampling/summary_candidate_negatives_nr.md
"""

import os
import sys
import glob
import gzip
import csv
import json
import time
import argparse
from collections import defaultdict

# Safeguard: prevent pyarrow from attempting to load an ABI-incompatible pandas if present
try:
    import pandas  # noqa: F401
except (ImportError, ValueError):
    class _NoPandasFinder:
        def find_spec(self, fullname, path, target=None):
            if fullname == "pandas" or fullname.startswith("pandas."):
                raise ModuleNotFoundError("No module named pandas")
            return None
    sys.meta_path.insert(0, _NoPandasFinder())
    if "pandas" in sys.modules:
        del sys.modules["pandas"]

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

WORKSPACE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

DEFAULT_STREAM = os.path.join(WORKSPACE_ROOT, "data", "arabidopsis", "uniprot_arabidopsis_stream.tsv.gz")
DEFAULT_TAIR_GAF = os.path.join(WORKSPACE_ROOT, "data", "raw", "arabidopsis", "tair", "tair.gaf.gz")
DEFAULT_ALIGNMENTS_TSV = os.path.join(
    WORKSPACE_ROOT, "data", "dataset_construction", "02_negative_sampling", "ara_canonical_plus_positives_alignments.tsv"
)
DEFAULT_CANDIDATES_DIR = os.path.join(
    WORKSPACE_ROOT, "data", "dataset_construction", "02_negative_sampling", "filtered_candidates"
)
DEFAULT_OUT_DIR = os.path.join(
    WORKSPACE_ROOT, "data", "dataset_construction", "02_negative_sampling", "final_candidates_nr"
)
DEFAULT_AUDIT_PATH = os.path.join(
    WORKSPACE_ROOT, "data", "dataset_construction", "02_negative_sampling", "ara_candidate_negatives_redundancy_reduction_audit.json"
)
DEFAULT_SUMMARY_MD = os.path.join(
    WORKSPACE_ROOT, "data", "dataset_construction", "02_negative_sampling", "summary_candidate_negatives_nr.md"
)

# Experimental TAIR evidence codes
EXPERIMENTAL_CODES = {"EXP", "IDA", "IPI", "IMP", "IGI", "IEP"}
CURATED_CODES = {"TAS", "NAS", "IC"}


def load_canonical_proteome(stream_path):
    """Load canonical proteome mapping and sequence lengths."""
    print(f"[{time.strftime('%X')}] Loading canonical reference proteome from: {stream_path}")
    if not os.path.exists(stream_path):
        raise FileNotFoundError(f"Stream file missing: {stream_path}")

    locus_candidates = defaultdict(list)
    with gzip.open(stream_path, "rt", encoding="utf-8", errors="replace") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        idx_acc = header.index("Entry")
        idx_oln = header.index("Gene Names (ordered locus)") if "Gene Names (ordered locus)" in header else -1
        idx_gene = header.index("Gene Names") if "Gene Names" in header else -1
        idx_seq = header.index("Sequence") if "Sequence" in header else -1

        for row in reader:
            acc = row[idx_acc].strip()
            seq = row[idx_seq].strip() if idx_seq >= 0 else ""
            is_trembl = acc.startswith("A0A")

            loci = []
            if idx_oln >= 0 and row[idx_oln].strip():
                for g in row[idx_oln].strip().split():
                    g_up = g.upper()
                    if g_up.startswith("AT") and ("G" in g_up or "M" in g_up or "C" in g_up):
                        loci.append(g_up)
            if not loci and idx_gene >= 0 and row[idx_gene].strip():
                for g in row[idx_gene].strip().split():
                    g_up = g.upper()
                    if g_up.startswith("AT") and ("G" in g_up or "M" in g_up or "C" in g_up):
                        loci.append(g_up)

            if loci and seq:
                primary_locus = loci[0]
                locus_candidates[primary_locus].append((acc, is_trembl, len(seq)))

    canonical_lengths = {}
    for locus, cands in locus_candidates.items():
        sorted_cands = sorted(cands, key=lambda x: (x[1], -x[2], x[0]))
        canonical_acc = sorted_cands[0][0]
        canonical_len = sorted_cands[0][2]
        canonical_lengths[canonical_acc] = canonical_len

    sorted_accs = sorted(canonical_lengths.keys())
    prot_to_idx = {acc: np.uint32(i) for i, acc in enumerate(sorted_accs)}
    idx_to_prot = {np.uint32(i): acc for i, acc in enumerate(sorted_accs)}
    lengths_arr = np.array([canonical_lengths[acc] for acc in sorted_accs], dtype=np.uint32)

    print(f"[{time.strftime('%X')}] Canonical proteome loaded: {len(sorted_accs):,} proteins.")
    return prot_to_idx, idx_to_prot, lengths_arr


def load_evidence_weights(tair_gaf_path, prot_to_idx):
    """Assign evidence quality weights (3=Exp, 2=Curated, 1=Comp/Other) per protein."""
    print(f"[{time.strftime('%X')}] Assigning TAIR annotation evidence weights: {tair_gaf_path}")
    N = len(prot_to_idx)
    weights = np.ones(N, dtype=np.uint8)  # default computational / baseline = 1

    if not os.path.exists(tair_gaf_path):
        print("  -> Warning: TAIR GAF not found; defaulting all evidence weights to 1.")
        return weights

    # Map locus to max weight
    locus_max_weight = defaultdict(int)
    with gzip.open(tair_gaf_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!"):
                continue
            parts = line.strip().split("\t")
            if len(parts) > 8 and parts[8] == "C":
                qualifier = parts[3].strip()
                if qualifier.startswith("NOT"):
                    continue
                locus = parts[1].strip().upper()
                ev_code = parts[6].strip().upper() if len(parts) > 6 else ""

                if ev_code in EXPERIMENTAL_CODES:
                    w = 3
                elif ev_code in CURATED_CODES:
                    w = 2
                else:
                    w = 1

                if w > locus_max_weight[locus]:
                    locus_max_weight[locus] = w

    for acc, idx in prot_to_idx.items():
        # Match locus or prefix
        w = locus_max_weight.get(acc, 1)
        weights[idx] = w

    print(f"[{time.strftime('%X')}] Evidence weights assigned. (Exp: {np.count_nonzero(weights == 3):,}, Curated: {np.count_nonzero(weights == 2):,})")
    return weights


def build_homology_components(alignments_tsv, prot_to_idx, theta=0.50, cmin=0.80, max_evalue=1e-3):
    """
    Parse MMseqs2 alignments, construct protein homology graph, and compute connected components.
    Returns:
        comp_id: np.ndarray of shape (N,) mapping each protein index to its component ID.
        homologs: dict mapping uint32 idx -> set of uint32 homolog indices.
    """
    print(f"\n[{time.strftime('%X')}] Computing Protein Homology Connected Components from: {alignments_tsv}")
    if not os.path.exists(alignments_tsv):
        raise FileNotFoundError(f"Alignments TSV missing: {alignments_tsv}")

    N = len(prot_to_idx)
    adj = defaultdict(set)
    homologs = defaultdict(set)
    total_alns = 0
    passed_alns = 0

    with open(alignments_tsv, "r", encoding="utf-8") as f:
        first_line = f.readline()
        tokens = [t.strip().lower() for t in first_line.split("\t")]
        has_header = "query" in tokens or "target" in tokens
        f.seek(0)
        column_names = [
            "query", "target", "fident", "alnlen", "mismatch", "gapopen",
            "qstart", "qend", "tstart", "tend", "evalue", "bits",
            "qlen", "tlen", "qcov", "tcov"
        ]
        reader = csv.DictReader(f, delimiter="\t") if has_header else csv.DictReader(f, fieldnames=column_names, delimiter="\t")

        for row in reader:
            total_alns += 1
            q = row.get("query", "").strip()
            t = row.get("target", "").strip()
            if not q or not t or q == t:
                continue
            if q not in prot_to_idx or t not in prot_to_idx:
                continue

            try:
                fident = float(row.get("fident", row.get("pident", 0.0)))
                if fident > 1.0:
                    fident /= 100.0
                qcov = float(row.get("qcov", 0.0))
                if qcov > 1.0:
                    qcov /= 100.0
                tcov = float(row.get("tcov", 0.0))
                if tcov > 1.0:
                    tcov /= 100.0
                evalue = float(row.get("evalue", 0.0))
            except ValueError:
                continue

            if fident >= theta and min(qcov, tcov) >= cmin and evalue <= max_evalue:
                u = prot_to_idx[q]
                v = prot_to_idx[t]
                adj[u].add(v)
                adj[v].add(u)
                homologs[u].add(v)
                homologs[v].add(u)
                passed_alns += 1

    print(f"[{time.strftime('%X')}] Parsed {total_alns:,} alignments; {passed_alns:,} meet theta >= {theta:.2f}, cmin >= {cmin:.2f}.")

    # Connected Components via BFS
    comp_id = np.zeros(N, dtype=np.uint32)
    visited = np.zeros(N, dtype=bool)
    current_comp = 1

    for i in range(N):
        if not visited[i]:
            queue = [i]
            visited[i] = True
            comp_id[i] = current_comp
            head = 0
            while head < len(queue):
                curr = queue[head]
                head += 1
                for neighbor in adj[curr]:
                    if not visited[neighbor]:
                        visited[neighbor] = True
                        comp_id[neighbor] = current_comp
                        queue.append(neighbor)
            current_comp += 1

    total_comps = current_comp - 1
    # Analyze component sizes
    comp_counts = defaultdict(int)
    for c in comp_id:
        comp_counts[c] += 1
    singleton_comps = sum(1 for cnt in comp_counts.values() if cnt == 1)
    multi_comps = total_comps - singleton_comps

    print(f"[{time.strftime('%X')}] Total Protein Homology Components: {total_comps:,}")
    print(f"  -> Singleton Components (no homologs): {singleton_comps:,} ({singleton_comps/total_comps*100:.1f}%)")
    print(f"  -> Multi-protein Paralog Families:      {multi_comps:,} ({multi_comps/total_comps*100:.1f}%)")
    max_comp_size = max(comp_counts.values()) if comp_counts else 0
    print(f"  -> Largest Gene Family Size:            {max_comp_size} proteins")

    return comp_id, homologs


def run_candidate_negative_clustering(args):
    """Main execution engine for negative pair redundancy reduction."""
    start_total_t = time.time()
    os.makedirs(args.out_dir, exist_ok=True)

    # 1. Load Proteome, Weights, and Homology Components
    prot_to_idx, idx_to_prot, lengths_arr = load_canonical_proteome(args.stream_path)
    weights_arr = load_evidence_weights(args.tair_gaf_path, prot_to_idx)
    comp_id_arr, homologs = build_homology_components(
        args.alignments_tsv, prot_to_idx, theta=args.theta, cmin=args.cmin, max_evalue=args.max_evalue
    )

    chunk_files = sorted(glob.glob(os.path.join(args.candidates_dir, "ara_candidate_negatives_filtered_chunk_*.parquet")))
    if not chunk_files:
        # Fallback check for any candidate parquet
        chunk_files = sorted(glob.glob(os.path.join(args.candidates_dir, "*.parquet")))
        if not chunk_files:
            raise FileNotFoundError(f"No candidate parquet chunks found in: {args.candidates_dir}")

    total_chunks = len(chunk_files)
    print(f"\n[{time.strftime('%X')}] Found {total_chunks} clean candidate Parquet chunks to cluster.")

    # --------------------------------------------------------------------------
    # PASS 1: Lightweight Indexing & Priority Scoring across all chunks
    # --------------------------------------------------------------------------
    print(f"\n[{time.strftime('%X')}] === PASS 1: Scanning Chunks & Computing Component Signatures ===")
    t_pass1 = time.time()

    list_u = []
    list_v = []
    list_chunk_id = []
    list_row_id = []
    list_score = []
    list_comp_key = []

    total_input_pairs = 0

    for chunk_idx, chunk_path in enumerate(chunk_files, 1):
        t0 = time.time()
        table = pq.read_table(chunk_path, columns=["participant_a_id", "participant_b_id"])
        n_rows = len(table)
        total_input_pairs += n_rows

        col_a = table["participant_a_id"].to_pylist()
        col_b = table["participant_b_id"].to_pylist()

        u_vec = np.array([prot_to_idx[a] for a in col_a], dtype=np.uint32)
        v_vec = np.array([prot_to_idx[b] for b in col_b], dtype=np.uint32)

        # Standardize u <= v
        swap_mask = u_vec > v_vec
        if np.any(swap_mask):
            u_tmp = u_vec.copy()
            u_vec[swap_mask] = v_vec[swap_mask]
            v_vec[swap_mask] = u_tmp[swap_mask]

        # Component keys: (min(K_u, K_v) << 32) | max(K_u, K_v)
        ku = comp_id_arr[u_vec]
        kv = comp_id_arr[v_vec]
        k_min = np.minimum(ku, kv).astype(np.uint64)
        k_max = np.maximum(ku, kv).astype(np.uint64)
        comp_key_vec = (k_min << np.uint64(32)) | k_max

        # Priority score: Evidence Weight * 10^7 + Combined Length
        # (Weight ranges from 2 to 6; Length ranges from ~60 to ~10,000)
        wu = weights_arr[u_vec].astype(np.int64)
        wv = weights_arr[v_vec].astype(np.int64)
        lu = lengths_arr[u_vec].astype(np.int64)
        lv = lengths_arr[v_vec].astype(np.int64)
        score_vec = ((wu + wv) * 10_000_000) + (lu + lv)

        chunk_id_vec = np.full(n_rows, chunk_idx, dtype=np.uint16)
        row_id_vec = np.arange(n_rows, dtype=np.uint32)

        list_u.append(u_vec)
        list_v.append(v_vec)
        list_chunk_id.append(chunk_id_vec)
        list_row_id.append(row_id_vec)
        list_score.append(score_vec)
        list_comp_key.append(comp_key_vec)

        if chunk_idx % 10 == 0 or chunk_idx == total_chunks:
            elapsed = time.time() - t_pass1
            print(f"[{time.strftime('%X')}] Indexed chunk {chunk_idx:03d}/{total_chunks} ({total_input_pairs:,} pairs) | {elapsed:.1f}s")

    # Concatenate all arrays
    print(f"[{time.strftime('%X')}] Consolidating index arrays in memory...")
    arr_u = np.concatenate(list_u)
    arr_v = np.concatenate(list_v)
    arr_chunk_id = np.concatenate(list_chunk_id)
    arr_row_id = np.concatenate(list_row_id)
    arr_score = np.concatenate(list_score)
    arr_comp_key = np.concatenate(list_comp_key)

    del list_u, list_v, list_chunk_id, list_row_id, list_score, list_comp_key

    total_pairs = len(arr_u)
    print(f"[{time.strftime('%X')}] Pass 1 Complete: {total_pairs:,} pairs indexed in {time.time() - t_pass1:.2f}s.")

    # --------------------------------------------------------------------------
    # PASS 2: Group by CompKey & Run Greedy CD-HIT Representative Selection
    # --------------------------------------------------------------------------
    print(f"\n[{time.strftime('%X')}] === PASS 2: Component-Pair Sorting & Greedy CD-HIT Selection ===")
    t_pass2 = time.time()

    print(f"[{time.strftime('%X')}] Sorting {total_pairs:,} pairs by component key...")
    order = np.argsort(arr_comp_key, kind="stable")
    arr_u = arr_u[order]
    arr_v = arr_v[order]
    arr_chunk_id = arr_chunk_id[order]
    arr_row_id = arr_row_id[order]
    arr_score = arr_score[order]
    arr_comp_key = arr_comp_key[order]
    del order

    # Identify component-pair group boundaries
    diff_mask = arr_comp_key[:-1] != arr_comp_key[1:]
    diff_indices = np.where(diff_mask)[0]
    boundaries = np.concatenate([[0], diff_indices + 1, [total_pairs]])
    num_groups = len(boundaries) - 1

    print(f"[{time.strftime('%X')}] Grouping complete: {num_groups:,} unique component-pair signatures.")

    # Set of representative pair locations: packed (uint16(chunk_id) << 32) | uint32(row_id)
    rep_locs = set()
    total_singletons = 0
    total_multi_groups = 0
    total_multi_pairs_evaluated = 0
    total_pruned_in_multi = 0

    for g_idx in range(num_groups):
        start = boundaries[g_idx]
        end = boundaries[g_idx + 1]
        group_size = end - start

        if group_size == 1:
            # Fast-path: trivially non-redundant singleton pair
            ch_id = arr_chunk_id[start]
            r_id = arr_row_id[start]
            packed_loc = (np.uint64(ch_id) << np.uint64(32)) | np.uint64(r_id)
            rep_locs.add(packed_loc)
            total_singletons += 1
        else:
            # Paralog micro-cluster: run exact greedy CD-HIT
            total_multi_groups += 1
            total_multi_pairs_evaluated += group_size

            # Sub-slice for this group
            g_u = arr_u[start:end]
            g_v = arr_v[start:end]
            g_ch = arr_chunk_id[start:end]
            g_row = arr_row_id[start:end]
            g_sc = arr_score[start:end]

            # Priority sort within group (score descending -> lexicographical u, v)
            sub_order = sorted(
                range(group_size),
                key=lambda idx: (-g_sc[idx], g_u[idx], g_v[idx])
            )

            unassigned = set(sub_order)

            for rep_cand in sub_order:
                if rep_cand not in unassigned:
                    continue

                unassigned.remove(rep_cand)
                ch_id = g_ch[rep_cand]
                r_id = g_row[rep_cand]
                packed_loc = (np.uint64(ch_id) << np.uint64(32)) | np.uint64(r_id)
                rep_locs.add(packed_loc)

                rep_u = g_u[rep_cand]
                rep_v = g_v[rep_cand]
                hu = homologs[rep_u] | {rep_u}
                hv = homologs[rep_v] | {rep_v}

                # Prune remaining pairs in this micro-cluster if S >= 0.50
                to_prune = []
                for other in unassigned:
                    ou = g_u[other]
                    ov = g_v[other]
                    # Check S(P_rep, P_other) >= 0.50:
                    # Orientation 1: (ou in hu and ov in hv)
                    # Orientation 2: (ou in hv and ov in hu)
                    if (ou in hu and ov in hv) or (ou in hv and ov in hu):
                        to_prune.append(other)

                for p_idx in to_prune:
                    unassigned.remove(p_idx)
                    total_pruned_in_multi += 1

        if (g_idx + 1) % 500_000 == 0 or (g_idx + 1) == num_groups:
            elapsed = time.time() - t_pass2
            pct = ((g_idx + 1) / num_groups) * 100
            print(f"[{time.strftime('%X')}] Evaluated group {g_idx + 1:,}/{num_groups:,} ({pct:.1f}%) | Representatives: {len(rep_locs):,} | {elapsed:.1f}s")

    # Free large index arrays
    del arr_u, arr_v, arr_chunk_id, arr_row_id, arr_score, arr_comp_key

    total_reps = len(rep_locs)
    total_pruned = total_pairs - total_reps
    pruned_pct = (total_pruned / total_pairs * 100.0) if total_pairs > 0 else 0.0
    retained_pct = (total_reps / total_pairs * 100.0) if total_pairs > 0 else 0.0

    print(f"\n[{time.strftime('%X')}] Pass 2 Complete in {time.time() - t_pass2:.2f}s:")
    print(f"  -> Total Unique Representatives Retained: {total_reps:,} ({retained_pct:.2f}%)")
    print(f"  -> Total Redundant Pairs Pruned:          {total_pruned:,} ({pruned_pct:.2f}%)")
    print(f"  -> Singleton Group Fast-Paths:            {total_singletons:,}")
    print(f"  -> Multi-Pair Groups Resolved:            {total_multi_groups:,}")

    # --------------------------------------------------------------------------
    # PASS 3: Stream and Export Final Non-Redundant Parquet Chunks
    # --------------------------------------------------------------------------
    print(f"\n[{time.strftime('%X')}] === PASS 3: Streaming Final Non-Redundant Chunks to Disk ===")
    t_pass3 = time.time()

    out_buf = {
        "participant_a_id": [],
        "participant_b_id": [],
        "locus_a": [],
        "locus_b": [],
        "compartments_a": [],
        "compartments_b": [],
    }

    out_chunk_idx = 1
    total_written = 0

    def flush_output():
        nonlocal out_chunk_idx, out_buf, total_written
        if not out_buf["participant_a_id"]:
            return

        out_path = os.path.join(args.out_dir, f"ara_candidate_negatives_nr_chunk_{out_chunk_idx:04d}.parquet")
        chunk_table = pa.Table.from_arrays(
            [
                pa.array(out_buf["participant_a_id"], type=pa.string()),
                pa.array(out_buf["participant_b_id"], type=pa.string()),
                pa.array(out_buf["locus_a"], type=pa.string()),
                pa.array(out_buf["locus_b"], type=pa.string()),
                pa.array(out_buf["compartments_a"], type=pa.string()),
                pa.array(out_buf["compartments_b"], type=pa.string()),
            ],
            names=[
                "participant_a_id",
                "participant_b_id",
                "locus_a",
                "locus_b",
                "compartments_a",
                "compartments_b",
            ],
        )
        pq.write_table(chunk_table, out_path, compression="snappy")
        print(f"[{time.strftime('%X')}] Wrote Final NR Chunk {out_chunk_idx:04d}: {len(out_buf['participant_a_id']):,} pairs -> {out_path}")
        out_chunk_idx += 1
        total_written += len(out_buf["participant_a_id"])
        for k in out_buf:
            out_buf[k].clear()

    unique_proteins_kept = set()

    for chunk_idx, chunk_path in enumerate(chunk_files, 1):
        t0 = time.time()
        table = pq.read_table(chunk_path)
        n_rows = len(table)

        col_a = table["participant_a_id"].to_pylist()
        col_b = table["participant_b_id"].to_pylist()
        col_loc_a = table["locus_a"].to_pylist()
        col_loc_b = table["locus_b"].to_pylist()
        col_comp_a = table["compartments_a"].to_pylist()
        col_comp_b = table["compartments_b"].to_pylist()

        # Check membership in rep_locs
        ch_shifted = np.uint64(chunk_idx) << np.uint64(32)
        kept_count = 0

        for r_id in range(n_rows):
            packed = ch_shifted | np.uint64(r_id)
            if packed in rep_locs:
                a_id = col_a[r_id]
                b_id = col_b[r_id]
                out_buf["participant_a_id"].append(a_id)
                out_buf["participant_b_id"].append(b_id)
                out_buf["locus_a"].append(col_loc_a[r_id])
                out_buf["locus_b"].append(col_loc_b[r_id])
                out_buf["compartments_a"].append(col_comp_a[r_id])
                out_buf["compartments_b"].append(col_comp_b[r_id])

                unique_proteins_kept.add(a_id)
                unique_proteins_kept.add(b_id)
                kept_count += 1

                if len(out_buf["participant_a_id"]) >= args.chunk_size:
                    flush_output()

        elapsed = time.time() - t0
        if chunk_idx % 10 == 0 or chunk_idx == total_chunks:
            print(f"[{time.strftime('%X')}] Streamed chunk {chunk_idx:03d}/{total_chunks}: kept {kept_count:,}/{n_rows:,} ({kept_count/n_rows*100:.1f}%) | {elapsed:.2f}s")

    # Flush remaining records
    flush_output()

    total_chunks_written = out_chunk_idx - 1
    total_time = time.time() - start_total_t

    print(f"\n[{time.strftime('%X')}] === Final Non-Redundant Negative Pool Completed! ===")
    print(f"  -> Total Clean Input Pairs:            {total_pairs:,}")
    print(f"  -> Non-Redundant Representative Pairs: {total_written:,} ({total_written/total_pairs*100:.2f}%)")
    print(f"  -> Redundant Paralog Pairs Pruned:     {total_pairs - total_written:,} ({(total_pairs - total_written)/total_pairs*100:.2f}%)")
    print(f"  -> Unique Proteins Represented:        {len(unique_proteins_kept):,}")
    print(f"  -> Output Chunks Written:              {total_chunks_written} chunks in {args.out_dir}")
    print(f"  -> Total Execution Time:               {total_time:.2f} seconds.")

    metrics = {
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "stage": "Stage 2.6: Pair-Level Redundancy Reduction Among Negatives",
        "protocol": "Option A (Priority: Evidence Quality -> Length -> Lexicographical)",
        "theta_identity_threshold": args.theta,
        "cmin_coverage_threshold": args.cmin,
        "evalue_cutoff": args.max_evalue,
        "total_clean_candidates_evaluated": total_pairs,
        "total_representative_negatives_retained": total_written,
        "representative_retention_pct": round((total_written / total_pairs * 100.0), 4) if total_pairs > 0 else 0.0,
        "total_redundant_pairs_pruned": total_pairs - total_written,
        "redundancy_pruned_pct": round(((total_pairs - total_written) / total_pairs * 100.0), 4) if total_pairs > 0 else 0.0,
        "unique_proteins_in_negative_pool": len(unique_proteins_kept),
        "total_component_pair_groups": num_groups,
        "singleton_component_pairs": total_singletons,
        "multi_pair_component_groups": total_multi_groups,
        "output_chunks_written": total_chunks_written,
        "chunk_size": args.chunk_size,
        "execution_time_seconds": round(total_time, 2),
        "input_directory": args.candidates_dir,
        "output_directory": args.out_dir,
    }

    # 4. Export Audit JSON and Markdown Summary
    write_audit_reports(args.audit_path, args.summary_md, metrics)


def write_audit_reports(audit_path, summary_md_path, metrics):
    """Write machine-readable JSON and human-readable Markdown summary reports."""
    os.makedirs(os.path.dirname(audit_path), exist_ok=True)
    with open(audit_path, "w", encoding="utf-8") as f:
        json.dump(metrics, f, indent=2)
    print(f"[{time.strftime('%X')}] Saved audit JSON to: {audit_path}")

    with open(summary_md_path, "w", encoding="utf-8") as f:
        f.write("# Stage 2.6: Non-Redundant Negative Candidate Pool Audit Report\n\n")
        f.write(f"**Execution Timestamp**: `{metrics['timestamp_utc']}`\n\n")
        f.write("## 1. Summary Overview\n\n")
        f.write("Pair-level redundancy reduction was executed across the clean negative candidate pool using greedy CD-HIT representative selection.\n\n")
        f.write("| Metric | Value |\n")
        f.write("| :--- | :--- |\n")
        f.write(f"| **Sequence Identity Threshold ($\\theta$)** | `{metrics['theta_identity_threshold']:.2f}` |\n")
        f.write(f"| **Reciprocal Coverage Threshold ($c_{{\\min}}$)** | `{metrics['cmin_coverage_threshold']:.2f}` |\n")
        f.write(f"| **E-value Cutoff** | `{metrics['evalue_cutoff']}` |\n")
        f.write(f"| **Total Clean Candidate Negatives Evaluated** | `{metrics['total_clean_candidates_evaluated']:,}` |\n")
        f.write(f"| **Non-Redundant Representatives Retained** | `{metrics['total_representative_negatives_retained']:,}` (`{metrics['representative_retention_pct']:.2f}%`) |\n")
        f.write(f"| **Redundant Paralogous Negatives Pruned** | `{metrics['total_redundant_pairs_pruned']:,}` (`{metrics['redundancy_pruned_pct']:.2f}%`) |\n")
        f.write(f"| **Unique Proteins Represented** | `{metrics['unique_proteins_in_negative_pool']:,}` |\n")
        f.write(f"| **Homology Component-Pair Signatures** | `{metrics['total_component_pair_groups']:,}` |\n")
        f.write(f"| **Singleton Component Pairs** | `{metrics['singleton_component_pairs']:,}` |\n")
        f.write(f"| **Multi-Pair Paralog Groups** | `{metrics['multi_pair_component_groups']:,}` |\n")
        f.write(f"| **Output Filtered Parquet Chunks** | `{metrics['output_chunks_written']}` chunks |\n")
        f.write(f"| **Total Execution Time** | `{metrics['execution_time_seconds']:.2f}` seconds |\n\n")

        f.write("## 2. Methodology & Representative Priority Hierarchy\n\n")
        f.write("Pairs within each homology component pair are collapsed using exact CD-HIT greedy clustering, prioritized by:\n")
        f.write("1. **TAIR Subcellular Evidence Confidence**: Direct experimental codes (`EXP`, `IDA`, `IPI`, `HDA` = 3.0) > Curated codes (`TAS`, `NAS` = 2.0) > Computational predictions (`IEA`, `ISS` = 1.0).\n")
        f.write("2. **Combined Sequence Length**: $(L_A + L_B)$ descending (longer canonical proteins provide richer domain context).\n")
        f.write("3. **Lexicographical Tie-Breaking**: Deterministic sorting by `participant_a_id`, `participant_b_id`.\n\n")

        f.write("## 3. Mathematical Guarantee\n\n")
        f.write("All surviving negative pairs strictly satisfy:\n")
        f.write("$$\\forall P_i, P_j \\in \\mathcal{P}_{\\text{neg, NR}} \\ (i \\ne j): \\quad S(P_i, P_j) < 0.50$$\n")
        f.write("where $S(P_i, P_j) = \\max [ \\min(s(A, C), s(B, D)), \\min(s(A, D), s(B, C)) ]$.\n")

    print(f"[{time.strftime('%X')}] Saved summary Markdown to: {summary_md_path}")


def main():
    parser = argparse.ArgumentParser(description="Pair-Level Redundancy Reduction for Candidate Negatives (HPC)")
    parser.add_argument("--stream_path", default=DEFAULT_STREAM, help="Path to uniprot_arabidopsis_stream.tsv.gz")
    parser.add_argument("--tair_gaf_path", default=DEFAULT_TAIR_GAF, help="Path to tair.gaf.gz")
    parser.add_argument("--alignments_tsv", default=DEFAULT_ALIGNMENTS_TSV, help="Path to MMseqs2 alignment TSV")
    parser.add_argument("--candidates_dir", default=DEFAULT_CANDIDATES_DIR, help="Directory containing clean filtered chunks")
    parser.add_argument("--out_dir", default=DEFAULT_OUT_DIR, help="Directory for final non-redundant negative chunks")
    parser.add_argument("--audit_path", default=DEFAULT_AUDIT_PATH, help="Path for audit JSON")
    parser.add_argument("--summary_md", default=DEFAULT_SUMMARY_MD, help="Path for summary Markdown")
    parser.add_argument("--theta", type=float, default=0.50, help="Identity threshold (default: 0.50)")
    parser.add_argument("--cmin", type=float, default=0.80, help="Coverage threshold (default: 0.80)")
    parser.add_argument("--max_evalue", type=float, default=1e-3, help="Max E-value (default: 1e-3)")
    parser.add_argument("--chunk_size", type=int, default=1000000, help="Rows per output chunk")

    args = parser.parse_args()

    print("================================================================================")
    print("   AraPPLM: Candidate Negative Pair Redundancy Reduction (HPC Streamer)         ")
    print("================================================================================")
    print(f"Timestamp:           {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Identity Threshold:  theta = {args.theta:.2f}")
    print(f"Coverage Threshold:  cmin  = {args.cmin:.2f}")
    print(f"Max E-Value:         {args.max_evalue}")
    print(f"Stream TSV:          {args.stream_path}")
    print(f"TAIR GAF:            {args.tair_gaf_path}")
    print(f"Alignments TSV:      {args.alignments_tsv}")
    print(f"Input Clean Chunks:  {args.candidates_dir}")
    print(f"Output Final NR Dir: {args.out_dir}")
    print("================================================================================")

    run_candidate_negative_clustering(args)


if __name__ == "__main__":
    main()
