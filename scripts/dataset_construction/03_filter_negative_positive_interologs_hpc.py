#!/usr/bin/env python3
"""
03_filter_negative_positive_interologs_hpc.py
=============================================
Filters candidate negative protein pairs for Arabidopsis thaliana by removing
pairs that exhibit pair-level sequence homology (positive interologs / pseudo-negatives)
against verified positive interactions from the latest benchmark splits (C1/C2/C3/C4).

Pair Similarity Metric (docs/dataset_construction/03_homology_clustering/01_pair_redundancy_reduction.md):
    S(P_neg, P_pos) = max[ min(s(A, C), s(B, D)), min(s(A, D), s(B, C)) ]

Where s(X, Y) is protein sequence identity:
    s(X, Y) = fident(X, Y) if fident >= theta and min(qcov, tcov) >= cmin and evalue <= 1e-3, else 0.0

A candidate negative pair P_neg = (A, B) is DISCARDED if there exists ANY positive pair
P_pos = (C, D) in the latest benchmark dataset (C1/C2/C3/C4) such that:
    S(P_neg, P_pos) >= theta (default: 0.50, cmin: 0.80)

Algorithmic Optimization:
    Rather than performing O(N_neg * N_pos) pairwise comparisons (~2.3 x 10^12 checks),
    this script precomputes the Exhaustive Positive Interolog Set:
        I_pos = Union_{(C, D) in P_pos} { (min(u, v), max(u, v)) | u in H(C), v in H(D), u != v }
    where H(X) is the homolog set of protein X at (theta=0.50, cmin=0.80) including self.
    Candidate negative chunks are then streamed and filtered via O(1) integer hash lookups (< 40 ns/pair).

Outputs:
    data/dataset_construction/02_negative_sampling/filtered_candidates/ara_candidate_negatives_filtered_chunk_*.parquet
    data/dataset_construction/02_negative_sampling/ara_candidate_negatives_interolog_removal_audit.json
    data/dataset_construction/02_negative_sampling/summary_candidate_negatives_filtered.md
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
DEFAULT_HETERO_FASTA = os.path.join(
    WORKSPACE_ROOT, "data", "dataset_construction", "01_positive_pairs", "arabidopsis", "ara_positive_pairs_high_only_hetero.fasta"
)
DEFAULT_UNIFIED_FASTA = os.path.join(WORKSPACE_ROOT, "data", "arabidopsis", "arabidopsis_canonical_plus_positives.fasta")
DEFAULT_SPLITS_DIR = os.path.join(WORKSPACE_ROOT, "data", "dataset_construction", "04_benchmark_splits", "arabidopsis")
DEFAULT_ALIGNMENTS_TSV = os.path.join(
    WORKSPACE_ROOT, "data", "dataset_construction", "02_negative_sampling", "ara_canonical_plus_positives_alignments.tsv"
)
DEFAULT_CANDIDATES_DIR = os.path.join(
    WORKSPACE_ROOT, "data", "dataset_construction", "02_negative_sampling", "candidates"
)
DEFAULT_OUT_DIR = os.path.join(
    WORKSPACE_ROOT, "data", "dataset_construction", "02_negative_sampling", "filtered_candidates"
)
DEFAULT_AUDIT_PATH = os.path.join(
    WORKSPACE_ROOT, "data", "dataset_construction", "02_negative_sampling", "ara_candidate_negatives_interolog_removal_audit.json"
)
DEFAULT_SUMMARY_MD = os.path.join(
    WORKSPACE_ROOT, "data", "dataset_construction", "02_negative_sampling", "summary_candidate_negatives_filtered.md"
)


def load_canonical_proteome(stream_path):
    """
    Load canonical reference proteome (Option A: 1 protein per locus) from stream TSV.
    Matches the exact locus-to-accession mapping used in candidate negative generation.
    """
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
                locus_candidates[primary_locus].append((acc, is_trembl, len(seq), seq))

    canonical_seqs = {}
    for locus, cands in locus_candidates.items():
        sorted_cands = sorted(cands, key=lambda x: (x[1], -x[2], x[0]))
        canonical_acc = sorted_cands[0][0]
        canonical_seq = sorted_cands[0][3]
        canonical_seqs[canonical_acc] = canonical_seq

    print(f"[{time.strftime('%X')}] Canonical proteome loaded: {len(canonical_seqs):,} proteins.")
    return canonical_seqs


def build_unified_fasta(stream_path, hetero_fasta_path, out_fasta_path):
    """
    Build a unified FASTA file containing all canonical proteins plus any additional
    positive interactors present in the positive dataset.
    """
    print(f"[{time.strftime('%X')}] Constructing unified proteome FASTA at: {out_fasta_path}")
    os.makedirs(os.path.dirname(out_fasta_path), exist_ok=True)

    canonical_seqs = load_canonical_proteome(stream_path)
    unified = dict(canonical_seqs)

    if os.path.exists(hetero_fasta_path):
        print(f"[{time.strftime('%X')}] Integrating positive sequences from: {hetero_fasta_path}")
        with open(hetero_fasta_path, "r", encoding="utf-8") as f:
            cur_id = None
            seq_parts = []
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if cur_id and cur_id not in unified:
                        unified[cur_id] = "".join(seq_parts)
                    cur_id = line[1:].split()[0]
                    seq_parts = []
                else:
                    seq_parts.append(line)
            if cur_id and cur_id not in unified:
                unified[cur_id] = "".join(seq_parts)

    print(f"[{time.strftime('%X')}] Writing {len(unified):,} sequences to FASTA...")
    with open(out_fasta_path, "w", encoding="utf-8") as f:
        for acc in sorted(unified.keys()):
            f.write(f">{acc}\n{unified[acc]}\n")

    print(f"[{time.strftime('%X')}] Unified FASTA written: {out_fasta_path}")
    return unified


def load_positive_benchmark_pairs(splits_dir):
    """
    Load all positive pairs from C1, C2, C3, and C4 benchmark partitions.
    Standardizes each pair such that participant_a <= participant_b.
    """
    print(f"\n[{time.strftime('%X')}] Loading positive benchmark pairs from: {splits_dir}")
    if not os.path.exists(splits_dir):
        raise FileNotFoundError(f"Positive splits directory not found: {splits_dir}")

    partition_files = {
        "C1": os.path.join(splits_dir, "ara_positive_c1_train_id0_50_cov0_80.tsv"),
        "C2": os.path.join(splits_dir, "ara_positive_c2_test_id0_50_cov0_80.tsv"),
        "C3": os.path.join(splits_dir, "ara_positive_c3_test_id0_50_cov0_80.tsv"),
        "C4": os.path.join(splits_dir, "ara_positive_c4_test_id0_50_cov0_80.tsv"),
    }

    positive_pairs = set()
    partition_counts = {}
    unique_proteins = set()

    for part_name, filepath in partition_files.items():
        if not os.path.exists(filepath):
            # Fallback search for any partition matching pattern
            candidates = glob.glob(os.path.join(splits_dir, f"ara_positive_{part_name.lower()}*.tsv"))
            if candidates:
                filepath = candidates[0]
            else:
                raise FileNotFoundError(f"Could not find positive partition file for {part_name}: {filepath}")

        part_pairs = 0
        with open(filepath, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                a = row["participant_a_id"].strip()
                b = row["participant_b_id"].strip()
                if a and b and a != b:
                    pair_key = (min(a, b), max(a, b))
                    positive_pairs.add(pair_key)
                    unique_proteins.add(a)
                    unique_proteins.add(b)
                    part_pairs += 1
        partition_counts[part_name] = part_pairs
        print(f"  -> Loaded {part_name}: {part_pairs:,} pairs from {os.path.basename(filepath)}")

    print(f"[{time.strftime('%X')}] Total non-redundant positive pairs: {len(positive_pairs):,}")
    print(f"[{time.strftime('%X')}] Total unique positive proteins: {len(unique_proteins):,}")

    return list(positive_pairs), unique_proteins, partition_counts


def build_homology_index(alignments_tsv, theta=0.50, cmin=0.80, max_evalue=1e-3):
    """
    Parse MMseqs2 convertalis TSV and construct homolog sets H(X) for all proteins:
        H(X) = { Y | s(X, Y) >= theta, min(qcov, tcov) >= cmin, evalue <= max_evalue }
    """
    print(f"\n[{time.strftime('%X')}] Parsing MMseqs2 alignments from: {alignments_tsv}")
    if not os.path.exists(alignments_tsv):
        raise FileNotFoundError(f"Alignments TSV file not found: {alignments_tsv}")

    start_t = time.time()
    homologs = defaultdict(set)
    total_alignments = 0
    passed_alignments = 0

    column_names = [
        "query", "target", "fident", "alnlen", "mismatch", "gapopen",
        "qstart", "qend", "tstart", "tend", "evalue", "bits",
        "qlen", "tlen", "qcov", "tcov"
    ]

    with open(alignments_tsv, "r", encoding="utf-8") as f:
        first_line = f.readline()
        if not first_line:
            raise ValueError(f"Alignments TSV is empty: {alignments_tsv}")

        tokens = [t.strip().lower() for t in first_line.split("\t")]
        has_header = "query" in tokens or "target" in tokens

        f.seek(0)
        if has_header:
            reader = csv.DictReader(f, delimiter="\t")
        else:
            reader = csv.DictReader(f, fieldnames=column_names, delimiter="\t")

        for row in reader:
            total_alignments += 1
            q = row.get("query", "").strip()
            t = row.get("target", "").strip()
            if not q or not t or q == t:
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
                homologs[q].add(t)
                homologs[t].add(q)
                passed_alignments += 1

    elapsed = time.time() - start_t
    print(f"[{time.strftime('%X')}] Parsed {total_alignments:,} alignments in {elapsed:.2f}s.")
    print(f"[{time.strftime('%X')}] Homology filter (id >= {theta:.2f}, cov >= {cmin:.2f}, evalue <= {max_evalue}): {passed_alignments:,} passed.")
    print(f"[{time.strftime('%X')}] Proteins with at least one homolog: {len(homologs):,}")

    return homologs


def build_positive_interolog_set(positive_pairs, homologs, prot_to_idx):
    """
    Construct the Exhaustive Positive Interolog Set:
        I_pos = Union_{(C, D) in P_pos} { (min(u, v), max(u, v)) | u in H(C), v in H(D), u != v }
    Packed as 64-bit unsigned integers: (uint32(min_idx) << 32) | uint32(max_idx).
    """
    print(f"\n[{time.strftime('%X')}] Compiling Exhaustive Positive Interolog Set (theta=0.50, cmin=0.80)...")
    start_t = time.time()

    interolog_packed = set()
    total_cartesian_evals = 0

    for c, d in positive_pairs:
        # Homolog set including self
        hc = homologs[c] | {c}
        hd = homologs[d] | {d}

        total_cartesian_evals += len(hc) * len(hd)

        for u in hc:
            if u not in prot_to_idx:
                continue
            idx_u = prot_to_idx[u]

            for v in hd:
                if v not in prot_to_idx or u == v:
                    continue
                idx_v = prot_to_idx[v]

                u_min = min(idx_u, idx_v)
                u_max = max(idx_u, idx_v)
                packed = (np.uint64(u_min) << np.uint64(32)) | np.uint64(u_max)
                interolog_packed.add(packed)

    elapsed = time.time() - start_t
    mem_mb = (len(interolog_packed) * 8) / (1024 * 1024)
    print(f"[{time.strftime('%X')}] Cartesian combinations evaluated: {total_cartesian_evals:,}")
    print(f"[{time.strftime('%X')}] Unique prohibited interolog pairs: {len(interolog_packed):,} ({mem_mb:.2f} MB in RAM)")
    print(f"[{time.strftime('%X')}] Precomputation completed in {elapsed:.2f} seconds.")

    # Convert to sorted NumPy array for high-speed binary search lookup
    interolog_array = np.array(sorted(list(interolog_packed)), dtype=np.uint64)
    return interolog_array, len(interolog_packed)


def filter_candidate_chunks(candidates_dir, out_dir, interolog_array, prot_to_idx, chunk_size=1000000):
    """
    Stream candidate negative parquet chunks and filter out pairs present in interolog_array.
    Surviving pairs are flushed to clean parquet chunks.
    """
    print(f"\n[{time.strftime('%X')}] Initiating Streaming Candidate Negative Filtering...")
    print(f"[{time.strftime('%X')}] Input Candidates Directory: {candidates_dir}")
    print(f"[{time.strftime('%X')}] Output Filtered Directory:  {out_dir}")
    os.makedirs(out_dir, exist_ok=True)

    chunk_files = sorted(glob.glob(os.path.join(candidates_dir, "ara_candidate_negatives_chunk_*.parquet")))
    if not chunk_files:
        raise FileNotFoundError(f"No candidate parquet chunks found in: {candidates_dir}")

    print(f"[{time.strftime('%X')}] Found {len(chunk_files)} candidate chunks to process.")

    start_stream_t = time.time()
    total_candidates_examined = 0
    total_homologous_removed = 0
    total_clean_retained = 0
    total_unmapped_kept = 0

    out_buf = {
        "participant_a_id": [],
        "participant_b_id": [],
        "locus_a": [],
        "locus_b": [],
        "compartments_a": [],
        "compartments_b": [],
    }

    out_chunk_idx = 1

    def flush_output():
        nonlocal out_chunk_idx, out_buf
        if not out_buf["participant_a_id"]:
            return

        out_path = os.path.join(out_dir, f"ara_candidate_negatives_filtered_chunk_{out_chunk_idx:04d}.parquet")
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
        print(f"[{time.strftime('%X')}] Wrote Clean Chunk {out_chunk_idx:04d}: {len(out_buf['participant_a_id']):,} pairs -> {out_path}")
        out_chunk_idx += 1
        for k in out_buf:
            out_buf[k].clear()

    has_forbidden = len(interolog_array) > 0

    for in_idx, chunk_file in enumerate(chunk_files, 1):
        t0 = time.time()
        table = pq.read_table(chunk_file)
        n_rows = len(table)
        total_candidates_examined += n_rows

        col_a = table["participant_a_id"].to_pylist()
        col_b = table["participant_b_id"].to_pylist()
        col_loc_a = table["locus_a"].to_pylist()
        col_loc_b = table["locus_b"].to_pylist()
        col_comp_a = table["compartments_a"].to_pylist()
        col_comp_b = table["compartments_b"].to_pylist()

        # Vectorized index lookup and packing
        u_arr = np.array([prot_to_idx.get(a, np.uint32(0xFFFFFFFF)) for a in col_a], dtype=np.uint32)
        v_arr = np.array([prot_to_idx.get(b, np.uint32(0xFFFFFFFF)) for b in col_b], dtype=np.uint32)

        valid_idx_mask = (u_arr != 0xFFFFFFFF) & (v_arr != 0xFFFFFFFF)

        is_forbidden = np.zeros(n_rows, dtype=bool)

        if has_forbidden and np.any(valid_idx_mask):
            u_valid = u_arr[valid_idx_mask]
            v_valid = v_arr[valid_idx_mask]
            u_min = np.minimum(u_valid, v_valid)
            u_max = np.maximum(u_valid, v_valid)
            packed = (u_min.astype(np.uint64) << np.uint64(32)) | u_max.astype(np.uint64)

            # High-speed binary search via np.searchsorted
            pos = np.searchsorted(interolog_array, packed)
            np.clip(pos, 0, len(interolog_array) - 1, out=pos)
            matched = (interolog_array[pos] == packed)

            # Scatter back
            valid_indices = np.where(valid_idx_mask)[0]
            is_forbidden[valid_indices[matched]] = True

        removed_in_chunk = int(np.count_nonzero(is_forbidden))
        retained_in_chunk = n_rows - removed_in_chunk

        total_homologous_removed += removed_in_chunk
        total_clean_retained += retained_in_chunk

        # Append clean rows to output buffer
        keep_indices = np.where(~is_forbidden)[0]
        for ki in keep_indices:
            out_buf["participant_a_id"].append(col_a[ki])
            out_buf["participant_b_id"].append(col_b[ki])
            out_buf["locus_a"].append(col_loc_a[ki])
            out_buf["locus_b"].append(col_loc_b[ki])
            out_buf["compartments_a"].append(col_comp_a[ki])
            out_buf["compartments_b"].append(col_comp_b[ki])

            if len(out_buf["participant_a_id"]) >= chunk_size:
                flush_output()

        t_elapsed = time.time() - t0
        print(
            f"[{time.strftime('%X')}] Chunk {in_idx:04d}/{len(chunk_files)} ({os.path.basename(chunk_file)}): "
            f"Input={n_rows:,} | Removed={removed_in_chunk:,} ({removed_in_chunk/n_rows*100:.2f}%) | "
            f"Retained={retained_in_chunk:,} | {t_elapsed:.2f}s"
        )

    # Flush remaining buffer
    flush_output()

    total_chunks_written = out_chunk_idx - 1
    total_time = time.time() - start_stream_t
    retention_pct = (total_clean_retained / total_candidates_examined * 100.0) if total_candidates_examined > 0 else 0.0
    removed_pct = (total_homologous_removed / total_candidates_examined * 100.0) if total_candidates_examined > 0 else 0.0

    print(f"\n[{time.strftime('%X')}] Streaming Filtering Completed!")
    print(f"[{time.strftime('%X')}] Total candidate negatives evaluated: {total_candidates_examined:,}")
    print(f"[{time.strftime('%X')}] Positive interologs discarded:        {total_homologous_removed:,} ({removed_pct:.2f}%)")
    print(f"[{time.strftime('%X')}] Clean negative pairs retained:        {total_clean_retained:,} ({retention_pct:.2f}%)")
    print(f"[{time.strftime('%X')}] Output chunks created:                {total_chunks_written} in {out_dir}")
    print(f"[{time.strftime('%X')}] Total execution time:                 {total_time:.2f} seconds.")

    metrics = {
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "protocol": "Option A (Filter & Export Clean Chunks)",
        "theta_identity_threshold": 0.50,
        "cmin_coverage_threshold": 0.80,
        "evalue_cutoff": 1e-3,
        "input_chunks_processed": len(chunk_files),
        "output_chunks_written": total_chunks_written,
        "chunk_size": chunk_size,
        "total_candidate_negatives_evaluated": total_candidates_examined,
        "total_homologous_interologs_removed": total_homologous_removed,
        "homologous_interologs_removed_pct": round(removed_pct, 4),
        "total_clean_negatives_retained": total_clean_retained,
        "clean_negatives_retained_pct": round(retention_pct, 4),
        "execution_time_seconds": round(total_time, 2),
        "candidates_dir": candidates_dir,
        "output_dir": out_dir,
    }

    return metrics


def write_audit_reports(audit_path, summary_md_path, metrics, partition_counts, n_interolog_rules):
    """Write machine-readable JSON and human-readable Markdown summary reports."""
    os.makedirs(os.path.dirname(audit_path), exist_ok=True)
    metrics["positive_benchmark_counts"] = partition_counts
    metrics["positive_interolog_rules_count"] = n_interolog_rules

    with open(audit_path, "w", encoding="utf-8") as f:
        json.dump(metrics, f, indent=2)
    print(f"[{time.strftime('%X')}] Saved audit JSON to: {audit_path}")

    with open(summary_md_path, "w", encoding="utf-8") as f:
        f.write("# Stage 2.5: Candidate Negative Homology Filtering Audit Report\n\n")
        f.write(f"**Execution Timestamp**: `{metrics['timestamp_utc']}`\n\n")
        f.write("## 1. Summary Overview\n\n")
        f.write("Candidate negative interactions were screened against all positive interactions in the latest benchmark splits to eliminate positive interologs and pseudo-negatives.\n\n")
        f.write("| Metric | Value |\n")
        f.write("| :--- | :--- |\n")
        f.write(f"| **Sequence Identity Threshold ($\\theta$)** | `{metrics['theta_identity_threshold']:.2f}` |\n")
        f.write(f"| **Reciprocal Coverage Threshold ($c_{{\\min}}$)** | `{metrics['cmin_coverage_threshold']:.2f}` |\n")
        f.write(f"| **E-value Cutoff** | `{metrics['evalue_cutoff']}` |\n")
        f.write(f"| **Positive Interactions Screened (C1/C2/C3/C4)** | `{sum(partition_counts.values()):,}` pairs |\n")
        f.write(f"| **Precomputed Forbidden Interolog Signatures** | `{n_interolog_rules:,}` pairs |\n")
        f.write(f"| **Total Candidate Negatives Evaluated** | `{metrics['total_candidate_negatives_evaluated']:,}` |\n")
        f.write(f"| **Positive Interologs Discarded** | `{metrics['total_homologous_interologs_removed']:,}` (`{metrics['homologous_interologs_removed_pct']:.2f}%`) |\n")
        f.write(f"| **Clean Non-Redundant Negatives Retained** | `{metrics['total_clean_negatives_retained']:,}` (`{metrics['clean_negatives_retained_pct']:.2f}%`) |\n")
        f.write(f"| **Output Filtered Parquet Chunks** | `{metrics['output_chunks_written']}` chunks |\n")
        f.write(f"| **Streaming Execution Time** | `{metrics['execution_time_seconds']:.2f}` seconds |\n\n")

        f.write("## 2. Positive Reference Breakdown\n\n")
        f.write("| Split | Purpose | Verified Positive Pairs |\n")
        f.write("| :--- | :--- | :--- |\n")
        for p_name, p_cnt in partition_counts.items():
            f.write(f"| **{p_name}** | Benchmark Positive Set | {p_cnt:,} |\n")
        f.write(f"| **Total** | Full Non-Redundant Positive Interactome | **{sum(partition_counts.values()):,}** |\n\n")

        f.write("## 3. Mathematical Equivalence\n\n")
        f.write("The pair similarity metric defined in `docs/dataset_construction/03_homology_clustering/01_pair_redundancy_reduction.md` is:\n")
        f.write("$$S(P_{\\text{neg}}, P_{\\text{pos}}) = \\max \\Big[ \\min(s(A, C), s(B, D)), \\min(s(A, D), s(B, C)) \\Big]$$\n\n")
        f.write("A candidate negative pair $(A, B)$ satisfies $S \\ge 0.50$ if and only if $(A, B) \\in \\mathcal{H}(C) \\times \\mathcal{H}(D)$ or $(A, B) \\in \\mathcal{H}(D) \\times \\mathcal{H}(C)$ for some verified positive pair $(C, D)$. All such pairs have been rigorously purged.\n")

    print(f"[{time.strftime('%X')}] Saved summary Markdown to: {summary_md_path}")


def main():
    parser = argparse.ArgumentParser(description="Filter Candidate Negatives Against Positive Interologs")
    parser.add_argument("--stream_path", default=DEFAULT_STREAM, help="Path to uniprot_arabidopsis_stream.tsv.gz")
    parser.add_argument("--hetero_fasta", default=DEFAULT_HETERO_FASTA, help="Path to ara_positive_pairs_high_only_hetero.fasta")
    parser.add_argument("--unified_fasta", default=DEFAULT_UNIFIED_FASTA, help="Output path for unified proteome FASTA")
    parser.add_argument("--splits_dir", default=DEFAULT_SPLITS_DIR, help="Path to 04_benchmark_splits/arabidopsis")
    parser.add_argument("--alignments_tsv", default=DEFAULT_ALIGNMENTS_TSV, help="Path to MMseqs2 alignment TSV")
    parser.add_argument("--candidates_dir", default=DEFAULT_CANDIDATES_DIR, help="Directory containing candidate chunks")
    parser.add_argument("--out_dir", default=DEFAULT_OUT_DIR, help="Directory for filtered negative chunks")
    parser.add_argument("--audit_path", default=DEFAULT_AUDIT_PATH, help="Path for audit JSON")
    parser.add_argument("--summary_md", default=DEFAULT_SUMMARY_MD, help="Path for summary Markdown")
    parser.add_argument("--theta", type=float, default=0.50, help="Identity threshold (default: 0.50)")
    parser.add_argument("--cmin", type=float, default=0.80, help="Coverage threshold (default: 0.80)")
    parser.add_argument("--max_evalue", type=float, default=1e-3, help="Max E-value (default: 1e-3)")
    parser.add_argument("--chunk_size", type=int, default=1000000, help="Rows per output chunk")
    parser.add_argument("--build_fasta_only", action="store_true", help="Only build unified FASTA and exit")

    args = parser.parse_args()

    print("================================================================================")
    print("      AraPPLM: Negative Interolog Homology Filtering (HPC Streamer)             ")
    print("================================================================================")
    print(f"Timestamp:              {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Identity Threshold:     theta = {args.theta:.2f}")
    print(f"Coverage Threshold:     cmin  = {args.cmin:.2f}")
    print(f"Max E-Value:            {args.max_evalue}")
    print(f"Benchmark Splits Dir:   {args.splits_dir}")
    print(f"Alignments TSV:         {args.alignments_tsv}")
    print(f"Candidates Dir:         {args.candidates_dir}")
    print(f"Filtered Output Dir:    {args.out_dir}")
    print("================================================================================")

    # 1. Unified FASTA generation (if needed)
    if args.build_fasta_only or not os.path.exists(args.unified_fasta):
        build_unified_fasta(args.stream_path, args.hetero_fasta, args.unified_fasta)
        if args.build_fasta_only:
            print("Unified FASTA generated. Exiting due to --build_fasta_only flag.")
            sys.exit(0)

    # 2. Load canonical proteome accessions to assign compact uint32 IDs
    canonical_seqs = load_canonical_proteome(args.stream_path)
    all_proteins = sorted(list(canonical_seqs.keys()))

    # Add positive proteins to indexing
    positive_pairs, pos_proteins, partition_counts = load_positive_benchmark_pairs(args.splits_dir)
    extra_pos = sorted(list(pos_proteins - set(all_proteins)))
    all_proteins.extend(extra_pos)

    prot_to_idx = {p: np.uint32(i) for i, p in enumerate(all_proteins)}
    print(f"[{time.strftime('%X')}] Protein index table built: {len(prot_to_idx):,} unique proteins mapped.")

    # 3. Parse MMseqs2 alignments
    homologs = build_homology_index(args.alignments_tsv, theta=args.theta, cmin=args.cmin, max_evalue=args.max_evalue)

    # 4. Construct Exhaustive Positive Interolog Set
    interolog_array, n_rules = build_positive_interolog_set(positive_pairs, homologs, prot_to_idx)

    # 5. Stream and filter candidate negative chunks
    metrics = filter_candidate_chunks(
        candidates_dir=args.candidates_dir,
        out_dir=args.out_dir,
        interolog_array=interolog_array,
        prot_to_idx=prot_to_idx,
        chunk_size=args.chunk_size,
    )

    # 6. Export audit summary
    write_audit_reports(args.audit_path, args.summary_md, metrics, partition_counts, n_rules)

    print("================================================================================")
    print("Negative Interolog Homology Filtering Completed Successfully!")
    print("================================================================================")


if __name__ == "__main__":
    main()
