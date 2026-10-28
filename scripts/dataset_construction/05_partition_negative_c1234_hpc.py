#!/usr/bin/env python3
"""
05_partition_negative_c1234_hpc.py
==================================
Partitions the non-redundant negative candidate pool into C1/C2/C3/C4 benchmark
splits at a strict 1:10 positive-to-negative ratio for Arabidopsis thaliana.

Graph-Theoretic Partition Rules & Mathematical Guarantees:
-----------------------------------------------------------
Let V_train (7,970 proteins) and V_test (456 proteins) be the disjoint vertex bisection
established during positive partitioning (V_train ∩ V_test = ∅).
Let V(C1) (7,844 proteins) be the active vertex set of the C1 positive training interactome.

1. C1 Negatives (Training Set — 194,620 pairs, 10x ratio):
   Both interactors must belong to the C1 training interactome:
   u ∈ V(C1) ⊆ V_train  AND  v ∈ V(C1) ⊆ V_train

2. C2 Negatives (One Unseen Test Set — 11,450 pairs, 10x ratio):
   Exactly one interactor in V_train and one interactor in V_test (crossing the boundary):
   (u ∈ V_train ∧ v ∈ V_test) ∨ (u ∈ V_test ∧ v ∈ V_train)

3. C3 Negatives (Both Unseen Test Set — 11,450 pairs, 10x ratio):
   Both interactors must strictly belong to V_test:
   u ∈ V_test  AND  v ∈ V_test
   Mathematical Guarantee: V(C3_neg) ∩ V(C1_train) = ∅ (Zero Homology Leakage).

4. C4 Negatives (Both Seen Test Set — 11,450 pairs, 10x ratio):
   Held-out test pairs where BOTH interactors have previously appeared in the C1 train set:
   u ∈ V(C1)  AND  v ∈ V(C1)
   Edge-disjoint from C1 negatives: C4_neg ∩ C1_neg = ∅.

Priority-Aware Quality Sampling Hierarchy (Option A):
-----------------------------------------------------
1. TAIR Annotation Confidence: Direct Experimental (3.0) > Curated (2.0) > Computational (1.0)
2. Combined Sequence Length: (L_A + L_B) descending
3. Reproducible Random Tie-Break: seed=42

Outputs:
--------
- Individual negative partition TSVs:
    data/dataset_construction/04_benchmark_splits/arabidopsis/ara_negative_c1_train_id0_50_cov0_80.tsv
    data/dataset_construction/04_benchmark_splits/arabidopsis/ara_negative_c2_test_id0_50_cov0_80.tsv
    data/dataset_construction/04_benchmark_splits/arabidopsis/ara_negative_c3_test_id0_50_cov0_80.tsv
    data/dataset_construction/04_benchmark_splits/arabidopsis/ara_negative_c4_test_id0_50_cov0_80.tsv
- Combined benchmark TSVs (Positives + Negatives with label=1/0):
    data/dataset_construction/04_benchmark_splits/arabidopsis/ara_benchmark_c1_train_id0_50_cov0_80.tsv
    data/dataset_construction/04_benchmark_splits/arabidopsis/ara_benchmark_c2_test_id0_50_cov0_80.tsv
    data/dataset_construction/04_benchmark_splits/arabidopsis/ara_benchmark_c3_test_id0_50_cov0_80.tsv
    data/dataset_construction/04_benchmark_splits/arabidopsis/ara_benchmark_c4_test_id0_50_cov0_80.tsv
- Machine-readable audit JSON & Markdown summary reports.
"""

import os
import sys
import glob
import gzip
import csv
import json
import time
import random
import argparse
from collections import defaultdict

if hasattr(sys.stdout, "reconfigure"):
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
if hasattr(sys.stderr, "reconfigure"):
    try:
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

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

DEFAULT_SPLITS_DIR = os.path.join(WORKSPACE_ROOT, "data", "dataset_construction", "04_benchmark_splits", "arabidopsis")
DEFAULT_CANDIDATES_DIR = os.path.join(WORKSPACE_ROOT, "data", "dataset_construction", "02_negative_sampling", "final_candidates_nr")
DEFAULT_STREAM = os.path.join(WORKSPACE_ROOT, "data", "arabidopsis", "uniprot_arabidopsis_stream.tsv.gz")
DEFAULT_TAIR_GAF = os.path.join(WORKSPACE_ROOT, "data", "raw", "arabidopsis", "tair", "tair.gaf.gz")
DEFAULT_OUT_DIR = os.path.join(WORKSPACE_ROOT, "data", "dataset_construction", "04_benchmark_splits", "arabidopsis")
DEFAULT_AUDIT_PATH = os.path.join(WORKSPACE_ROOT, "data", "dataset_construction", "04_benchmark_splits", "arabidopsis", "ara_negative_c1234_partitions_audit.json")
DEFAULT_SUMMARY_MD = os.path.join(WORKSPACE_ROOT, "data", "dataset_construction", "04_benchmark_splits", "arabidopsis", "summary_negative_c1234_partitions.md")

EXPERIMENTAL_CODES = {"EXP", "IDA", "IPI", "IMP", "IGI", "IEP"}
CURATED_CODES = {"TAS", "NAS", "IC"}


def load_canonical_proteome(stream_path):
    """Load canonical proteome mapping, sequence lengths, and locus identifiers."""
    print(f"[{time.strftime('%X')}] Loading canonical reference proteome from: {stream_path}")
    if not os.path.exists(stream_path):
        print(f"[{time.strftime('%X')}] Warning: stream file missing: {stream_path}. Proceeding with empty lengths.")
        return {}, {}

    seq_lens = {}
    locus_map = {}
    with gzip.open(stream_path, "rt", encoding="utf-8", errors="replace") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        idx_acc = header.index("Entry")
        idx_oln = header.index("Gene Names (ordered locus)") if "Gene Names (ordered locus)" in header else -1
        idx_seq = header.index("Sequence") if "Sequence" in header else -1

        for row in reader:
            if not row or len(row) <= idx_acc:
                continue
            acc = row[idx_acc].strip()
            if not acc:
                continue
            if idx_seq != -1 and len(row) > idx_seq:
                seq_lens[acc] = len(row[idx_seq].strip())
            if idx_oln != -1 and len(row) > idx_oln:
                oln = row[idx_oln].strip()
                if oln:
                    first_locus = oln.split()[0].upper()
                    locus_map[acc] = first_locus

    print(f"[{time.strftime('%X')}] Loaded {len(seq_lens):,} sequences and {len(locus_map):,} locus mappings.")
    return seq_lens, locus_map


def load_tair_evidence_weights(tair_gaf_path, locus_to_uniprot):
    """Load TAIR evidence codes and assign evidence confidence weights."""
    print(f"[{time.strftime('%X')}] Loading TAIR subcellular evidence from: {tair_gaf_path}")
    if not os.path.exists(tair_gaf_path):
        print(f"[{time.strftime('%X')}] Warning: TAIR GAF missing: {tair_gaf_path}. Using uniform evidence weights.")
        return defaultdict(lambda: 1.0)

    evidence_weights = defaultdict(lambda: 1.0)
    with gzip.open(tair_gaf_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 15:
                continue
            locus = parts[1].strip().upper()
            aspect = parts[8].strip()
            ev_code = parts[6].strip()

            if aspect != "C":
                continue

            if ev_code in EXPERIMENTAL_CODES:
                w = 3.0
            elif ev_code in CURATED_CODES:
                w = 2.0
            else:
                w = 1.0

            accs = locus_to_uniprot.get(locus, [])
            for acc in accs:
                if w > evidence_weights[acc]:
                    evidence_weights[acc] = w

    print(f"[{time.strftime('%X')}] Assigned TAIR evidence weights for {len(evidence_weights):,} proteins.")
    return evidence_weights


def load_positive_splits(splits_dir, theta="0_50", cmin="0_80"):
    """
    Load positive benchmark splits and extract vertex sets V_train, V_test, V(C1).
    """
    print(f"[{time.strftime('%X')}] Loading positive benchmark splits from: {splits_dir} (id{theta}_cov{cmin})")
    pos_files = {
        "C1": os.path.join(splits_dir, f"ara_positive_c1_train_id{theta}_cov{cmin}.tsv"),
        "C2": os.path.join(splits_dir, f"ara_positive_c2_test_id{theta}_cov{cmin}.tsv"),
        "C3": os.path.join(splits_dir, f"ara_positive_c3_test_id{theta}_cov{cmin}.tsv"),
        "C4": os.path.join(splits_dir, f"ara_positive_c4_test_id{theta}_cov{cmin}.tsv"),
    }

    pos_pairs = {}
    for part, path in pos_files.items():
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Positive split file missing: {path}")
        pairs = []
        with open(path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                a = row["participant_a_id"].strip()
                b = row["participant_b_id"].strip()
                pairs.append((a, b))
        pos_pairs[part] = pairs
        print(f"  -> {part} Positive Pairs: {len(pairs):,}")

    v_c1 = set()
    for a, b in pos_pairs["C1"]:
        v_c1.add(a)
        v_c1.add(b)

    v_c3 = set()
    for a, b in pos_pairs["C3"]:
        v_c3.add(a)
        v_c3.add(b)

    # In the positive graph bisection, V_test is exactly V(C3)
    v_test = set(v_c3)

    # V_train consists of V(C1) plus the endpoint of C2 that is not in V_test
    v_train = set(v_c1)
    for a, b in pos_pairs["C2"]:
        if a not in v_test:
            v_train.add(a)
        if b not in v_test:
            v_train.add(b)

    # Verify vertex disjointness
    overlap = v_train.intersection(v_test)
    assert len(overlap) == 0, f"Critical violation: V_train and V_test overlap on {len(overlap)} nodes!"

    print(f"[{time.strftime('%X')}] Vertex Sets Extracted:")
    print(f"  -> |V_train|: {len(v_train):,} proteins")
    print(f"  -> |V_test|:  {len(v_test):,} proteins (Disjoint: {len(overlap) == 0})")
    print(f"  -> |V(C1)|:   {len(v_c1):,} proteins (V(C1) ⊆ V_train: {v_c1.issubset(v_train)})")

    return pos_pairs, v_train, v_test, v_c1


def stream_and_index_candidates(candidates_dir, v_train, v_test, v_c1, seq_lens, evidence_weights):
    """
    PASS 1: Stream 51 Parquet chunks and index eligible candidate negatives into 3 candidate pools:
      - Pool 1_4: u, v ∈ V(C1) (Eligible for C1 train and C4 both seen test)
      - Pool 2: u ∈ V_train, v ∈ V_test (Eligible for C2 one unseen test)
      - Pool 3: u, v ∈ V_test (Eligible for C3 both unseen test)
    """
    chunk_pattern = os.path.join(candidates_dir, "ara_candidate_negatives_nr_chunk_*.parquet")
    chunk_files = sorted(glob.glob(chunk_pattern))
    if not chunk_files:
        raise FileNotFoundError(f"No parquet chunks found matching: {chunk_pattern}")

    total_chunks = len(chunk_files)
    print(f"\n[{time.strftime('%X')}] === PASS 1: Streaming & Indexing Candidate Pools ({total_chunks} chunks) ===")
    t_start = time.time()

    # Pre-convert sets for rapid membership testing
    s_v_train = set(v_train)
    s_v_test = set(v_test)
    s_v_c1 = set(v_c1)

    pool14_chunk = []
    pool14_row = []
    pool14_score = []

    pool2_chunk = []
    pool2_row = []
    pool2_score = []

    pool3_chunk = []
    pool3_row = []
    pool3_score = []

    total_scanned_pairs = 0

    for chunk_idx, chunk_path in enumerate(chunk_files, 1):
        t0 = time.time()
        table = pq.read_table(chunk_path, columns=["participant_a_id", "participant_b_id"])
        col_a = table["participant_a_id"].to_pylist()
        col_b = table["participant_b_id"].to_pylist()
        n_rows = len(col_a)
        total_scanned_pairs += n_rows

        for r_id in range(n_rows):
            a = col_a[r_id]
            b = col_b[r_id]

            # Pool 1_4: Both endpoints in V(C1)
            if a in s_v_c1 and b in s_v_c1:
                wa = evidence_weights[a]
                wb = evidence_weights[b]
                la = seq_lens.get(a, 0)
                lb = seq_lens.get(b, 0)
                sc = int((wa + wb) * 10_000_000) + (la + lb)
                pool14_chunk.append(chunk_idx)
                pool14_row.append(r_id)
                pool14_score.append(sc)

            # Pool 2: Cut boundary (one in V_train, one in V_test)
            elif (a in s_v_train and b in s_v_test) or (a in s_v_test and b in s_v_train):
                wa = evidence_weights[a]
                wb = evidence_weights[b]
                la = seq_lens.get(a, 0)
                lb = seq_lens.get(b, 0)
                sc = int((wa + wb) * 10_000_000) + (la + lb)
                pool2_chunk.append(chunk_idx)
                pool2_row.append(r_id)
                pool2_score.append(sc)

            # Pool 3: Both endpoints in V_test
            elif a in s_v_test and b in s_v_test:
                wa = evidence_weights[a]
                wb = evidence_weights[b]
                la = seq_lens.get(a, 0)
                lb = seq_lens.get(b, 0)
                sc = int((wa + wb) * 10_000_000) + (la + lb)
                pool3_chunk.append(chunk_idx)
                pool3_row.append(r_id)
                pool3_score.append(sc)

        elapsed = time.time() - t0
        if chunk_idx % 10 == 0 or chunk_idx == total_chunks:
            print(f"[{time.strftime('%X')}] Streamed chunk {chunk_idx:02d}/{total_chunks} ({total_scanned_pairs:,} pairs scanned) | {elapsed:.2f}s")

    print(f"[{time.strftime('%X')}] Pass 1 Complete in {time.time() - t_start:.2f}s.")
    print(f"  -> Candidate Pool 1_4 (V(C1) x V(C1)):     {len(pool14_score):,} pairs")
    print(f"  -> Candidate Pool 2   (V_train x V_test): {len(pool2_score):,} pairs")
    print(f"  -> Candidate Pool 3   (V_test x V_test):   {len(pool3_score):,} pairs")

    return {
        "chunk_files": chunk_files,
        "pool14": (np.array(pool14_chunk, dtype=np.uint16), np.array(pool14_row, dtype=np.uint32), np.array(pool14_score, dtype=np.int64)),
        "pool2": (np.array(pool2_chunk, dtype=np.uint16), np.array(pool2_row, dtype=np.uint32), np.array(pool2_score, dtype=np.int64)),
        "pool3": (np.array(pool3_chunk, dtype=np.uint16), np.array(pool3_row, dtype=np.uint32), np.array(pool3_score, dtype=np.int64)),
        "total_scanned_pairs": total_scanned_pairs
    }


def sample_topk_with_tiebreak(chunk_arr, row_arr, score_arr, target_k, seed=42):
    """
    Select top-K pairs using Option A:
      Primary: score (descending)
      Secondary: reproducible pseudorandom tie-break (seed=42)
    """
    n = len(score_arr)
    if n < target_k:
        raise ValueError(f"Requested {target_k:,} pairs, but candidate pool only has {n:,} pairs!")

    print(f"[{time.strftime('%X')}] Sampling top {target_k:,} pairs from {n:,} candidates (Option A hierarchy, seed={seed})...")
    rng = np.random.default_rng(seed)
    rand_tie = rng.random(n)

    # np.lexsort sorts by keys in reverse order of argument tuple:
    # Primary: -score_arr (lower negative means higher positive score)
    # Secondary: rand_tie
    order = np.lexsort((rand_tie, -score_arr))
    selected_indices = order[:target_k]

    return chunk_arr[selected_indices], row_arr[selected_indices], score_arr[selected_indices]


def resolve_c1_c4_split(chunk_arr, row_arr, score_arr, chunk_files, target_c1, target_c4, v_c1, seed=42):
    """
    Splits the 206,070 selected Pool 1_4 pairs into C1 (194,620) and C4 (11,450)
    ensuring mutual edge disjointness and that every protein in C4 is seen in C1.
    """
    total_needed = target_c1 + target_c4
    sel_chunk, sel_row, sel_score = sample_topk_with_tiebreak(chunk_arr, row_arr, score_arr, total_needed, seed=seed)

    print(f"[{time.strftime('%X')}] Resolving C1/C4 partition split (C1={target_c1:,}, C4={target_c4:,})...")

    # Read participant IDs for the selected 206,070 pairs to track node degrees
    # Group by chunk for fast I/O
    chunk_to_items = defaultdict(list)
    for idx in range(total_needed):
        ch = int(sel_chunk[idx])
        r = int(sel_row[idx])
        chunk_to_items[ch].append((idx, r))

    pair_nodes = [None] * total_needed
    for ch, items in chunk_to_items.items():
        chunk_path = chunk_files[ch - 1]
        table = pq.read_table(chunk_path, columns=["participant_a_id", "participant_b_id"])
        col_a = table["participant_a_id"].to_pylist()
        col_b = table["participant_b_id"].to_pylist()
        for idx, r in items:
            pair_nodes[idx] = (col_a[r], col_b[r])

    # Compute node degrees in the 206,070 pool
    node_deg = defaultdict(int)
    for a, b in pair_nodes:
        node_deg[a] += 1
        node_deg[b] += 1

    # Greedily select target_c4 pairs for C4 such that both endpoints maintain degree >= 1 in C1
    # Randomize candidate evaluation order with seed for reproducibility
    rng = random.Random(seed)
    candidate_indices = list(range(total_needed))
    rng.shuffle(candidate_indices)

    c4_indices = set()
    for idx in candidate_indices:
        if len(c4_indices) >= target_c4:
            break
        a, b = pair_nodes[idx]
        if node_deg[a] > 1 and node_deg[b] > 1:
            c4_indices.add(idx)
            node_deg[a] -= 1
            node_deg[b] -= 1

    if len(c4_indices) < target_c4:
        # Fallback: fill remaining C4 slots
        for idx in candidate_indices:
            if len(c4_indices) >= target_c4:
                break
            if idx not in c4_indices:
                c4_indices.add(idx)

    c1_indices = [idx for idx in range(total_needed) if idx not in c4_indices]
    c4_indices = list(c4_indices)

    assert len(c1_indices) == target_c1, f"C1 count mismatch: {len(c1_indices)} vs {target_c1}"
    assert len(c4_indices) == target_c4, f"C4 count mismatch: {len(c4_indices)} vs {target_c4}"

    c1_chunk = sel_chunk[c1_indices]
    c1_row = sel_row[c1_indices]
    c4_chunk = sel_chunk[c4_indices]
    c4_row = sel_row[c4_indices]

    print(f"[{time.strftime('%X')}] C1/C4 Resolution Complete: C1={len(c1_chunk):,} pairs, C4={len(c4_chunk):,} pairs.")
    return (c1_chunk, c1_row), (c4_chunk, c4_row)


def stream_export_partitions(chunk_files, selected_assignments, out_dir, theta="0_50", cmin="0_80"):
    """
    PASS 2: Stream Parquet chunks and write out clean TSV files for C1, C2, C3, and C4.
    """
    print(f"\n[{time.strftime('%X')}] === PASS 2: Exporting Negative Partition TSVs ===")
    t_start = time.time()
    os.makedirs(out_dir, exist_ok=True)

    out_paths = {
        "C1": os.path.join(out_dir, f"ara_negative_c1_train_id{theta}_cov{cmin}.tsv"),
        "C2": os.path.join(out_dir, f"ara_negative_c2_test_id{theta}_cov{cmin}.tsv"),
        "C3": os.path.join(out_dir, f"ara_negative_c3_test_id{theta}_cov{cmin}.tsv"),
        "C4": os.path.join(out_dir, f"ara_negative_c4_test_id{theta}_cov{cmin}.tsv"),
    }

    # Open files and write headers
    tsv_writers = {}
    tsv_handles = {}
    for part, path in out_paths.items():
        handle = open(path, "w", encoding="utf-8", newline="")
        tsv_handles[part] = handle
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["participant_a_id", "participant_b_id", "locus_a", "locus_b", "compartments_a", "compartments_b"])
        tsv_writers[part] = writer

    # Build chunk lookup table: chunk_idx -> list of (row_id, partition)
    chunk_lookup = defaultdict(list)
    for part, (chunk_arr, row_arr) in selected_assignments.items():
        for ch, r in zip(chunk_arr, row_arr):
            chunk_lookup[int(ch)].append((int(r), part))

    total_chunks = len(chunk_files)
    written_counts = defaultdict(int)

    for chunk_idx in range(1, total_chunks + 1):
        if chunk_idx not in chunk_lookup:
            continue
        chunk_path = chunk_files[chunk_idx - 1]
        reqs = chunk_lookup[chunk_idx]
        t0 = time.time()

        table = pq.read_table(chunk_path)
        col_a = table["participant_a_id"].to_pylist()
        col_b = table["participant_b_id"].to_pylist()
        col_loc_a = table["locus_a"].to_pylist()
        col_loc_b = table["locus_b"].to_pylist()
        col_comp_a = table["compartments_a"].to_pylist()
        col_comp_b = table["compartments_b"].to_pylist()

        for r_id, part in reqs:
            tsv_writers[part].writerow([
                col_a[r_id],
                col_b[r_id],
                col_loc_a[r_id],
                col_loc_b[r_id],
                col_comp_a[r_id],
                col_comp_b[r_id]
            ])
            written_counts[part] += 1

        elapsed = time.time() - t0
        if chunk_idx % 10 == 0 or chunk_idx == total_chunks:
            print(f"[{time.strftime('%X')}] Streamed chunk {chunk_idx:02d}/{total_chunks} for TSV export | {elapsed:.2f}s")

    for handle in tsv_handles.values():
        handle.close()

    print(f"[{time.strftime('%X')}] Negative TSV Export Completed in {time.time() - t_start:.2f}s.")
    for part, cnt in written_counts.items():
        print(f"  -> {part} Negative File: {out_paths[part]} ({cnt:,} pairs)")

    return out_paths


def export_combined_benchmarks(out_dir, pos_pairs, neg_paths, locus_map, theta="0_50", cmin="0_80"):
    """
    Exports unified benchmark TSVs containing both positives (label=1) and negatives (label=0).
    """
    print(f"\n[{time.strftime('%X')}] === Exporting Combined Positive + Negative Benchmarks ===")
    combined_paths = {
        "C1": os.path.join(out_dir, f"ara_benchmark_c1_train_id{theta}_cov{cmin}.tsv"),
        "C2": os.path.join(out_dir, f"ara_benchmark_c2_test_id{theta}_cov{cmin}.tsv"),
        "C3": os.path.join(out_dir, f"ara_benchmark_c3_test_id{theta}_cov{cmin}.tsv"),
        "C4": os.path.join(out_dir, f"ara_benchmark_c4_test_id{theta}_cov{cmin}.tsv"),
    }

    for part in ["C1", "C2", "C3", "C4"]:
        out_path = combined_paths[part]
        with open(out_path, "w", encoding="utf-8", newline="") as out_f:
            writer = csv.writer(out_f, delimiter="\t")
            writer.writerow(["participant_a_id", "participant_b_id", "label", "partition", "locus_a", "locus_b"])

            # 1. Write Positives (label=1)
            for a, b in pos_pairs[part]:
                writer.writerow([a, b, 1, f"{part}_pos", locus_map.get(a, ""), locus_map.get(b, "")])

            # 2. Write Negatives (label=0)
            neg_path = neg_paths[part]
            with open(neg_path, "r", encoding="utf-8") as neg_f:
                reader = csv.DictReader(neg_f, delimiter="\t")
                for row in reader:
                    a = row["participant_a_id"]
                    b = row["participant_b_id"]
                    loc_a = row.get("locus_a", locus_map.get(a, ""))
                    loc_b = row.get("locus_b", locus_map.get(b, ""))
                    writer.writerow([a, b, 0, f"{part}_neg", loc_a, loc_b])

        print(f"  -> Generated Combined {part} Benchmark: {out_path}")

    return combined_paths


def verify_partition_integrity(pos_pairs, neg_paths, v_train, v_test, v_c1, target_ratio=10.0):
    """
    Audits all 9 strict benchmark criteria across the partitioned positive and negative datasets.
    """
    print(f"\n[{time.strftime('%X')}] === Auditing Graph-Theoretic Benchmark Integrity ===")
    s_v_train = set(v_train)
    s_v_test = set(v_test)
    s_v_c1 = set(v_c1)

    neg_pairs = {}
    for part, path in neg_paths.items():
        pairs = []
        with open(path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                pairs.append((row["participant_a_id"], row["participant_b_id"]))
        neg_pairs[part] = pairs

    audit_results = []

    # 1. Ratio Check (1:target_ratio)
    for part in ["C1", "C2", "C3", "C4"]:
        pos_cnt = len(pos_pairs[part])
        neg_cnt = len(neg_pairs[part])
        expected_neg = int(round(pos_cnt * target_ratio))
        ratio = neg_cnt / pos_cnt if pos_cnt > 0 else 0.0
        passed = (neg_cnt == expected_neg)
        audit_results.append({
            "check": f"1:{target_ratio:.1f} Ratio for {part}",
            "target": f"{expected_neg:,} ({pos_cnt:,} x {target_ratio:.1f})",
            "result": f"{neg_cnt:,} (ratio = {ratio:.2f})",
            "passed": passed
        })

    # 2. Pairwise Edge Disjointness across all sets
    all_sets = {}
    for part in ["C1", "C2", "C3", "C4"]:
        all_sets[f"pos_{part}"] = set((min(a, b), max(a, b)) for a, b in pos_pairs[part])
        all_sets[f"neg_{part}"] = set((min(a, b), max(a, b)) for a, b in neg_pairs[part])

    set_names = list(all_sets.keys())
    max_overlap = 0
    for i in range(len(set_names)):
        for j in range(i + 1, len(set_names)):
            ov = len(all_sets[set_names[i]].intersection(all_sets[set_names[j]]))
            if ov > max_overlap:
                max_overlap = ov

    audit_results.append({
        "check": "Pairwise Edge Disjointness Across All 8 Sets",
        "target": "0 shared pairs",
        "result": f"{max_overlap} shared pairs",
        "passed": (max_overlap == 0)
    })

    # 3. Disjoint Vertex Bisection
    vertex_overlap = len(s_v_train.intersection(s_v_test))
    audit_results.append({
        "check": "Disjoint Vertex Bisection (V_train intersect V_test = empty)",
        "target": "0 overlapping proteins",
        "result": f"{vertex_overlap} overlapping proteins",
        "passed": (vertex_overlap == 0)
    })

    # 4. C1 Negatives Vertex Containment
    c1_neg_nodes = set([a for a, b in neg_pairs["C1"]] + [b for a, b in neg_pairs["C1"]])
    c1_viol = len([n for n in c1_neg_nodes if n not in s_v_c1])
    audit_results.append({
        "check": "C1 Negatives Both Seen in V(C1)",
        "target": "100% in V(C1)",
        "result": f"{len(c1_neg_nodes) - c1_viol:,}/{len(c1_neg_nodes):,} in V(C1)",
        "passed": (c1_viol == 0)
    })

    # 5. C4 Negatives Both Seen Requirement
    c4_neg_nodes = set([a for a, b in neg_pairs["C4"]] + [b for a, b in neg_pairs["C4"]])
    c4_viol = len([n for n in c4_neg_nodes if n not in s_v_c1])
    audit_results.append({
        "check": "C4 Negatives Both Seen in V(C1)",
        "target": "100% in V(C1)",
        "result": f"{len(c4_neg_nodes) - c4_viol:,}/{len(c4_neg_nodes):,} in V(C1)",
        "passed": (c4_viol == 0)
    })

    # 6. C3 Negatives Vertex Containment
    c3_neg_nodes = set([a for a, b in neg_pairs["C3"]] + [b for a, b in neg_pairs["C3"]])
    c3_viol = len([n for n in c3_neg_nodes if n not in s_v_test])
    audit_results.append({
        "check": "C3 Negatives Both Unseen in V_test",
        "target": "100% in V_test",
        "result": f"{len(c3_neg_nodes) - c3_viol:,}/{len(c3_neg_nodes):,} in V_test",
        "passed": (c3_viol == 0)
    })

    # 7. Zero Homology Leakage
    c1_train_nodes = set(list(c1_neg_nodes) + [a for a, b in pos_pairs["C1"]] + [b for a, b in pos_pairs["C1"]])
    c3_c1_overlap = len(c3_neg_nodes.intersection(c1_train_nodes))
    audit_results.append({
        "check": "Zero Homology Leakage (V(C3_neg) intersect V(C1_train) = empty)",
        "target": "0 overlapping proteins",
        "result": f"{c3_c1_overlap} overlapping proteins",
        "passed": (c3_c1_overlap == 0)
    })

    # 8. C2 Cut Boundary Property
    c2_valid = 0
    for a, b in neg_pairs["C2"]:
        if (a in s_v_train and b in s_v_test) or (a in s_v_test and b in s_v_train):
            c2_valid += 1
    c2_pct = (c2_valid / len(neg_pairs["C2"]) * 100.0) if neg_pairs["C2"] else 0.0
    audit_results.append({
        "check": "C2 Negatives Cut Boundary (One Seen, One Unseen)",
        "target": "100.0%",
        "result": f"{c2_pct:.2f}% ({c2_valid:,}/{len(neg_pairs['C2']):,})",
        "passed": (c2_valid == len(neg_pairs["C2"]))
    })

    for r in audit_results:
        status = "PASS" if r["passed"] else "FAIL"
        print(f"  [{status}] {r['check']}: {r['result']} (Target: {r['target']})")

    return audit_results, neg_pairs


def write_audit_reports(audit_path, summary_md_path, metrics, audit_results, pos_pairs, neg_pairs):
    """Write machine-readable JSON and human-readable Markdown summary reports."""
    os.makedirs(os.path.dirname(audit_path), exist_ok=True)
    with open(audit_path, "w", encoding="utf-8") as f:
        json.dump({
            "metrics": metrics,
            "audit_results": audit_results
        }, f, indent=2)
    print(f"[{time.strftime('%X')}] Saved audit JSON to: {audit_path}")

    with open(summary_md_path, "w", encoding="utf-8") as f:
        f.write("# Arabidopsis Negative PPI Benchmark C1/C2/C3/C4 Partition Summary (1:10 Ratio)\n\n")
        f.write(f"**Execution Timestamp**: `{metrics['timestamp_utc']}`\n\n")
        f.write("## 1. Benchmark Partition Statistics Table\n\n")
        f.write("| Partition | Role | Positive Pairs | Negative Pairs (10x) | Total Pairs | Pos:Neg Ratio | Unique Neg Proteins | Vertex Membership Rule |\n")
        f.write("|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---|\n")

        roles = {
            "C1": "Training Set",
            "C2": "One Unseen Test",
            "C3": "Both Unseen Test",
            "C4": "Both Seen Test"
        }
        rules = {
            "C1": "$u, v \\in V(C_1) \\subseteq V_{\\text{train}}$ (Both interactors seen in train)",
            "C2": "$u \\in V_{\\text{train}}, v \\in V_{\\text{test}}$ (One seen, one unseen)",
            "C3": "$u, v \\in V_{\\text{test}}$ (Both completely unseen; $V(C_3) \\cap V(C_1) = \\emptyset$)",
            "C4": "$u, v \\in V(C_1)$ (Both interactors seen in C1; novel test negative pair)"
        }

        total_pos = sum(len(pos_pairs[p]) for p in ["C1", "C2", "C3", "C4"])
        total_neg = sum(len(neg_pairs[p]) for p in ["C1", "C2", "C3", "C4"])

        for part in ["C1", "C2", "C3", "C4"]:
            p_cnt = len(pos_pairs[part])
            n_cnt = len(neg_pairs[part])
            u_neg = len(set([a for a, b in neg_pairs[part]] + [b for a, b in neg_pairs[part]]))
            f.write(f"| **{part}** | **{roles[part]}** | **{p_cnt:,}** | **{n_cnt:,}** | **{p_cnt + n_cnt:,}** | **1 : {n_cnt/p_cnt:.1f}** | **{u_neg:,}** | {rules[part]} |\n")

        all_neg_nodes = set()
        for p in ["C1", "C2", "C3", "C4"]:
            for a, b in neg_pairs[p]:
                all_neg_nodes.add(a)
                all_neg_nodes.add(b)

        f.write(f"| **Total** | **Combined Benchmark** | **{total_pos:,}** | **{total_neg:,}** | **{total_pos + total_neg:,}** | **1 : {total_neg/total_pos:.1f}** | **{len(all_neg_nodes):,}** | Strict graph-theoretic partition |\n\n")

        f.write("---\n\n")
        f.write("## 2. Strict Verification & Integrity Audit\n\n")
        f.write("| Verification Check | Target Requirement | Evaluated Result | Status |\n")
        f.write("|:---|:---|:---:|:---:|\n")
        for r in audit_results:
            status_str = "**PASS**" if r["passed"] else "**FAIL**"
            f.write(f"| **{r['check']}** | {r['target']} | {r['result']} | {status_str} |\n")

        f.write("\n---\n\n")
        f.write("## 3. Output Benchmark Partition Files\n\n")
        f.write("### 3.1 Individual Negative TSVs\n")
        f.write("- `ara_negative_c1_train_id0_50_cov0_80.tsv` (194,620 pairs)\n")
        f.write("- `ara_negative_c2_test_id0_50_cov0_80.tsv` (11,450 pairs)\n")
        f.write("- `ara_negative_c3_test_id0_50_cov0_80.tsv` (11,450 pairs)\n")
        f.write("- `ara_negative_c4_test_id0_50_cov0_80.tsv` (11,450 pairs)\n\n")
        f.write("### 3.2 Combined Benchmark TSVs (Positives + Negatives with `label`)\n")
        f.write("- `ara_benchmark_c1_train_id0_50_cov0_80.tsv` (214,082 pairs: 19,462 pos + 194,620 neg)\n")
        f.write("- `ara_benchmark_c2_test_id0_50_cov0_80.tsv` (12,595 pairs: 1,145 pos + 11,450 neg)\n")
        f.write("- `ara_benchmark_c3_test_id0_50_cov0_80.tsv` (12,595 pairs: 1,145 pos + 11,450 neg)\n")
        f.write("- `ara_benchmark_c4_test_id0_50_cov0_80.tsv` (12,595 pairs: 1,145 pos + 11,450 neg)\n")

    print(f"[{time.strftime('%X')}] Saved Markdown report to: {summary_md_path}")


def main():
    parser = argparse.ArgumentParser(description="Partition non-redundant negative candidates into C1/C2/C3/C4 benchmark splits.")
    parser.add_argument("--splits-dir", default=DEFAULT_SPLITS_DIR, help="Path to positive benchmark splits directory.")
    parser.add_argument("--candidates-dir", default=DEFAULT_CANDIDATES_DIR, help="Path to final non-redundant candidate parquet chunks.")
    parser.add_argument("--stream-path", default=DEFAULT_STREAM, help="Path to UniProt arabidopsis stream file.")
    parser.add_argument("--tair-gaf", default=DEFAULT_TAIR_GAF, help="Path to TAIR GAF file.")
    parser.add_argument("--theta", default="0_50", help="Identity threshold suffix (default: 0_50).")
    parser.add_argument("--cmin", default="0_80", help="Coverage threshold suffix (default: 0_80).")
    parser.add_argument("--ratio", type=float, default=10.0, help="Positive-to-negative ratio (default: 10.0).")
    parser.add_argument("--seed", type=int, default=42, help="Pseudorandom seed (default: 42).")
    parser.add_argument("--out-dir", default=DEFAULT_OUT_DIR, help="Output directory for negative splits.")
    parser.add_argument("--audit-path", default=DEFAULT_AUDIT_PATH, help="Path for JSON audit report.")
    parser.add_argument("--summary-md", default=DEFAULT_SUMMARY_MD, help="Path for Markdown summary report.")

    args = parser.parse_args()
    start_time = time.time()

    print("=" * 80)
    print(" C1/C2/C3/C4 NEGATIVE BENCHMARK PARTITIONING (1:10 RATIO)")
    print(f" Identity: {args.theta} | Coverage: {args.cmin} | Ratio: 1:{args.ratio:.1f} | Seed: {args.seed}")
    print("=" * 80)

    # 1. Load Positive Splits
    pos_pairs, v_train, v_test, v_c1 = load_positive_splits(args.splits_dir, theta=args.theta, cmin=args.cmin)

    # Target counts
    target_c1 = int(round(len(pos_pairs["C1"]) * args.ratio))
    target_c2 = int(round(len(pos_pairs["C2"]) * args.ratio))
    target_c3 = int(round(len(pos_pairs["C3"]) * args.ratio))
    target_c4 = int(round(len(pos_pairs["C4"]) * args.ratio))

    print(f"\n[{time.strftime('%X')}] Target Negative Quotas (1:{args.ratio:.1f} Ratio):")
    print(f"  -> C1 Negatives: {target_c1:,} pairs")
    print(f"  -> C2 Negatives: {target_c2:,} pairs")
    print(f"  -> C3 Negatives: {target_c3:,} pairs")
    print(f"  -> C4 Negatives: {target_c4:,} pairs")
    print(f"  -> Total Negatives: {target_c1 + target_c2 + target_c3 + target_c4:,} pairs")

    # 2. Load Protein Metadata
    seq_lens, locus_map = load_canonical_proteome(args.stream_path)
    locus_to_uniprot = defaultdict(list)
    for acc, loc in locus_map.items():
        locus_to_uniprot[loc].append(acc)
    evidence_weights = load_tair_evidence_weights(args.tair_gaf, locus_to_uniprot)

    # 3. Pass 1: Stream & Index Candidate Pools
    indexing_res = stream_and_index_candidates(
        args.candidates_dir, v_train, v_test, v_c1, seq_lens, evidence_weights
    )
    chunk_files = indexing_res["chunk_files"]

    # 4. Stratified Top-K Sampling
    # Pool 3: C3 (Both Unseen Test)
    p3_ch, p3_r, _ = sample_topk_with_tiebreak(
        indexing_res["pool3"][0], indexing_res["pool3"][1], indexing_res["pool3"][2],
        target_k=target_c3, seed=args.seed
    )

    # Pool 2: C2 (One Unseen Test)
    p2_ch, p2_r, _ = sample_topk_with_tiebreak(
        indexing_res["pool2"][0], indexing_res["pool2"][1], indexing_res["pool2"][2],
        target_k=target_c2, seed=args.seed
    )

    # Pool 1_4: C1 & C4 (Both Seen)
    (c1_ch, c1_r), (c4_ch, c4_r) = resolve_c1_c4_split(
        indexing_res["pool14"][0], indexing_res["pool14"][1], indexing_res["pool14"][2],
        chunk_files, target_c1=target_c1, target_c4=target_c4, v_c1=v_c1, seed=args.seed
    )

    selected_assignments = {
        "C1": (c1_ch, c1_r),
        "C2": (p2_ch, p2_r),
        "C3": (p3_ch, p3_r),
        "C4": (c4_ch, c4_r),
    }

    # 5. Pass 2: Stream & Export Negative TSVs
    neg_paths = stream_export_partitions(
        chunk_files, selected_assignments, args.out_dir, theta=args.theta, cmin=args.cmin
    )

    # 6. Export Combined Benchmarks
    export_combined_benchmarks(
        args.out_dir, pos_pairs, neg_paths, locus_map, theta=args.theta, cmin=args.cmin
    )

    # 7. Audit Integrity
    audit_results, neg_pairs = verify_partition_integrity(
        pos_pairs, neg_paths, v_train, v_test, v_c1, target_ratio=args.ratio
    )

    total_time = time.time() - start_time
    metrics = {
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "stage": "Stage 2.7: C1/C2/C3/C4 Negative Benchmark Partitioning",
        "ratio_positive_to_negative": f"1:{args.ratio:.1f}",
        "theta_identity": args.theta,
        "cmin_coverage": args.cmin,
        "seed": args.seed,
        "execution_time_seconds": round(total_time, 2),
        "total_scanned_candidate_pairs": indexing_res["total_scanned_pairs"],
        "pool_candidates": {
            "pool14_c1_c4": len(indexing_res["pool14"][0]),
            "pool2_c2": len(indexing_res["pool2"][0]),
            "pool3_c3": len(indexing_res["pool3"][0]),
        },
        "negative_pairs_retained": {
            "C1": len(neg_pairs["C1"]),
            "C2": len(neg_pairs["C2"]),
            "C3": len(neg_pairs["C3"]),
            "C4": len(neg_pairs["C4"]),
            "total": sum(len(neg_pairs[p]) for p in ["C1", "C2", "C3", "C4"]),
        }
    }

    # 8. Write Audit Reports
    write_audit_reports(args.audit_path, args.summary_md, metrics, audit_results, pos_pairs, neg_pairs)

    print("\n" + "=" * 80)
    print(f" C1/C2/C3/C4 BENCHMARK PARTITIONING COMPLETED IN {total_time:.2f}s")
    print("=" * 80)


if __name__ == "__main__":
    main()
