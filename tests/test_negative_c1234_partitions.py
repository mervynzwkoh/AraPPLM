#!/usr/bin/env python3
"""
test_negative_c1234_partitions.py
=================================
Automated test suite verifying the C1/C2/C3/C4 negative benchmark partitioning logic,
Option A quality prioritization hierarchy, exact 1:10 quotas, vertex containment,
and zero homology leakage guarantees.
"""

import os
import sys
import csv
import json
import tempfile
import unittest

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

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, REPO_ROOT)

import importlib.util

spec = importlib.util.spec_from_file_location(
    "partition_mod",
    os.path.join(REPO_ROOT, "scripts", "dataset_construction", "05_partition_negative_c1234_hpc.py")
)
partition_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(partition_mod)


class TestNegativeC1234Partitions(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_synthetic_negative_partitioning_pipeline(self):
        """
        Verify the end-to-end partitioning logic on a synthetic candidate dataset:
        Ensures exact quotas, edge disjointness, vertex rules, and zero leakage.
        """
        # Define synthetic vertex sets:
        # Train: T1..T10
        # Test:  U1..U6
        # C1 positive nodes: T1..T8 (subset of V_train)
        v_train = {f"T{i}" for i in range(1, 11)}
        v_test = {f"U{i}" for i in range(1, 7)}
        v_c1 = {f"T{i}" for i in range(1, 9)}

        pos_pairs = {
            "C1": [("T1", "T2"), ("T2", "T3"), ("T3", "T4"), ("T4", "T5")],  # 4 pos -> 8 neg (2x ratio for test)
            "C2": [("T1", "U1"), ("T2", "U2")],                              # 2 pos -> 4 neg
            "C3": [("U1", "U2"), ("U3", "U4")],                              # 2 pos -> 4 neg
            "C4": [("T5", "T6"), ("T7", "T8")],                              # 2 pos -> 4 neg
        }

        # Sequence lengths and weights for Option A
        seq_lens = {f"T{i}": 100 * i for i in range(1, 11)}
        seq_lens.update({f"U{i}": 50 * i for i in range(1, 7)})
        evidence_weights = {f"T{i}": 2.0 for i in range(1, 11)}
        evidence_weights["T1"] = 3.0
        evidence_weights["T2"] = 3.0
        evidence_weights.update({f"U{i}": 1.0 for i in range(1, 7)})

        # Create synthetic candidates across 2 Parquet chunks
        # 1. Pool 1_4 (C1 & C4): pairs from V(C1) x V(C1)
        p14_candidates = [
            ("T1", "T3"), ("T1", "T4"), ("T2", "T4"), ("T2", "T5"),
            ("T3", "T5"), ("T3", "T6"), ("T4", "T6"), ("T4", "T7"),
            ("T5", "T7"), ("T5", "T8"), ("T6", "T8"), ("T1", "T5"),
            ("T2", "T6"), ("T3", "T7")
        ]

        # 2. Pool 2 (C2): pairs from V_train x V_test
        p2_candidates = [
            ("T1", "U3"), ("T2", "U4"), ("T3", "U1"), ("T4", "U2"),
            ("T5", "U5"), ("T6", "U6"), ("T7", "U1")
        ]

        # 3. Pool 3 (C3): pairs from V_test x V_test
        p3_candidates = [
            ("U1", "U3"), ("U1", "U4"), ("U2", "U5"), ("U2", "U6"),
            ("U3", "U5"), ("U4", "U6"), ("U5", "U6")
        ]

        # Ignore pool: external proteins
        p_ext_candidates = [("T1", "X1"), ("U1", "X2"), ("X1", "X2")]

        all_candidates = p14_candidates + p2_candidates + p3_candidates + p_ext_candidates
        mid = len(all_candidates) // 2
        chunk1 = all_candidates[:mid]
        chunk2 = all_candidates[mid:]

        candidates_dir = os.path.join(self.temp_dir.name, "candidates")
        os.makedirs(candidates_dir, exist_ok=True)

        for c_idx, pairs in enumerate([chunk1, chunk2], 1):
            schema = pa.schema([
                ("participant_a_id", pa.string()),
                ("participant_b_id", pa.string()),
                ("locus_a", pa.string()),
                ("locus_b", pa.string()),
                ("compartments_a", pa.string()),
                ("compartments_b", pa.string()),
            ])
            tbl = pa.Table.from_arrays([
                pa.array([p[0] for p in pairs]),
                pa.array([p[1] for p in pairs]),
                pa.array([f"LOC_{p[0]}" for p in pairs]),
                pa.array([f"LOC_{p[1]}" for p in pairs]),
                pa.array(["Cytosol" for _ in pairs]),
                pa.array(["Nucleus" for _ in pairs]),
            ], schema=schema)
            pq.write_table(tbl, os.path.join(candidates_dir, f"ara_candidate_negatives_nr_chunk_{c_idx:03d}.parquet"))

        # Run Pass 1: Stream and Index
        idx_res = partition_mod.stream_and_index_candidates(
            candidates_dir, v_train, v_test, v_c1, seq_lens, evidence_weights
        )

        self.assertEqual(len(idx_res["pool14"][0]), len(p14_candidates))
        self.assertEqual(len(idx_res["pool2"][0]), len(p2_candidates))
        self.assertEqual(len(idx_res["pool3"][0]), len(p3_candidates))

        # Sample targets:
        # C1 target: 8, C4 target: 4, C2 target: 4, C3 target: 4
        target_c1 = 8
        target_c4 = 4
        target_c2 = 4
        target_c3 = 4

        # C3 Top-K
        p3_ch, p3_r, _ = partition_mod.sample_topk_with_tiebreak(
            idx_res["pool3"][0], idx_res["pool3"][1], idx_res["pool3"][2],
            target_k=target_c3, seed=42
        )
        self.assertEqual(len(p3_ch), target_c3)

        # C2 Top-K
        p2_ch, p2_r, _ = partition_mod.sample_topk_with_tiebreak(
            idx_res["pool2"][0], idx_res["pool2"][1], idx_res["pool2"][2],
            target_k=target_c2, seed=42
        )
        self.assertEqual(len(p2_ch), target_c2)

        # C1/C4 Split
        (c1_ch, c1_r), (c4_ch, c4_r) = partition_mod.resolve_c1_c4_split(
            idx_res["pool14"][0], idx_res["pool14"][1], idx_res["pool14"][2],
            idx_res["chunk_files"], target_c1=target_c1, target_c4=target_c4, v_c1=v_c1, seed=42
        )
        self.assertEqual(len(c1_ch), target_c1)
        self.assertEqual(len(c4_ch), target_c4)

        # Export Negative TSVs
        out_dir = os.path.join(self.temp_dir.name, "out_splits")
        selected_assignments = {
            "C1": (c1_ch, c1_r),
            "C2": (p2_ch, p2_r),
            "C3": (p3_ch, p3_r),
            "C4": (c4_ch, c4_r),
        }
        neg_paths = partition_mod.stream_export_partitions(
            idx_res["chunk_files"], selected_assignments, out_dir, theta="0_50", cmin="0_80"
        )

        # Export Combined Benchmarks
        locus_map = {f"T{i}": f"LOC_T{i}" for i in range(1, 11)}
        locus_map.update({f"U{i}": f"LOC_U{i}" for i in range(1, 7)})
        combined_paths = partition_mod.export_combined_benchmarks(
            out_dir, pos_pairs, neg_paths, locus_map, theta="0_50", cmin="0_80"
        )

        # Verify Integrity (with target_ratio=2.0 for synthetic test)
        audit_results, neg_pairs = partition_mod.verify_partition_integrity(
            pos_pairs, neg_paths, v_train, v_test, v_c1, target_ratio=2.0
        )

        for check in audit_results:
            self.assertTrue(check["passed"], f"Verification check failed: {check['check']} - {check['result']}")

        # Assert combined files have exact line counts (pos + neg + header)
        for part, p_cnt, n_cnt in [("C1", 4, 8), ("C2", 2, 4), ("C3", 2, 4), ("C4", 2, 4)]:
            with open(combined_paths[part], "r", encoding="utf-8") as f:
                lines = f.readlines()
                self.assertEqual(len(lines), p_cnt + n_cnt + 1)  # +1 for header

    def test_sample_topk_capping_behavior(self):
        """
        Verify that sample_topk_with_tiebreak properly caps when allow_cap=True
        and raises ValueError when allow_cap=False.
        """
        chunk_arr = np.array([1, 1, 1], dtype=np.uint16)
        row_arr = np.array([0, 1, 2], dtype=np.uint32)
        score_arr = np.array([100, 200, 150], dtype=np.int64)

        # 1. When allow_cap=True and requested > available
        sel_c, sel_r, sel_s = partition_mod.sample_topk_with_tiebreak(
            chunk_arr, row_arr, score_arr, target_k=10, seed=42, allow_cap=True, pool_name="TestPool"
        )
        self.assertEqual(len(sel_c), 3)
        self.assertEqual(list(sel_s), [200, 150, 100])

        # 2. When allow_cap=False and requested > available -> raises ValueError
        with self.assertRaises(ValueError):
            partition_mod.sample_topk_with_tiebreak(
                chunk_arr, row_arr, score_arr, target_k=10, seed=42, allow_cap=False
            )

    def test_disk_benchmark_partitions_if_present(self):
        """
        If the actual negative partition files are generated in data/dataset_construction/04_benchmark_splits/,
        verifies that they satisfy all 9 strict integrity rules and exact counts.
        """
        splits_dir = os.path.join(REPO_ROOT, "data", "dataset_construction", "04_benchmark_splits", "arabidopsis")
        neg_c1 = os.path.join(splits_dir, "ara_negative_c1_train_id0_50_cov0_80.tsv")
        if not os.path.isfile(neg_c1):
            self.skipTest("Negative partition files not yet generated on disk (run pipeline first).")

        pos_pairs, v_train, v_test, v_c1 = partition_mod.load_positive_splits(splits_dir, theta="0_50", cmin="0_80")
        neg_paths = {
            "C1": neg_c1,
            "C2": os.path.join(splits_dir, "ara_negative_c2_test_id0_50_cov0_80.tsv"),
            "C3": os.path.join(splits_dir, "ara_negative_c3_test_id0_50_cov0_80.tsv"),
            "C4": os.path.join(splits_dir, "ara_negative_c4_test_id0_50_cov0_80.tsv"),
        }

        audit_results, _ = partition_mod.verify_partition_integrity(pos_pairs, neg_paths, v_train, v_test, v_c1)
        for check in audit_results:
            self.assertTrue(check["passed"], f"Real dataset check failed: {check['check']} - {check['result']}")


if __name__ == "__main__":
    unittest.main()
