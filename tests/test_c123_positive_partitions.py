#!/usr/bin/env python3
"""
test_c123_positive_partitions.py
=================================
Automated test suite verifying the integrity, disjointness, and C1/C2/C3/C4 criteria
of the partitioned positive PPI benchmark dataset for id0_50_cov0_80 using unittest.
"""

import os
import csv
import unittest

BASE_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "data", "dataset_construction", "04_benchmark_splits", "arabidopsis")
)
ORIG_PATH = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "data", "dataset_construction", "03_homology_clustering", "arabidopsis", "ara_positive_pairs_nr_id0_50_cov0_80.tsv")
)

THETA = "0_50"
CMIN = "0_80"
EXPECTED_TOTAL = 22897
EXPECTED_C1 = 19462
EXPECTED_C2 = 1145
EXPECTED_C3 = 1145
EXPECTED_C4 = 1145


def load_tsv_pairs(path):
    assert os.path.isfile(path), f"File missing: {path}"
    pairs = []
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        assert "participant_a_id" in reader.fieldnames
        assert "participant_b_id" in reader.fieldnames
        assert "cluster_id" in reader.fieldnames
        for row in reader:
            pairs.append((row["participant_a_id"].strip(), row["participant_b_id"].strip()))
    return pairs


class TestPositiveC1234Partitions(unittest.TestCase):

    def test_partition_files_exist_and_valid(self):
        """Verify all 4 partition files exist and have expected non-zero sizes."""
        for part in ["c1_train", "c2_test", "c3_test", "c4_test"]:
            fname = f"ara_positive_{part}_id{THETA}_cov{CMIN}.tsv"
            fpath = os.path.join(BASE_DIR, fname)
            self.assertTrue(os.path.isfile(fpath), f"Missing partition file: {fname}")
            self.assertGreater(os.path.getsize(fpath), 1000, f"File unexpectedly small: {fname}")

    def test_partition_counts_and_proportions(self):
        """Verify exact 85% C1, 5% C2, 5% C3, 5% C4 pair counts."""
        c1_path = os.path.join(BASE_DIR, f"ara_positive_c1_train_id{THETA}_cov{CMIN}.tsv")
        c2_path = os.path.join(BASE_DIR, f"ara_positive_c2_test_id{THETA}_cov{CMIN}.tsv")
        c3_path = os.path.join(BASE_DIR, f"ara_positive_c3_test_id{THETA}_cov{CMIN}.tsv")
        c4_path = os.path.join(BASE_DIR, f"ara_positive_c4_test_id{THETA}_cov{CMIN}.tsv")

        c1_pairs = load_tsv_pairs(c1_path)
        c2_pairs = load_tsv_pairs(c2_path)
        c3_pairs = load_tsv_pairs(c3_path)
        c4_pairs = load_tsv_pairs(c4_path)

        self.assertEqual(len(c1_pairs), EXPECTED_C1, f"C1 pairs count mismatch: {len(c1_pairs)} vs {EXPECTED_C1}")
        self.assertEqual(len(c2_pairs), EXPECTED_C2, f"C2 pairs count mismatch: {len(c2_pairs)} vs {EXPECTED_C2}")
        self.assertEqual(len(c3_pairs), EXPECTED_C3, f"C3 pairs count mismatch: {len(c3_pairs)} vs {EXPECTED_C3}")
        self.assertEqual(len(c4_pairs), EXPECTED_C4, f"C4 pairs count mismatch: {len(c4_pairs)} vs {EXPECTED_C4}")

        total = len(c1_pairs) + len(c2_pairs) + len(c3_pairs) + len(c4_pairs)
        self.assertEqual(total, EXPECTED_TOTAL, f"Total pairs mismatch: {total} vs {EXPECTED_TOTAL}")

        pct_c1 = len(c1_pairs) / total * 100
        pct_c2 = len(c2_pairs) / total * 100
        pct_c3 = len(c3_pairs) / total * 100
        pct_c4 = len(c4_pairs) / total * 100

        self.assertAlmostEqual(pct_c1, 85.0, delta=0.01)
        self.assertAlmostEqual(pct_c2, 5.0, delta=0.01)
        self.assertAlmostEqual(pct_c3, 5.0, delta=0.01)
        self.assertAlmostEqual(pct_c4, 5.0, delta=0.01)

    def test_edge_disjointness(self):
        """Verify that all partitions are mutually pairwise edge-disjoint."""
        c1_pairs = load_tsv_pairs(os.path.join(BASE_DIR, f"ara_positive_c1_train_id{THETA}_cov{CMIN}.tsv"))
        c2_pairs = load_tsv_pairs(os.path.join(BASE_DIR, f"ara_positive_c2_test_id{THETA}_cov{CMIN}.tsv"))
        c3_pairs = load_tsv_pairs(os.path.join(BASE_DIR, f"ara_positive_c3_test_id{THETA}_cov{CMIN}.tsv"))
        c4_pairs = load_tsv_pairs(os.path.join(BASE_DIR, f"ara_positive_c4_test_id{THETA}_cov{CMIN}.tsv"))

        set_c1 = set((min(a, b), max(a, b)) for a, b in c1_pairs)
        set_c2 = set((min(a, b), max(a, b)) for a, b in c2_pairs)
        set_c3 = set((min(a, b), max(a, b)) for a, b in c3_pairs)
        set_c4 = set((min(a, b), max(a, b)) for a, b in c4_pairs)

        self.assertEqual(len(set_c1.intersection(set_c2)), 0, "C1 and C2 share pairs!")
        self.assertEqual(len(set_c1.intersection(set_c3)), 0, "C1 and C3 share pairs!")
        self.assertEqual(len(set_c1.intersection(set_c4)), 0, "C1 and C4 share pairs!")
        self.assertEqual(len(set_c2.intersection(set_c3)), 0, "C2 and C3 share pairs!")
        self.assertEqual(len(set_c2.intersection(set_c4)), 0, "C2 and C4 share pairs!")
        self.assertEqual(len(set_c3.intersection(set_c4)), 0, "C3 and C4 share pairs!")

    def test_vertex_membership_and_disjointness_criteria(self):
        """Verify vertex set containment, C4 seen criteria, and zero C1/C3 overlap."""
        c1_pairs = load_tsv_pairs(os.path.join(BASE_DIR, f"ara_positive_c1_train_id{THETA}_cov{CMIN}.tsv"))
        c2_pairs = load_tsv_pairs(os.path.join(BASE_DIR, f"ara_positive_c2_test_id{THETA}_cov{CMIN}.tsv"))
        c3_pairs = load_tsv_pairs(os.path.join(BASE_DIR, f"ara_positive_c3_test_id{THETA}_cov{CMIN}.tsv"))
        c4_pairs = load_tsv_pairs(os.path.join(BASE_DIR, f"ara_positive_c4_test_id{THETA}_cov{CMIN}.tsv"))

        c1_prots = set([p[0] for p in c1_pairs] + [p[1] for p in c1_pairs])
        c3_prots = set([p[0] for p in c3_pairs] + [p[1] for p in c3_pairs])
        c4_prots = set([p[0] for p in c4_pairs] + [p[1] for p in c4_pairs])

        # 1. Zero overlap between C1 and C3
        self.assertEqual(len(c1_prots.intersection(c3_prots)), 0, "C1 and C3 share proteins!")

        # 2. Zero overlap between C4 and C3
        self.assertEqual(len(c4_prots.intersection(c3_prots)), 0, "C4 and C3 share proteins!")

        # 3. C4 Both Seen Requirement: every protein in C4 must have appeared in C1 train set
        self.assertTrue(c4_prots.issubset(c1_prots), "C4 contains proteins that did not appear in C1!")

        # 4. C2 boundary cut
        for a, b in c2_pairs:
            # exactly one should be in unseen (c3_prots or V_test)
            in_c3 = (a in c3_prots) or (b in c3_prots)
            self.assertTrue(in_c3, f"C2 pair ({a}, {b}) has no unseen interactor!")

    def test_matches_original_dataset(self):
        """Verify partition union matches the original non-redundant positive input exactly."""
        self.assertTrue(os.path.isfile(ORIG_PATH), f"Original TSV missing: {ORIG_PATH}")
        orig_pairs = load_tsv_pairs(ORIG_PATH)
        orig_set = set((min(a, b), max(a, b)) for a, b in orig_pairs)

        c1_pairs = load_tsv_pairs(os.path.join(BASE_DIR, f"ara_positive_c1_train_id{THETA}_cov{CMIN}.tsv"))
        c2_pairs = load_tsv_pairs(os.path.join(BASE_DIR, f"ara_positive_c2_test_id{THETA}_cov{CMIN}.tsv"))
        c3_pairs = load_tsv_pairs(os.path.join(BASE_DIR, f"ara_positive_c3_test_id{THETA}_cov{CMIN}.tsv"))
        c4_pairs = load_tsv_pairs(os.path.join(BASE_DIR, f"ara_positive_c4_test_id{THETA}_cov{CMIN}.tsv"))

        union_set = set((min(a, b), max(a, b)) for a, b in (c1_pairs + c2_pairs + c3_pairs + c4_pairs))
        self.assertEqual(len(union_set), len(orig_set))
        self.assertEqual(union_set, orig_set, "Union of partitions does not match original non-redundant set!")


if __name__ == "__main__":
    unittest.main()
