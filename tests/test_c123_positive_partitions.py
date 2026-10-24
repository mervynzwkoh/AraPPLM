#!/usr/bin/env python3
"""
test_c123_positive_partitions.py
=================================
Automated test suite verifying the integrity, disjointness, and C1/C2/C3 criteria
of all partitioned positive PPI datasets across all 8 parameter sweeps using unittest.
"""

import os
import csv
import unittest

SWEEPS = [
    ("0_40", "0_50", 19939),
    ("0_40", "0_70", 20798),
    ("0_40", "0_80", 21480),
    ("0_40", "0_90", 22546),
    ("0_70", "0_50", 26170),
    ("0_70", "0_70", 26307),
    ("0_70", "0_80", 26413),
    ("0_70", "0_90", 26807),
]

BASE_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "data", "dataset_construction", "04_benchmark_splits", "arabidopsis")
)
ORIG_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "data", "dataset_construction", "03_homology_clustering", "arabidopsis")
)


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


class TestPositiveC123Partitions(unittest.TestCase):

    def test_partition_files_exist_and_valid(self):
        """Verify all 24 partition files exist and have non-zero content."""
        for theta, cmin, expected_total in SWEEPS:
            with self.subTest(sweep=f"id{theta}_cov{cmin}"):
                for part in ["c1_train", "c2_test", "c3_test"]:
                    fname = f"ara_positive_{part}_id{theta}_cov{cmin}.tsv"
                    fpath = os.path.join(BASE_DIR, fname)
                    self.assertTrue(os.path.isfile(fpath), f"Missing partition file: {fname}")
                    self.assertGreater(os.path.getsize(fpath), 1000, f"File unexpectedly small: {fname}")

    def test_partition_proportions_and_disjointness(self):
        """Verify exact 80/10/10 edge proportions and zero C1/C3 protein overlap."""
        for theta, cmin, expected_total in SWEEPS:
            with self.subTest(sweep=f"id{theta}_cov{cmin}"):
                c1_path = os.path.join(BASE_DIR, f"ara_positive_c1_train_id{theta}_cov{cmin}.tsv")
                c2_path = os.path.join(BASE_DIR, f"ara_positive_c2_test_id{theta}_cov{cmin}.tsv")
                c3_path = os.path.join(BASE_DIR, f"ara_positive_c3_test_id{theta}_cov{cmin}.tsv")

                c1_pairs = load_tsv_pairs(c1_path)
                c2_pairs = load_tsv_pairs(c2_path)
                c3_pairs = load_tsv_pairs(c3_path)

                # 1. Edge Count Exhaustiveness
                total = len(c1_pairs) + len(c2_pairs) + len(c3_pairs)
                self.assertEqual(total, expected_total, f"Expected {expected_total} pairs, got {total}")

                # 2. Proportion targets: 80% / 10% / 10%
                pct_c1 = len(c1_pairs) / total * 100
                pct_c2 = len(c2_pairs) / total * 100
                pct_c3 = len(c3_pairs) / total * 100
                self.assertAlmostEqual(pct_c1, 80.0, delta=0.1, msg=f"C1 percentage {pct_c1:.2f}% out of bounds")
                self.assertAlmostEqual(pct_c2, 10.0, delta=0.1, msg=f"C2 percentage {pct_c2:.2f}% out of bounds")
                self.assertAlmostEqual(pct_c3, 10.0, delta=0.1, msg=f"C3 percentage {pct_c3:.2f}% out of bounds")

                # 3. Canonical Edge Disjointness
                set_c1 = set((min(a, b), max(a, b)) for a, b in c1_pairs)
                set_c2 = set((min(a, b), max(a, b)) for a, b in c2_pairs)
                set_c3 = set((min(a, b), max(a, b)) for a, b in c3_pairs)

                self.assertEqual(len(set_c1.intersection(set_c2)), 0, "C1 and C2 share pairs!")
                self.assertEqual(len(set_c1.intersection(set_c3)), 0, "C1 and C3 share pairs!")
                self.assertEqual(len(set_c2.intersection(set_c3)), 0, "C2 and C3 share pairs!")

                # 4. Protein Disjointness: Zero overlap between C1 and C3
                c1_prots = set([p[0] for p in c1_pairs] + [p[1] for p in c1_pairs])
                c3_prots = set([p[0] for p in c3_pairs] + [p[1] for p in c3_pairs])
                overlap = c1_prots.intersection(c3_prots)
                self.assertEqual(len(overlap), 0, f"Found {len(overlap)} overlapping proteins between C1 and C3!")

    def test_matching_with_original_positive_dataset(self):
        """Verify partition union matches the original non-redundant positive input."""
        for theta, cmin, expected_total in SWEEPS:
            with self.subTest(sweep=f"id{theta}_cov{cmin}"):
                orig_path = os.path.join(ORIG_DIR, f"ara_positive_pairs_nr_id{theta}_cov{cmin}.tsv")
                orig_pairs = load_tsv_pairs(orig_path)
                self.assertEqual(len(orig_pairs), expected_total)

                orig_set = set((min(a, b), max(a, b)) for a, b in orig_pairs)

                c1_path = os.path.join(BASE_DIR, f"ara_positive_c1_train_id{theta}_cov{cmin}.tsv")
                c2_path = os.path.join(BASE_DIR, f"ara_positive_c2_test_id{theta}_cov{cmin}.tsv")
                c3_path = os.path.join(BASE_DIR, f"ara_positive_c3_test_id{theta}_cov{cmin}.tsv")

                c1_pairs = load_tsv_pairs(c1_path)
                c2_pairs = load_tsv_pairs(c2_path)
                c3_pairs = load_tsv_pairs(c3_path)

                combined_set = set((min(a, b), max(a, b)) for a, b in (c1_pairs + c2_pairs + c3_pairs))
                self.assertEqual(combined_set, orig_set, "Combined partition pairs do not match original dataset!")


if __name__ == "__main__":
    unittest.main()
