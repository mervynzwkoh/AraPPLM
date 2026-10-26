#!/usr/bin/env python3
"""
test_negative_interolog_filter.py
=================================
Automated test suite verifying the pair-level homology filtering logic,
mathematical equivalence, and Parquet streaming filter for candidate negative pairs.
"""

import os
import sys
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

# Add repo root and scripts directory to path
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, REPO_ROOT)
sys.path.insert(0, os.path.join(REPO_ROOT, "scripts", "dataset_construction"))

import importlib.util

spec = importlib.util.spec_from_file_location(
    "filter_module",
    os.path.join(REPO_ROOT, "scripts", "dataset_construction", "03_filter_negative_positive_interologs_hpc.py")
)
filter_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(filter_mod)


class TestNegativeInterologFilter(unittest.TestCase):

    def test_pair_similarity_metric(self):
        """
        Verify exact evaluation of the pair similarity metric:
        S(P_neg, P_pos) = max[ min(s(A, C), s(B, D)), min(s(A, D), s(B, C)) ]
        """
        # Case 1: Orientation 1 match (A~C: 0.85, B~D: 0.65) -> S = min(0.85, 0.65) = 0.65
        s_ac, s_bd = 0.85, 0.65
        s_ad, s_bc = 0.10, 0.15
        s_pair = max(min(s_ac, s_bd), min(s_ad, s_bc))
        self.assertAlmostEqual(s_pair, 0.65)
        self.assertGreaterEqual(s_pair, 0.50)  # Filtered!

        # Case 2: Orientation 2 match (A~D: 0.90, B~C: 0.55) -> S = min(0.90, 0.55) = 0.55
        s_ac, s_bd = 0.20, 0.10
        s_ad, s_bc = 0.90, 0.55
        s_pair = max(min(s_ac, s_bd), min(s_ad, s_bc))
        self.assertAlmostEqual(s_pair, 0.55)
        self.assertGreaterEqual(s_pair, 0.50)  # Filtered!

        # Case 3: One-sided homology only (A~C: 0.95, but B~D: 0.0, B~C: 0.0) -> S = 0.0
        s_ac, s_bd = 0.95, 0.0
        s_ad, s_bc = 0.0, 0.0
        s_pair = max(min(s_ac, s_bd), min(s_ad, s_bc))
        self.assertAlmostEqual(s_pair, 0.0)
        self.assertLess(s_pair, 0.50)  # Retained!

    def test_homology_filtering_thresholds(self):
        """Verify strict application of theta=0.50 and cmin=0.80."""
        # Simulated MMseqs2 alignment rows
        raw_rows = [
            # High id, high cov -> PASS
            {"query": "P1", "target": "P2", "fident": "0.60", "qcov": "0.85", "tcov": "0.90", "evalue": "1e-10"},
            # High id, LOW cov -> FAIL (cmin < 0.80)
            {"query": "P1", "target": "P3", "fident": "0.75", "qcov": "0.60", "tcov": "0.85", "evalue": "1e-10"},
            # LOW id, high cov -> FAIL (theta < 0.50)
            {"query": "P1", "target": "P4", "fident": "0.45", "qcov": "0.90", "tcov": "0.95", "evalue": "1e-10"},
            # High id, high cov, HIGH evalue -> FAIL
            {"query": "P1", "target": "P5", "fident": "0.80", "qcov": "0.85", "tcov": "0.85", "evalue": "0.05"},
        ]

        with tempfile.NamedTemporaryFile("w", delete=False, suffix=".tsv", encoding="utf-8") as f:
            f.write("query\ttarget\tfident\talnlen\tmismatch\tgapopen\tqstart\tqend\ttstart\ttend\tevalue\tbits\tqlen\ttlen\tqcov\ttcov\n")
            for r in raw_rows:
                f.write(f"{r['query']}\t{r['target']}\t{r['fident']}\t100\t0\t0\t1\t100\t1\t100\t{r['evalue']}\t200\t100\t100\t{r['qcov']}\t{r['tcov']}\n")
            temp_aln = f.name

        try:
            homologs = filter_mod.build_homology_index(temp_aln, theta=0.50, cmin=0.80, max_evalue=1e-3)
            # P1 should ONLY have P2 as homolog
            self.assertIn("P2", homologs["P1"])
            self.assertNotIn("P3", homologs["P1"])  # Failed coverage
            self.assertNotIn("P4", homologs["P1"])  # Failed identity
            self.assertNotIn("P5", homologs["P1"])  # Failed evalue
        finally:
            if os.path.exists(temp_aln):
                os.remove(temp_aln)

    def test_interolog_set_generation(self):
        """
        Verify exhaustive interolog set compilation:
        If positive pair is (C, D) and H(C) = {C, C_paralog}, H(D) = {D, D_paralog},
        all 4 combinations must be forbidden.
        """
        positive_pairs = [("PROT_C", "PROT_D")]
        homologs = {
            "PROT_C": {"PROT_C_HOM"},
            "PROT_D": {"PROT_D_HOM"},
        }
        all_prots = ["PROT_C", "PROT_D", "PROT_C_HOM", "PROT_D_HOM", "PROT_UNRELATED_1", "PROT_UNRELATED_2"]
        prot_to_idx = {p: np.uint32(i) for i, p in enumerate(all_prots)}

        interolog_array, count = filter_mod.build_positive_interolog_set(positive_pairs, homologs, prot_to_idx)

        # 4 forbidden pairs: (C, D), (C, D_HOM), (C_HOM, D), (C_HOM, D_HOM)
        self.assertEqual(count, 4)

        def is_forbidden(p1, p2):
            u, v = prot_to_idx[p1], prot_to_idx[p2]
            packed = (np.uint64(min(u, v)) << np.uint64(32)) | np.uint64(max(u, v))
            pos = np.searchsorted(interolog_array, packed)
            return pos < len(interolog_array) and interolog_array[pos] == packed

        self.assertTrue(is_forbidden("PROT_C", "PROT_D"))
        self.assertTrue(is_forbidden("PROT_C", "PROT_D_HOM"))
        self.assertTrue(is_forbidden("PROT_C_HOM", "PROT_D"))
        self.assertTrue(is_forbidden("PROT_C_HOM", "PROT_D_HOM"))

        # Unrelated pairs must NOT be forbidden
        self.assertFalse(is_forbidden("PROT_C", "PROT_UNRELATED_1"))
        self.assertFalse(is_forbidden("PROT_C_HOM", "PROT_UNRELATED_1"))
        self.assertFalse(is_forbidden("PROT_UNRELATED_1", "PROT_UNRELATED_2"))

    def test_filter_candidate_chunks_end_to_end(self):
        """Verify streaming Parquet filtering on a synthetic chunk."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cand_dir = os.path.join(tmpdir, "candidates")
            out_dir = os.path.join(tmpdir, "filtered")
            os.makedirs(cand_dir, exist_ok=True)

            # Create synthetic candidate negative parquet chunk
            # 5 pairs:
            # 1. ("POS_C", "POS_D") -> EXACT POSITIVE -> MUST DISCARD
            # 2. ("POS_C_HOM", "POS_D_HOM") -> PARALOG INTEROLOG -> MUST DISCARD
            # 3. ("POS_C", "CLEAN_1") -> ONE-SIDED -> MUST RETAIN
            # 4. ("CLEAN_1", "CLEAN_2") -> UNRELATED -> MUST RETAIN
            # 5. ("CLEAN_2", "POS_D") -> ONE-SIDED -> MUST RETAIN

            chunk_data = {
                "participant_a_id": ["POS_C", "POS_C_HOM", "POS_C", "CLEAN_1", "CLEAN_2"],
                "participant_b_id": ["POS_D", "POS_D_HOM", "CLEAN_1", "CLEAN_2", "POS_D"],
                "locus_a": ["AT1G01", "AT1G02", "AT1G01", "AT2G01", "AT2G02"],
                "locus_b": ["AT3G01", "AT3G02", "AT2G01", "AT2G02", "AT3G01"],
                "compartments_a": ["nucleus", "nucleus", "nucleus", "chloroplast", "chloroplast"],
                "compartments_b": ["cytosol", "cytosol", "chloroplast", "mitochondrion", "cytosol"],
            }

            table = pa.Table.from_arrays(
                [pa.array(chunk_data[col], type=pa.string()) for col in chunk_data],
                names=list(chunk_data.keys()),
            )
            chunk_path = os.path.join(cand_dir, "ara_candidate_negatives_chunk_0001.parquet")
            pq.write_table(table, chunk_path)

            all_prots = ["POS_C", "POS_D", "POS_C_HOM", "POS_D_HOM", "CLEAN_1", "CLEAN_2"]
            prot_to_idx = {p: np.uint32(i) for i, p in enumerate(all_prots)}

            positive_pairs = [("POS_C", "POS_D")]
            homologs = {
                "POS_C": {"POS_C_HOM"},
                "POS_D": {"POS_D_HOM"},
            }

            interolog_array, count = filter_mod.build_positive_interolog_set(positive_pairs, homologs, prot_to_idx)

            metrics = filter_mod.filter_candidate_chunks(
                candidates_dir=cand_dir,
                out_dir=out_dir,
                interolog_array=interolog_array,
                prot_to_idx=prot_to_idx,
                chunk_size=100,
            )

            self.assertEqual(metrics["total_candidate_negatives_evaluated"], 5)
            self.assertEqual(metrics["total_homologous_interologs_removed"], 2)
            self.assertEqual(metrics["total_clean_negatives_retained"], 3)

            # Check output parquet content
            out_files = glob_files = [os.path.join(out_dir, f) for f in os.listdir(out_dir) if f.endswith(".parquet")]
            self.assertEqual(len(out_files), 1)

            out_table = pq.read_table(out_files[0])
            self.assertEqual(len(out_table), 3)

            out_pairs = list(zip(out_table["participant_a_id"].to_pylist(), out_table["participant_b_id"].to_pylist()))
            self.assertIn(("POS_C", "CLEAN_1"), out_pairs)
            self.assertIn(("CLEAN_1", "CLEAN_2"), out_pairs)
            self.assertIn(("CLEAN_2", "POS_D"), out_pairs)
            self.assertNotIn(("POS_C", "POS_D"), out_pairs)
            self.assertNotIn(("POS_C_HOM", "POS_D_HOM"), out_pairs)


if __name__ == "__main__":
    unittest.main()
