#!/usr/bin/env python3
"""
test_negative_redundancy_clustering.py
======================================
Automated test suite verifying the pair-level redundancy reduction logic,
component decoupling theorem, priority hierarchy, and Parquet streaming for negatives.
"""

import os
import sys
import gzip
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

# Add repo root to path
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, REPO_ROOT)

import importlib.util

spec = importlib.util.spec_from_file_location(
    "cluster_mod",
    os.path.join(REPO_ROOT, "scripts", "dataset_construction", "04_cluster_candidate_negatives_mmseqs_hpc.py")
)
cluster_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cluster_mod)


class TestNegativeRedundancyClustering(unittest.TestCase):

    def test_homology_components_and_orthogonality(self):
        """
        Verify that connected components partition proteins such that
        cross-component pairs have strictly S(P_1, P_2) == 0.0.
        """
        raw_rows = [
            # A1 and A2 are homologs (Component 1)
            {"query": "A1", "target": "A2", "fident": "0.70", "qcov": "0.90", "tcov": "0.90", "evalue": "1e-15"},
            # B1 and B2 are homologs (Component 2)
            {"query": "B1", "target": "B2", "fident": "0.65", "qcov": "0.85", "tcov": "0.85", "evalue": "1e-10"},
            # C1 and D1 are singletons (Components 3 and 4)
        ]

        with tempfile.NamedTemporaryFile("w", delete=False, suffix=".tsv", encoding="utf-8") as f:
            f.write("query\ttarget\tfident\talnlen\tmismatch\tgapopen\tqstart\tqend\ttstart\ttend\tevalue\tbits\tqlen\ttlen\tqcov\ttcov\n")
            for r in raw_rows:
                f.write(f"{r['query']}\t{r['target']}\t{r['fident']}\t100\t0\t0\t1\t100\t1\t100\t{r['evalue']}\t200\t100\t100\t{r['qcov']}\t{r['tcov']}\n")
            temp_aln = f.name

        try:
            all_prots = ["A1", "A2", "B1", "B2", "C1", "D1"]
            prot_to_idx = {p: np.uint32(i) for i, p in enumerate(all_prots)}

            comp_id_arr, homologs = cluster_mod.build_homology_components(
                temp_aln, prot_to_idx, theta=0.50, cmin=0.80, max_evalue=1e-3
            )

            # A1 and A2 in same component
            self.assertEqual(comp_id_arr[prot_to_idx["A1"]], comp_id_arr[prot_to_idx["A2"]])
            # B1 and B2 in same component
            self.assertEqual(comp_id_arr[prot_to_idx["B1"]], comp_id_arr[prot_to_idx["B2"]])
            # Different components
            self.assertNotEqual(comp_id_arr[prot_to_idx["A1"]], comp_id_arr[prot_to_idx["B1"]])
            self.assertNotEqual(comp_id_arr[prot_to_idx["C1"]], comp_id_arr[prot_to_idx["D1"]])

            # Orthogonality: P1 = (A1, B1), P2 = (C1, D1)
            # CompKeys differ -> S(P1, P2) == 0.0
            k1 = (min(comp_id_arr[prot_to_idx["A1"]], comp_id_arr[prot_to_idx["B1"]]), max(comp_id_arr[prot_to_idx["A1"]], comp_id_arr[prot_to_idx["B1"]]))
            k2 = (min(comp_id_arr[prot_to_idx["C1"]], comp_id_arr[prot_to_idx["D1"]]), max(comp_id_arr[prot_to_idx["C1"]], comp_id_arr[prot_to_idx["D1"]]))
            self.assertNotEqual(k1, k2)

        finally:
            if os.path.exists(temp_aln):
                os.remove(temp_aln)

    def test_end_to_end_mock_clustering(self):
        """
        Verify end-to-end clustering on synthetic negative chunks:
        Pairs in the same paralog micro-cluster are reduced to the top representative.
        Pairs in singleton components are retained.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            cand_dir = os.path.join(tmpdir, "filtered_candidates")
            final_dir = os.path.join(tmpdir, "final_candidates_nr")
            os.makedirs(cand_dir, exist_ok=True)

            # Synthetic alignments: A1 ~ A2 (Component 1), B1 ~ B2 (Component 2), S1, S2 (Singletons)
            aln_rows = [
                {"query": "A1", "target": "A2", "fident": "0.75", "qcov": "0.90", "tcov": "0.90", "evalue": "1e-20"},
                {"query": "B1", "target": "B2", "fident": "0.80", "qcov": "0.85", "tcov": "0.85", "evalue": "1e-20"},
            ]
            aln_tsv = os.path.join(tmpdir, "alignments.tsv")
            with open(aln_tsv, "w", encoding="utf-8") as f:
                f.write("query\ttarget\tfident\talnlen\tmismatch\tgapopen\tqstart\tqend\ttstart\ttend\tevalue\tbits\tqlen\ttlen\tqcov\ttcov\n")
                for r in aln_rows:
                    f.write(f"{r['query']}\t{r['target']}\t{r['fident']}\t100\t0\t0\t1\t100\t1\t100\t{r['evalue']}\t200\t100\t100\t{r['qcov']}\t{r['tcov']}\n")

            # Synthetic canonical stream TSV
            stream_tsv_gz = os.path.join(tmpdir, "stream.tsv.gz")
            with gzip.open(stream_tsv_gz, "wt", encoding="utf-8") as f:
                f.write("Entry\tEntry Name\tGene Names\tGene Names (primary)\tGene Names (ordered locus)\tGene Names (ORF)\tSequence\n")
                f.write("A1\tA1_ARATH\t\t\tAT1G01010\t\t" + "M" * 500 + "\n")  # len 500
                f.write("A2\tA2_ARATH\t\t\tAT1G01020\t\t" + "M" * 300 + "\n")  # len 300
                f.write("B1\tB1_ARATH\t\t\tAT2G01010\t\t" + "M" * 600 + "\n")  # len 600
                f.write("B2\tB2_ARATH\t\t\tAT2G01020\t\t" + "M" * 200 + "\n")  # len 200
                f.write("S1\tS1_ARATH\t\t\tAT3G01010\t\t" + "M" * 400 + "\n")  # singleton
                f.write("S2\tS2_ARATH\t\t\tAT4G01010\t\t" + "M" * 450 + "\n")  # singleton

            # Synthetic TAIR GAF
            tair_gaf_gz = os.path.join(tmpdir, "tair.gaf.gz")
            with gzip.open(tair_gaf_gz, "wt", encoding="utf-8") as f:
                # A1 and B1 have experimental evidence (EXP = weight 3)
                f.write("TAIR\tAT1G01010\tA1\t\tGO:0005634\t\tEXP\t\tC\t\t\t\t\n")
                f.write("TAIR\tAT2G01010\tB1\t\tGO:0005737\t\tEXP\t\tC\t\t\t\t\n")
                # A2 and B2 have computational evidence (IEA = weight 1)
                f.write("TAIR\tAT1G01020\tA2\t\tGO:0005634\t\tIEA\t\tC\t\t\t\t\n")
                f.write("TAIR\tAT2G01020\tB2\t\tGO:0005737\t\tIEA\t\tC\t\t\t\t\n")

            # Synthetic Candidate Chunks:
            # Group 1: 4 paralogous pairs connecting Comp(A) and Comp(B):
            #   (A1, B1) -> Weight = 3+3=6, Len = 1100 -> HIGHEST PRIORITY -> MUST BE REP!
            #   (A1, B2) -> Weight = 3+1=4, Len = 700 -> PRUNED
            #   (A2, B1) -> Weight = 1+3=4, Len = 900 -> PRUNED
            #   (A2, B2) -> Weight = 1+1=2, Len = 500 -> PRUNED
            # Group 2: Singleton pair (S1, S2) -> Weight = 1+1=2, Len = 850 -> MUST BE RETAINED!
            chunk_data = {
                "participant_a_id": ["A1", "A1", "A2", "A2", "S1"],
                "participant_b_id": ["B1", "B2", "B1", "B2", "S2"],
                "locus_a": ["AT1G01010", "AT1G01010", "AT1G01020", "AT1G01020", "AT3G01010"],
                "locus_b": ["AT2G01010", "AT2G01020", "AT2G01010", "AT2G01020", "AT4G01010"],
                "compartments_a": ["nucleus", "nucleus", "nucleus", "nucleus", "chloroplast"],
                "compartments_b": ["cytosol", "cytosol", "cytosol", "cytosol", "mitochondrion"],
            }
            chunk_table = pa.Table.from_arrays(
                [pa.array(chunk_data[k], type=pa.string()) for k in chunk_data],
                names=list(chunk_data.keys()),
            )
            pq.write_table(chunk_table, os.path.join(cand_dir, "ara_candidate_negatives_filtered_chunk_0001.parquet"))

            audit_path = os.path.join(tmpdir, "audit.json")
            summary_md = os.path.join(tmpdir, "summary.md")

            # Mock args
            class MockArgs:
                pass

            args = MockArgs()
            args.stream_path = stream_tsv_gz
            args.tair_gaf_path = tair_gaf_gz
            args.alignments_tsv = aln_tsv
            args.candidates_dir = cand_dir
            args.out_dir = final_dir
            args.audit_path = audit_path
            args.summary_md = summary_md
            args.theta = 0.50
            args.cmin = 0.80
            args.max_evalue = 1e-3
            args.chunk_size = 100

            cluster_mod.run_candidate_negative_clustering(args)

            # Verification:
            # Out of 5 input pairs:
            # 2 representatives retained: (A1, B1) and (S1, S2)
            # 3 pairs pruned: (A1, B2), (A2, B1), (A2, B2)
            out_files = [os.path.join(final_dir, f) for f in os.listdir(final_dir) if f.endswith(".parquet")]
            self.assertEqual(len(out_files), 1)

            out_table = pq.read_table(out_files[0])
            self.assertEqual(len(out_table), 2)

            reps = list(zip(out_table["participant_a_id"].to_pylist(), out_table["participant_b_id"].to_pylist()))
            self.assertIn(("A1", "B1"), reps)
            self.assertIn(("S1", "S2"), reps)
            self.assertNotIn(("A1", "B2"), reps)
            self.assertNotIn(("A2", "B1"), reps)
            self.assertNotIn(("A2", "B2"), reps)

            # Check audit JSON
            with open(audit_path, "r", encoding="utf-8") as f:
                audit = json.load(f)
            self.assertEqual(audit["total_clean_candidates_evaluated"], 5)
            self.assertEqual(audit["total_representative_negatives_retained"], 2)
            self.assertEqual(audit["total_redundant_pairs_pruned"], 3)
            self.assertEqual(audit["unique_proteins_in_negative_pool"], 4)


if __name__ == "__main__":
    unittest.main()
