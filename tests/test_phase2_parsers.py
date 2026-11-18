"""
Unit tests for Parsers and Assay Harmonizers.
Tests IntAct, BioGRID, STRING, and XL-MS parsers against representative frozen records.
"""
import os
import sys
import unittest

WORKSPACE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if WORKSPACE_ROOT not in sys.path:
    sys.path.insert(0, WORKSPACE_ROOT)

from scripts.dataset_construction.parsing_and_harmonization.harmonize_assays import AssayHarmonizer
from scripts.dataset_construction.parsing_and_harmonization.parse_intact import parse_intact_source, INTACT_SOURCES
from scripts.dataset_construction.parsing_and_harmonization.parse_biogrid import parse_biogrid_source, BIOGRID_SOURCES
from scripts.dataset_construction.parsing_and_harmonization.parse_string import parse_string_source, STRING_SOURCES
from scripts.dataset_construction.parsing_and_harmonization.parse_xlms import parse_xlms_records, XLMS_SOURCE

class TestPhase2Parsers(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.harmonizer = AssayHarmonizer()

    # -------------------------------------------------------------------------
    # 1. ASSAY HARMONIZER TESTS
    # -------------------------------------------------------------------------
    def test_harmonizer_y2h(self):
        res = self.harmonizer.harmonize_psi_mi('psi-mi:"MI:0018"(two hybrid)')
        self.assertEqual(res["assay_family"], "Y2H")
        self.assertTrue(res["supports_direct_binding"])
        self.assertTrue(res["supports_physical_association"])
        self.assertEqual(res["native_in_planta"], "FALSE")
        self.assertFalse(res["needs_manual_review"])

    def test_harmonizer_coip(self):
        res = self.harmonizer.harmonize_psi_mi('psi-mi:"MI:0006"(anti bait coimmunoprecipitation)')
        self.assertEqual(res["assay_family"], "coIP")
        self.assertFalse(res["supports_direct_binding"])
        self.assertTrue(res["supports_physical_association"])
        self.assertEqual(res["native_in_planta"], "TRUE")

    def test_harmonizer_biogrid_genetic(self):
        res = self.harmonizer.harmonize_biogrid_system("Synthetic Lethality", "genetic")
        self.assertEqual(res["assay_family"], "genetic")
        self.assertTrue(res["assay_is_genetic"])
        self.assertFalse(res["supports_direct_binding"])

    def test_harmonizer_bifc(self):
        res = self.harmonizer.harmonize_psi_mi('psi-mi:"MI:0809"(bimolecular fluorescence complementation)')
        self.assertEqual(res["assay_family"], "BiFC")
        self.assertTrue(res["supports_proximity_only"])

    # -------------------------------------------------------------------------
    # 2. INTACT PARSER TESTS
    # -------------------------------------------------------------------------
    def test_intact_rice_4530(self):
        src = next((s for s in INTACT_SOURCES if s["source_id"] == "SRC_B2_INTACT_RICE_4530"), None)
        if not src or not os.path.exists(os.path.join(WORKSPACE_ROOT, src["file_path"])):
            self.skipTest(f"Raw source file not found: {src['file_path'] if src else 'N/A'}")
        recs, rej = parse_intact_source(src, self.harmonizer)
        self.assertEqual(len(recs), 378)
        self.assertEqual(len(rej), 0)
        rec = recs[0]
        self.assertEqual(rec["source_id"], "SRC_B2_INTACT_RICE_4530")
        self.assertEqual(rec["participant_a_taxid"], 4530)
        self.assertTrue(rec["eligible_for_training"])
        self.assertFalse(rec["is_external_holdout"])
        self.assertTrue(rec["evidence_id"].startswith("EV_INTACT_RICE4530_"))
        self.assertIn("raw_fields_json", rec)
        self.assertGreater(len(rec["raw_fields_json"]), 10)

    def test_intact_soybean_holdout(self):
        src = next((s for s in INTACT_SOURCES if s["source_id"] == "SRC_C3_INTACT_SOYBEAN"), None)
        if not src or not os.path.exists(os.path.join(WORKSPACE_ROOT, src["file_path"])):
            self.skipTest(f"Raw source file not found: {src['file_path'] if src else 'N/A'}")
        recs, rej = parse_intact_source(src, self.harmonizer)
        self.assertEqual(len(recs), 21)
        self.assertEqual(len(rej), 0)
        for r in recs:
            self.assertTrue(r["is_external_holdout"])
            self.assertFalse(r["eligible_for_training"])

    # -------------------------------------------------------------------------
    # 3. BIOGRID PARSER TESTS
    # -------------------------------------------------------------------------
    def test_biogrid_maize_holdout(self):
        src = next((s for s in BIOGRID_SOURCES if s["source_id"] == "SRC_C1_BIOGRID_MAIZE"), None)
        if not src or not os.path.exists(os.path.join(WORKSPACE_ROOT, src["file_path"])):
            self.skipTest(f"Raw source file not found: {src['file_path'] if src else 'N/A'}")
        recs, rej = parse_biogrid_source(src, self.harmonizer)
        self.assertEqual(len(recs), 18)
        self.assertEqual(len(rej), 0)
        for r in recs:
            self.assertTrue(r["is_external_holdout"])
            self.assertFalse(r["eligible_for_training"])
            self.assertEqual(r["source_resource"], "BioGRID")

    def test_biogrid_tomato_holdout(self):
        src = next((s for s in BIOGRID_SOURCES if s["source_id"] == "SRC_C2_BIOGRID_TOMATO"), None)
        if not src or not os.path.exists(os.path.join(WORKSPACE_ROOT, src["file_path"])):
            self.skipTest(f"Raw source file not found: {src['file_path'] if src else 'N/A'}")
        recs, rej = parse_biogrid_source(src, self.harmonizer)
        self.assertEqual(len(recs), 141)
        self.assertEqual(len(rej), 0)

    # -------------------------------------------------------------------------
    # 4. STRING PARSER TESTS
    # -------------------------------------------------------------------------
    def test_string_parser_sample(self):
        src = STRING_SOURCES[0]
        if not os.path.exists(os.path.join(WORKSPACE_ROOT, src["file_path"])):
            self.skipTest(f"Raw source file not found: {src['file_path']}")
        recs, rej = parse_string_source(src, max_rows=50)
        self.assertEqual(len(recs), 50)
        self.assertEqual(len(rej), 0)
        rec = recs[0]
        self.assertEqual(rec["source_resource"], "STRING")
        self.assertEqual(rec["native_experimental_status"], "unresolved")
        self.assertIsNotNone(rec["string_combined_score"])
        self.assertEqual(rec["participant_a_id_namespace"], "string")
        self.assertTrue(rec["eligible_for_training"])

    def test_string_maize_holdout(self):
        src = next((s for s in STRING_SOURCES if s["source_id"] == "SRC_C1_STRING_MAIZE_PHYSICAL"), None)
        if not src or not os.path.exists(os.path.join(WORKSPACE_ROOT, src["file_path"])):
            self.skipTest(f"Raw source file not found: {src['file_path'] if src else 'N/A'}")
        recs, rej = parse_string_source(src, max_rows=10)
        self.assertEqual(len(recs), 10)
        for r in recs:
            self.assertTrue(r["is_external_holdout"])
            self.assertFalse(r["eligible_for_training"])

    # -------------------------------------------------------------------------
    # 5. XL-MS PARSER TESTS
    # -------------------------------------------------------------------------
    def test_xlms_parser(self):
        if not os.path.exists(os.path.join(WORKSPACE_ROOT, XLMS_SOURCE["file_path"])):
            self.skipTest(f"Raw source file not found: {XLMS_SOURCE['file_path']}")
        recs, rej, stats = parse_xlms_records()
        self.assertEqual(len(recs), 390526)
        self.assertEqual(len(rej), 0)
        rec = recs[0]
        self.assertEqual(rec["source_id"], "SRC_A4_XLMS_2026_PLINK_SEARCH")
        self.assertTrue(rec["is_temporal_holdout"])
        self.assertFalse(rec["eligible_for_training"])
        self.assertEqual(rec["assay_family"], "XL_MS")
        self.assertEqual(rec["native_in_planta"], "TRUE")
        self.assertEqual(rec["pmid"], "39133827")
        self.assertEqual(rec["doi"], "10.1073/pnas.2519615123")
        self.assertEqual(rec["publication_year"], "2026")
        self.assertEqual(rec["xlms_fdr"], "0.05")
        self.assertIn(rec["xlms_interprotein_flag"], [True, False])
        self.assertIn(rec["xlms_intraprotein_flag"], [True, False])

if __name__ == "__main__":
    unittest.main()
