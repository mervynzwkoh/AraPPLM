"""
Unit tests for Phase 2 Parsers and Harmonizers.
Tests IntAct, BioGRID, STRING, and XL-MS parsers against representative frozen records.
"""
import os
import sys
import pytest

WORKSPACE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if WORKSPACE_ROOT not in sys.path:
    sys.path.insert(0, WORKSPACE_ROOT)

from scripts.phase2.harmonize_assays import AssayHarmonizer
from scripts.phase2.parse_intact import parse_intact_source, INTACT_SOURCES
from scripts.phase2.parse_biogrid import parse_biogrid_source, BIOGRID_SOURCES
from scripts.phase2.parse_string import parse_string_source, STRING_SOURCES
from scripts.phase2.parse_xlms import parse_xlms_records, XLMS_SOURCE

@pytest.fixture(scope="module")
def harmonizer():
    return AssayHarmonizer()

# -----------------------------------------------------------------------------
# 1. ASSAY HARMONIZER TESTS
# -----------------------------------------------------------------------------
def test_harmonizer_y2h(harmonizer):
    res = harmonizer.harmonize_psi_mi('psi-mi:"MI:0018"(two hybrid)')
    assert res["assay_family"] == "Y2H"
    assert res["supports_direct_binding"] is True
    assert res["supports_physical_association"] is True
    assert res["native_in_planta"] == "FALSE"
    assert res["needs_manual_review"] is False

def test_harmonizer_coip(harmonizer):
    res = harmonizer.harmonize_psi_mi('psi-mi:"MI:0006"(anti bait coimmunoprecipitation)')
    assert res["assay_family"] == "coIP"
    assert res["supports_direct_binding"] is False
    assert res["supports_physical_association"] is True
    assert res["native_in_planta"] == "TRUE"

def test_harmonizer_biogrid_genetic(harmonizer):
    res = harmonizer.harmonize_biogrid_system("Synthetic Lethality", "genetic")
    assert res["assay_family"] == "genetic"
    assert res["assay_is_genetic"] is True
    assert res["supports_direct_binding"] is False

def test_harmonizer_bifc(harmonizer):
    res = harmonizer.harmonize_psi_mi('psi-mi:"MI:0809"(bimolecular fluorescence complementation)')
    assert res["assay_family"] == "BiFC"
    assert res["supports_proximity_only"] is True

# -----------------------------------------------------------------------------
# 2. INTACT PARSER TESTS
# -----------------------------------------------------------------------------
def test_intact_rice_4530(harmonizer):
    src = next(s for s in INTACT_SOURCES if s["source_id"] == "SRC_B2_INTACT_RICE_4530")
    recs, rej = parse_intact_source(src, harmonizer)
    assert len(recs) == 378
    assert len(rej) == 0
    rec = recs[0]
    assert rec["source_id"] == "SRC_B2_INTACT_RICE_4530"
    assert rec["participant_a_taxid"] == 4530
    assert rec["eligible_for_training"] is True
    assert rec["is_external_holdout"] is False
    assert rec["evidence_id"].startswith("EV_INTACT_RICE4530_")
    assert "raw_fields_json" in rec
    assert len(rec["raw_fields_json"]) > 10

def test_intact_soybean_holdout(harmonizer):
    src = next(s for s in INTACT_SOURCES if s["source_id"] == "SRC_C3_INTACT_SOYBEAN")
    recs, rej = parse_intact_source(src, harmonizer)
    assert len(recs) == 21
    assert len(rej) == 0
    for r in recs:
        assert r["is_external_holdout"] is True
        assert r["eligible_for_training"] is False

# -----------------------------------------------------------------------------
# 3. BIOGRID PARSER TESTS
# -----------------------------------------------------------------------------
def test_biogrid_maize_holdout(harmonizer):
    src = next(s for s in BIOGRID_SOURCES if s["source_id"] == "SRC_C1_BIOGRID_MAIZE")
    recs, rej = parse_biogrid_source(src, harmonizer)
    assert len(recs) == 18
    assert len(rej) == 0
    for r in recs:
        assert r["is_external_holdout"] is True
        assert r["eligible_for_training"] is False
        assert r["source_resource"] == "BioGRID"

def test_biogrid_tomato_holdout(harmonizer):
    src = next(s for s in BIOGRID_SOURCES if s["source_id"] == "SRC_C2_BIOGRID_TOMATO")
    recs, rej = parse_biogrid_source(src, harmonizer)
    assert len(recs) == 141
    assert len(rej) == 0

# -----------------------------------------------------------------------------
# 4. STRING PARSER TESTS
# -----------------------------------------------------------------------------
def test_string_parser_sample():
    src = STRING_SOURCES[0]
    recs, rej = parse_string_source(src, max_rows=50)
    assert len(recs) == 50
    assert len(rej) == 0
    rec = recs[0]
    assert rec["source_resource"] == "STRING"
    assert rec["native_experimental_status"] == "unresolved"
    assert rec["string_combined_score"] is not None
    assert rec["participant_a_id_namespace"] == "string"
    assert rec["eligible_for_training"] is True

def test_string_maize_holdout():
    src = next(s for s in STRING_SOURCES if s["source_id"] == "SRC_C1_STRING_MAIZE_PHYSICAL")
    recs, rej = parse_string_source(src, max_rows=10)
    assert len(recs) == 10
    for r in recs:
        assert r["is_external_holdout"] is True
        assert r["eligible_for_training"] is False

# -----------------------------------------------------------------------------
# 5. XL-MS PARSER TESTS
# -----------------------------------------------------------------------------
def test_xlms_parser():
    recs, rej, stats = parse_xlms_records()
    assert len(recs) == 390526
    assert len(rej) == 0
    rec = recs[0]
    assert rec["source_id"] == "SRC_A4_XLMS_2026_PLINK_SEARCH"
    assert rec["is_temporal_holdout"] is True
    assert rec["eligible_for_training"] is False
    assert rec["assay_family"] == "XL_MS"
    assert rec["native_in_planta"] == "TRUE"
    assert rec["pmid"] == "39133827"
    assert rec["doi"] == "10.1073/pnas.2519615123"
    assert rec["publication_year"] == "2026"
    assert rec["xlms_fdr"] == "0.05"
    assert rec["xlms_interprotein_flag"] in [True, False]
    assert rec["xlms_intraprotein_flag"] in [True, False]
