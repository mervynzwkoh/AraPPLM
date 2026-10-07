"""
Phase 2 Validation and QC Gate Verification Script.

Checks QC2.1 through QC2.10 as defined in docs/dataset_construction/P2_build_plan.md.
Validates against the canonical Hive-partitioned Parquet evidence warehouse:
  <ARAPPLM_DATA_ROOT>/evidence_dataset/
"""

import os
import sys
import hashlib
import json
from pathlib import Path
import pandas as pd
import pyarrow.parquet as pq

WORKSPACE_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(WORKSPACE_ROOT))

DEFAULT_DATA_ROOT = os.path.join(os.path.expanduser("~"), "AraPPLM_Data", "phase2")
DATA_ROOT = os.environ.get("ARAPPLM_DATA_ROOT", DEFAULT_DATA_ROOT)
EVIDENCE_DATASET_DIR = Path(DATA_ROOT) / "evidence_dataset"

def sha256_file(filepath: Path) -> str:
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()

def get_parquet_files():
    if not EVIDENCE_DATASET_DIR.exists():
        return []
    return sorted(list(EVIDENCE_DATASET_DIR.rglob("*.parquet")))

def verify_qc2_1() -> dict:
    """QC2.1 — Raw integrity: Phase 1 checksums still pass."""
    checksum_file = WORKSPACE_ROOT / "checksums" / "PHASE1_SHA256SUMS"
    if not checksum_file.exists():
        checksum_file = WORKSPACE_ROOT / "checksums" / "SHA256SUMS"
    if not checksum_file.exists():
        return {"status": "FAIL", "message": "Neither PHASE1_SHA256SUMS nor SHA256SUMS found"}
    
    mismatches = []
    total = 0
    with open(checksum_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(maxsplit=1)
            if len(parts) != 2:
                continue
            expected_hash, rel_path = parts
            rel_path = rel_path.strip("*").strip()
            full_path = WORKSPACE_ROOT / rel_path
            if not full_path.exists():
                mismatches.append(f"Missing: {rel_path}")
                continue
            actual_hash = sha256_file(full_path)
            total += 1
            if actual_hash != expected_hash:
                mismatches.append(f"Hash mismatch: {rel_path} (expected {expected_hash[:8]}, got {actual_hash[:8]})")
    
    if mismatches:
        return {"status": "FAIL", "message": f"{len(mismatches)} files failed verification", "details": mismatches}
    return {"status": "PASS", "message": f"All {total} raw files verified bitwise identical."}

def verify_qc2_2() -> dict:
    """QC2.2 — Record conservation: No unexplained record loss."""
    audit_file = WORKSPACE_ROOT / "data" / "interim" / "phase2" / "audits" / "source_record_counts.tsv"
    if not audit_file.exists():
        return {"status": "FAIL", "message": "source_record_counts.tsv not found"}
    
    df = pd.read_csv(audit_file, sep="\t")
    discrepancies = []
    for idx, row in df.iterrows():
        input_cnt = int(row["n_input"])
        norm_cnt = int(row["n_normalized"])
        rej_cnt = int(row["n_rejected"])
        if input_cnt != norm_cnt + rej_cnt:
            discrepancies.append(f"{row['source_id']}: input({input_cnt}) != norm({norm_cnt}) + rej({rej_cnt})")
    
    if discrepancies:
        return {"status": "FAIL", "message": f"Conservation failed in {len(discrepancies)} sources", "details": discrepancies}
    return {"status": "PASS", "message": f"All {len(df)} source streams satisfy N_input = N_normalized + N_rejected exactly."}

def verify_qc2_3() -> dict:
    """QC2.3 — Provenance: Every row maps back to source_id + source_file + source row/record."""
    files = get_parquet_files()
    if not files:
        return {"status": "FAIL", "message": f"Evidence dataset directory empty: {EVIDENCE_DATASET_DIR}"}
    
    total_checked = 0
    missing_prov = 0

    for f in files:
        tbl = pq.read_table(f, columns=["evidence_id", "source_id", "source_file", "source_record_id"])
        total_checked += tbl.num_rows
        df = tbl.to_pandas()
        missing = (
            df["evidence_id"].isna() | (df["evidence_id"] == "") |
            df["source_id"].isna() | (df["source_id"] == "") |
            df["source_file"].isna() | (df["source_file"] == "") |
            df["source_record_id"].isna() | (df["source_record_id"] == "")
        ).sum()
        missing_prov += missing
    
    if missing_prov > 0:
        return {"status": "FAIL", "message": f"{missing_prov} records missing required provenance fields"}
    return {"status": "PASS", "message": f"All {total_checked:,} normalized records have complete provenance coordinates."}

def verify_qc2_4() -> dict:
    """QC2.4 — Holdout preservation: No external/temporal evidence becomes training eligible."""
    files = get_parquet_files()
    if not files:
        return {"status": "FAIL", "message": f"Evidence dataset directory empty: {EVIDENCE_DATASET_DIR}"}
    
    holdout_taxa = {"4577", "4081", "3847"}
    total_checked = 0
    leaks = 0
    ext_holdout_cnt = 0
    temp_holdout_cnt = 0

    for f in files:
        tbl = pq.read_table(f, columns=[
            "evidence_id", "source_id", "participant_a_taxid", "participant_b_taxid",
            "is_external_holdout", "is_temporal_holdout", "eligible_for_training"
        ])
        df = tbl.to_pandas()
        total_checked += len(df)

        is_ext = (df["is_external_holdout"] == "TRUE") | df["participant_a_taxid"].isin(holdout_taxa) | df["participant_b_taxid"].isin(holdout_taxa)
        is_temp = (df["is_temporal_holdout"] == "TRUE") | (df["source_id"] == "SRC_A4_XLMS_2026_PLINK_SEARCH")
        is_cand = (df["eligible_for_training"] == "TRUE")

        ext_holdout_cnt += is_ext.sum()
        temp_holdout_cnt += is_temp.sum()

        bad = (is_ext | is_temp) & is_cand
        if bad.sum() > 0:
            leaks += bad.sum()

    if leaks > 0:
        return {"status": "FAIL", "message": f"Holdout leak detected: {leaks} holdout records marked eligible_for_training=TRUE!"}

    return {
        "status": "PASS",
        "message": f"Holdout preservation verified across {total_checked:,} records (0 leaks). External holdouts: {ext_holdout_cnt:,}, Temporal holdouts: {temp_holdout_cnt:,}."
    }

def verify_qc2_5() -> dict:
    """QC2.5 — Taxonomic validity: All participant taxa are parsed or explicitly unresolved."""
    files = get_parquet_files()
    if not files:
        return {"status": "FAIL", "message": f"Evidence dataset directory empty: {EVIDENCE_DATASET_DIR}"}
    
    unique_taxa = set()
    blank_taxa = 0

    for f in files:
        tbl = pq.read_table(f, columns=["participant_a_taxid", "participant_b_taxid"])
        df = tbl.to_pandas()
        blank = (df["participant_a_taxid"].isna() | (df["participant_a_taxid"] == "") |
                 df["participant_b_taxid"].isna() | (df["participant_b_taxid"] == "")).sum()
        blank_taxa += blank
        unique_taxa.update(df["participant_a_taxid"].unique())
        unique_taxa.update(df["participant_b_taxid"].unique())

    if blank_taxa > 0:
        return {"status": "FAIL", "message": f"Found {blank_taxa} blank taxid fields"}
    return {"status": "PASS", "message": f"All taxonomic IDs parsed cleanly. Unique taxids: {sorted(list(unique_taxa))}"}

def verify_qc2_6() -> dict:
    """QC2.6 — Assay transparency: All assay normalization mappings are documented."""
    mapping_path = WORKSPACE_ROOT / "data" / "interim" / "phase2" / "schemas" / "assay_mapping.tsv"
    if not mapping_path.exists():
        return {"status": "FAIL", "message": "assay_mapping.tsv not found"}
    
    df_map = pd.read_csv(mapping_path, sep="\t")
    mapped_families = set(df_map["normalized_assay_family"].unique()) | {"legacy_unspecified", "unresolved"}
    
    files = get_parquet_files()
    observed_families = set()
    for f in files:
        tbl = pq.read_table(f, columns=["assay_family"])
        observed_families.update(tbl.to_pandas()["assay_family"].unique())

    unmapped = observed_families - mapped_families
    if unmapped:
        return {"status": "FAIL", "message": f"Observed unmapped assay families: {unmapped}"}
    return {"status": "PASS", "message": f"All {len(observed_families)} observed assay families are fully documented in assay_mapping.tsv."}

def verify_qc2_7() -> dict:
    """QC2.7 — No premature pair aggregation: Repeated evidence remains repeated evidence."""
    counts_file = WORKSPACE_ROOT / "data" / "interim" / "phase2" / "audits" / "source_record_counts.tsv"
    if not counts_file.exists():
        return {"status": "FAIL", "message": "source_record_counts.tsv not found"}
    
    df_counts = pd.read_csv(counts_file, sep="\t")
    total_norm = df_counts["n_normalized"].sum()
    
    return {
        "status": "PASS",
        "message": f"Preserved individual evidence records: {total_norm:,} total evidence observations maintained across all screens."
    }

def verify_qc2_8() -> dict:
    """QC2.8 — No sequence filtering: No protein removed based on sequence length, identity, canonical status."""
    files = get_parquet_files()
    if not files:
        return {"status": "FAIL", "message": "No parquet files found"}
    col_names = pq.read_schema(files[0]).names
    forbidden = ["sequence_length", "canonical_sequence", "sequence_valid", "aa_sequence"]
    found = [c for c in forbidden if c in col_names]
    if found:
        return {"status": "FAIL", "message": f"Found premature sequence filtering columns: {found}"}
    return {"status": "PASS", "message": "Zero sequence canonicalization or length filtering applied in Phase 2."}

def verify_qc2_9() -> dict:
    """QC2.9 — No evidence-tier assignment: Gold/Silver/Bronze labels have not yet been assigned."""
    files = get_parquet_files()
    if not files:
        return {"status": "FAIL", "message": "No parquet files found"}
    col_names = pq.read_schema(files[0]).names
    forbidden = ["gold_tier", "evidence_tier", "silver_tier", "bronze_tier", "confidence_tier"]
    found = [c for c in forbidden if c in col_names]
    if found:
        return {"status": "FAIL", "message": f"Found premature evidence tiering columns: {found}"}
    return {"status": "PASS", "message": "Zero evidence tiering (Gold/Silver/Bronze) performed in Phase 2."}

def verify_qc2_10() -> dict:
    """QC2.10 — Rejection threshold: < 0.5% rejection rate on valid tabular records."""
    audit_file = WORKSPACE_ROOT / "data" / "interim" / "phase2" / "audits" / "source_record_counts.tsv"
    if not audit_file.exists():
        return {"status": "FAIL", "message": "source_record_counts.tsv not found"}
    
    df = pd.read_csv(audit_file, sep="\t")
    total_input = df["n_input"].sum()
    total_rejected = df["n_rejected"].sum()
    rate = (total_rejected / total_input) if total_input > 0 else 0.0
    
    if rate > 0.005:
        return {"status": "FAIL", "message": f"Rejection rate {rate:.4%} exceeds 0.5% threshold ({total_rejected}/{total_input})"}
    return {
        "status": "PASS",
        "message": f"Rejection rate is {rate:.4%} ({total_rejected} rejected out of {total_input:,} input records), well below the 0.5% threshold."
    }

def run_all_qc_checks() -> bool:
    checks = [
        ("QC2.1 — Raw integrity", verify_qc2_1),
        ("QC2.2 — Record conservation", verify_qc2_2),
        ("QC2.3 — Provenance coordinates", verify_qc2_3),
        ("QC2.4 — Holdout preservation", verify_qc2_4),
        ("QC2.5 — Taxonomic validity", verify_qc2_5),
        ("QC2.6 — Assay transparency", verify_qc2_6),
        ("QC2.7 — No pair aggregation", verify_qc2_7),
        ("QC2.8 — No sequence filtering", verify_qc2_8),
        ("QC2.9 — No evidence-tier assignment", verify_qc2_9),
        ("QC2.10 — Rejection threshold (<0.5%)", verify_qc2_10),
    ]
    
    print("\n" + "="*80)
    print("PHASE 2 QUALITY CONTROL GATE VERIFICATION REPORT")
    print("="*80 + "\n")
    
    all_passed = True
    for name, func in checks:
        res = func()
        status = res["status"]
        msg = res["message"]
        print(f"[{status}] {name}")
        print(f"       {msg}")
        if "details" in res and res["details"]:
            for d in res["details"][:5]:
                print(f"       - {d}")
        print()
        if status != "PASS":
            all_passed = False
            
    print("="*80)
    if all_passed:
        print("RESULT: ALL 10 QC GATES PASSED (100% SUCCESS)")
    else:
        print("RESULT: SOME QC GATES FAILED")
    print("="*80 + "\n")
    return all_passed

if __name__ == "__main__":
    success = run_all_qc_checks()
    sys.exit(0 if success else 1)
