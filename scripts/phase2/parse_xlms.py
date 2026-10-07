"""
Phase 2 Parser: 2026 PhoX XL-MS (Trinh et al., 2026)
Parses Total_XL_plink3-2_v3.csv (390,526 records) from PRIDE PXD038379.
Preserves peptide sequences, linker, score, evalue, FDR.
Flags:
  xlms_interprotein_flag (Protein_Type == 'Inter-Protein')
  xlms_intraprotein_flag (Protein_Type == 'Intra-Protein')
  eligible_for_training = FALSE (strict temporal holdout)
  candidate_temporal_holdout = TRUE
Uses streaming generator parse_xlms_stream(batch_size=50000) for low-memory execution.
"""
import os
import sys
import csv
import json
import hashlib

WORKSPACE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if WORKSPACE_ROOT not in sys.path:
    sys.path.insert(0, WORKSPACE_ROOT)

XLMS_SOURCE = {
    "source_id": "SRC_A4_XLMS_2026_PLINK_SEARCH",
    "species_code": "ARA",
    "file_path": "data/raw/arabidopsis/xlms_2026/Total_XL_plink3-2_v3.csv",
    "is_external_holdout": False,
    "is_temporal_holdout": True,
    "eligible_for_training": False
}

def parse_xlms_peptides(peptide_field: str):
    if not peptide_field or peptide_field == "-":
        return "NA", "NA"
    if ")-" in peptide_field:
        parts = peptide_field.split(")-", 1)
        pep_a = parts[0] + ")"
        pep_b = parts[1]
        return pep_a.strip(), pep_b.strip()
    return peptide_field.strip(), "NA"

def parse_xlms_stream(batch_size=50000):
    full_path = os.path.join(WORKSPACE_ROOT, XLMS_SOURCE["file_path"])
    if not os.path.exists(full_path):
        raise FileNotFoundError(f"Missing raw file: {full_path}")

    batch_records = []
    batch_rejected = []
    row_number = 0

    with open(full_path, "r", encoding="utf-8", errors="replace") as f:
        reader = csv.reader(f)
        header = next(reader, None)
        row_number += 1

        for cols in reader:
            row_number += 1
            if not cols or len(cols) == 0:
                continue

            if len(cols) < 32:
                batch_rejected.append({
                    "source_id": XLMS_SOURCE["source_id"],
                    "source_file": XLMS_SOURCE["file_path"],
                    "source_row": row_number,
                    "reason": f"malformed_column_count:{len(cols)}_expected_32",
                    "raw_record": ",".join(cols)
                })
                continue

            order_id = cols[0].strip()
            prot_a = cols[29].strip()
            prot_b = cols[30].strip()
            prot_type = cols[14].strip()
            linker = cols[6].strip()
            score_str = cols[10].strip()
            evalue_str = cols[9].strip()
            pep_raw = cols[4].strip()

            if not prot_a and not prot_b:
                batch_rejected.append({
                    "source_id": XLMS_SOURCE["source_id"],
                    "source_file": XLMS_SOURCE["file_path"],
                    "source_row": row_number,
                    "reason": "missing_both_proteins",
                    "raw_record": ",".join(cols)
                })
                continue

            if not prot_b:
                prot_b = prot_a

            is_inter = (prot_type.lower() == "inter-protein")
            is_intra = (prot_type.lower() == "intra-protein")

            pep_a, pep_b = parse_xlms_peptides(pep_raw)

            try:
                score_val = float(score_str) if score_str else None
            except ValueError:
                score_val = None

            try:
                evalue_val = float(evalue_str) if evalue_str else None
            except ValueError:
                evalue_val = None

            evidence_id = f"EV_XLMS_PHOX2026_{row_number:08d}"

            pair_tuple = sorted([f"uniprot_or_raw:{prot_a}", f"uniprot_or_raw:{prot_b}"])
            pair_key = f"{pair_tuple[0]}|{pair_tuple[1]}"

            fp_str = f"{pair_key}|XLMS_PhoX_2026|{order_id}|{pep_raw}"
            evidence_fp = hashlib.sha256(fp_str.encode("utf-8")).hexdigest()

            raw_dict = {
                "Order": order_id,
                "Spectrum": cols[1].strip() if len(cols) > 1 else "",
                "Peptide": pep_raw,
                "Linker": linker,
                "Score": score_str,
                "E-value": evalue_str,
                "Protein_Type": prot_type,
                "Proteins1": prot_a,
                "Proteins2": prot_b
            }
            raw_json = json.dumps(raw_dict, ensure_ascii=False)

            rec = {
                "evidence_id": evidence_id,
                "source_id": XLMS_SOURCE["source_id"],
                "source_resource": "PRIDE / pLink 3.2",
                "source_release": "Trinh et al., 2026",
                "source_file": XLMS_SOURCE["file_path"],
                "source_row_number": row_number,
                "source_record_id": order_id,
                "retrieval_date": "2026-09-21",
                "participant_a_original_id": prot_a,
                "participant_a_id_namespace": "uniprot_or_raw",
                "participant_a_id_value": prot_a,
                "participant_a_original_name": prot_a,
                "participant_a_aliases_raw": "NA",
                "participant_a_taxid": 3702,
                "participant_a_species": "Arabidopsis thaliana",
                "participant_a_biological_role": "unspecified",
                "participant_a_experimental_role": "unspecified",
                "participant_b_original_id": prot_b,
                "participant_b_id_namespace": "uniprot_or_raw",
                "participant_b_id_value": prot_b,
                "participant_b_original_name": prot_b,
                "participant_b_aliases_raw": "NA",
                "participant_b_taxid": 3702,
                "participant_b_species": "Arabidopsis thaliana",
                "participant_b_biological_role": "unspecified",
                "participant_b_experimental_role": "unspecified",
                "same_species_pair": True,
                "plant_plant_pair": True,
                "host_pathogen_pair": False,
                "taxon_ambiguous": False,
                "unordered_pair_key_provisional": pair_key,
                "interaction_semantics": "proximity",
                "interaction_type_raw": "cross-linking",
                "interaction_type_psi_mi": "MI:0030",
                "interaction_type_psi_mi_name": "cross-linking",
                "detection_method_raw": "cross-linking study",
                "detection_method_psi_mi": "MI:0030",
                "detection_method_psi_mi_name": "cross-linking study",
                "assay_family": "XL_MS",
                "assay_supports_direct_binding": True,
                "assay_supports_physical_association": True,
                "assay_supports_proximity_only": False,
                "assay_is_genetic": False,
                "assay_is_computational": False,
                "native_in_planta": "TRUE",
                "publication_id_raw": "Trinh et al. (2026)",
                "pmid": "39133827",
                "doi": "10.1073/pnas.2519615123",
                "publication_year": "2026",
                "publication_source": "PRIDE_PXD066234",
                "throughput_raw": "High Throughput Cross-linking MS",
                "throughput_class": "high_throughput",
                "screen_id": "PhoX_XLMS_Arabidopsis_2026",
                "study_id": "PXD066234",
                "intact_miscore": None,
                "biogrid_score": None,
                "string_combined_score": None,
                "string_physical_combined_score": None,
                "string_experimental_score": None,
                "string_database_score": None,
                "string_textmining_score": None,
                "native_experimental_status": "native_in_vivo",
                "source_confidence_raw": f"score:{score_str};evalue:{evalue_str}",
                "source_confidence_name": "pLink Score and E-value",
                "xlms_interprotein_flag": is_inter,
                "xlms_intraprotein_flag": is_intra,
                "xlms_peptide_a": pep_a,
                "xlms_peptide_b": pep_b,
                "xlms_linker": linker,
                "xlms_score": score_val,
                "xlms_evalue": evalue_val,
                "xlms_fdr": "0.05",
                "is_external_holdout": False,
                "is_temporal_holdout": True,
                "eligible_for_training": False,
                "evidence_fingerprint": evidence_fp,
                "possible_cross_database_duplicate": False,
                "needs_manual_review": False,
                "review_reason": "NA",
                "raw_fields_json": raw_json
            }
            batch_records.append(rec)

            if len(batch_records) >= batch_size:
                yield batch_records, batch_rejected
                batch_records = []
                batch_rejected = []

        if batch_records or batch_rejected:
            yield batch_records, batch_rejected

def parse_xlms_records():
    all_recs = []
    all_rej = []
    for b_recs, b_rej in parse_xlms_stream(batch_size=50000):
        all_recs.extend(b_recs)
        all_rej.extend(b_rej)
    inter_count = sum(1 for r in all_recs if r["xlms_interprotein_flag"])
    intra_count = sum(1 for r in all_recs if r["xlms_intraprotein_flag"])
    stats = {
        "total_records": len(all_recs),
        "inter_protein_records": inter_count,
        "intra_protein_records": intra_count,
        "rejected_records": len(all_rej)
    }
    return all_recs, all_rej, stats

if __name__ == "__main__":
    recs, rej, stats = parse_xlms_records()
    print("Total XL-MS normalized records:", len(recs))
    print("Total XL-MS rejected records:", len(rej))
    print("Stats:", stats)
