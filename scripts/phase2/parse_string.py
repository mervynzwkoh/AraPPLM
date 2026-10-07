"""
Phase 2 Parser: STRING Physical Links Detailed v12.0
Parses physical links for Arabidopsis (3702), Rice (39947), Maize (4577), Tomato (4081), and Soybean (3847).
Preserves individual score channels (experimental, database, textmining, combined_score).
Records native_experimental_status = 'unresolved' in strict accordance with P2_build_plan.md §18.
Omits raw_fields_json for STRING since all raw columns are preserved verbatim in dedicated columns
and the source file is checksummed. Uses bounded chunk streaming for rock-solid memory safety.
"""
import os
import sys
import gzip
import hashlib

WORKSPACE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if WORKSPACE_ROOT not in sys.path:
    sys.path.insert(0, WORKSPACE_ROOT)

STRING_SOURCES = [
    {
        "source_id": "SRC_A3_STRING_ARA_PHYSICAL",
        "file_path": "data/raw/arabidopsis/string/3702.protein.physical.links.detailed.v12.0.txt.gz",
        "species_code": "ARA",
        "species_default": "Arabidopsis thaliana",
        "taxid_default": 3702,
        "is_external_holdout": False,
        "is_temporal_holdout": False,
        "eligible_for_training": True
    },
    {
        "source_id": "SRC_B4_STRING_RICE_PHYSICAL",
        "file_path": "data/raw/rice/string/39947.protein.physical.links.detailed.v12.0.txt.gz",
        "species_code": "RICE",
        "species_default": "Oryza sativa Japonica",
        "taxid_default": 39947,
        "is_external_holdout": False,
        "is_temporal_holdout": False,
        "eligible_for_training": True
    },
    {
        "source_id": "SRC_C1_STRING_MAIZE_PHYSICAL",
        "file_path": "data/raw/external_holdout/maize/4577.protein.physical.links.detailed.v12.0.txt.gz",
        "species_code": "MAIZE",
        "species_default": "Zea mays",
        "taxid_default": 4577,
        "is_external_holdout": True,
        "is_temporal_holdout": False,
        "eligible_for_training": False
    },
    {
        "source_id": "SRC_C2_STRING_TOMATO_PHYSICAL",
        "file_path": "data/raw/external_holdout/tomato/4081.protein.physical.links.detailed.v12.0.txt.gz",
        "species_code": "TOMATO",
        "species_default": "Solanum lycopersicum",
        "taxid_default": 4081,
        "is_external_holdout": True,
        "is_temporal_holdout": False,
        "eligible_for_training": False
    },
    {
        "source_id": "SRC_C3_STRING_SOYBEAN_PHYSICAL",
        "file_path": "data/raw/external_holdout/soybean/3847.protein.physical.links.detailed.v12.0.txt.gz",
        "species_code": "SOYBEAN",
        "species_default": "Glycine max",
        "taxid_default": 3847,
        "is_external_holdout": True,
        "is_temporal_holdout": False,
        "eligible_for_training": False
    }
]

def parse_string_stream(src_config, batch_size=100000, max_rows=None):
    """
    Generator yielding (batch_records, batch_rejected) chunks of size batch_size.
    Omits raw_fields_json to eliminate serialization overhead while preserving full provenance.
    """
    full_path = os.path.join(WORKSPACE_ROOT, src_config["file_path"])
    if not os.path.exists(full_path):
        raise FileNotFoundError(f"Missing raw file: {full_path}")

    batch_records = []
    batch_rejected = []
    row_number = 0
    yielded_count = 0

    source_id = src_config["source_id"]
    source_file = src_config["file_path"]
    species_code = src_config["species_code"]
    species_default = src_config["species_default"]
    taxid_default = src_config["taxid_default"]
    is_external = src_config["is_external_holdout"]
    is_temporal = src_config["is_temporal_holdout"]
    is_eligible = src_config["eligible_for_training"]

    with gzip.open(full_path, "rt", encoding="utf-8", errors="replace") as f:
        header_line = f.readline()
        row_number += 1

        for line in f:
            row_number += 1
            if max_rows and yielded_count + len(batch_records) >= max_rows:
                break

            parts = line.strip().split()
            if len(parts) < 6:
                batch_rejected.append({
                    "source_id": source_id,
                    "source_file": source_file,
                    "source_row": row_number,
                    "reason": f"malformed_column_count:{len(parts)}_expected_6",
                    "raw_record": line.strip()
                })
                continue

            prot1, prot2, exp_sc, db_sc, tm_sc, comb_sc = parts[:6]

            # Fast id splitting
            val1 = prot1.split(".", 1)[1] if "." in prot1 else prot1
            val2 = prot2.split(".", 1)[1] if "." in prot2 else prot2
            tax1 = int(prot1.split(".", 1)[0]) if "." in prot1 and prot1.split(".", 1)[0].isdigit() else taxid_default
            tax2 = int(prot2.split(".", 1)[0]) if "." in prot2 and prot2.split(".", 1)[0].isdigit() else taxid_default

            evidence_id = f"EV_STRING_{species_code}_{row_number:08d}"

            # Canonical unordered pair key
            if val1 <= val2:
                pair_key = f"string:{val1}|string:{val2}"
            else:
                pair_key = f"string:{val2}|string:{val1}"

            # Evidence fingerprint
            evidence_fp = hashlib.sha256(f"{pair_key}|STRING_v12.0|{source_id}".encode("utf-8")).hexdigest()

            rec = {
                "evidence_id": evidence_id,
                "source_id": source_id,
                "source_resource": "STRING",
                "source_release": "v12.0",
                "source_file": source_file,
                "source_row_number": row_number,
                "source_record_id": f"{prot1}_{prot2}",
                "retrieval_date": "2026-09-21",
                "participant_a_original_id": prot1,
                "participant_a_id_namespace": "string",
                "participant_a_id_value": val1,
                "participant_a_original_name": val1,
                "participant_a_aliases_raw": "NA",
                "participant_a_taxid": tax1,
                "participant_a_species": species_default,
                "participant_a_biological_role": "NA",
                "participant_a_experimental_role": "NA",
                "participant_b_original_id": prot2,
                "participant_b_id_namespace": "string",
                "participant_b_id_value": val2,
                "participant_b_original_name": val2,
                "participant_b_aliases_raw": "NA",
                "participant_b_taxid": tax2,
                "participant_b_species": species_default,
                "participant_b_biological_role": "NA",
                "participant_b_experimental_role": "NA",
                "same_species_pair": True,
                "plant_plant_pair": True,
                "host_pathogen_pair": False,
                "taxon_ambiguous": False,
                "unordered_pair_key_provisional": pair_key,
                "interaction_semantics": "computational",
                "interaction_type_raw": "physical_association_channel",
                "interaction_type_psi_mi": "NA",
                "interaction_type_psi_mi_name": "NA",
                "detection_method_raw": "STRING multi-channel physical prediction / evidence integration",
                "detection_method_psi_mi": "NA",
                "detection_method_psi_mi_name": "NA",
                "assay_family": "computational",
                "assay_supports_direct_binding": False,
                "assay_supports_physical_association": True,
                "assay_supports_proximity_only": False,
                "assay_is_genetic": False,
                "assay_is_computational": True,
                "native_in_planta": "unresolved",
                "publication_id_raw": "Szklarczyk et al. (2023)",
                "pmid": "36370105",
                "doi": "10.1093/nar/gkac1000",
                "publication_year": "2023",
                "publication_source": "STRING_v12",
                "throughput_raw": "Computational Integration",
                "throughput_class": "high_throughput",
                "screen_id": "STRING_v12_physical",
                "study_id": "STRING_v12",
                "intact_miscore": None,
                "biogrid_score": None,
                "string_combined_score": int(comb_sc),
                "string_physical_combined_score": int(comb_sc),
                "string_experimental_score": int(exp_sc),
                "string_database_score": int(db_sc),
                "string_textmining_score": int(tm_sc),
                "native_experimental_status": "unresolved",
                "source_confidence_raw": f"combined:{comb_sc};exp:{exp_sc};db:{db_sc};tm:{tm_sc}",
                "source_confidence_name": "STRING Physical Scores",
                "xlms_interprotein_flag": False,
                "xlms_intraprotein_flag": False,
                "xlms_peptide_a": "NA",
                "xlms_peptide_b": "NA",
                "xlms_linker": "NA",
                "xlms_score": None,
                "xlms_evalue": None,
                "xlms_fdr": None,
                "is_external_holdout": is_external,
                "is_temporal_holdout": is_temporal,
                "eligible_for_training": is_eligible,
                "evidence_fingerprint": evidence_fp,
                "possible_cross_database_duplicate": False,
                "needs_manual_review": False,
                "review_reason": "NA",
                "raw_fields_json": "OMITTED_STRING_FIELDS_PRESERVED_IN_COLUMNS"
            }
            batch_records.append(rec)

            if len(batch_records) >= batch_size:
                yield batch_records, batch_rejected
                yielded_count += len(batch_records)
                batch_records = []
                batch_rejected = []

        if batch_records or batch_rejected:
            yield batch_records, batch_rejected

def parse_string_source(src_config, max_rows=None):
    """
    Convenience method retained for tests.
    """
    all_recs = []
    all_rej = []
    for batch_recs, batch_rej in parse_string_stream(src_config, batch_size=50000, max_rows=max_rows):
        all_recs.extend(batch_recs)
        all_rej.extend(batch_rej)
        if max_rows and len(all_recs) >= max_rows:
            all_recs = all_recs[:max_rows]
            break
    return all_recs, all_rej
