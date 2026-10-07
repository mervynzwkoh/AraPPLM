"""
Phase 2 Parser: IntAct MITAB 2.7
Parses PSI-MI TAB 2.7 format for Arabidopsis, Rice, Maize, Tomato, and Soybean.
Preserves all 42 original columns in raw_fields_json.
"""
import os
import sys
import re
import json
import hashlib

WORKSPACE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if WORKSPACE_ROOT not in sys.path:
    sys.path.insert(0, WORKSPACE_ROOT)

from scripts.phase2.harmonize_assays import AssayHarmonizer

INTACT_SOURCES = [
    {
        "source_id": "SRC_A1_INTACT_ARA",
        "species_code": "ARA",
        "file_path": "data/raw/arabidopsis/intact/intact_arabidopsis_taxid3702_release252.mitab27.txt",
        "is_external_holdout": False,
        "is_temporal_holdout": False,
        "eligible_for_training": True,
        "species_default": "Arabidopsis thaliana",
        "taxid_default": 3702,
        "output_subdir": "arabidopsis"
    },
    {
        "source_id": "SRC_B2_INTACT_RICE_39947",
        "species_code": "RICE39947",
        "file_path": "data/raw/rice/intact/intact_rice_taxid39947_release252.mitab27.txt",
        "is_external_holdout": False,
        "is_temporal_holdout": False,
        "eligible_for_training": True,
        "species_default": "Oryza sativa Japonica",
        "taxid_default": 39947,
        "output_subdir": "rice"
    },
    {
        "source_id": "SRC_B2_INTACT_RICE_4530",
        "species_code": "RICE4530",
        "file_path": "data/raw/rice/intact/intact_rice_taxid4530_release252.mitab27.txt",
        "is_external_holdout": False,
        "is_temporal_holdout": False,
        "eligible_for_training": True,
        "species_default": "Oryza sativa",
        "taxid_default": 4530,
        "output_subdir": "rice"
    },
    {
        "source_id": "SRC_C1_INTACT_MAIZE",
        "species_code": "MAIZE",
        "file_path": "data/raw/external_holdout/maize/intact_maize_taxid4577_release252.mitab27.txt",
        "is_external_holdout": True,
        "is_temporal_holdout": False,
        "eligible_for_training": False,
        "species_default": "Zea mays",
        "taxid_default": 4577,
        "output_subdir": "external_holdout/maize"
    },
    {
        "source_id": "SRC_C2_INTACT_TOMATO",
        "species_code": "TOMATO",
        "file_path": "data/raw/external_holdout/tomato/intact_tomato_taxid4081_release252.mitab27.txt",
        "is_external_holdout": True,
        "is_temporal_holdout": False,
        "eligible_for_training": False,
        "species_default": "Solanum lycopersicum",
        "taxid_default": 4081,
        "output_subdir": "external_holdout/tomato"
    },
    {
        "source_id": "SRC_C3_INTACT_SOYBEAN",
        "species_code": "SOYBEAN",
        "file_path": "data/raw/external_holdout/soybean/intact_soybean_taxid3847_release252.mitab27.txt",
        "is_external_holdout": True,
        "is_temporal_holdout": False,
        "eligible_for_training": False,
        "species_default": "Glycine max",
        "taxid_default": 3847,
        "output_subdir": "external_holdout/soybean"
    }
]

MITAB_COLUMN_NAMES = [
    "id_a", "id_b", "alt_ids_a", "alt_ids_b", "aliases_a", "aliases_b",
    "detection_method", "first_author", "publication_id", "taxid_a", "taxid_b",
    "interaction_type", "source_db", "interaction_id", "confidence_values",
    "expansion_method", "biological_role_a", "biological_role_b",
    "experimental_role_a", "experimental_role_b", "type_a", "type_b",
    "xref_a", "xref_b", "interaction_xref", "annotation_a", "annotation_b",
    "interaction_annotation", "host_organism", "interaction_parameters",
    "creation_date", "update_date", "checksum_a", "checksum_b",
    "interaction_checksum", "negative", "feature_a", "feature_b",
    "stoichiometry_a", "stoichiometry_b", "id_method_a", "id_method_b"
]

VIRIDIPLANTAE_TAXIDS = {3702, 39947, 4530, 4577, 4081, 3847, 3694, 3811, 4565, 3711, 3704, 3827}

def parse_id(id_str):
    if not id_str or id_str == "-":
        return "NA", "NA"
    if ":" in id_str:
        ns, val = id_str.split(":", 1)
        return ns.strip(), val.strip()
    return "unresolved", id_str.strip()

def parse_taxid(taxid_str, default_taxid):
    if not taxid_str or taxid_str == "-":
        return default_taxid
    m = re.search(r'taxid:(\d+)', taxid_str)
    if m:
        return int(m.group(1))
    return default_taxid

def parse_symbol(alias_str, default_id):
    if not alias_str or alias_str == "-":
        return default_id
    # Look for psi-mi display_short or gene name
    for item in alias_str.split("|"):
        if "(gene name)" in item or "(display_short)" in item:
            val = item.split("(")[0]
            if ":" in val:
                val = val.split(":", 1)[1]
            return val.strip()
    # Fallback to first alias
    first = alias_str.split("|")[0].split("(")[0]
    if ":" in first:
        first = first.split(":", 1)[1]
    return first.strip()

def parse_intact_source(src_config, harmonizer):
    full_path = os.path.join(WORKSPACE_ROOT, src_config["file_path"])
    records = []
    rejected = []

    if not os.path.exists(full_path):
        raise FileNotFoundError(f"Missing raw file: {full_path}")

    row_number = 0
    with open(full_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            row_number += 1
            raw_line = line.strip()
            if not raw_line:
                continue

            # Header row
            if raw_line.startswith("#"):
                continue

            cols = line.rstrip("\r\n").split("\t")
            if len(cols) < 42:
                rejected.append({
                    "source_id": src_config["source_id"],
                    "source_file": src_config["file_path"],
                    "source_row": row_number,
                    "reason": f"malformed_column_count:{len(cols)}_expected_42",
                    "raw_record": raw_line
                })
                continue

            id_a_raw = cols[0].strip()
            id_b_raw = cols[1].strip()

            if (not id_a_raw or id_a_raw == "-") and (not id_b_raw or id_b_raw == "-"):
                rejected.append({
                    "source_id": src_config["source_id"],
                    "source_file": src_config["file_path"],
                    "source_row": row_number,
                    "reason": "missing_both_participants",
                    "raw_record": raw_line
                })
                continue

            ns_a, val_a = parse_id(id_a_raw)
            ns_b, val_b = parse_id(id_b_raw)

            name_a = parse_symbol(cols[4], val_a)
            name_b = parse_symbol(cols[5], val_b)

            tax_a = parse_taxid(cols[9], src_config["taxid_default"])
            tax_b = parse_taxid(cols[10], src_config["taxid_default"])

            # Interaction ID e.g. intact:EBI-26995264
            int_id_col = cols[13].strip()
            source_rec_id = "NA"
            if int_id_col and int_id_col != "-":
                for item in int_id_col.split("|"):
                    if item.startswith("intact:"):
                        source_rec_id = item.replace("intact:", "")
                        break
                if source_rec_id == "NA":
                    source_rec_id = int_id_col.split("|")[0]

            # Evidence ID (deterministic)
            evidence_id = f"EV_INTACT_{src_config['species_code']}_{row_number:08d}"

            # Assay harmonization
            det_method_raw = cols[6].strip()
            assay_info = harmonizer.harmonize_psi_mi(det_method_raw)

            # Interaction type & semantics
            type_raw = cols[11].strip()
            type_mi = "NA"
            type_name = "NA"
            m_type = re.search(r'psi-mi:"(MI:\d{4})"\(([^)]+)\)', type_raw)
            if m_type:
                type_mi = m_type.group(1)
                type_name = m_type.group(2)

            semantics = harmonizer.infer_interaction_semantics(assay_info, raw_type_str=type_raw)

            # Publication metadata
            pub_col = cols[8].strip()
            pmid = "NA"
            doi = "NA"
            if pub_col and pub_col != "-":
                for item in pub_col.split("|"):
                    if item.startswith("pubmed:"):
                        pmid = item.replace("pubmed:", "").strip()
                    elif item.startswith("doi:"):
                        doi = item.replace("doi:", "").strip()

            # Publication year from author string e.g. "Alexander et al. (2021)"
            author_col = cols[7].strip()
            pub_year = "NA"
            m_yr = re.search(r'\((\d{4})\)', author_col)
            if m_yr:
                pub_year = m_yr.group(1)

            # Confidence (MIscore)
            conf_col = cols[14].strip()
            miscore = None
            if "intact-miscore:" in conf_col:
                m_sc = re.search(r'intact-miscore:([0-9.]+)', conf_col)
                if m_sc:
                    try:
                        miscore = float(m_sc.group(1))
                    except ValueError:
                        miscore = None

            # Roles
            bio_role_a = cols[16].strip() if cols[16] != "-" else "NA"
            bio_role_b = cols[17].strip() if cols[17] != "-" else "NA"
            exp_role_a = cols[18].strip() if cols[18] != "-" else "NA"
            exp_role_b = cols[19].strip() if cols[19] != "-" else "NA"

            # Flags
            same_sp = (tax_a == tax_b)
            plant_a = (tax_a in VIRIDIPLANTAE_TAXIDS)
            plant_b = (tax_b in VIRIDIPLANTAE_TAXIDS)
            plant_plant = (plant_a and plant_b)
            host_pathogen = (plant_a != plant_b)
            taxon_ambig = False

            # Provisional unordered pair key
            pair_tuple = sorted([f"{ns_a}:{val_a}", f"{ns_b}:{val_b}"])
            pair_key = f"{pair_tuple[0]}|{pair_tuple[1]}"

            # Evidence fingerprint (participant A, B, pmid, assay, source)
            fp_str = f"{pair_key}|{pmid}|{assay_info['assay_family']}|{src_config['source_id']}"
            evidence_fp = hashlib.sha256(fp_str.encode("utf-8")).hexdigest()

            # Raw fields dictionary
            raw_dict = {MITAB_COLUMN_NAMES[i]: cols[i] for i in range(len(cols))}
            raw_json = json.dumps(raw_dict, ensure_ascii=False)

            # Throughput
            throughput_class = "unknown"
            if "high throughput" in raw_line.lower() or "two hybrid array" in det_method_raw.lower() or "pooling" in det_method_raw.lower():
                throughput_class = "high_throughput"
            elif "low throughput" in raw_line.lower():
                throughput_class = "low_throughput"

            rec = {
                "evidence_id": evidence_id,
                "source_id": src_config["source_id"],
                "source_resource": "IntAct",
                "source_release": "Release 252",
                "source_file": src_config["file_path"],
                "source_row_number": row_number,
                "source_record_id": source_rec_id,
                "retrieval_date": "2026-09-21",
                "participant_a_original_id": id_a_raw,
                "participant_a_id_namespace": ns_a,
                "participant_a_id_value": val_a,
                "participant_a_original_name": name_a,
                "participant_a_aliases_raw": cols[4],
                "participant_a_taxid": tax_a,
                "participant_a_species": src_config["species_default"] if same_sp else f"taxid_{tax_a}",
                "participant_a_biological_role": bio_role_a,
                "participant_a_experimental_role": exp_role_a,
                "participant_b_original_id": id_b_raw,
                "participant_b_id_namespace": ns_b,
                "participant_b_id_value": val_b,
                "participant_b_original_name": name_b,
                "participant_b_aliases_raw": cols[5],
                "participant_b_taxid": tax_b,
                "participant_b_species": src_config["species_default"] if same_sp else f"taxid_{tax_b}",
                "participant_b_biological_role": bio_role_b,
                "participant_b_experimental_role": exp_role_b,
                "same_species_pair": same_sp,
                "plant_plant_pair": plant_plant,
                "host_pathogen_pair": host_pathogen,
                "taxon_ambiguous": taxon_ambig,
                "unordered_pair_key_provisional": pair_key,
                "interaction_semantics": semantics,
                "interaction_type_raw": type_raw,
                "interaction_type_psi_mi": type_mi,
                "interaction_type_psi_mi_name": type_name,
                "detection_method_raw": assay_info["detection_method_raw"],
                "detection_method_psi_mi": assay_info["detection_method_psi_mi"],
                "detection_method_psi_mi_name": assay_info["detection_method_psi_mi_name"],
                "assay_family": assay_info["assay_family"],
                "assay_supports_direct_binding": assay_info["supports_direct_binding"],
                "assay_supports_physical_association": assay_info["supports_physical_association"],
                "assay_supports_proximity_only": assay_info["supports_proximity_only"],
                "assay_is_genetic": assay_info["assay_is_genetic"],
                "assay_is_computational": assay_info["assay_is_computational"],
                "native_in_planta": assay_info["native_in_planta"],
                "publication_id_raw": pub_col,
                "pmid": pmid,
                "doi": doi,
                "publication_year": pub_year,
                "publication_source": "IntAct_MITAB",
                "throughput_raw": "NA",
                "throughput_class": throughput_class,
                "screen_id": cols[13].split("|")[0],
                "study_id": cols[8].split("|")[0],
                "intact_miscore": miscore,
                "biogrid_score": None,
                "string_combined_score": None,
                "string_physical_combined_score": None,
                "string_experimental_score": None,
                "string_database_score": None,
                "string_textmining_score": None,
                "native_experimental_status": "not_applicable",
                "source_confidence_raw": conf_col,
                "source_confidence_name": "intact-miscore" if miscore is not None else "none",
                "xlms_interprotein_flag": False,
                "xlms_intraprotein_flag": False,
                "xlms_peptide_a": "NA",
                "xlms_peptide_b": "NA",
                "xlms_linker": "NA",
                "xlms_score": None,
                "xlms_evalue": None,
                "xlms_fdr": None,
                "is_external_holdout": True if (tax_a in {4577, 4081, 3847} or tax_b in {4577, 4081, 3847}) else src_config["is_external_holdout"],
                "is_temporal_holdout": src_config["is_temporal_holdout"],
                "eligible_for_training": False if (tax_a in {4577, 4081, 3847} or tax_b in {4577, 4081, 3847}) else src_config["eligible_for_training"],
                "evidence_fingerprint": evidence_fp,
                "possible_cross_database_duplicate": False,
                "needs_manual_review": assay_info["needs_manual_review"],
                "review_reason": assay_info["review_reason"],
                "raw_fields_json": raw_json
            }
            records.append(rec)

    return records, rejected

def run_intact_parsing():
    harmonizer = AssayHarmonizer()
    all_intact_records = []
    all_intact_rejected = []
    stats = {}

    for src in INTACT_SOURCES:
        print(f"Parsing IntAct source: {src['source_id']} ({src['file_path']})...")
        recs, rej = parse_intact_source(src, harmonizer)
        all_intact_records.extend(recs)
        all_intact_rejected.extend(rej)
        stats[src["source_id"]] = {
            "n_input": len(recs) + len(rej),
            "n_normalized": len(recs),
            "n_rejected": len(rej)
        }
        print(f"  -> Normalized: {len(recs)}, Rejected: {len(rej)}")

    return all_intact_records, all_intact_rejected, stats

if __name__ == "__main__":
    recs, rej, stats = run_intact_parsing()
    print("Total IntAct normalized records:", len(recs))
    print("Total IntAct rejected records:", len(rej))
    print("Stats:", stats)
