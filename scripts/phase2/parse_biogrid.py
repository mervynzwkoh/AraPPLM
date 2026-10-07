"""
Phase 2 Parser: BioGRID TAB 3.0
Parses BioGRID TAB 3.0 format for Arabidopsis, Rice, Maize, Tomato, and Soybean.
Preserves both physical and genetic interactions with appropriate semantic tagging.
Preserves all 37 original columns in raw_fields_json.
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

BIOGRID_SOURCES = [
    {
        "source_id": "SRC_A2_BIOGRID_ARA",
        "species_code": "ARA",
        "file_path": "data/raw/arabidopsis/biogrid/BIOGRID-ORGANISM-Arabidopsis_thaliana_Columbia-5.0.261.tab3.txt",
        "is_external_holdout": False,
        "is_temporal_holdout": False,
        "eligible_for_training": True,
        "species_default": "Arabidopsis thaliana (Columbia)",
        "taxid_default": 3702,
        "output_subdir": "arabidopsis"
    },
    {
        "source_id": "SRC_B3_BIOGRID_RICE",
        "species_code": "RICE",
        "file_path": "data/raw/rice/biogrid/BIOGRID-ORGANISM-Oryza_sativa_Japonica-5.0.261.tab3.txt",
        "is_external_holdout": False,
        "is_temporal_holdout": False,
        "eligible_for_training": True,
        "species_default": "Oryza sativa Japonica",
        "taxid_default": 39947,
        "output_subdir": "rice"
    },
    {
        "source_id": "SRC_C1_BIOGRID_MAIZE",
        "species_code": "MAIZE",
        "file_path": "data/raw/external_holdout/maize/BIOGRID-ORGANISM-Zea_mays-5.0.261.tab3.txt",
        "is_external_holdout": True,
        "is_temporal_holdout": False,
        "eligible_for_training": False,
        "species_default": "Zea mays",
        "taxid_default": 4577,
        "output_subdir": "external_holdout/maize"
    },
    {
        "source_id": "SRC_C2_BIOGRID_TOMATO",
        "species_code": "TOMATO",
        "file_path": "data/raw/external_holdout/tomato/BIOGRID-ORGANISM-Solanum_lycopersicum-5.0.261.tab3.txt",
        "is_external_holdout": True,
        "is_temporal_holdout": False,
        "eligible_for_training": False,
        "species_default": "Solanum lycopersicum",
        "taxid_default": 4081,
        "output_subdir": "external_holdout/tomato"
    },
    {
        "source_id": "SRC_C3_BIOGRID_SOYBEAN",
        "species_code": "SOYBEAN",
        "file_path": "data/raw/external_holdout/soybean/BIOGRID-ORGANISM-Glycine_max-5.0.261.tab3.txt",
        "is_external_holdout": True,
        "is_temporal_holdout": False,
        "eligible_for_training": False,
        "species_default": "Glycine max",
        "taxid_default": 3847,
        "output_subdir": "external_holdout/soybean"
    }
]

BIOGRID_COLUMN_NAMES = [
    "biogrid_interaction_id", "entrez_gene_a", "entrez_gene_b",
    "biogrid_id_a", "biogrid_id_b", "systematic_name_a", "systematic_name_b",
    "official_symbol_a", "official_symbol_b", "synonyms_a", "synonyms_b",
    "experimental_system", "experimental_system_type", "author", "pub_source",
    "organism_id_a", "organism_id_b", "throughput", "score", "modification",
    "qualifications", "tags", "source_database", "swissprot_accessions_a",
    "trembl_accessions_a", "refseq_accessions_a", "swissprot_accessions_b",
    "trembl_accessions_b", "refseq_accessions_b", "ontology_term_categories",
    "ontology_term_qualifier_ids", "ontology_term_names", "additional_qualifications",
    "organism_name_a", "organism_name_b", "organism_a", "organism_b"
]

VIRIDIPLANTAE_TAXIDS = {3702, 39947, 4530, 4577, 4081, 3847, 3694, 3811, 4565, 3711, 3704, 3827}

def parse_biogrid_participant_id(cols, idx_sym, idx_sys, idx_sp, idx_tr, idx_entrez, idx_bg):
    # Priority for id_value: Swiss-Prot accession > Systematic Name > Official Symbol > Entrez Gene
    sp = cols[idx_sp].strip()
    if sp and sp != "-":
        acc = sp.split("|")[0].strip()
        return "uniprotkb", acc

    tr = cols[idx_tr].strip()
    if tr and tr != "-":
        acc = tr.split("|")[0].strip()
        return "uniprotkb", acc

    sys_name = cols[idx_sys].strip()
    if sys_name and sys_name != "-":
        return "systematic_name", sys_name

    sym = cols[idx_sym].strip()
    if sym and sym != "-":
        return "gene_symbol", sym

    entrez = cols[idx_entrez].strip()
    if entrez and entrez != "-":
        return "entrezgene", entrez

    bg = cols[idx_bg].strip()
    return "biogrid", bg

def parse_biogrid_source(src_config, harmonizer):
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
            if row_number == 1 and ("#BioGRID" in raw_line or "BIOGRID ID" in raw_line.upper() or raw_line.startswith("#")):
                continue

            cols = line.rstrip("\r\n").split("\t")
            if len(cols) < 37:
                rejected.append({
                    "source_id": src_config["source_id"],
                    "source_file": src_config["file_path"],
                    "source_row": row_number,
                    "reason": f"malformed_column_count:{len(cols)}_expected_37",
                    "raw_record": raw_line
                })
                continue

            biogrid_id = cols[0].strip()
            exp_system = cols[11].strip()
            exp_type = cols[12].strip()

            # Participant A & B
            ns_a, val_a = parse_biogrid_participant_id(cols, 7, 5, 23, 24, 1, 3)
            ns_b, val_b = parse_biogrid_participant_id(cols, 8, 6, 26, 27, 2, 4)

            name_a = cols[7].strip() if cols[7] != "-" else val_a
            name_b = cols[8].strip() if cols[8] != "-" else val_b

            try:
                tax_a = int(cols[15].strip()) if cols[15].strip() not in ["-", ""] else src_config["taxid_default"]
            except ValueError:
                tax_a = src_config["taxid_default"]

            try:
                tax_b = int(cols[16].strip()) if cols[16].strip() not in ["-", ""] else src_config["taxid_default"]
            except ValueError:
                tax_b = src_config["taxid_default"]

            sp_name_a = cols[33].strip() if cols[33] != "-" else src_config["species_default"]
            sp_name_b = cols[34].strip() if cols[34] != "-" else src_config["species_default"]

            # Deterministic evidence ID
            evidence_id = f"EV_BIOGRID_{src_config['species_code']}_{row_number:08d}"

            # Assay harmonization
            assay_info = harmonizer.harmonize_biogrid_system(exp_system, exp_type)
            is_genetic = (exp_type.lower() == "genetic")
            semantics = harmonizer.infer_interaction_semantics(assay_info, raw_type_str=exp_type, is_genetic=is_genetic)

            # Publication metadata
            pub_source = cols[14].strip()
            pmid = "NA"
            if "PUBMED:" in pub_source.upper():
                pmid = re.sub(r'(?i)PUBMED:', '', pub_source).strip()

            author_str = cols[13].strip()
            pub_year = "NA"
            m_yr = re.search(r'\((\d{4})\)', author_str)
            if m_yr:
                pub_year = m_yr.group(1)

            # Throughput
            tp_raw = cols[17].strip()
            tp_class = "unknown"
            if "high throughput" in tp_raw.lower():
                tp_class = "high_throughput"
            elif "low throughput" in tp_raw.lower():
                tp_class = "low_throughput"

            # Score
            score_raw = cols[18].strip()
            biogrid_score = score_raw if (score_raw and score_raw != "-") else None

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

            # Evidence fingerprint
            fp_str = f"{pair_key}|{pmid}|{assay_info['assay_family']}|{src_config['source_id']}"
            evidence_fp = hashlib.sha256(fp_str.encode("utf-8")).hexdigest()

            # Raw fields dictionary
            raw_dict = {BIOGRID_COLUMN_NAMES[i]: cols[i] for i in range(len(cols))}
            raw_json = json.dumps(raw_dict, ensure_ascii=False)

            rec = {
                "evidence_id": evidence_id,
                "source_id": src_config["source_id"],
                "source_resource": "BioGRID",
                "source_release": "5.0.261",
                "source_file": src_config["file_path"],
                "source_row_number": row_number,
                "source_record_id": biogrid_id,
                "retrieval_date": "2026-09-21",
                "participant_a_original_id": f"{ns_a}:{val_a}",
                "participant_a_id_namespace": ns_a,
                "participant_a_id_value": val_a,
                "participant_a_original_name": name_a,
                "participant_a_aliases_raw": cols[9],
                "participant_a_taxid": tax_a,
                "participant_a_species": sp_name_a,
                "participant_a_biological_role": "NA",
                "participant_a_experimental_role": "NA",
                "participant_b_original_id": f"{ns_b}:{val_b}",
                "participant_b_id_namespace": ns_b,
                "participant_b_id_value": val_b,
                "participant_b_original_name": name_b,
                "participant_b_aliases_raw": cols[10],
                "participant_b_taxid": tax_b,
                "participant_b_species": sp_name_b,
                "participant_b_biological_role": "NA",
                "participant_b_experimental_role": "NA",
                "same_species_pair": same_sp,
                "plant_plant_pair": plant_plant,
                "host_pathogen_pair": host_pathogen,
                "taxon_ambiguous": taxon_ambig,
                "unordered_pair_key_provisional": pair_key,
                "interaction_semantics": semantics,
                "interaction_type_raw": exp_type,
                "interaction_type_psi_mi": "NA",
                "interaction_type_psi_mi_name": exp_type,
                "detection_method_raw": exp_system,
                "detection_method_psi_mi": assay_info["detection_method_psi_mi"],
                "detection_method_psi_mi_name": assay_info["detection_method_psi_mi_name"],
                "assay_family": assay_info["assay_family"],
                "assay_supports_direct_binding": assay_info["supports_direct_binding"],
                "assay_supports_physical_association": assay_info["supports_physical_association"],
                "assay_supports_proximity_only": assay_info["supports_proximity_only"],
                "assay_is_genetic": is_genetic,
                "assay_is_computational": False,
                "native_in_planta": assay_info["native_in_planta"],
                "publication_id_raw": pub_source,
                "pmid": pmid,
                "doi": "NA",
                "publication_year": pub_year,
                "publication_source": "BioGRID_TAB3",
                "throughput_raw": tp_raw,
                "throughput_class": tp_class,
                "screen_id": "NA",
                "study_id": author_str,
                "intact_miscore": None,
                "biogrid_score": biogrid_score,
                "string_combined_score": None,
                "string_physical_combined_score": None,
                "string_experimental_score": None,
                "string_database_score": None,
                "string_textmining_score": None,
                "native_experimental_status": "not_applicable",
                "source_confidence_raw": score_raw if biogrid_score else "NA",
                "source_confidence_name": "BioGRID Score" if biogrid_score else "none",
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

def run_biogrid_parsing():
    harmonizer = AssayHarmonizer()
    all_biogrid_records = []
    all_biogrid_rejected = []
    stats = {}

    for src in BIOGRID_SOURCES:
        print(f"Parsing BioGRID source: {src['source_id']} ({src['file_path']})...")
        recs, rej = parse_biogrid_source(src, harmonizer)
        all_biogrid_records.extend(recs)
        all_biogrid_rejected.extend(rej)
        stats[src["source_id"]] = {
            "n_input": len(recs) + len(rej),
            "n_normalized": len(recs),
            "n_rejected": len(rej)
        }
        print(f"  -> Normalized: {len(recs)}, Rejected: {len(rej)}")

    return all_biogrid_records, all_biogrid_rejected, stats

if __name__ == "__main__":
    recs, rej, stats = run_biogrid_parsing()
    print("Total BioGRID normalized records:", len(recs))
    print("Total BioGRID rejected records:", len(rej))
    print("Stats:", stats)
