"""
Collate Positive Sample Dataset for Rice (Oryza sativa)
======================================================
Parses raw interaction data for Rice from:
1. IntAct (taxid 39947 and 4530) - retains physical association interaction type
2. BioGRID (taxid 39947) - excludes proximity and genetic interaction types
3. STRING (taxid 39947, physical links full) - retains direct experimental interactions (experiments > 0),
   strictly excluding transferred interologs (experiments == 0)

Retains ONLY genuine database provenance fields (eliminating internal guideline
and benchmark-split artifacts), deduplicates into unique positive pairs, and
exports clean provenance tables and benchmark datasets.
"""

import os
import sys
import json
import gzip
import argparse
import collections
from datetime import datetime, timezone
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

WORKSPACE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if WORKSPACE_ROOT not in sys.path:
    sys.path.insert(0, WORKSPACE_ROOT)

from scripts.phase2.harmonize_assays import AssayHarmonizer
from scripts.phase2.parse_intact import INTACT_SOURCES, parse_intact_source
from scripts.phase2.parse_biogrid import BIOGRID_SOURCES, parse_biogrid_source


def infer_molecule_type(r, side="a"):
    """Infer participant molecule type directly from database fields."""
    raw_str = r.get("raw_fields_json", "{}")
    try:
        raw_dict = json.loads(raw_str)
    except Exception:
        raw_dict = {}

    if r.get("source_resource") in ("BioGRID", "STRING"):
        return "protein"

    t = raw_dict.get(f"type_{side}", "")
    if "MI:0326" in t or "protein" in t.lower():
        return "protein"
    elif "MI:0320" in t or "ribonucleic" in t.lower():
        return "ribonucleic_acid"
    elif "MI:0319" in t or "deoxyribonucleic" in t.lower():
        return "deoxyribonucleic_acid"
    return "protein"


def extract_first_author(r):
    """Extract primary author citation from raw database fields."""
    try:
        raw_dict = json.loads(r.get("raw_fields_json", "{}"))
        return raw_dict.get("first_author", raw_dict.get("author", "NA"))
    except Exception:
        return "NA"


def build_rice_positives(output_dir=None, interim_audit=None):
    if output_dir is None:
        output_dir = os.path.join(WORKSPACE_ROOT, "data", "processed", "rice")
    os.makedirs(output_dir, exist_ok=True)

    harmonizer = AssayHarmonizer()

    # ---------------------------------------------------------
    # 1. Ingest Raw IntAct Sources
    # ---------------------------------------------------------
    intact_configs = [s for s in INTACT_SOURCES if s["species_code"] in ("RICE39947", "RICE4530")]
    raw_intact_records = []
    raw_intact_counts = {}
    for cfg in intact_configs:
        recs, rejs = parse_intact_source(cfg, harmonizer)
        raw_intact_counts[cfg["source_id"]] = {
            "file_path": cfg["file_path"],
            "parsed_records": len(recs),
            "rejected_records": len(rejs),
            "species_code": cfg["species_code"],
            "taxid": cfg["taxid_default"]
        }
        raw_intact_records.extend(recs)

    # Filter IntAct: physical association interaction type
    intact_positive = [r for r in raw_intact_records if r["interaction_semantics"] == "physical_association"]
    intact_excluded = [r for r in raw_intact_records if r["interaction_semantics"] != "physical_association"]

    # ---------------------------------------------------------
    # 2. Ingest Raw BioGRID Source
    # ---------------------------------------------------------
    bg_configs = [s for s in BIOGRID_SOURCES if s["species_code"] == "RICE"]
    raw_bg_records = []
    raw_bg_counts = {}
    for cfg in bg_configs:
        recs, rejs = parse_biogrid_source(cfg, harmonizer)
        raw_bg_counts[cfg["source_id"]] = {
            "file_path": cfg["file_path"],
            "parsed_records": len(recs),
            "rejected_records": len(rejs),
            "species_code": cfg["species_code"],
            "taxid": cfg["taxid_default"]
        }
        raw_bg_records.extend(recs)

    # Filter BioGRID: exclude proximity and genetic interaction types
    bg_positive = [r for r in raw_bg_records if r["interaction_semantics"] not in ("proximity", "genetic")]
    bg_excluded = [r for r in raw_bg_records if r["interaction_semantics"] in ("proximity", "genetic")]

    # ---------------------------------------------------------
    # 3. Ingest STRING Direct Physical Interactions (Full Decoupled File)
    # ---------------------------------------------------------
    string_file_path = os.path.join(WORKSPACE_ROOT, "data", "raw", "rice", "string", "39947.protein.physical.links.full.v12.0.txt.gz")
    aliases_path = os.path.join(WORKSPACE_ROOT, "data", "raw", "rice", "string", "39947.protein.aliases.v12.0.txt.gz")

    aliases_dict = collections.defaultdict(lambda: {"gene_name": None, "locus_tag": None, "aliases": set()})
    if os.path.exists(aliases_path):
        with gzip.open(aliases_path, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    sid, alias, src = parts[0], parts[1], parts[2]
                    ac = sid.replace("39947.", "")
                    aliases_dict[ac]["aliases"].add(alias)
                    if src == "UniProt_GN_Name" and not aliases_dict[ac]["gene_name"]:
                        aliases_dict[ac]["gene_name"] = alias
                    elif src == "UniProt_GN_OrderedLocusNames" and not aliases_dict[ac]["locus_tag"]:
                        aliases_dict[ac]["locus_tag"] = alias

    string_positive_evidence = []
    string_total_raw = 0
    string_excluded_interologs = 0
    string_excluded_neither = 0

    if os.path.exists(string_file_path):
        with gzip.open(string_file_path, "rt", encoding="utf-8", errors="replace") as f:
            header = f.readline().strip().split()
            # columns: protein1 protein2 homology experiments experiments_transferred database database_transferred textmining textmining_transferred combined_score
            for line_no, line in enumerate(f, start=2):
                parts = line.strip().split()
                if len(parts) < 10:
                    continue
                string_total_raw += 1
                p1_full, p2_full = parts[0], parts[1]
                exp = int(parts[3])
                exp_tr = int(parts[4])
                db = int(parts[5])
                db_tr = int(parts[6])
                tm = int(parts[7])
                tm_tr = int(parts[8])
                comb = int(parts[9])

                if exp > 0:
                    p_a = p1_full.replace("39947.", "")
                    p_b = p2_full.replace("39947.", "")
                    pair_key = f"{min(p_a, p_b)}|{max(p_a, p_b)}"

                    alias_a = aliases_dict[p_a]
                    alias_b = aliases_dict[p_b]
                    name_a = alias_a["gene_name"] or alias_a["locus_tag"] or p_a
                    name_b = alias_b["gene_name"] or alias_b["locus_tag"] or p_b
                    aliases_a_str = "|".join(sorted(alias_a["aliases"])) if alias_a["aliases"] else "NA"
                    aliases_b_str = "|".join(sorted(alias_b["aliases"])) if alias_b["aliases"] else "NA"

                    rec = {
                        "source_resource": "STRING",
                        "source_release": "v12.0",
                        "source_file": "data/raw/rice/string/39947.protein.physical.links.full.v12.0.txt.gz",
                        "source_row_number": line_no,
                        "source_record_id": f"{p1_full}-{p2_full}",
                        "canonical_pair_key": pair_key,
                        "participant_a_id": p_a,
                        "participant_a_namespace": "uniprotkb",
                        "participant_a_original_id": p1_full,
                        "participant_a_name": name_a,
                        "participant_a_aliases": aliases_a_str,
                        "participant_a_taxid": 39947,
                        "participant_a_species": "Oryza sativa Japonica Group",
                        "participant_a_biological_role": "unspecified",
                        "participant_a_experimental_role": "unspecified",
                        "participant_a_molecule_type": "protein",
                        "participant_b_id": p_b,
                        "participant_b_namespace": "uniprotkb",
                        "participant_b_original_id": p2_full,
                        "participant_b_name": name_b,
                        "participant_b_aliases": aliases_b_str,
                        "participant_b_taxid": 39947,
                        "participant_b_species": "Oryza sativa Japonica Group",
                        "participant_b_biological_role": "unspecified",
                        "participant_b_experimental_role": "unspecified",
                        "participant_b_molecule_type": "protein",
                        "is_protein_protein": "TRUE",
                        "same_species_pair": "TRUE",
                        "interaction_type_raw": "physical association",
                        "interaction_type_psi_mi": "MI:0915",
                        "interaction_type_name": "physical association",
                        "detection_method_raw": f"direct experimental assay (STRING experiments score {exp})",
                        "detection_method_psi_mi": "MI:0045",
                        "detection_method_name": "experimental interaction detection",
                        "first_author": "Szklarczyk et al. (2023)",
                        "publication_id_raw": "PMID:36370105",
                        "pmid": "36370105",
                        "publication_year": "2023",
                        "throughput_raw": "high_throughput",
                        "source_confidence_raw": f"experiments:{exp};combined:{comb};exp_transferred:{exp_tr};database:{db};textmining:{tm}",
                        "intact_miscore": "NA"
                    }
                    string_positive_evidence.append(rec)
                else:
                    if exp_tr > 0:
                        string_excluded_interologs += 1
                    else:
                        string_excluded_neither += 1

    # ---------------------------------------------------------
    # 4. Clean Evidence Records: Retain ONLY Database Provenance
    # ---------------------------------------------------------
    clean_evidence_records = []
    
    # Process IntAct and BioGRID
    for r in (intact_positive + bg_positive):
        mol_a = infer_molecule_type(r, "a")
        mol_b = infer_molecule_type(r, "b")
        p_a = str(r["participant_a_id_value"]).strip()
        p_b = str(r["participant_b_id_value"]).strip()
        pair_key = f"{min(p_a, p_b)}|{max(p_a, p_b)}"

        author = extract_first_author(r)

        rec = {
            "source_resource": r["source_resource"],
            "source_release": r["source_release"],
            "source_file": r["source_file"],
            "source_row_number": r["source_row_number"],
            "source_record_id": r["source_record_id"],
            "canonical_pair_key": pair_key,
            "participant_a_id": p_a,
            "participant_a_namespace": r["participant_a_id_namespace"],
            "participant_a_original_id": r["participant_a_original_id"],
            "participant_a_name": r["participant_a_original_name"],
            "participant_a_aliases": r["participant_a_aliases_raw"],
            "participant_a_taxid": r["participant_a_taxid"],
            "participant_a_species": r["participant_a_species"],
            "participant_a_biological_role": r["participant_a_biological_role"],
            "participant_a_experimental_role": r["participant_a_experimental_role"],
            "participant_a_molecule_type": mol_a,
            "participant_b_id": p_b,
            "participant_b_namespace": r["participant_b_id_namespace"],
            "participant_b_original_id": r["participant_b_original_id"],
            "participant_b_name": r["participant_b_original_name"],
            "participant_b_aliases": r["participant_b_aliases_raw"],
            "participant_b_taxid": r["participant_b_taxid"],
            "participant_b_species": r["participant_b_species"],
            "participant_b_biological_role": r["participant_b_biological_role"],
            "participant_b_experimental_role": r["participant_b_experimental_role"],
            "participant_b_molecule_type": mol_b,
            "is_protein_protein": "TRUE" if (mol_a == "protein" and mol_b == "protein") else "FALSE",
            "same_species_pair": "TRUE" if r["same_species_pair"] else "FALSE",
            "interaction_type_raw": r["interaction_type_raw"],
            "interaction_type_psi_mi": r["interaction_type_psi_mi"],
            "interaction_type_name": r["interaction_type_psi_mi_name"],
            "detection_method_raw": r["detection_method_raw"],
            "detection_method_psi_mi": r["detection_method_psi_mi"],
            "detection_method_name": r["detection_method_psi_mi_name"],
            "first_author": author,
            "publication_id_raw": r["publication_id_raw"],
            "pmid": r["pmid"],
            "publication_year": r["publication_year"],
            "throughput_raw": r["throughput_raw"],
            "source_confidence_raw": r["source_confidence_raw"],
            "intact_miscore": r["intact_miscore"]
        }
        clean_evidence_records.append(rec)

    # Append STRING direct evidence records
    clean_evidence_records.extend(string_positive_evidence)

    df_evidence = pd.DataFrame(clean_evidence_records)

    evidence_tsv_path = os.path.join(output_dir, "rice_positive_evidence.tsv")
    evidence_parquet_path = os.path.join(output_dir, "rice_positive_evidence.parquet")
    df_evidence.to_csv(evidence_tsv_path, sep="\t", index=False)
    
    arrow_schema = pa.schema([(c, pa.string()) for c in df_evidence.columns])
    table_evidence = pa.Table.from_pandas(df_evidence.astype(str), schema=arrow_schema)
    pq.write_table(table_evidence, evidence_parquet_path, compression="snappy")

    # ---------------------------------------------------------
    # 5. Deduplicate into Unique Positive Pairs (Pure DB Provenance)
    # ---------------------------------------------------------
    pair_evidence_map = collections.defaultdict(list)
    for r in clean_evidence_records:
        pair_tuple = tuple(r["canonical_pair_key"].split("|"))
        pair_evidence_map[pair_tuple].append(r)

    sorted_pair_keys = sorted(pair_evidence_map.keys(), key=lambda p: (p[0], p[1]))

    dedup_rows = []
    for idx, pair_key in enumerate(sorted_pair_keys, start=1):
        ev_list = pair_evidence_map[pair_key]
        p_a, p_b = pair_key
        is_homo = (p_a == p_b)

        sources = sorted(list(set(r["source_resource"] for r in ev_list)))
        int_types = sorted(list(set(r["interaction_type_name"] for r in ev_list if r["interaction_type_name"] != "NA")))
        det_methods = sorted(list(set(r["detection_method_name"] for r in ev_list if r["detection_method_name"] != "NA")))
        authors = sorted(list(set(r["first_author"] for r in ev_list if r["first_author"] != "NA")))
        
        pmids = sorted(list(set(
            str(r["pmid"]) for r in ev_list 
            if r["pmid"] and str(r["pmid"]) not in ("NA", "unassigned", "None", "")
        )))
        
        years = sorted(list(set(
            str(r["publication_year"]) for r in ev_list 
            if r["publication_year"] and str(r["publication_year"]) not in ("NA", "None", "")
        )))

        tp = sorted(list(set(str(r["throughput_raw"]) for r in ev_list if str(r["throughput_raw"]) not in ("NA", "None", ""))))
        src_rec_ids = sorted(list(set(str(r["source_record_id"]) for r in ev_list if str(r["source_record_id"]) != "NA")))

        # Names / Symbols from DB
        names_a = set(r["participant_a_name"] for r in ev_list if r["participant_a_id"] == p_a)
        names_a.update(r["participant_b_name"] for r in ev_list if r["participant_b_id"] == p_a)
        names_a.discard("NA")

        names_b = set(r["participant_a_name"] for r in ev_list if r["participant_a_id"] == p_b)
        names_b.update(r["participant_b_name"] for r in ev_list if r["participant_b_id"] == p_b)
        names_b.discard("NA")

        # Molecule types from DB
        types_a = set(r["participant_a_molecule_type"] for r in ev_list if r["participant_a_id"] == p_a)
        types_a.update(r["participant_b_molecule_type"] for r in ev_list if r["participant_b_id"] == p_a)

        types_b = set(r["participant_a_molecule_type"] for r in ev_list if r["participant_a_id"] == p_b)
        types_b.update(r["participant_b_molecule_type"] for r in ev_list if r["participant_b_id"] == p_b)

        type_a_str = "|".join(sorted(types_a)) if types_a else "protein"
        type_b_str = "|".join(sorted(types_b)) if types_b else "protein"

        # IntAct MIscore
        miscores = []
        for r in ev_list:
            if r["intact_miscore"] not in ("NA", "None", None):
                try:
                    miscores.append(float(r["intact_miscore"]))
                except ValueError:
                    pass
        max_miscore = f"{max(miscores):.4f}" if miscores else "NA"

        dedup_rows.append({
            "pair_id": f"RICE_POS_{idx:04d}",
            "canonical_pair_key": f"{p_a}|{p_b}",
            "participant_a_id": p_a,
            "participant_b_id": p_b,
            "participant_a_name": ";".join(sorted(names_a)) if names_a else p_a,
            "participant_b_name": ";".join(sorted(names_b)) if names_b else p_b,
            "participant_a_type": type_a_str,
            "participant_b_type": type_b_str,
            "is_protein_protein_pair": "TRUE" if (type_a_str == "protein" and type_b_str == "protein") else "FALSE",
            "is_homomeric": "TRUE" if is_homo else "FALSE",
            "pair_type": "homomer" if is_homo else "heteromer",
            "evidence_count": len(ev_list),
            "sources": "|".join(sources),
            "interaction_types": "|".join(int_types) if int_types else "physical association",
            "detection_methods": "|".join(det_methods),
            "first_authors": "|".join(authors) if authors else "NA",
            "pmid_count": len(pmids),
            "pmids": "|".join(pmids) if pmids else "NA",
            "publication_years": "|".join(years) if years else "NA",
            "throughput_classes": "|".join(tp) if tp else "unassigned",
            "max_intact_miscore": max_miscore,
            "source_record_ids": "|".join(src_rec_ids)
        })

    df_dedup = pd.DataFrame(dedup_rows)

    dedup_tsv_path = os.path.join(output_dir, "rice_positive_pairs_deduplicated.tsv")
    dedup_parquet_path = os.path.join(output_dir, "rice_positive_pairs_deduplicated.parquet")
    df_dedup.to_csv(dedup_tsv_path, sep="\t", index=False)
    
    arrow_schema_dedup = pa.schema([(c, pa.string()) for c in df_dedup.columns])
    table_dedup = pa.Table.from_pandas(df_dedup.astype(str), schema=arrow_schema_dedup)
    pq.write_table(table_dedup, dedup_parquet_path, compression="snappy")

    # ---------------------------------------------------------
    # 6. Benchmark Dataset Generation (Deferred)
    # ---------------------------------------------------------
    # Benchmark dataset files are deferred pending criteria finalization.
    # Master evidence and deduplicated pairs tables are the primary deliverables.
    df_prot_only = df_dedup[df_dedup["is_protein_protein_pair"] == "TRUE"]

    # ---------------------------------------------------------
    # 7. Audit & Provenance Verification
    # ---------------------------------------------------------
    unique_proteins = sorted(list(set(df_dedup["participant_a_id"]).union(set(df_dedup["participant_b_id"]))))
    unique_pure_proteins = sorted(list(set(df_prot_only["participant_a_id"]).union(set(df_prot_only["participant_b_id"]))))

    hetero_pairs = df_dedup[df_dedup["pair_type"] == "heteromer"]
    homo_pairs = df_dedup[df_dedup["pair_type"] == "homomer"]
    evidence_dist = collections.Counter(df_dedup["evidence_count"])

    audit_summary = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "species": "Oryza sativa (Rice)",
        "provenance_standard": "Database-Only Provenance (Internal Guidelines / Holdout Flags Stripped)",
        "filtering_rules": {
            "IntAct": "Retain interaction_semantics == 'physical_association'",
            "BioGRID": "Exclude proximity and genetic interactions",
            "STRING": "Retain direct experiments > 0 (excluding transferred interologs where experiments == 0)"
        },
        "evidence_table_columns_count": len(df_evidence.columns),
        "dedup_table_columns_count": len(df_dedup.columns),
        "filtering_reconciliation": {
            "intact_raw_records": len(raw_intact_records),
            "intact_positive_retained": len(intact_positive),
            "biogrid_raw_records": len(raw_bg_records),
            "biogrid_positive_retained": len(bg_positive),
            "biogrid_excluded": len(bg_excluded),
            "string_raw_physical_records": string_total_raw,
            "string_positive_retained": len(string_positive_evidence),
            "string_excluded_transferred_interologs": string_excluded_interologs,
            "string_excluded_neither_experimental": string_excluded_neither,
            "total_positive_evidence_records": len(clean_evidence_records)
        },
        "deduplicated_pairs_summary": {
            "total_unique_pairs": len(df_dedup),
            "heteromeric_pairs": len(hetero_pairs),
            "homomeric_pairs": len(homo_pairs),
            "pure_protein_protein_pairs": len(df_prot_only),
            "nucleic_acid_associated_pairs": len(df_dedup) - len(df_prot_only),
            "unique_proteins_all": len(unique_proteins),
            "unique_proteins_pure_protein_pairs": len(unique_pure_proteins),
            "evidence_count_distribution": {str(k): v for k, v in sorted(evidence_dist.items())}
        },
        "literature_provenance": {
            "total_unique_pmids": len(set(
                p for p_str in df_dedup["pmids"] if p_str != "NA"
                for p in p_str.split("|")
            )),
            "intact_pmids_count": len(set(r["pmid"] for r in intact_positive if r["pmid"] not in ("NA", "", None))),
            "biogrid_pmids_count": len(set(r["pmid"] for r in bg_positive if r["pmid"] not in ("NA", "", None))),
            "string_pmids_count": 1,
            "unique_authors_count": len(set(df_evidence["first_author"]))
        },
        "output_files": {
            "evidence_tsv": evidence_tsv_path,
            "evidence_parquet": evidence_parquet_path,
            "dedup_pairs_tsv": dedup_tsv_path,
            "dedup_pairs_parquet": dedup_parquet_path
        }
    }

    # ---------------------------------------------------------
    # 8. Compute HIPPIE Quality Scores & Intermediate Provenance
    # ---------------------------------------------------------
    from scripts.data_preparation.compute_hippie_scores import run_hippie_pipeline
    hippie_stats = run_hippie_pipeline(output_dir)
    audit_summary["hippie_scoring_summary"] = hippie_stats
    audit_summary["output_files"]["pair_hippie_provenance_tsv"] = os.path.join(output_dir, "rice_positive_hippie_provenance.tsv")
    audit_summary["output_files"]["pair_hippie_provenance_parquet"] = os.path.join(output_dir, "rice_positive_hippie_provenance.parquet")
    audit_summary["output_files"]["evidence_hippie_provenance_tsv"] = os.path.join(output_dir, "rice_positive_evidence_hippie_provenance.tsv")
    audit_summary["output_files"]["evidence_hippie_provenance_parquet"] = os.path.join(output_dir, "rice_positive_evidence_hippie_provenance.parquet")

    audit_path = os.path.join(output_dir, "rice_positive_audit.json")
    with open(audit_path, "w", encoding="utf-8") as f:
        json.dump(audit_summary, f, indent=2)

    if interim_audit:
        os.makedirs(os.path.dirname(os.path.abspath(interim_audit)), exist_ok=True)
        with open(interim_audit, "w", encoding="utf-8") as f:
            json.dump(audit_summary, f, indent=2)

    print(f"CLEAN DATABASE PROVENANCE COLLATE SUCCESSFUL:")
    print(f"  - Evidence Table Columns:       {len(df_evidence.columns) + 1} (including record_hippie_score)")
    print(f"  - Deduplicated Table Columns:   {len(df_dedup.columns) + 2} (including hippie_score & level)")
    print(f"  - Collated Positive Records:    {len(clean_evidence_records):,}")
    print(f"  - Deduplicated Unique Pairs:    {len(df_dedup):,}")
    print(f"    * Pure Protein-Protein Pairs: {len(df_prot_only):,}")
    print(f"    * Heteromeric Pairs:          {len(hetero_pairs):,}")
    print(f"    * Homomeric Pairs:            {len(homo_pairs):,}")
    print(f"    * High Confidence HIPPIE (>=0.72): {hippie_stats['high_confidence_count']:,} ({hippie_stats['high_confidence_count']/len(df_dedup)*100:.1f}%)")

    return audit_summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Collate positive Rice PPI dataset with clean database-only provenance")
    parser.add_argument("--output-dir", default=os.path.join(WORKSPACE_ROOT, "data", "processed", "rice"),
                        help="Directory to write processed positive dataset files")
    parser.add_argument("--interim-audit", default=os.path.join(WORKSPACE_ROOT, "data", "interim", "phase2", "audits", "rice_positive_collate_audit.json"),
                        help="Optional interim audit JSON path")
    args = parser.parse_args()

    build_rice_positives(output_dir=args.output_dir, interim_audit=args.interim_audit)
