"""
Collate Positive Sample Datasets for External Holdout Species
============================================================
Parses, filters, and standardizes physical protein-protein interaction evidence for:
1. Maize (Zea mays, taxid 4577)
2. Tomato (Solanum lycopersicum, taxid 4081)
3. Soybean (Glycine max, taxid 3847)

Sources:
- IntAct (Release 252) - retains physical association / direct binary interactions
- BioGRID (Release 5.0.261) - retains physical interactions (excludes genetic and proximity)

Exports standardized datasets matching Arabidopsis and Rice standards:
- <species>_positive_evidence.tsv / .parquet (with record_hippie_score)
- <species>_positive_pairs_deduplicated.tsv / .parquet (with hippie_score & level)
- <species>_positive_hippie_provenance.tsv / .parquet (pair-level HIPPIE audit)
- <species>_positive_evidence_hippie_provenance.tsv / .parquet (evidence-level HIPPIE audit)
- <species>_positive_audit.json (provenance audit)
"""

import os
import sys
import json
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
from scripts.data_preparation.compute_hippie_scores import run_hippie_pipeline

EXTERNAL_SPECIES_CONFIG = {
    "maize": {
        "species_code": "MAIZE",
        "scientific_name": "Zea mays",
        "common_name": "Maize",
        "taxid": 4577,
        "prefix": "maize",
        "pair_id_prefix": "MAIZE_POS"
    },
    "tomato": {
        "species_code": "TOMATO",
        "scientific_name": "Solanum lycopersicum",
        "common_name": "Tomato",
        "taxid": 4081,
        "prefix": "tomato",
        "pair_id_prefix": "TOMATO_POS"
    },
    "soybean": {
        "species_code": "SOYBEAN",
        "scientific_name": "Glycine max",
        "common_name": "Soybean",
        "taxid": 3847,
        "prefix": "soybean",
        "pair_id_prefix": "SOYBEAN_POS"
    }
}


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


def build_species_positives(species_key, output_base=None):
    cfg_meta = EXTERNAL_SPECIES_CONFIG[species_key]
    sp_code = cfg_meta["species_code"]
    prefix = cfg_meta["prefix"]
    taxid = cfg_meta["taxid"]
    common_name = cfg_meta["common_name"]
    sci_name = cfg_meta["scientific_name"]
    pair_prefix = cfg_meta["pair_id_prefix"]

    if output_base is None:
        output_dir = os.path.join(WORKSPACE_ROOT, "data", "processed", prefix)
    else:
        output_dir = os.path.join(output_base, prefix)
    os.makedirs(output_dir, exist_ok=True)

    harmonizer = AssayHarmonizer()

    print("=" * 70)
    print(f"COLLATING POSITIVE PPI DATASET: {common_name.upper()} ({sci_name}, TaxID: {taxid})")
    print(f"Target Output Directory: {output_dir}")
    print("=" * 70)

    # ---------------------------------------------------------
    # 1. Ingest Raw IntAct Sources
    # ---------------------------------------------------------
    intact_cfg = [s for s in INTACT_SOURCES if s["species_code"] == sp_code][0]
    raw_intact_records, raw_intact_rejections = parse_intact_source(intact_cfg, harmonizer)
    print(f"Parsed IntAct: {len(raw_intact_records):,} records ({len(raw_intact_rejections):,} rejected)")

    # Filter IntAct: physical interactions only (physical_association, direct_binary)
    # Exclude proximity / colocalization
    intact_positive = []
    intact_excluded = []
    for r in raw_intact_records:
        itype = str(r["interaction_semantics"])
        raw_t = str(r["interaction_type_raw"])
        if itype in ("physical_association", "direct_binary") or "MI:0915" in raw_t or "MI:0407" in raw_t or "MI:0217" in raw_t:
            if not ("MI:0403" in raw_t or "MI:2364" in raw_t or itype == "proximity"):
                intact_positive.append(r)
            else:
                intact_excluded.append(r)
        else:
            intact_excluded.append(r)
    print(f"  -> Retained {len(intact_positive):,} physical IntAct records ({len(intact_excluded):,} excluded).")

    # ---------------------------------------------------------
    # 2. Ingest Raw BioGRID Sources
    # ---------------------------------------------------------
    bg_cfg = [s for s in BIOGRID_SOURCES if s["species_code"] == sp_code][0]
    raw_bg_records, raw_bg_rejections = parse_biogrid_source(bg_cfg, harmonizer)
    print(f"Parsed BioGRID: {len(raw_bg_records):,} records ({len(raw_bg_rejections):,} rejected)")

    # Filter BioGRID: physical interactions only (exclude genetic and proximity/PCA)
    bg_positive = []
    bg_excluded = []
    for r in raw_bg_records:
        itype = str(r["interaction_semantics"])
        sys_raw = str(r.get("detection_method_raw", "")).lower()
        if itype in ("physical_association", "direct_binary") and sys_raw != "pca" and itype != "genetic":
            bg_positive.append(r)
        else:
            bg_excluded.append(r)
    print(f"  -> Retained {len(bg_positive):,} physical BioGRID records ({len(bg_excluded):,} excluded).")

    # ---------------------------------------------------------
    # 3. Clean Evidence Records: Standardize Provenance Fields
    # ---------------------------------------------------------
    clean_evidence_records = []
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

    df_evidence = pd.DataFrame(clean_evidence_records)

    evidence_tsv_path = os.path.join(output_dir, f"{prefix}_positive_evidence.tsv")
    evidence_parquet_path = os.path.join(output_dir, f"{prefix}_positive_evidence.parquet")
    df_evidence.to_csv(evidence_tsv_path, sep="\t", index=False)

    arrow_schema = pa.schema([(c, pa.string()) for c in df_evidence.columns])
    table_evidence = pa.Table.from_pandas(df_evidence.astype(str), schema=arrow_schema)
    pq.write_table(table_evidence, evidence_parquet_path, compression="snappy")
    print(f"Exported preliminary evidence table: {evidence_tsv_path} ({len(df_evidence):,} rows)")

    # ---------------------------------------------------------
    # 4. Deduplicate into Canonical Positive Pairs
    # ---------------------------------------------------------
    pair_evidence_map = collections.defaultdict(list)
    for r in clean_evidence_records:
        pair_tuple = tuple(r["canonical_pair_key"].split("|"))
        pair_evidence_map[pair_tuple].append(r)

    sorted_pair_keys = sorted(pair_evidence_map.keys(), key=lambda p: (p[0], p[1]))
    id_padding = 6 if len(sorted_pair_keys) >= 10000 else 5

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

        names_a = set(r["participant_a_name"] for r in ev_list if r["participant_a_id"] == p_a)
        names_a.update(r["participant_b_name"] for r in ev_list if r["participant_b_id"] == p_a)
        names_a.discard("NA")

        names_b = set(r["participant_a_name"] for r in ev_list if r["participant_a_id"] == p_b)
        names_b.update(r["participant_b_name"] for r in ev_list if r["participant_b_id"] == p_b)
        names_b.discard("NA")

        types_a = set(r["participant_a_molecule_type"] for r in ev_list if r["participant_a_id"] == p_a)
        types_a.update(r["participant_b_molecule_type"] for r in ev_list if r["participant_b_id"] == p_a)

        types_b = set(r["participant_a_molecule_type"] for r in ev_list if r["participant_a_id"] == p_b)
        types_b.update(r["participant_b_molecule_type"] for r in ev_list if r["participant_b_id"] == p_b)

        type_a_str = "|".join(sorted(types_a)) if types_a else "protein"
        type_b_str = "|".join(sorted(types_b)) if types_b else "protein"

        miscores = []
        for r in ev_list:
            if r["intact_miscore"] not in ("NA", "None", None):
                try:
                    miscores.append(float(r["intact_miscore"]))
                except (ValueError, TypeError):
                    pass
        max_mi = f"{max(miscores):.4f}" if miscores else "NA"

        is_prot_prot = "TRUE" if (type_a_str == "protein" and type_b_str == "protein") else "FALSE"

        dedup_rows.append({
            "pair_id": f"{pair_prefix}_{idx:0{id_padding}d}",
            "canonical_pair_key": f"{p_a}|{p_b}",
            "participant_a_id": p_a,
            "participant_b_id": p_b,
            "participant_a_name": "|".join(sorted(names_a)) if names_a else p_a,
            "participant_b_name": "|".join(sorted(names_b)) if names_b else p_b,
            "participant_a_type": type_a_str,
            "participant_b_type": type_b_str,
            "is_protein_protein_pair": is_prot_prot,
            "is_homomeric": "TRUE" if is_homo else "FALSE",
            "pair_type": "homomer" if is_homo else "heteromer",
            "evidence_count": len(ev_list),
            "sources": "|".join(sources),
            "interaction_types": "|".join(int_types) if int_types else "physical association",
            "detection_methods": "|".join(det_methods) if det_methods else "unspecified",
            "first_authors": "|".join(authors) if authors else "NA",
            "pmid_count": len(pmids),
            "pmids": "|".join(pmids) if pmids else "NA",
            "publication_years": "|".join(years) if years else "NA",
            "throughput_classes": "|".join(tp) if tp else "unspecified",
            "max_intact_miscore": max_mi,
            "source_record_ids": "|".join(src_rec_ids)
        })

    df_dedup = pd.DataFrame(dedup_rows)

    dedup_tsv_path = os.path.join(output_dir, f"{prefix}_positive_pairs_deduplicated.tsv")
    dedup_parquet_path = os.path.join(output_dir, f"{prefix}_positive_pairs_deduplicated.parquet")
    df_dedup.to_csv(dedup_tsv_path, sep="\t", index=False)

    arrow_schema_dedup = pa.schema([(c, pa.string()) for c in df_dedup.columns])
    table_dedup = pa.Table.from_pandas(df_dedup.astype(str), schema=arrow_schema_dedup)
    pq.write_table(table_dedup, dedup_parquet_path, compression="snappy")
    print(f"Exported preliminary dedup pairs: {dedup_tsv_path} ({len(df_dedup):,} pairs)")

    # ---------------------------------------------------------
    # 5. Compute HIPPIE Quality Scores & Intermediate Provenance
    # ---------------------------------------------------------
    print("Computing HIPPIE scores for", common_name, "...")
    hippie_stats = run_hippie_pipeline(data_dir=output_dir, species=prefix)

    # ---------------------------------------------------------
    # 6. Build and Export Comprehensive Audit JSON
    # ---------------------------------------------------------
    df_prot_only = df_dedup[df_dedup["is_protein_protein_pair"] == "TRUE"]
    unique_proteins = sorted(list(set(df_dedup["participant_a_id"]).union(set(df_dedup["participant_b_id"]))))
    unique_pure_proteins = sorted(list(set(df_prot_only["participant_a_id"]).union(set(df_prot_only["participant_b_id"]))))

    hetero_pairs = df_dedup[df_dedup["pair_type"] == "heteromer"]
    homo_pairs = df_dedup[df_dedup["pair_type"] == "homomer"]
    evidence_dist = collections.Counter(df_dedup["evidence_count"])

    audit_summary = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "species_code": sp_code,
        "species_common": common_name,
        "species_scientific": sci_name,
        "taxonomy_id": taxid,
        "provenance_standard": "Database-Only Provenance (Internal Guidelines / Holdout Flags Stripped)",
        "filtering_rules": {
            "IntAct": "Retain interaction_semantics in ('physical_association', 'direct_binary'); exclude colocalization",
            "BioGRID": "Retain physical systems; exclude proximity (PCA) and genetic interactions"
        },
        "filtering_reconciliation": {
            "intact_raw_records": len(raw_intact_records),
            "intact_positive_retained": len(intact_positive),
            "intact_excluded": len(intact_excluded),
            "biogrid_raw_records": len(raw_bg_records),
            "biogrid_positive_retained": len(bg_positive),
            "biogrid_excluded": len(bg_excluded),
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
            "unique_authors_count": len(set(df_evidence["first_author"]))
        },
        "hippie_scoring_summary": hippie_stats,
        "output_files": {
            "evidence_tsv": evidence_tsv_path,
            "evidence_parquet": evidence_parquet_path,
            "dedup_pairs_tsv": dedup_tsv_path,
            "dedup_pairs_parquet": dedup_parquet_path,
            "pair_hippie_provenance_tsv": os.path.join(output_dir, f"{prefix}_positive_hippie_provenance.tsv"),
            "pair_hippie_provenance_parquet": os.path.join(output_dir, f"{prefix}_positive_hippie_provenance.parquet"),
            "evidence_hippie_provenance_tsv": os.path.join(output_dir, f"{prefix}_positive_evidence_hippie_provenance.tsv"),
            "evidence_hippie_provenance_parquet": os.path.join(output_dir, f"{prefix}_positive_evidence_hippie_provenance.parquet")
        }
    }

    audit_path = os.path.join(output_dir, f"{prefix}_positive_audit.json")
    with open(audit_path, "w", encoding="utf-8") as f:
        json.dump(audit_summary, f, indent=2)

    print(f"\n[DONE] Finished {common_name}:")
    print(f"  Evidence records: {len(clean_evidence_records):,}")
    print(f"  Unique pairs:     {len(df_dedup):,}")
    print(f"  Audit file:       {audit_path}")

    return audit_summary


def main():
    parser = argparse.ArgumentParser(description="Collate positive datasets with HIPPIE scores for external species")
    parser.add_argument("--species", default="all", choices=["all", "maize", "tomato", "soybean"],
                        help="Target species (default: all)")
    parser.add_argument("--output-base", default=None, help="Base directory for output (default: data/processed)")
    args = parser.parse_args()

    targets = ["maize", "tomato", "soybean"] if args.species == "all" else [args.species]

    summaries = {}
    for sp in targets:
        summaries[sp] = build_species_positives(sp, output_base=args.output_base)

    print("\n" + "=" * 70)
    print("ALL TARGET EXTERNAL SPECIES PROCESSED SUCCESSFULLY!")
    print("=" * 70)
    for sp, s in summaries.items():
        pairs = s["deduplicated_pairs_summary"]["total_unique_pairs"]
        evs = s["filtering_reconciliation"]["total_positive_evidence_records"]
        hippie = s["hippie_scoring_summary"]
        print(f"{sp.upper():8s}: {evs:,} evidence rows -> {pairs:,} pairs | Mean HIPPIE: {hippie['mean_score']:.4f} ({hippie['high_confidence_count']} high conf)")


if __name__ == "__main__":
    main()
