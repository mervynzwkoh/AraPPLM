"""
Analyze External Holdout Species Breakdown from Raw Datasets and Parquet Warehouse
=================================================================================
Analyzes interactome records for Zea mays (Maize, taxid 4577),
Solanum lycopersicum (Tomato, taxid 4081), and Glycine max (Soybean, taxid 3847).

Computes:
1. Total raw and normalized records available per species and per data source.
2. Filtered physical PPIs (experimental physical, STRING physical links, STRING experimental > 0).
3. Breakdown by interaction semantics and assay families.
4. Deduplication metrics (unique canonical pairs, homomers vs. heteromers, unique proteins).
5. Exports comprehensive audit JSON to data/interim/phase2/audits/external_source_breakdown.json.
"""

import os
import sys
import glob
import json
import gzip
import collections
import pandas as pd
import pyarrow.parquet as pq

WORKSPACE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DEFAULT_DATA_ROOT = os.path.join(os.path.expanduser("~"), "AraPPLM_Data", "phase2")
DATA_ROOT = os.environ.get("ARAPPLM_DATA_ROOT", DEFAULT_DATA_ROOT)
AUDIT_DIR = os.path.join(WORKSPACE_ROOT, "data", "interim", "phase2", "audits")

SPECIES_METADATA = {
    "MAIZE": {
        "scientific_name": "Zea mays",
        "common_name": "Maize",
        "taxonomy_id": 4577,
        "raw_dir": os.path.join(WORKSPACE_ROOT, "data", "raw", "external_holdout", "maize"),
        "intact_file": "intact_maize_taxid4577_release252.mitab27.txt",
        "biogrid_file": "BIOGRID-ORGANISM-Zea_mays-5.0.261.tab3.txt",
        "string_file": "4577.protein.physical.links.detailed.v12.0.txt.gz"
    },
    "TOMATO": {
        "scientific_name": "Solanum lycopersicum",
        "common_name": "Tomato",
        "taxonomy_id": 4081,
        "raw_dir": os.path.join(WORKSPACE_ROOT, "data", "raw", "external_holdout", "tomato"),
        "intact_file": "intact_tomato_taxid4081_release252.mitab27.txt",
        "biogrid_file": "BIOGRID-ORGANISM-Solanum_lycopersicum-5.0.261.tab3.txt",
        "string_file": "4081.protein.physical.links.detailed.v12.0.txt.gz"
    },
    "SOYBEAN": {
        "scientific_name": "Glycine max",
        "common_name": "Soybean",
        "taxonomy_id": 3847,
        "raw_dir": os.path.join(WORKSPACE_ROOT, "data", "raw", "external_holdout", "soybean"),
        "intact_file": "intact_soybean_taxid3847_release252.mitab27.txt",
        "biogrid_file": "BIOGRID-ORGANISM-Glycine_max-5.0.261.tab3.txt",
        "string_file": "3847.protein.physical.links.detailed.v12.0.txt.gz"
    }
}


def analyze_species(sp_code, sp_meta):
    print(f"\nAnalyzing {sp_meta['common_name']} ({sp_meta['scientific_name']}, TaxID: {sp_meta['taxonomy_id']})...")
    res_dict = {
        "species_code": sp_code,
        "scientific_name": sp_meta["scientific_name"],
        "common_name": sp_meta["common_name"],
        "taxonomy_id": sp_meta["taxonomy_id"],
        "sources": {}
    }

    # 1. IntAct MITAB 2.7
    intact_path = os.path.join(sp_meta["raw_dir"], sp_meta["intact_file"])
    with open(intact_path, "r", encoding="utf-8", errors="ignore") as f:
        intact_lines = [l.strip().split("\t") for l in f if l.strip()]

    intact_total = len(intact_lines)
    intact_phys_records = 0
    intact_excluded_records = 0
    intact_by_type = collections.Counter()
    intact_by_method = collections.Counter()
    intact_phys_pairs = set()
    intact_phys_homo = set()
    intact_phys_hetero = set()
    intact_phys_proteins = set()

    for r in intact_lines:
        itype = r[11] if len(r) > 11 else "unknown"
        imethod = r[6] if len(r) > 6 else "unknown"
        intact_by_type[itype] += 1
        intact_by_method[imethod] += 1

        p1 = r[0].split(":")[-1]
        p2 = r[1].split(":")[-1]
        pair = tuple(sorted([p1, p2]))

        # Physical filter
        if any(term in itype for term in ["MI:0915", "MI:0407", "MI:0914", "MI:0217"]) and not any(term in itype for term in ["MI:0403", "MI:2364"]):
            intact_phys_records += 1
            intact_phys_pairs.add(pair)
            intact_phys_proteins.add(p1)
            intact_phys_proteins.add(p2)
            if p1 == p2:
                intact_phys_homo.add(pair)
            else:
                intact_phys_hetero.add(pair)
        else:
            intact_excluded_records += 1

    res_dict["sources"]["INTACT"] = {
        "source_resource": "IntAct",
        "file_name": sp_meta["intact_file"],
        "total_records": intact_total,
        "physical_records": intact_phys_records,
        "excluded_records": intact_excluded_records,
        "unique_physical_pairs_total": len(intact_phys_pairs),
        "unique_physical_pairs_hetero": len(intact_phys_hetero),
        "unique_physical_pairs_homo": len(intact_phys_homo),
        "unique_physical_proteins": len(intact_phys_proteins),
        "interaction_types": dict(intact_by_type),
        "detection_methods": dict(intact_by_method.most_common(10))
    }

    # 2. BioGRID TAB 3.0
    biogrid_path = os.path.join(sp_meta["raw_dir"], sp_meta["biogrid_file"])
    with open(biogrid_path, "r", encoding="utf-8", errors="ignore") as f:
        bg_lines = [l.strip().split("\t") for l in f if l.strip() and not l.startswith("#")]

    bg_total = len(bg_lines)
    bg_phys_records = 0
    bg_excluded_records = 0
    bg_by_sys_type = collections.Counter()
    bg_by_sys_name = collections.Counter()
    bg_phys_pairs = set()
    bg_phys_homo = set()
    bg_phys_hetero = set()
    bg_phys_proteins = set()

    for r in bg_lines:
        sys_type = r[12] if len(r) > 12 else "unknown"
        sys_name = r[11] if len(r) > 11 else "unknown"
        bg_by_sys_type[sys_type] += 1
        bg_by_sys_name[sys_name] += 1

        p1 = r[1]
        p2 = r[2]
        pair = tuple(sorted([p1, p2]))

        # Physical filter
        if sys_type == "physical" and sys_name != "PCA":
            bg_phys_records += 1
            bg_phys_pairs.add(pair)
            bg_phys_proteins.add(p1)
            bg_phys_proteins.add(p2)
            if p1 == p2:
                bg_phys_homo.add(pair)
            else:
                bg_phys_hetero.add(pair)
        else:
            bg_excluded_records += 1

    res_dict["sources"]["BIOGRID"] = {
        "source_resource": "BioGRID",
        "file_name": sp_meta["biogrid_file"],
        "total_records": bg_total,
        "physical_records": bg_phys_records,
        "excluded_records": bg_excluded_records,
        "unique_physical_pairs_total": len(bg_phys_pairs),
        "unique_physical_pairs_hetero": len(bg_phys_hetero),
        "unique_physical_pairs_homo": len(bg_phys_homo),
        "unique_physical_proteins": len(bg_phys_proteins),
        "system_types": dict(bg_by_sys_type),
        "experimental_systems": dict(bg_by_sys_name.most_common(10))
    }

    # 3. STRING Physical Links Detailed v12.0
    string_path = os.path.join(sp_meta["raw_dir"], sp_meta["string_file"])
    str_total = 0
    str_exp_gt_0 = 0
    str_exp_0 = 0
    with gzip.open(string_path, "rt", encoding="utf-8") as f:
        _ = f.readline()
        for line in f:
            str_total += 1
            parts = line.split()
            exp = int(parts[2])
            if exp > 0:
                str_exp_gt_0 += 1
            else:
                str_exp_0 += 1

    res_dict["sources"]["STRING"] = {
        "source_resource": "STRING Physical Links v12.0",
        "file_name": sp_meta["string_file"],
        "total_records": str_total,
        "unique_canonical_pairs": str_total // 2,
        "experimental_gt_0_records": str_exp_gt_0,
        "experimental_gt_0_canonical_pairs": str_exp_gt_0 // 2,
        "transferred_or_computational_only_records": str_exp_0,
        "transferred_or_computational_canonical_pairs": str_exp_0 // 2
    }

    # 4. Summary Aggregations
    curated_exp_total = intact_total + bg_total
    curated_exp_physical = intact_phys_records + bg_phys_records
    curated_exp_excluded = intact_excluded_records + bg_excluded_records

    res_dict["summary"] = {
        "total_records_all_sources": intact_total + bg_total + str_total,
        "curated_experimental_total": curated_exp_total,
        "curated_experimental_physical": curated_exp_physical,
        "curated_experimental_excluded": curated_exp_excluded,
        "string_physical_links_all": str_total,
        "string_physical_links_experimental_gt_0": str_exp_gt_0,
        "string_physical_links_transferred_only": str_exp_0,
        # Total Physical PPI Records under different thresholds/definitions:
        "total_physical_ppis_curated_experimental_only": curated_exp_physical,
        "total_physical_ppis_curated_plus_string_experimental": curated_exp_physical + str_exp_gt_0,
        "total_physical_ppis_curated_plus_all_string_physical": curated_exp_physical + str_total
    }

    return res_dict


def main():
    os.makedirs(AUDIT_DIR, exist_ok=True)
    all_species_results = {}
    grand_totals = {
        "total_records_all_sources": 0,
        "curated_experimental_total": 0,
        "curated_experimental_physical": 0,
        "curated_experimental_excluded": 0,
        "string_physical_links_all": 0,
        "string_physical_links_experimental_gt_0": 0,
        "string_physical_links_transferred_only": 0,
        "total_physical_ppis_curated_experimental_only": 0,
        "total_physical_ppis_curated_plus_string_experimental": 0,
        "total_physical_ppis_curated_plus_all_string_physical": 0
    }

    for sp_code, sp_meta in SPECIES_METADATA.items():
        res = analyze_species(sp_code, sp_meta)
        all_species_results[sp_code] = res
        for k in grand_totals:
            grand_totals[k] += res["summary"][k]

    audit_output = {
        "generated_timestamp_utc": pd.Timestamp.now("UTC").isoformat(),
        "species_breakdown": all_species_results,
        "grand_totals": grand_totals
    }

    output_path = os.path.join(AUDIT_DIR, "external_source_breakdown.json")
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(audit_output, f, indent=2)

    print(f"\nWrote full audit to: {output_path}")


if __name__ == "__main__":
    main()
