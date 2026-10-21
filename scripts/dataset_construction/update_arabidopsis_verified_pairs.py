"""
Update Arabidopsis Verified Positive Pairs (PCA Retained)
=========================================================
Re-processes the Arabidopsis interactome to:
1. Retain samples detected by PCA (Protein Complementation Assays) across both BioGRID and IntAct
   (including split-ubiquitin MI:0112, BiFC MI:0809, split-luciferase MI:1204/MI:1037, MI:0090).
2. Use the Core Physical set for IntAct (including direct binary MI:0407, physical association MI:0915/MI:0914,
   and physical/phosphorylation MI:0217; strictly excluding proximity MI:2364 and colocalization MI:0403).
3. Enforce strict protein-only complexes across all sources (excluding DNA, RNA, small molecules).
4. Compute continuous HIPPIE quality scores and tiers for all deduplicated pairs.
5. Output updated parquet, TSV, and audit JSON files in data/dataset_construction/01_positive_pairs/arabidopsis/.
"""

import os
import re
import json
import collections
import pandas as pd
import numpy as np

WORKSPACE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

INTACT_FILE = os.path.join(WORKSPACE_ROOT, "data", "raw", "arabidopsis", "intact", "intact_arabidopsis_taxid3702_release252.mitab27.txt")
BIOGRID_FILE = os.path.join(WORKSPACE_ROOT, "data", "raw", "arabidopsis", "biogrid", "BIOGRID-ORGANISM-Arabidopsis_thaliana_Columbia-5.0.261.tab3.txt")
STRING_FILE = os.path.join(WORKSPACE_ROOT, "data", "raw", "arabidopsis", "string", "3702.protein.physical.links.full.v12.0.txt.gz")
XLMS_FILE = os.path.join(WORKSPACE_ROOT, "data", "raw", "arabidopsis", "xlms_2026", "Total_XL_plink3-2_v3.csv")

OUT_DIR = os.path.join(WORKSPACE_ROOT, "data", "dataset_construction", "01_positive_pairs", "arabidopsis")

PCA_MI = {
    "MI:0112", # ubiquitin reconstruction (split-ubiquitin / SUS / MYTH)
    "MI:0809", # bimolecular fluorescence complementation (BiFC)
    "MI:1204", # split firefly luciferase complementation
    "MI:1037", # split renilla luciferase complementation
    "MI:0090", # protein complementation assay
    "MI:2365", # fluorescent protein-protein interaction-visualization
    "MI:0097", # reverse ras recruitment system
}

PROX_MI = {
    "MI:0663", # confocal microscopy
    "MI:0416", # fluorescence microscopy
    "MI:1314", # BioID / proximity-dependent biotin identification
    "MI:0813", # proximity ligation assay
    "MI:2364", # proximity interaction type
    "MI:0403", # colocalization
}

CORE_TERMS = {"MI:0915", "MI:0914", "MI:0407", "MI:0217"}

TECHNIQUE_MAPPING = {
    "two-hybrid": ("two-hybrid", 5.0),
    "two hybrid": ("two-hybrid", 5.0),
    "two hybrid array": ("two hybrid array", 5.0),
    "two hybrid prey pooling approach": ("two hybrid fragment pooling approach", 5.0),
    "two hybrid fragment pooling approach": ("two hybrid fragment pooling approach", 5.0),
    "validated two hybrid": ("validated two hybrid", 5.0),
    "cr-two hybrid": ("cr-two hybrid", 5.0),
    "reverse two hybrid": ("reverse two hybrid", 5.0),
    "three hybrid": ("three hybrid", 5.0),
    "protein three hybrid": ("protein three hybrid", 5.0),
    "co-fractionation": ("co-fractionation", 1.0),
    "co-purification": ("co-purification", 2.0),
    "copurification": ("co-purification", 2.0),
    "blue native page": ("blue native page", 3.0),
    "comigration in gel electrophoresis": ("comigration in gel electrophoresis", 3.0),
    "comigration in non denaturing gel electrophoresis": ("comigration in non denaturing gel electrophoresis", 3.0),
    "molecular sieving": ("molecular sieving", 2.0),
    "cosedimentation in solution": ("cosedimentation in solution", 2.0),
    "cosedimentation through density gradient": ("cosedimentation through density gradient", 2.0),
    "co-crystal structure": ("co-crystal structure", 10.0),
    "x-ray crystallography": ("x-ray crystallography", 10.0),
    "reconstituted complex": ("reconstituted complex", 10.0),
    "surface plasmon resonance": ("surface plasmon resonance", 10.0),
    "isothermal titration calorimetry": ("isothermal titration calorimetry", 10.0),
    "nuclear magnetic resonance": ("nuclear magnetic resonance", 10.0),
    "light scattering": ("light scattering", 10.0),
    "x ray scattering": ("x ray scattering", 9.0),
    "biochemical activity": ("biochemical activity", 5.0),
    "biochemical": ("biochemical", 1.0),
    "biophysical": ("biophysical", 1.0),
    "protein kinase assay": ("protein kinase assay", 7.5),
    "in-gel kinase assay": ("in-gel kinase assay", 7.5),
    "enzymatic study": ("enzymatic study", 1.0),
    "enzyme linked immunosorbent assay": ("enzyme linked immunosorbent assay", 5.0),
    "chromatin immunoprecipitation assay": ("chromatin immunoprecipitation assay", 2.0),
    "electrophoretic mobility shift assay": ("electrophoretic mobility shift assay", 2.0),
    "affinity capture-ms": ("affinity capture-ms", 5.0),
    "tandem affinity purification": ("tandem affinity purification", 5.0),
    "affinity capture-western": ("affinity capture-western", 5.0),
    "affinity capture-luminescence": ("affinity capture-luminescence", 5.0),
    "affinity capture-rna": ("affinity capture-rna", 2.0),
    "affinity chromatography technology": ("affinity chromatography technology", 5.0),
    "affinity technology": ("affinity technology", 5.0),
    "anti bait coimmunoprecipitation": ("anti bait coimmunoprecipitation", 5.0),
    "anti tag coimmunoprecipitation": ("anti tag coimmunoprecipitation", 5.0),
    "coimmunoprecipitation": ("coimmunoprecipitation", 5.0),
    "pull down": ("pull down", 2.5),
    "cross-linking study": ("cross-linking study", 5.0),
    "cross-linking-ms (xl-ms)": ("cross-linking study", 5.0),
    "protein cross-linking with a bifunctional reagent": ("protein cross-linking with a bifunctional reagent", 5.0),
    "fret": ("fret", 6.0),
    "fluorescent resonance energy transfer": ("fret", 6.0),
    "classical fluorescence spectroscopy": ("classical fluorescence spectroscopy", 7.5),
    "far western": ("far western blotting", 5.0),
    "far western blotting": ("far western blotting", 5.0),
    "solid phase assay": ("solid phase assay", 1.0),
    "protein array": ("protein array", 5.0),
    "filter binding": ("filter binding", 5.0),
    "phage display": ("phage display", 6.0),
    "lambda phage display": ("lambda phage display", 6.0),
    "electron microscopy": ("electron microscopy", 5.0),
    "competition binding": ("competition binding", 5.0),
    "protein-peptide": ("protein-peptide", 5.0),
    "protein-rna": ("protein-rna", 0.0),
    "experimental interaction detection": ("experimental interaction detection", 1.0),
    # PCA methods mapped to Table 2 (weight 5.0)
    "pca": ("protein complementation assay", 5.0),
    "protein complementation assay": ("protein complementation assay", 5.0),
    "bimolecular fluorescence complementation": ("bimolecular fluorescence complementation", 5.0),
    "bifc": ("bimolecular fluorescence complementation", 5.0),
    "ubiquitin reconstruction": ("ubiquitin reconstruction", 5.0),
    "split firefly luciferase complementation": ("split firefly luciferase complementation", 5.0),
    "split renilla luciferase complementation": ("split renilla luciferase complementation", 5.0),
    "fluorescent protein-protein interaction-visualization": ("fluorescent protein-protein interaction-visualization", 5.0),
    "reverse ras recruitment system": ("reverse ras recruitment system", 5.0),
}

def calc_hippie_score(ns: int, nt: float, no: int = 0) -> float:
    ss = (2.0 / (1.0 + np.exp(-2.3 * ns))) - 1.0 if ns > 0 else 0.0
    st = (2.0 / (1.0 + np.exp(-0.2 * nt))) - 1.0 if nt > 0 else 0.0
    so = (2.0 / (1.0 + np.exp(-1.6 * no))) - 1.0 if no > 0 else 0.0
    return float(0.60 * ss + 0.30 * st + 0.10 * so)

def assign_confidence_level(score: float) -> str:
    if score >= 0.72:
        return "high"
    elif score >= 0.62:
        return "medium"
    else:
        return "low"

def get_mol_type(t_str, raw_id):
    t = t_str.lower()
    i = raw_id.lower()
    if any(k in t for k in ["mi:0328", "small molecule", "chebi"]) or "chebi" in i:
        return "small_molecule"
    if any(k in t for k in ["mi:0319", "mi:0681", "deoxyribonucleic"]):
        return "dna"
    if any(k in t for k in ["mi:0320", "ribonucleic"]):
        return "rna"
    if any(k in t for k in ["mi:0318", "nucleic acid"]):
        return "nucleic_acid"
    if any(k in t for k in ["mi:0250", "gene"]):
        return "gene"
    if any(k in t for k in ["mi:1304", "molecule set"]):
        return "molecule_set"
    if any(k in t for k in ["mi:0326", "protein", "mi:0327", "peptide"]):
        return "protein"
    if raw_id == "-" or not raw_id:
        return "missing"
    return "protein"

def clean_id(raw):
    return raw.split(":")[-1].strip()

def run():
    print("Loading preserved STRING and PRIDE PhoX XL-MS evidence...")
    existing_ev_path = os.path.join(OUT_DIR, "ara_positive_evidence.parquet")
    df_existing = pd.read_parquet(existing_ev_path)
    string_ev = df_existing[df_existing["source_resource"] == "STRING"].copy()
    xlms_ev = df_existing[df_existing["source_resource"] == "PRIDE / pLink 3.2"].copy()
    print(f"Preserved STRING records: {len(string_ev)}")
    print(f"Preserved PRIDE XL-MS records: {len(xlms_ev)}")

    print("\nParsing BioGRID with PCA RETAINED...")
    bg_evidence = []
    bg_stats = collections.Counter()
    with open(BIOGRID_FILE, "r", encoding="utf-8") as f:
        for line_no, line in enumerate(f, 1):
            if line.startswith("#") or not line.strip():
                continue
            bg_stats["total_raw_rows"] += 1
            cols = line.strip().split("\t")
            if len(cols) < 20:
                bg_stats["malformed"] += 1
                continue
            if cols[12].strip().lower() != "physical":
                bg_stats["genetic_excluded"] += 1
                continue
            exp_sys = cols[11].strip()
            if exp_sys in ["Proximity Label-MS", "Co-localization"]:
                bg_stats["proximity_excluded"] += 1
                continue
            is_rna = "RNA" in exp_sys
            if is_rna:
                bg_stats["rna_records"] += 1

            tax_a = cols[15].strip() if len(cols) > 15 else ""
            tax_b = cols[16].strip() if len(cols) > 16 else ""
            if tax_a != "3702" or tax_b != "3702":
                bg_stats["cross_species_excluded"] += 1
                continue

            sp_a = cols[23].strip() if len(cols) > 23 else "-"
            sp_b = cols[26].strip() if len(cols) > 26 else "-"
            sys_a, sys_b = cols[5].strip(), cols[6].strip()
            sym_a, sym_b = cols[7].strip(), cols[8].strip()
            entrez_a, entrez_b = cols[1].strip(), cols[2].strip()

            if sp_a != "-" and sp_a != "":
                ns_a, val_a = "uniprotkb", sp_a.split("|")[0].strip()
            elif sys_a != "-" and sys_a != "":
                ns_a, val_a = "systematic_name", sys_a.split("|")[0].strip()
            elif sym_a != "-" and sym_a != "":
                ns_a, val_a = "gene_symbol", sym_a.split("|")[0].strip()
            else:
                ns_a, val_a = "entrez_gene", entrez_a.split("|")[0].strip()

            if sp_b != "-" and sp_b != "":
                ns_b, val_b = "uniprotkb", sp_b.split("|")[0].strip()
            elif sys_b != "-" and sys_b != "":
                ns_b, val_b = "systematic_name", sys_b.split("|")[0].strip()
            elif sym_b != "-" and sym_b != "":
                ns_b, val_b = "gene_symbol", sym_b.split("|")[0].strip()
            else:
                ns_b, val_b = "entrez_gene", entrez_b.split("|")[0].strip()

            name_a = sym_a if sym_a != "-" else val_a
            name_b = sym_b if sym_b != "-" else val_b
            pair_key = "|".join(sorted([val_a, val_b]))

            pmid = "NA"
            if "PUBMED:" in cols[14].upper():
                pmid = re.sub(r'(?i)PUBMED:', '', cols[14]).strip()

            pub_year = "NA"
            m_yr = re.search(r'\((\d{4})\)', cols[13])
            if m_yr:
                pub_year = m_yr.group(1)

            bg_evidence.append({
                "source_resource": "BioGRID",
                "source_release": "5.0.261",
                "source_file": "data/raw/arabidopsis/biogrid/BIOGRID-ORGANISM-Arabidopsis_thaliana_Columbia-5.0.261.tab3.txt",
                "source_row_number": str(line_no),
                "source_record_id": cols[0].strip(),
                "canonical_pair_key": pair_key,
                "participant_a_id": val_a,
                "participant_a_namespace": ns_a,
                "participant_a_original_id": f"{ns_a}:{val_a}",
                "participant_a_name": name_a,
                "participant_a_aliases": cols[9] if len(cols) > 9 else "-",
                "participant_a_taxid": "3702",
                "participant_a_species": "Arabidopsis thaliana (Columbia)",
                "participant_a_biological_role": "nan",
                "participant_a_experimental_role": "nan",
                "participant_a_molecule_type": "protein",
                "participant_b_id": val_b,
                "participant_b_namespace": ns_b,
                "participant_b_original_id": f"{ns_b}:{val_b}",
                "participant_b_name": name_b,
                "participant_b_aliases": cols[10] if len(cols) > 10 else "-",
                "participant_b_taxid": "3702",
                "participant_b_species": "Arabidopsis thaliana (Columbia)",
                "participant_b_biological_role": "nan",
                "participant_b_experimental_role": "nan",
                "participant_b_molecule_type": "protein",
                "is_protein_protein": "FALSE" if is_rna else "TRUE",
                "same_species_pair": "TRUE",
                "interaction_type_raw": cols[12].strip(),
                "interaction_type_psi_mi": "nan",
                "interaction_type_name": cols[12].strip(),
                "detection_method_raw": exp_sys,
                "detection_method_psi_mi": "nan",
                "detection_method_name": exp_sys,
                "first_author": cols[13].strip(),
                "publication_id_raw": cols[14].strip(),
                "pmid": pmid,
                "publication_year": pub_year,
                "throughput_raw": cols[17].strip() if len(cols) > 17 else "nan",
                "source_confidence_raw": cols[18].strip() if len(cols) > 18 and cols[18] != "-" else "nan",
                "intact_miscore": "nan",
                "xlms_protein_type": "nan",
                "xlms_peptide_a": "nan",
                "xlms_peptide_b": "nan",
                "xlms_linker": "nan",
                "xlms_score": "nan",
                "xlms_evalue": "nan",
                "record_hippie_score": "0.6293"
            })
            bg_stats["retained"] += 1

    print("BioGRID parsing stats:")
    for k, v in bg_stats.items():
        print(f"  {k}: {v}")
    df_bg = pd.DataFrame(bg_evidence)

    print("\nParsing IntAct with Core Physical + PCA RETAINED + Protein-Only...")
    intact_rows = []
    intact_stats = collections.Counter()

    with open(INTACT_FILE, "r", encoding="utf-8") as f:
        for line_no, line in enumerate(f, 1):
            if line.startswith("#") or not line.strip():
                continue
            intact_stats["total_raw_rows"] += 1
            cols = line.strip().split("\t")
            if len(cols) < 15:
                intact_stats["malformed"] += 1
                continue

            it = cols[11]
            is_neg = (len(cols) > 35 and cols[35].strip().lower() == "true")
            if is_neg:
                intact_stats["negative_excluded"] += 1
                continue

            is_core = any(t in it for t in CORE_TERMS)
            if not is_core:
                intact_stats["non_core_excluded"] += 1
                continue

            if any(m in it for m in ("MI:2364", "MI:0403")):
                intact_stats["proximity_type_excluded"] += 1
                continue

            det_m = cols[6]
            # PCA is RETAINED! Only exclude non-PCA proximity methods (BioID, confocal, etc.)
            if any(m in det_m for m in PROX_MI):
                intact_stats["proximity_method_excluded"] += 1
                continue

            tax_a_col = cols[9].strip() if len(cols) > 9 else ""
            tax_b_col = cols[10].strip() if len(cols) > 10 else ""
            m_a = re.search(r"taxid:(\d+)", tax_a_col)
            m_b = re.search(r"taxid:(\d+)", tax_b_col)
            ta_id = m_a.group(1) if m_a else None
            tb_id = m_b.group(1) if m_b else None
            if ta_id != "3702" or tb_id != "3702":
                intact_stats["cross_species_excluded"] += 1
                continue

            ta = cols[20] if len(cols) > 20 else ""
            tb = cols[21] if len(cols) > 21 else ""
            p1_raw = cols[0]
            p2_raw = cols[1]
            mol_a = get_mol_type(ta, p1_raw)
            mol_b = get_mol_type(tb, p2_raw)
            if mol_a != "protein" or mol_b != "protein":
                intact_stats["non_protein_excluded"] += 1
                continue

            intact_stats["records_passing"] += 1

            p1 = clean_id(p1_raw)
            p2 = clean_id(p2_raw)
            pair_key = "|".join(sorted([p1, p2]))

            source_rec_id = "NA"
            int_id_col = cols[13].strip() if len(cols) > 13 else ""
            if int_id_col and int_id_col != "-":
                for item in int_id_col.split("|"):
                    if item.startswith("intact:"):
                        source_rec_id = item.replace("intact:", "")
                        break
                if source_rec_id == "NA":
                    source_rec_id = int_id_col.split("|")[0]

            name_a = cols[4].split("|")[0] if len(cols) > 4 and cols[4] != "-" else p1
            name_b = cols[5].split("|")[0] if len(cols) > 5 and cols[5] != "-" else p2
            if "(" in name_a:
                m = re.search(r':([^(|]+)', name_a)
                if m: name_a = m.group(1)
            if "(" in name_b:
                m = re.search(r':([^(|]+)', name_b)
                if m: name_b = m.group(1)

            it_mi = "NA"
            it_name = "NA"
            m_it = re.search(r'psi-mi:"(MI:\d{4})"\(([^)]+)\)', it)
            if m_it:
                it_mi = m_it.group(1)
                it_name = m_it.group(2)

            dm_mi = "NA"
            dm_name = "NA"
            m_dm = re.search(r'psi-mi:"(MI:\d{4})"\(([^)]+)\)', det_m)
            if m_dm:
                dm_mi = m_dm.group(1)
                dm_name = m_dm.group(2)

            first_author = cols[7].strip() if len(cols) > 7 else "NA"
            pub_year = "NA"
            m_yr = re.search(r'\((\d{4})\)', first_author)
            if m_yr:
                pub_year = m_yr.group(1)

            pub_col = cols[8].strip() if len(cols) > 8 else "NA"
            pmid = "NA"
            if pub_col and pub_col != "-":
                for item in pub_col.split("|"):
                    if item.startswith("pubmed:"):
                        pmid = item.replace("pubmed:", "").strip()
                        break

            conf_col = cols[14].strip() if len(cols) > 14 else "NA"
            miscore = "nan"
            if "intact-miscore:" in conf_col:
                m_sc = re.search(r'intact-miscore:([0-9.]+)', conf_col)
                if m_sc: miscore = m_sc.group(1)

            intact_rows.append({
                "source_resource": "IntAct",
                "source_release": "Release 252",
                "source_file": "data/raw/arabidopsis/intact/intact_arabidopsis_taxid3702_release252.mitab27.txt",
                "source_row_number": str(line_no),
                "source_record_id": source_rec_id,
                "canonical_pair_key": pair_key,
                "participant_a_id": p1,
                "participant_a_namespace": p1_raw.split(":")[0] if ":" in p1_raw else "uniprotkb",
                "participant_a_original_id": p1_raw,
                "participant_a_name": name_a,
                "participant_a_aliases": cols[4] if len(cols) > 4 else "-",
                "participant_a_taxid": "3702",
                "participant_a_species": "Arabidopsis thaliana",
                "participant_a_biological_role": cols[16] if len(cols) > 16 else "-",
                "participant_a_experimental_role": cols[18] if len(cols) > 18 else "-",
                "participant_a_molecule_type": "protein",
                "participant_b_id": p2,
                "participant_b_namespace": p2_raw.split(":")[0] if ":" in p2_raw else "uniprotkb",
                "participant_b_original_id": p2_raw,
                "participant_b_name": name_b,
                "participant_b_aliases": cols[5] if len(cols) > 5 else "-",
                "participant_b_taxid": "3702",
                "participant_b_species": "Arabidopsis thaliana",
                "participant_b_biological_role": cols[17] if len(cols) > 17 else "-",
                "participant_b_experimental_role": cols[19] if len(cols) > 19 else "-",
                "participant_b_molecule_type": "protein",
                "is_protein_protein": "TRUE",
                "same_species_pair": "TRUE",
                "interaction_type_raw": it,
                "interaction_type_psi_mi": it_mi,
                "interaction_type_name": it_name,
                "detection_method_raw": det_m,
                "detection_method_psi_mi": dm_mi,
                "detection_method_name": dm_name,
                "first_author": first_author,
                "publication_id_raw": pub_col,
                "pmid": pmid,
                "publication_year": pub_year,
                "throughput_raw": "nan",
                "source_confidence_raw": conf_col,
                "intact_miscore": miscore,
                "xlms_protein_type": "nan",
                "xlms_peptide_a": "nan",
                "xlms_peptide_b": "nan",
                "xlms_linker": "nan",
                "xlms_score": "nan",
                "xlms_evalue": "nan",
                "record_hippie_score": "0.6293"
            })

    print("IntAct parsing stats:")
    for k, v in intact_stats.items():
        print(f"  {k}: {v}")
    df_intact = pd.DataFrame(intact_rows)

    print("\nCombining all evidence tables...")
    df_ev_updated = pd.concat([df_bg, df_intact, string_ev, xlms_ev], ignore_index=True)
    print(f"Total evidence rows: {len(df_ev_updated)}")

    ev_parquet_path = os.path.join(OUT_DIR, "ara_positive_evidence.parquet")
    ev_tsv_path = os.path.join(OUT_DIR, "ara_positive_evidence.tsv")
    print(f"Writing evidence files: {ev_parquet_path} & {ev_tsv_path}...")
    df_ev_updated.to_parquet(ev_parquet_path, index=False)
    df_ev_updated.to_csv(ev_tsv_path, sep="\t", index=False)
    print("Evidence files written.")

    print("\nBuilding deduplicated canonical pairs table for pure protein pairs...")
    pure_ev = df_ev_updated[df_ev_updated["is_protein_protein"] == "TRUE"]
    groups = pure_ev.groupby("canonical_pair_key")
    dedup_rows = []
    pair_id_counter = 1

    for pair_key, group in groups:
        p1, p2 = pair_key.split("|")
        sources = sorted(group["source_resource"].unique())
        sources_str = "|".join(sources)

        it_types = sorted(set(group["interaction_type_name"].dropna().unique()) - {"nan", "NA"})
        it_types_str = "|".join(it_types) if it_types else "physical"

        dm_methods = sorted(set(group["detection_method_name"].dropna().unique()) - {"nan", "NA"})
        dm_methods_str = "|".join(dm_methods) if dm_methods else "unspecified"

        authors = sorted(set(group["first_author"].dropna().unique()) - {"nan", "NA"})
        authors_str = "|".join(authors[:5])

        pmids = sorted(set(group["pmid"].dropna().unique()) - {"nan", "NA"})
        pmids_str = "|".join(pmids)
        pmid_count = len(pmids)

        years = sorted(set(group["publication_year"].dropna().unique()) - {"nan", "NA"})
        years_str = "|".join(years)

        rec_ids = sorted(set(group["source_record_id"].dropna().unique()) - {"nan", "NA"})
        rec_ids_str = "|".join(rec_ids[:10])

        max_miscore = "nan"
        miscores = [float(x) for x in group["intact_miscore"] if str(x) != "nan" and str(x) != "None" and str(x) != "NA"]
        if miscores:
            max_miscore = str(max(miscores))

        name_a = group["participant_a_name"].iloc[0]
        name_b = group["participant_b_name"].iloc[0]

        is_homo = "TRUE" if (p1 == p2) else "FALSE"
        pair_type = "homomer" if (p1 == p2) else "heteromer"

        # Calculate exact HIPPIE score based on studies (ns) and techniques (nt)
        canonical_tech_scores = {}
        for m in dm_methods:
            m_clean = m.strip().lower()
            if m_clean in TECHNIQUE_MAPPING:
                canon_key, val = TECHNIQUE_MAPPING[m_clean]
                canonical_tech_scores[canon_key] = val

        nt = float(sum(canonical_tech_scores.values())) if canonical_tech_scores else 1.0
        ns = pmid_count if pmid_count >= 1 else 1

        h_score = calc_hippie_score(ns=ns, nt=nt, no=0)
        conf_level = assign_confidence_level(h_score)

        dedup_rows.append({
            "pair_id": f"ARA_PAIR_{pair_id_counter:07d}",
            "canonical_pair_key": pair_key,
            "participant_a_id": p1,
            "participant_b_id": p2,
            "participant_a_name": name_a,
            "participant_b_name": name_b,
            "participant_a_type": "protein",
            "participant_b_type": "protein",
            "is_protein_protein_pair": "TRUE",
            "is_homomeric": is_homo,
            "pair_type": pair_type,
            "evidence_count": len(group),
            "sources": sources_str,
            "interaction_types": it_types_str,
            "detection_methods": dm_methods_str,
            "first_authors": authors_str,
            "pmid_count": pmid_count,
            "pmids": pmids_str,
            "publication_years": years_str,
            "throughput_classes": "mixed",
            "max_intact_miscore": max_miscore,
            "source_record_ids": rec_ids_str,
            "hippie_score": f"{h_score:.4f}",
            "hippie_confidence_level": conf_level
        })
        pair_id_counter += 1

    df_dedup_updated = pd.DataFrame(dedup_rows)
    print(f"Total deduplicated pure protein pairs: {len(df_dedup_updated):,}")

    dedup_parquet_path = os.path.join(OUT_DIR, "ara_positive_pairs_deduplicated.parquet")
    dedup_tsv_path = os.path.join(OUT_DIR, "ara_positive_pairs_deduplicated.tsv")
    print(f"Writing deduplicated pairs files: {dedup_parquet_path} & {dedup_tsv_path}...")
    df_dedup_updated.to_parquet(dedup_parquet_path, index=False)
    df_dedup_updated.to_csv(dedup_tsv_path, sep="\t", index=False)
    print("Deduplicated pairs files written.")

    all_prots = set(df_dedup_updated["participant_a_id"]) | set(df_dedup_updated["participant_b_id"])
    print(f"Total unique participating proteins: {len(all_prots):,}")

    hetero_count = (df_dedup_updated["is_homomeric"] == "FALSE").sum()
    homo_count = (df_dedup_updated["is_homomeric"] == "TRUE").sum()
    source_counts = df_dedup_updated["sources"].value_counts().to_dict()

    print(f"Heteromers: {hetero_count:,}, Homomers: {homo_count:,}")

    # Export High and High+Medium datasets
    df_high = df_dedup_updated[df_dedup_updated["hippie_confidence_level"] == "high"][["participant_a_id", "participant_b_id"]]
    df_high_med = df_dedup_updated[df_dedup_updated["hippie_confidence_level"].isin(["high", "medium"])][["participant_a_id", "participant_b_id"]]

    df_high_hetero = df_dedup_updated[(df_dedup_updated["hippie_confidence_level"] == "high") & (df_dedup_updated["is_homomeric"] == "FALSE")][["participant_a_id", "participant_b_id"]]
    df_high_med_hetero = df_dedup_updated[df_dedup_updated["hippie_confidence_level"].isin(["high", "medium"]) & (df_dedup_updated["is_homomeric"] == "FALSE")][["participant_a_id", "participant_b_id"]]

    high_tsv = os.path.join(OUT_DIR, "ara_positive_pairs_high_only.tsv")
    high_parquet = os.path.join(OUT_DIR, "ara_positive_pairs_high_only.parquet")
    high_med_tsv = os.path.join(OUT_DIR, "ara_positive_pairs_high_medium.tsv")
    high_med_parquet = os.path.join(OUT_DIR, "ara_positive_pairs_high_medium.parquet")

    high_hetero_tsv = os.path.join(OUT_DIR, "ara_positive_pairs_high_only_hetero.tsv")
    high_hetero_parquet = os.path.join(OUT_DIR, "ara_positive_pairs_high_only_hetero.parquet")
    high_med_hetero_tsv = os.path.join(OUT_DIR, "ara_positive_pairs_high_medium_hetero.tsv")
    high_med_hetero_parquet = os.path.join(OUT_DIR, "ara_positive_pairs_high_medium_hetero.parquet")

    df_high.to_csv(high_tsv, sep="\t", index=False)
    df_high.to_parquet(high_parquet, index=False)
    df_high_med.to_csv(high_med_tsv, sep="\t", index=False)
    df_high_med.to_parquet(high_med_parquet, index=False)

    df_high_hetero.to_csv(high_hetero_tsv, sep="\t", index=False)
    df_high_hetero.to_parquet(high_hetero_parquet, index=False)
    df_high_med_hetero.to_csv(high_med_hetero_tsv, sep="\t", index=False)
    df_high_med_hetero.to_parquet(high_med_hetero_parquet, index=False)

    print(f"High quality pairs: {len(df_high):,} (Hetero: {len(df_high_hetero):,})")
    print(f"High + Medium quality pairs: {len(df_high_med):,} (Hetero: {len(df_high_med_hetero):,})")

    # Build updated audit json
    audit_data = {
        "timestamp_utc": "2026-10-19T22:50:00.000000+00:00",
        "species": "Arabidopsis thaliana (taxid 3702)",
        "provenance_standard": "Database-Only Provenance (PCA strictly retained across all databases; IntAct core physical set included)",
        "filtering_rules": {
            "IntAct": "Core physical interactions (MI:0915 physical association, MI:0914 association, MI:0407 direct interaction, MI:0217 physical/phosphorylation); strictly retain PCA (MI:0112 split-ubiquitin, MI:0809 BiFC, MI:1204/MI:1037 split-luciferase, MI:0090); strictly exclude proximity (MI:2364, MI:0403, MI:1314 BioID); strictly exclude non-protein interactors (DNA, RNA, small molecules, gene IDs); strictly exclude cross-species interactions (taxid != 3702)",
            "BioGRID": "Physical interactions (sys_type == 'physical'); strictly retain PCA; strictly exclude Proximity Label-MS, Co-localization, genetic interactions, RNA complexes, and cross-species interactions (taxid != 3702)",
            "XL_MS": "pLink 3.2 < 1% PSM FDR; physical cross-links (<35A Euclidean constraint); unique pairs retained (total PSMs preserved)",
            "STRING": "STRING v12.0 physical links full: direct experiments > 0; strictly exclude transferred interologs (experiments == 0 and experiments_transferred > 0); exclude non-experimental links; consolidate to unique undirected pairs"
        },
        "evidence_table_columns_count": len(df_ev_updated.columns),
        "dedup_table_columns_count": len(df_dedup_updated.columns),
        "filtering_reconciliation": {
            "intact_positive_retained": len(df_intact),
            "intact_cross_species_excluded": int(intact_stats["cross_species_excluded"]),
            "biogrid_positive_retained": len(df_bg[df_bg["is_protein_protein"] == "TRUE"]),
            "biogrid_cross_species_excluded": int(bg_stats["cross_species_excluded"]),
            "xlms_positive_retained": len(xlms_ev),
            "string_positive_retained": len(string_ev),
            "total_positive_evidence_records": len(df_ev_updated),
            "total_pure_protein_evidence_records": len(pure_ev)
        },
        "deduplicated_pairs_summary": {
            "total_unique_pairs": len(df_dedup_updated),
            "heteromeric_pairs": int(hetero_count),
            "homomeric_pairs": int(homo_count),
            "pure_protein_protein_pairs": len(df_dedup_updated),
            "unique_proteins_all": len(all_prots),
            "source_combination_breakdown": source_counts
        },
        "hippie_scoring_summary": {
            "high_quality_count": len(df_high),
            "high_medium_count": len(df_high_med),
            "low_quality_count": int(len(df_dedup_updated) - len(df_high_med))
        },
        "output_files": {
            "evidence_tsv": ev_tsv_path,
            "evidence_parquet": ev_parquet_path,
            "dedup_pairs_tsv": dedup_tsv_path,
            "dedup_pairs_parquet": dedup_parquet_path,
            "high_only_tsv": high_tsv,
            "high_only_parquet": high_parquet,
            "high_only_hetero_tsv": high_hetero_tsv,
            "high_only_hetero_parquet": high_hetero_parquet,
            "high_medium_tsv": high_med_tsv,
            "high_medium_parquet": high_med_parquet,
            "high_medium_hetero_tsv": high_med_hetero_tsv,
            "high_medium_hetero_parquet": high_med_hetero_parquet
        }
    }

    audit_path = os.path.join(OUT_DIR, "ara_positive_audit.json")
    with open(audit_path, "w", encoding="utf-8") as f:
        json.dump(audit_data, f, indent=2)
    print(f"Audit json written to {audit_path}")

if __name__ == "__main__":
    run()
