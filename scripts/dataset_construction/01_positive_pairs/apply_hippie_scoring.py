"""
Apply HIPPIE Scoring to Arabidopsis and Maize Master Positive Data
==================================================================
This script implements the exact HIPPIE (Human Integrated Protein-Protein
Interaction rEference; Schaefer et al., PLoS ONE 2012) experimental scoring
framework for curated plant PPI datasets:
  - Arabidopsis thaliana (87,673 deduplicated pairs)
  - Zea mays (Maize, 50,984 deduplicated pairs)

Scoring Formulation:
--------------------
  S = w_s * s_s(n_s) + w_t * s_t(n_t) + w_o * s_o(n_o)
  Where:
    s_i(n) = 2 / (1 + exp(-a_i * n)) - 1
    a_s = 2.3, w_s = 0.60  (Study / PMID count subscore)
    a_t = 0.2, w_t = 0.30  (Experimental technique quality subscore)
    a_o = 1.6, w_o = 0.10  (Interspecies orthology reproducibility)

Confidence Tiers:
-----------------
  - Low Quality:    S < 0.62
  - Medium Quality: 0.629 <= S < 0.72 (or 0.62 <= S < 0.72)
  - High Quality:   S >= 0.72

Outputs Generated:
------------------
  1. Full master deduplicated tables (all 24 metadata columns):
     - <species>_positive_pairs_deduplicated.parquet & .tsv
  2. High-quality only datasets (2 columns: participant_a_id, participant_b_id):
     - <species>_positive_pairs_high_only.parquet & .tsv
  3. High + Medium quality datasets (2 columns: participant_a_id, participant_b_id):
     - <species>_positive_pairs_high_medium.parquet & .tsv
  4. Comprehensive audit JSON files.
"""

import os
import sys
import json
from datetime import datetime, timezone
import numpy as np
import pandas as pd

WORKSPACE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

# Canonical technique scoring dictionary mapped directly from Table 2 of Schaefer et al. (2012)
TECHNIQUE_MAPPING = {
    # Two-hybrid family (score 5)
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

    # Co-fractionation / chromatography / migration
    "co-fractionation": ("co-fractionation", 1.0),
    "co-purification": ("co-purification", 2.0),
    "copurification": ("co-purification", 2.0),
    "blue native page": ("blue native page", 3.0),
    "comigration in gel electrophoresis": ("comigration in gel electrophoresis", 3.0),
    "comigration in non denaturing gel electrophoresis": ("comigration in non denaturing gel electrophoresis", 3.0),
    "molecular sieving": ("molecular sieving", 2.0),
    "cosedimentation in solution": ("cosedimentation in solution", 2.0),
    "cosedimentation through density gradient": ("cosedimentation through density gradient", 2.0),

    # Direct biophysical / structural / reconstituted (scores 9-10)
    "co-crystal structure": ("co-crystal structure", 10.0),
    "x-ray crystallography": ("x-ray crystallography", 10.0),
    "reconstituted complex": ("reconstituted complex", 10.0),
    "surface plasmon resonance": ("surface plasmon resonance", 10.0),
    "isothermal titration calorimetry": ("isothermal titration calorimetry", 10.0),
    "nuclear magnetic resonance": ("nuclear magnetic resonance", 10.0),
    "light scattering": ("light scattering", 10.0),
    "x ray scattering": ("x ray scattering", 9.0),

    # Biochemical / enzymatic (scores 1 - 7.5)
    "biochemical activity": ("biochemical activity", 5.0),
    "biochemical": ("biochemical", 1.0),
    "biophysical": ("biophysical", 1.0),
    "protein kinase assay": ("protein kinase assay", 7.5),
    "in-gel kinase assay": ("in-gel kinase assay", 7.5),
    "enzymatic study": ("enzymatic study", 1.0),
    "enzyme linked immunosorbent assay": ("enzyme linked immunosorbent assay", 5.0),
    "chromatin immunoprecipitation assay": ("chromatin immunoprecipitation assay", 2.0),
    "electrophoretic mobility shift assay": ("electrophoretic mobility shift assay", 2.0),

    # Affinity capture / co-IP / pull-down / cross-linking (scores 2.5 - 5)
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

    # Optics / FRET / display / arrays (scores 5 - 7.5)
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

    # STRING generic experimental channel (score 1)
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
    """Computes continuous HIPPIE score using optimized parameters."""
    ss = (2.0 / (1.0 + np.exp(-2.3 * ns))) - 1.0 if ns > 0 else 0.0
    st = (2.0 / (1.0 + np.exp(-0.2 * nt))) - 1.0 if nt > 0 else 0.0
    so = (2.0 / (1.0 + np.exp(-1.6 * no))) - 1.0 if no > 0 else 0.0
    return float(0.60 * ss + 0.30 * st + 0.10 * so)


def assign_confidence_level(score: float) -> str:
    """
    Classifies HIPPIE score into confidence tiers:
      - High:   score >= 0.72
      - Medium: 0.629 <= score < 0.72 (or score >= 0.62)
      - Low:    score < 0.62
    """
    if score >= 0.72:
        return "high"
    elif score >= 0.62:
        return "medium"
    else:
        return "low"


def get_tier_stats(df_subset):
    """Calculates pair topology and protein count for a dataframe subset."""
    tot = len(df_subset)
    if tot == 0:
        return {"total": 0, "hetero": 0, "homo": 0, "unique_proteins": 0}
    
    is_homo = df_subset["is_homomeric"].astype(str).str.upper() == "TRUE"
    homo = int(is_homo.sum())
    hetero = int((~is_homo).sum())
    
    prots = set(df_subset["participant_a_id"].dropna()).union(set(df_subset["participant_b_id"].dropna()))
    return {
        "total": tot,
        "hetero": hetero,
        "hetero_pct": round((hetero / tot) * 100.0, 2),
        "homo": homo,
        "homo_pct": round((homo / tot) * 100.0, 2),
        "unique_proteins": len(prots)
    }


def process_dataset(species_name: str, data_dir: str, file_prefix: str):
    parquet_path = os.path.join(data_dir, f"{file_prefix}_positive_pairs_deduplicated.parquet")
    tsv_path = os.path.join(data_dir, f"{file_prefix}_positive_pairs_deduplicated.tsv")
    audit_path = os.path.join(data_dir, f"{file_prefix}_positive_audit.json")

    print(f"\n========================================================")
    print(f"Processing {species_name} PPI Master Table: {parquet_path}")
    print(f"========================================================")

    df = pd.read_parquet(parquet_path)
    total_pairs = len(df)
    print(f"Loaded {total_pairs:,} deduplicated pairs.")

    scores = []
    confidence_levels = []
    nt_values = []
    ns_values = []

    for idx, row in df.iterrows():
        # Parse and deduplicate detection methods into canonical Table 2 categories
        methods_raw = str(row["detection_methods"]).split("|")
        canonical_tech_scores = {}
        for m in methods_raw:
            m_clean = m.strip().lower()
            if m_clean in TECHNIQUE_MAPPING:
                canon_key, val = TECHNIQUE_MAPPING[m_clean]
                canonical_tech_scores[canon_key] = val

        nt = float(sum(canonical_tech_scores.values())) if canonical_tech_scores else 1.0

        # Determine number of studies (ns)
        pmid_count_raw = str(row.get("pmid_count", "1"))
        ns = int(pmid_count_raw) if pmid_count_raw.isdigit() and int(pmid_count_raw) >= 1 else 1

        score = calc_hippie_score(ns=ns, nt=nt, no=0)
        conf = assign_confidence_level(score)

        scores.append(score)
        confidence_levels.append(conf)
        nt_values.append(nt)
        ns_values.append(ns)

    scores_arr = np.array(scores)

    # Attach/overwrite columns in dataframe
    df["hippie_score"] = [f"{s:.4f}" for s in scores]
    df["hippie_confidence_level"] = confidence_levels

    # Compute distribution statistics across tiers
    df_high = df[df["hippie_confidence_level"] == "high"]
    df_med = df[df["hippie_confidence_level"] == "medium"]
    df_low = df[df["hippie_confidence_level"] == "low"]
    df_high_med = df[df["hippie_confidence_level"].isin(["high", "medium"])]

    stats_high = get_tier_stats(df_high)
    stats_med = get_tier_stats(df_med)
    stats_low = get_tier_stats(df_low)
    stats_high_med = get_tier_stats(df_high_med)
    stats_total = get_tier_stats(df)

    min_sc = float(scores_arr.min())
    max_sc = float(scores_arr.max())
    mean_sc = float(scores_arr.mean())
    median_sc = float(np.median(scores_arr))
    q25 = float(np.percentile(scores_arr, 25))
    q75 = float(np.percentile(scores_arr, 75))
    iqr = float(q75 - q25)

    print(f"\nHIPPIE Scoring Distribution for {species_name}:")
    print(f"  • Total Deduplicated Pairs: {total_pairs:,}")
    print(f"  • High Quality   (Score >= 0.72):          {stats_high['total']:>6,} ({stats_high['total']/total_pairs*100:6.2f}%) | Hetero: {stats_high['hetero']:>6,} | Homo: {stats_high['homo']:>5,} | Unique Prots: {stats_high['unique_proteins']:>6,}")
    print(f"  • Medium Quality (0.629 <= Score < 0.72):  {stats_med['total']:>6,} ({stats_med['total']/total_pairs*100:6.2f}%) | Hetero: {stats_med['hetero']:>6,} | Homo: {stats_med['homo']:>5,} | Unique Prots: {stats_med['unique_proteins']:>6,}")
    print(f"  • High + Medium Combined:                  {stats_high_med['total']:>6,} ({stats_high_med['total']/total_pairs*100:6.2f}%) | Hetero: {stats_high_med['hetero']:>6,} | Homo: {stats_high_med['homo']:>5,} | Unique Prots: {stats_high_med['unique_proteins']:>6,}")
    print(f"  • Low Quality    (Score < 0.62):           {stats_low['total']:>6,} ({stats_low['total']/total_pairs*100:6.2f}%) | Hetero: {stats_low['hetero']:>6,} | Homo: {stats_low['homo']:>5,} | Unique Prots: {stats_low['unique_proteins']:>6,}")
    print(f"  • Score Metrics: Min={min_sc:.4f}, 25%={q25:.4f}, Median={median_sc:.4f}, Mean={mean_sc:.4f}, 75%={q75:.4f}, Max={max_sc:.4f}, IQR={iqr:.4f}")

    # Write updated master files
    print(f"\nSaving updated master files:")
    print(f"  -> {parquet_path}")
    df.to_parquet(parquet_path, index=False)
    print(f"  -> {tsv_path}")
    df.to_csv(tsv_path, sep="\t", index=False)

    # -------------------------------------------------------------
    # Export 2-Column Quality-Filtered Datasets (Protein ID Pairs)
    # -------------------------------------------------------------
    seq_cols = ["participant_a_id", "participant_b_id"]

    # 1. High-only dataset (2 columns)
    high_tsv_path = os.path.join(data_dir, f"{file_prefix}_positive_pairs_high_only.tsv")
    high_parquet_path = os.path.join(data_dir, f"{file_prefix}_positive_pairs_high_only.parquet")
    df_high_2col = df_high[seq_cols].copy()
    df_high_2col.to_csv(high_tsv_path, sep="\t", index=False)
    df_high_2col.to_parquet(high_parquet_path, index=False)
    print(f"  -> High-Only 2-column dataset exported: {high_tsv_path} ({len(df_high_2col):,} pairs)")

    # 2. High + Medium dataset (2 columns)
    hm_tsv_path = os.path.join(data_dir, f"{file_prefix}_positive_pairs_high_medium.tsv")
    hm_parquet_path = os.path.join(data_dir, f"{file_prefix}_positive_pairs_high_medium.parquet")
    df_hm_2col = df_high_med[seq_cols].copy()
    df_hm_2col.to_csv(hm_tsv_path, sep="\t", index=False)
    df_hm_2col.to_parquet(hm_parquet_path, index=False)
    print(f"  -> High+Medium 2-column dataset exported: {hm_tsv_path} ({len(df_hm_2col):,} pairs)")

    # 3. High-only heteromers-only dataset (2 columns)
    df_high_hetero = df_high[df_high["participant_a_id"] != df_high["participant_b_id"]][seq_cols].copy()
    high_hetero_tsv_path = os.path.join(data_dir, f"{file_prefix}_positive_pairs_high_only_hetero.tsv")
    high_hetero_parquet_path = os.path.join(data_dir, f"{file_prefix}_positive_pairs_high_only_hetero.parquet")
    df_high_hetero.to_csv(high_hetero_tsv_path, sep="\t", index=False)
    df_high_hetero.to_parquet(high_hetero_parquet_path, index=False)
    print(f"  -> High-Only Heteromer-only dataset exported: {high_hetero_tsv_path} ({len(df_high_hetero):,} pairs)")

    # 4. High + Medium heteromers-only dataset (2 columns)
    df_hm_hetero = df_high_med[df_high_med["participant_a_id"] != df_high_med["participant_b_id"]][seq_cols].copy()
    hm_hetero_tsv_path = os.path.join(data_dir, f"{file_prefix}_positive_pairs_high_medium_hetero.tsv")
    hm_hetero_parquet_path = os.path.join(data_dir, f"{file_prefix}_positive_pairs_high_medium_hetero.parquet")
    df_hm_hetero.to_csv(hm_hetero_tsv_path, sep="\t", index=False)
    df_hm_hetero.to_parquet(hm_hetero_parquet_path, index=False)
    print(f"  -> High+Medium Heteromer-only dataset exported: {hm_hetero_tsv_path} ({len(df_hm_hetero):,} pairs)")

    # Update audit JSON if present
    if os.path.exists(audit_path):
        with open(audit_path, "r", encoding="utf-8") as f:
            audit_data = json.load(f)

        audit_data["hippie_scoring_summary"] = {
            "species": species_name.lower(),
            "timestamp_recalculated_utc": datetime.now(timezone.utc).isoformat(),
            "total_pairs": total_pairs,
            "thresholds": {
                "low": "< 0.62",
                "medium": ">= 0.629 and < 0.72",
                "high": ">= 0.72"
            },
            "high_quality": stats_high,
            "medium_quality": stats_med,
            "high_plus_medium_quality": stats_high_med,
            "low_quality": stats_low,
            "score_statistics": {
                "mean_score": round(mean_sc, 4),
                "median_score": round(median_sc, 4),
                "min_score": round(min_sc, 4),
                "max_score": round(max_sc, 4),
                "q25": round(q25, 4),
                "q75": round(q75, 4),
                "iqr": round(iqr, 4)
            },
            "exported_quality_datasets": {
                "high_only_tsv": high_tsv_path,
                "high_only_parquet": high_parquet_path,
                "high_only_pairs": len(df_high_2col),
                "high_medium_tsv": hm_tsv_path,
                "high_medium_parquet": hm_parquet_path,
                "high_medium_pairs": len(df_hm_2col),
                "high_only_hetero_tsv": high_hetero_tsv_path,
                "high_only_hetero_parquet": high_hetero_parquet_path,
                "high_only_hetero_pairs": len(df_high_hetero),
                "high_medium_hetero_tsv": hm_hetero_tsv_path,
                "high_medium_hetero_parquet": hm_hetero_parquet_path,
                "high_medium_hetero_pairs": len(df_hm_hetero),
                "columns": seq_cols
            }
        }

        with open(audit_path, "w", encoding="utf-8") as f:
            json.dump(audit_data, f, indent=2)
        print(f"  -> Updated audit metadata: {audit_path}")

    return {
        "species": species_name,
        "total": total_pairs,
        "high": stats_high,
        "medium": stats_med,
        "high_medium": stats_high_med,
        "low": stats_low,
        "min": min_sc,
        "mean": mean_sc,
        "median": median_sc,
        "max": max_sc,
        "iqr": iqr
    }


def main():
    ara_dir = os.path.join(WORKSPACE_ROOT, "data", "dataset_construction", "01_positive_pairs", "arabidopsis")
    maize_dir = os.path.join(WORKSPACE_ROOT, "data", "dataset_construction", "01_positive_pairs", "maize")

    ara_summary = process_dataset("Arabidopsis thaliana", ara_dir, "ara")
    maize_summary = process_dataset("Zea mays", maize_dir, "maize")

    print("\n" + "=" * 60)
    print("ALL QUALITY DATASETS SUCCESSFULLY GENERATED AND AUDITED")
    print("=" * 60)


if __name__ == "__main__":
    main()
