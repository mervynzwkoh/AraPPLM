#!/usr/bin/env python3
"""
Phase 2: Cluster Quality Enrichment & Plant Threshold Derivation
===============================================================
Analyzes unsupervised clusters derived from PPLM paired embeddings, evaluates
statistical enrichment for experimental gold standards (PhoX XL-MS, multi-PMID,
biophysical binding), and derives empirical plant-calibrated thresholds.

Key Outputs:
  - data/processed/arabidopsis/plant_calibrated_thresholds.json
  - data/processed/arabidopsis/cluster_enrichment_summary.tsv
  - results/figures/ara_pplm_umap_clusters.png
  - results/figures/ara_pplm_umap_hippie.png
  - results/figures/ara_pplm_umap_gold_standards.png
  - results/figures/ara_cluster_quality_enrichment.png

Usage:
  python scripts/data_preparation/calibrate_plant_thresholds.py \
      --input data/processed/arabidopsis/ara_embeddings_clustered.parquet \
      --output_dir data/processed/arabidopsis \
      --figures_dir results/figures
"""

import os
import sys
import json
import argparse
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

# Matplotlib in headless/HPC mode
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

WORKSPACE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def parse_args():
    parser = argparse.ArgumentParser(description="Calibrate plant PPI thresholds from PPLM embedding clusters")
    parser.add_argument("--input", default="data/processed/arabidopsis/ara_embeddings_clustered.parquet",
                        help="Path to clustered embeddings table (.parquet or .tsv)")
    parser.add_argument("--output_dir", default="data/processed/arabidopsis",
                        help="Directory to save calibration json and tsv summary")
    parser.add_argument("--figures_dir", default="results/figures",
                        help="Directory to save generated figures")
    return parser.parse_args()


def load_dataset(input_path):
    print(f"Loading clustered dataset from: {input_path}")
    if input_path.endswith(".parquet") and os.path.exists(input_path):
        try:
            return pd.read_parquet(input_path)
        except Exception:
            pass
    tsv_path = input_path.replace(".parquet", ".tsv")
    if os.path.exists(tsv_path):
        return pd.read_csv(tsv_path, sep="\t")
    raise FileNotFoundError(f"Neither {input_path} nor {tsv_path} found.")


def analyze_clusters(df, cluster_col="pca_cluster_kmeans"):
    """
    Computes provenance quality statistics and hypergeometric enrichment
    for each cluster.
    """
    total_N = len(df)
    
    # Define gold standard flag: PhoX XL-MS OR >= 2 PMIDs
    is_xlms = df["sources"].astype(str).str.contains("PRIDE|pLink|XL-MS", case=False, na=False)
    pmid_cnt = pd.to_numeric(df["pmid_count"], errors="coerce").fillna(0)
    is_multipmid = pmid_cnt >= 2
    is_gold = is_xlms | is_multipmid
    total_gold = is_gold.sum()
    
    df["_is_gold"] = is_gold
    df["_hippie_float"] = pd.to_numeric(df["hippie_score"], errors="coerce")
    
    miscore_series = pd.to_numeric(df["max_intact_miscore"], errors="coerce")
    df["_miscore_float"] = miscore_series

    cluster_stats = []
    unique_clusters = sorted([c for c in df[cluster_col].unique() if c != -1])
    
    # Also evaluate noise cluster if HDBSCAN
    if -1 in df[cluster_col].unique():
        unique_clusters.append(-1)

    for c in unique_clusters:
        sub = df[df[cluster_col] == c]
        k = len(sub)
        if k == 0:
            continue
            
        k_gold = sub["_is_gold"].sum()
        k_xlms = (sub["sources"].astype(str).str.contains("PRIDE|pLink|XL-MS", case=False, na=False)).sum()
        k_multipmid = (pd.to_numeric(sub["pmid_count"], errors="coerce") >= 2).sum()
        
        # Hypergeometric p-value (survival function: P(X >= k_gold))
        p_val = hypergeom.sf(k_gold - 1, total_N, total_gold, k)
        
        h_scores = sub["_hippie_float"].dropna()
        m_scores = sub["_miscore_float"].dropna()
        
        cluster_stats.append({
            "cluster_id": c,
            "total_pairs": k,
            "pct_of_dataset": np.round(k / total_N * 100, 2),
            "gold_pairs_count": int(k_gold),
            "gold_pairs_pct": np.round(k_gold / k * 100, 2),
            "xlms_pairs_count": int(k_xlms),
            "multipmid_pairs_count": int(k_multipmid),
            "hypergeometric_p_value": float(p_val),
            "is_gold_enriched": bool(p_val < 0.01 and (k_gold / k) > (total_gold / total_N)),
            "hippie_mean": np.round(float(h_scores.mean()), 4) if len(h_scores) > 0 else None,
            "hippie_median": np.round(float(h_scores.median()), 4) if len(h_scores) > 0 else None,
            "hippie_q25": np.round(float(h_scores.quantile(0.25)), 4) if len(h_scores) > 0 else None,
            "hippie_q75": np.round(float(h_scores.quantile(0.75)), 4) if len(h_scores) > 0 else None,
            "miscore_median": np.round(float(m_scores.median()), 4) if len(m_scores) > 0 else None,
            "miscore_q25": np.round(float(m_scores.quantile(0.25)), 4) if len(m_scores) > 0 else None,
        })
        
    return pd.DataFrame(cluster_stats)


def plot_diagnostics(df, cluster_stats, figures_dir):
    os.makedirs(figures_dir, exist_ok=True)
    print(f"\nGenerating diagnostic figures in: {figures_dir}")
    
    # 1. UMAP by K-Means Clusters
    plt.figure(figsize=(9, 7))
    plt.scatter(df["umap_x"], df["umap_y"], c=df["pca_cluster_kmeans"], cmap="tab20", s=4, alpha=0.6)
    plt.title("PPLM Paired Embeddings UMAP — K-Means Clusters", fontsize=14, fontweight="bold")
    plt.xlabel("UMAP 1")
    plt.ylabel("UMAP 2")
    cbar = plt.colorbar()
    cbar.set_label("Cluster ID")
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, "ara_pplm_umap_clusters.png"), dpi=300)
    plt.close()
    
    # 2. UMAP by HIPPIE Score
    plt.figure(figsize=(9, 7))
    sc = plt.scatter(df["umap_x"], df["umap_y"], c=df["_hippie_float"], cmap="viridis", s=4, alpha=0.7)
    plt.title("PPLM Paired Embeddings UMAP — HIPPIE Quality Score", fontsize=14, fontweight="bold")
    plt.xlabel("UMAP 1")
    plt.ylabel("UMAP 2")
    cbar = plt.colorbar(sc)
    cbar.set_label("HIPPIE Score")
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, "ara_pplm_umap_hippie.png"), dpi=300)
    plt.close()

    # 3. UMAP Gold Standards Highlight
    plt.figure(figsize=(9, 7))
    # Background
    plt.scatter(df.loc[~df["_is_gold"], "umap_x"], df.loc[~df["_is_gold"], "umap_y"],
                c="lightgrey", s=3, alpha=0.4, label="Standard Evidence")
    # Gold
    plt.scatter(df.loc[df["_is_gold"], "umap_x"], df.loc[df["_is_gold"], "umap_y"],
                c="crimson", s=6, alpha=0.8, label="Gold Standard (PhoX XL-MS / Multi-PMID)")
    plt.title("PPLM Paired Embeddings UMAP — Gold Standard Evidence", fontsize=14, fontweight="bold")
    plt.xlabel("UMAP 1")
    plt.ylabel("UMAP 2")
    plt.legend(markerscale=3, loc="upper right")
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, "ara_pplm_umap_gold_standards.png"), dpi=300)
    plt.close()

    # 4. Cluster Quality Enrichment Bar Chart
    plt.figure(figsize=(10, 5))
    valid_stats = cluster_stats[cluster_stats["cluster_id"] != -1]
    colors = ["forestgreen" if row["is_gold_enriched"] else "steelblue" for _, row in valid_stats.iterrows()]
    plt.bar(valid_stats["cluster_id"].astype(str), valid_stats["gold_pairs_pct"], color=colors, edgecolor="black", alpha=0.8)
    plt.axhline(df["_is_gold"].mean() * 100, color="red", linestyle="--", label=f"Global Average ({df['_is_gold'].mean()*100:.1f}%)")
    plt.title("Proportion of Gold-Standard Evidence Across Clusters", fontsize=13, fontweight="bold")
    plt.xlabel("K-Means Cluster ID")
    plt.ylabel("% Gold-Standard Pairs")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, "ara_cluster_quality_enrichment.png"), dpi=300)
    plt.close()


def main():
    args = parse_args()
    print("=" * 70)
    print("CALIBRATING PLANT PPI CONFIDENCE THRESHOLDS FROM PPLM CLUSTERS")
    print("=" * 70)

    # 1. Load dataset
    input_path = os.path.join(WORKSPACE_ROOT, args.input) if not os.path.isabs(args.input) else args.input
    df = load_dataset(input_path)

    # 2. Analyze K-Means clusters
    print("\n1. Performing quality enrichment analysis on K-Means clusters...")
    cluster_stats = analyze_clusters(df, cluster_col="pca_cluster_kmeans")
    
    # 3. Display summary
    print("\nCluster Enrichment Summary:")
    print(cluster_stats[["cluster_id", "total_pairs", "gold_pairs_pct", "hypergeometric_p_value", "is_gold_enriched", "hippie_median", "hippie_q25"]])

    # 4. Save cluster summary table
    output_dir = os.path.join(WORKSPACE_ROOT, args.output_dir) if not os.path.isabs(args.output_dir) else args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    summary_tsv_path = os.path.join(output_dir, "cluster_enrichment_summary.tsv")
    cluster_stats.to_csv(summary_tsv_path, sep="\t", index=False)
    print(f"\nSaved cluster enrichment summary table: {summary_tsv_path}")

    # 5. Derive Plant-Calibrated Thresholds
    print("\n2. Deriving plant-calibrated thresholds...")
    enriched_clusters = cluster_stats[cluster_stats["is_gold_enriched"]]
    
    if len(enriched_clusters) > 0:
        # Enriched clusters identified
        calibrated_hippie = float(enriched_clusters["hippie_q25"].min())
        calibrated_miscore = float(enriched_clusters["miscore_q25"].dropna().min()) if enriched_clusters["miscore_q25"].dropna().size > 0 else 0.40
        enriched_cluster_ids = enriched_clusters["cluster_id"].tolist()
        print(f"   Found {len(enriched_clusters)} gold-standard enriched clusters: {enriched_cluster_ids}")
    else:
        # Fallback to top quartile of clusters
        top_clusters = cluster_stats.sort_values(by="gold_pairs_pct", ascending=False).head(3)
        calibrated_hippie = float(top_clusters["hippie_q25"].min())
        calibrated_miscore = float(top_clusters["miscore_q25"].dropna().min()) if top_clusters["miscore_q25"].dropna().size > 0 else 0.40
        enriched_cluster_ids = top_clusters["cluster_id"].tolist()
        print(f"   No formally enriched clusters via p < 0.01; using top 3 quality clusters: {enriched_cluster_ids}")

    # Compute qualified counts under calibrated rules
    h_pass = df["_hippie_float"] >= calibrated_hippie
    m_pass = (df["_miscore_float"] >= calibrated_miscore) & df["_miscore_float"].notna()
    calibrated_positives = df[h_pass | m_pass | df["_is_gold"]]

    calibration_audit = {
        "calibrated_hippie_threshold": np.round(calibrated_hippie, 4),
        "calibrated_intact_miscore_gate": np.round(calibrated_miscore, 4),
        "enriched_cluster_ids": enriched_cluster_ids,
        "total_deduplicated_pairs": len(df),
        "calibrated_high_confidence_pairs": len(calibrated_positives),
        "calibrated_high_confidence_pct": np.round(len(calibrated_positives) / len(df) * 100, 2),
        "breakdown": {
            "hippie_qualifying": int(h_pass.sum()),
            "intact_miscore_qualifying": int(m_pass.sum()),
            "gold_evidence_qualifying": int(df["_is_gold"].sum())
        },
        "comparison_with_legacy_human_cutoffs": {
            "legacy_hippie_072_count": int((df["_hippie_float"] >= 0.72).sum()),
            "legacy_miscore_045_count": int((df["_miscore_float"] >= 0.45).sum())
        }
    }

    calib_json_path = os.path.join(output_dir, "plant_calibrated_thresholds.json")
    with open(calib_json_path, "w") as f:
        json.dump(calibration_audit, f, indent=2)
    print(f"Saved plant calibrated thresholds to: {calib_json_path}")
    print(json.dumps(calibration_audit, indent=2))

    # 6. Plot diagnostics
    figures_dir = os.path.join(WORKSPACE_ROOT, args.figures_dir) if not os.path.isabs(args.figures_dir) else args.figures_dir
    plot_diagnostics(df, cluster_stats, figures_dir)
    print("\nPhase 2 Threshold Calibration Successfully Completed.")


if __name__ == "__main__":
    main()
