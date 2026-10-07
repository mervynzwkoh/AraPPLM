#!/usr/bin/env python3
"""
Phase 2: Dimensionality Reduction & Unsupervised Clustering of PPLM Embeddings
==============================================================================
Performs PCA (top 50 components), K-Means, HDBSCAN, and 2D UMAP projection
on the 2,560-dimensional PPLM paired embeddings extracted in Phase 1.

Key Outputs:
  - data/processed/arabidopsis/ara_embeddings_pca50.npy (N x 50 float32 matrix)
  - data/processed/arabidopsis/ara_embeddings_clustered.parquet (aligned metadata + cluster labels + UMAP coordinates)
  - data/processed/arabidopsis/ara_embeddings_clustered.tsv (TSV backup)

Usage:
  python scripts/data_preparation/cluster_pair_embeddings.py \
      --embeddings data/processed/arabidopsis/ara_positive_embeddings_2560d.npy \
      --meta data/processed/arabidopsis/ara_positive_embeddings_meta.parquet \
      --output data/processed/arabidopsis/ara_embeddings_clustered.parquet \
      --pca_components 50 \
      --kmeans_k 15 \
      --hdbscan_min_cluster_size 50
"""

import os
import sys
import time
import argparse
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

# Try importing HDBSCAN from sklearn or standalone hdbscan
try:
    from sklearn.cluster import HDBSCAN
    HAS_HDBSCAN = True
except ImportError:
    try:
        import hdbscan
        HDBSCAN = hdbscan.HDBSCAN
        HAS_HDBSCAN = True
    except ImportError:
        HAS_HDBSCAN = False

# Try importing UMAP
try:
    import umap
    HAS_UMAP = True
except ImportError:
    HAS_UMAP = False
    from sklearn.manifold import TSNE

WORKSPACE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def parse_args():
    parser = argparse.ArgumentParser(description="Cluster PPLM paired embeddings and compute UMAP projection")
    parser.add_argument("--embeddings", default="data/processed/arabidopsis/ara_positive_embeddings_2560d.npy",
                        help="Path to .npy embeddings array")
    parser.add_argument("--meta", default="data/processed/arabidopsis/ara_positive_embeddings_meta.parquet",
                        help="Path to aligned metadata (.parquet or .tsv)")
    parser.add_argument("--pca_output", default="data/processed/arabidopsis/ara_embeddings_pca50.npy",
                        help="Path to save PCA 50 matrix")
    parser.add_argument("--output", default="data/processed/arabidopsis/ara_embeddings_clustered.parquet",
                        help="Output path for clustered table (.parquet)")
    parser.add_argument("--pca_components", type=int, default=50,
                        help="Number of PCA components (default: 50)")
    parser.add_argument("--kmeans_k", type=int, default=15,
                        help="Number of K-Means clusters (default: 15)")
    parser.add_argument("--hdbscan_min_cluster_size", type=int, default=50,
                        help="HDBSCAN min_cluster_size (default: 50)")
    parser.add_argument("--hdbscan_min_samples", type=int, default=15,
                        help="HDBSCAN min_samples (default: 15)")
    parser.add_argument("--n_neighbors", type=int, default=30,
                        help="UMAP n_neighbors (default: 30)")
    parser.add_argument("--min_dist", type=float, default=0.1,
                        help="UMAP min_dist (default: 0.1)")
    parser.add_argument("--random_state", type=int, default=42,
                        help="Random seed for reproducibility")
    return parser.parse_args()


def load_inputs(emb_path, meta_path):
    print(f"Loading metadata from: {meta_path}")
    if not os.path.exists(meta_path):
        tsv_path = meta_path.replace(".parquet", ".tsv")
        if os.path.exists(tsv_path):
            df_meta = pd.read_csv(tsv_path, sep="\t")
        else:
            raise FileNotFoundError(f"Metadata file not found: {meta_path}")
    else:
        try:
            df_meta = pd.read_parquet(meta_path)
        except Exception:
            tsv_path = meta_path.replace(".parquet", ".tsv")
            df_meta = pd.read_csv(tsv_path, sep="\t")

    num_meta = len(df_meta)
    print(f"Loaded metadata rows: {num_meta:,}")

    print(f"Loading embeddings from: {emb_path}")
    if not os.path.exists(emb_path):
        raise FileNotFoundError(f"Embeddings file not found: {emb_path}")

    # Support both standard .npy files and raw np.memmap binary arrays
    emb_data = None
    try:
        emb_data = np.load(emb_path, mmap_mode="r")
        num_samples, dim = emb_data.shape
        print(f"Loaded embeddings via np.load: {num_samples:,} samples x {dim:,} dimensions.")
    except Exception as e:
        print(f"np.load was unable to parse header ({e}). Falling back to raw binary np.memmap...")
        file_bytes = os.path.getsize(emb_path)
        dim = 2560
        num_samples = file_bytes // (dim * 4)  # 4 bytes per float32
        emb_data = np.memmap(emb_path, dtype="float32", mode="r", shape=(num_samples, dim))
        print(f"Loaded embeddings via raw np.memmap: {num_samples:,} samples x {dim:,} dimensions.")

    if len(df_meta) != num_samples:
        print(f"[Warning] Metadata count ({len(df_meta)}) != embeddings count ({num_samples}). Truncating to min.")
        min_len = min(len(df_meta), num_samples)
        df_meta = df_meta.iloc[:min_len].copy()
        emb_data = emb_data[:min_len]

    return emb_data, df_meta


def run_pca(embeddings, n_components=50, random_state=42):
    print(f"\n1. Fitting PCA (n_components={n_components})...")
    start = time.time()
    pca = PCA(n_components=n_components, random_state=random_state)
    pca_emb = pca.fit_transform(embeddings)
    elapsed = time.time() - start
    
    total_var = np.sum(pca.explained_variance_ratio_) * 100
    print(f"   PCA completed in {elapsed:.1f}s.")
    print(f"   Total variance explained by top {n_components} PCs: {total_var:.2f}%")
    print(f"   Top 5 individual PC variance ratios: {np.round(pca.explained_variance_ratio_[:5]*100, 2)}")
    return pca_emb, pca


def run_kmeans(pca_emb, k=15, random_state=42):
    print(f"\n2. Fitting K-Means (k={k})...")
    start = time.time()
    kmeans = KMeans(n_clusters=k, random_state=random_state, n_init=10)
    labels = kmeans.fit_predict(pca_emb)
    elapsed = time.time() - start
    print(f"   K-Means completed in {elapsed:.1f}s.")
    return labels


def run_hdbscan(pca_emb, min_cluster_size=50, min_samples=15):
    if not HAS_HDBSCAN:
        print("\n[Warning] HDBSCAN not installed. Skipping HDBSCAN clustering.")
        return np.full(len(pca_emb), -1, dtype=int)

    print(f"\n3. Fitting HDBSCAN (min_cluster_size={min_cluster_size}, min_samples={min_samples})...")
    start = time.time()
    try:
        clusterer = HDBSCAN(min_cluster_size=min_cluster_size, min_samples=min_samples)
        labels = clusterer.fit_predict(pca_emb)
    except Exception as e:
        print(f"   HDBSCAN failed ({e}). Defaulting to noise label -1.")
        labels = np.full(len(pca_emb), -1, dtype=int)
        
    elapsed = time.time() - start
    num_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    noise_count = np.sum(labels == -1)
    noise_pct = noise_count / len(labels) * 100
    print(f"   HDBSCAN completed in {elapsed:.1f}s.")
    print(f"   Identified {num_clusters} dense clusters. Noise points: {noise_count:,} ({noise_pct:.1f}%)")
    return labels


def run_umap_projection(pca_emb, n_neighbors=30, min_dist=0.1, random_state=42):
    start = time.time()
    if HAS_UMAP:
        print(f"\n4. Computing UMAP 2D projection (n_neighbors={n_neighbors}, min_dist={min_dist}, metric='cosine')...")
        reducer = umap.UMAP(
            n_components=2,
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            metric="cosine",
            random_state=random_state,
            n_jobs=-1
        )
        umap_coords = reducer.fit_transform(pca_emb)
        print(f"   UMAP projection completed in {time.time() - start:.1f}s.")
    else:
        print("\n[Warning] umap-learn not installed. Falling back to t-SNE 2D projection...")
        tsne = TSNE(n_components=2, random_state=random_state, n_jobs=-1)
        umap_coords = tsne.fit_transform(pca_emb)
        print(f"   t-SNE projection completed in {time.time() - start:.1f}s.")
        
    return umap_coords


def main():
    args = parse_args()
    print("=" * 70)
    print("PHASE 2: PPLM EMBEDDING CLUSTERING & 2D PROJECTION")
    print("=" * 70)

    # 1. Load inputs
    emb_path = os.path.join(WORKSPACE_ROOT, args.embeddings) if not os.path.isabs(args.embeddings) else args.embeddings
    meta_path = os.path.join(WORKSPACE_ROOT, args.meta) if not os.path.isabs(args.meta) else args.meta
    embeddings, df_meta = load_inputs(emb_path, meta_path)

    # 2. PCA
    pca_emb, _ = run_pca(embeddings, n_components=args.pca_components, random_state=args.random_state)
    
    # Save PCA array
    pca_out = os.path.join(WORKSPACE_ROOT, args.pca_output) if not os.path.isabs(args.pca_output) else args.pca_output
    os.makedirs(os.path.dirname(pca_out), exist_ok=True)
    np.save(pca_out, pca_emb)
    print(f"   Saved PCA matrix to: {pca_out}")

    # 3. K-Means
    kmeans_labels = run_kmeans(pca_emb, k=args.kmeans_k, random_state=args.random_state)

    # 4. HDBSCAN
    hdbscan_labels = run_hdbscan(
        pca_emb,
        min_cluster_size=args.hdbscan_min_cluster_size,
        min_samples=args.hdbscan_min_samples
    )

    # 5. UMAP 2D
    umap_coords = run_umap_projection(
        pca_emb,
        n_neighbors=args.n_neighbors,
        min_dist=args.min_dist,
        random_state=args.random_state
    )

    # 6. Assemble output table
    print("\n5. Assembling clustered dataset...")
    df_result = df_meta.copy()
    df_result["pca_cluster_kmeans"] = kmeans_labels
    df_result["pca_cluster_hdbscan"] = hdbscan_labels
    df_result["umap_x"] = np.round(umap_coords[:, 0], 5)
    df_result["umap_y"] = np.round(umap_coords[:, 1], 5)

    # 7. Export outputs
    output_path = os.path.join(WORKSPACE_ROOT, args.output) if not os.path.isabs(args.output) else args.output
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    try:
        df_result.to_parquet(output_path, compression="snappy", index=False)
        print(f"   Successfully exported Parquet to: {output_path}")
    except Exception as e:
        print(f"   Parquet write failed ({e}). Writing TSV only.")

    tsv_out = output_path.replace(".parquet", ".tsv")
    df_result.to_csv(tsv_out, sep="\t", index=False)
    print(f"   Successfully exported TSV to: {tsv_out}")

    print("\n" + "=" * 70)
    print("CLUSTERING COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
