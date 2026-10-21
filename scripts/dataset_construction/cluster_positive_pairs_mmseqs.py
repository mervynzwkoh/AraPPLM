#!/usr/bin/env python3
"""
cluster_positive_pairs_mmseqs.py
================================
Accelerated Pair-Level Redundancy Reduction for Arabidopsis Protein-Protein Interactions
using MMseqs2 alignments and sparse inverted-index candidate evaluation.

Similarity metric between pair P_i = (A, B) and P_j = (C, D):
    S(P_i, P_j) = max[ min(s(A,C), s(B,D)), min(s(A,D), s(B,C)) ]

where s(X, Y) is protein sequence similarity (identity with coverage filter).

Usage:
    # Single threshold evaluation:
    python3 cluster_positive_pairs_mmseqs.py \\
        --pairs_tsv data/dataset_construction/01_positive_pairs/arabidopsis/ara_positive_pairs_high_only_hetero.tsv \\
        --alignments_tsv data/dataset_construction/03_homology_clustering/arabidopsis/ara_prot_all_vs_all.tsv \\
        --theta 0.40 --cmin 0.50

    # Full multi-parameter grid sweep (theta in [0.40, 0.70] x cmin in [0.40, 0.50, 0.70]):
    python3 cluster_positive_pairs_mmseqs.py \\
        --pairs_tsv data/dataset_construction/01_positive_pairs/arabidopsis/ara_positive_pairs_high_only_hetero.tsv \\
        --alignments_tsv data/dataset_construction/03_homology_clustering/arabidopsis/ara_prot_all_vs_all.tsv \\
        --sweep
"""

import os
import sys
import csv
import time
import argparse
from collections import defaultdict


def load_fasta_lengths(fasta_path):
    """Load protein sequence lengths from FASTA file."""
    if not fasta_path or not os.path.isfile(fasta_path):
        return {}
    lengths = {}
    current_id = None
    current_len = 0
    with open(fasta_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    lengths[current_id] = current_len
                current_id = line[1:].split()[0]
                current_len = 0
            else:
                current_len += len(line)
        if current_id is not None:
            lengths[current_id] = current_len
    return lengths


def load_metadata_scores(metadata_tsv):
    """Load HIPPIE score and evidence count from deduplicated TSV if available."""
    metadata = {}
    if not metadata_tsv or not os.path.isfile(metadata_tsv):
        return metadata
    
    with open(metadata_tsv, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            a = row.get("participant_a_id", "").strip()
            b = row.get("participant_b_id", "").strip()
            if not a or not b:
                continue
            pair_key = (min(a, b), max(a, b))
            try:
                score = float(row.get("hippie_score", 0.0))
            except ValueError:
                score = 0.0
            try:
                evidence = int(row.get("evidence_count", 1))
            except ValueError:
                evidence = 1
            metadata[pair_key] = (score, evidence)
    return metadata


def load_positive_pairs(pairs_tsv, metadata_map, fasta_lengths):
    """
    Load pairs from TSV and attach metadata (score, evidence_count, total_length).
    Standardizes each pair such that participant_a_id <= participant_b_id.
    """
    pairs = []
    seen = set()
    with open(pairs_tsv, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            a = row.get("participant_a_id", "").strip()
            b = row.get("participant_b_id", "").strip()
            if not a or not b:
                continue
            pair_key = (min(a, b), max(a, b))
            if pair_key in seen:
                continue
            seen.add(pair_key)
            
            # Lookup metadata
            if pair_key in metadata_map:
                score, evidence = metadata_map[pair_key]
            else:
                try:
                    score = float(row.get("score", row.get("hippie_score", 0.72)))
                except ValueError:
                    score = 0.72
                try:
                    evidence = int(row.get("evidence_count", 1))
                except ValueError:
                    evidence = 1
            
            len_a = fasta_lengths.get(pair_key[0], 0)
            len_b = fasta_lengths.get(pair_key[1], 0)
            total_len = len_a + len_b

            pairs.append({
                "pair_id": len(pairs) + 1,
                "a": pair_key[0],
                "b": pair_key[1],
                "score": score,
                "evidence_count": evidence,
                "total_len": total_len,
            })
    return pairs


def load_raw_alignments(alignments_tsv):
    """
    Parse MMseqs2 convertalis TSV once into memory.
    Expected columns: query, target, fident, ..., evalue, ..., qcov, tcov
    """
    print(f"Loading MMseqs2 alignments from: {alignments_tsv}")
    start_t = time.time()
    alignments = []
    
    with open(alignments_tsv, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            q = row.get("query", "").strip()
            t = row.get("target", "").strip()
            if not q or not t or q == t:
                continue
            
            # Handle fident or pident
            if "fident" in row:
                try:
                    fident = float(row["fident"])
                except ValueError:
                    continue
                if fident > 1.0:
                    fident /= 100.0
            elif "pident" in row:
                try:
                    fident = float(row["pident"]) / 100.0
                except ValueError:
                    continue
            else:
                continue
            
            try:
                qcov = float(row.get("qcov", 0.0))
                if qcov > 1.0:
                    qcov /= 100.0
                tcov = float(row.get("tcov", 0.0))
                if tcov > 1.0:
                    tcov /= 100.0
                evalue = float(row.get("evalue", 0.0))
            except ValueError:
                continue
            
            alignments.append((q, t, fident, qcov, tcov, evalue))
            
    print(f"Loaded {len(alignments):,} raw alignments in {time.time() - start_t:.2f}s")
    return alignments


def build_homology_index(raw_alignments, theta, cmin, max_evalue=1e-3):
    """
    Build single-protein similarity lookup and homology sets for a given (theta, cmin).
    """
    sim = {}
    homologs = defaultdict(set)
    
    for q, t, fident, qcov, tcov, evalue in raw_alignments:
        if fident >= theta and min(qcov, tcov) >= cmin and evalue <= max_evalue:
            pair = (min(q, t), max(q, t))
            # Keep highest fident if multiple HSPs exist
            if pair not in sim or fident > sim[pair]:
                sim[pair] = fident
            homologs[q].add(t)
            homologs[t].add(q)
            
    return sim, homologs


def run_pair_clustering(pairs, raw_alignments, theta, cmin, max_evalue=1e-3):
    """
    Perform accelerated pair clustering for a given (theta, cmin).
    """
    t0 = time.time()
    n_pairs = len(pairs)
    
    # 1. Build homology lookup
    sim, homologs = build_homology_index(raw_alignments, theta, cmin, max_evalue)
    
    def get_prot_sim(p1, p2):
        if p1 == p2:
            return 1.0
        return sim.get((min(p1, p2), max(p1, p2)), 0.0)

    # 2. Build inverted index: protein -> set of pair indices
    prot_to_pairs = defaultdict(set)
    for idx, p in enumerate(pairs):
        prot_to_pairs[p["a"]].add(idx)
        prot_to_pairs[p["b"]].add(idx)

    # 3. Find candidate redundant pairs and compute S(P_i, P_j)
    redundant_neighbors = defaultdict(list)
    candidate_checks = 0
    redundant_edges_count = 0
    
    for i, p_i in enumerate(pairs):
        a_i, b_i = p_i["a"], p_i["b"]
        
        # Homolog sets including self
        h_a = homologs[a_i] | {a_i}
        h_b = homologs[b_i] | {b_i}
        
        # Set of pairs touching H(a_i)
        pairs_touching_a = set()
        for x in h_a:
            pairs_touching_a.update(prot_to_pairs[x])
            
        # Set of pairs touching H(b_i)
        pairs_touching_b = set()
        for y in h_b:
            pairs_touching_b.update(prot_to_pairs[y])
            
        # Candidate pairs must be in intersection
        candidates = pairs_touching_a & pairs_touching_b
        
        for j in candidates:
            if j <= i:
                continue
            candidate_checks += 1
            p_j = pairs[j]
            c_j, d_j = p_j["a"], p_j["b"]
            
            # Orientation 1: (A~C, B~D)
            s_ac = get_prot_sim(a_i, c_j)
            s_bd = get_prot_sim(b_i, d_j)
            s1 = min(s_ac, s_bd)
            
            # Orientation 2: (A~D, B~C)
            s_ad = get_prot_sim(a_i, d_j)
            s_bc = get_prot_sim(b_i, c_j)
            s2 = min(s_ad, s_bc)
            
            s_pair = max(s1, s2)
            if s_pair >= theta:
                redundant_neighbors[i].append((j, s_pair))
                redundant_neighbors[j].append((i, s_pair))
                redundant_edges_count += 1

    # 4. Greedy Representative Selection (CD-HIT style)
    # Priority: HIPPIE score desc -> evidence_count desc -> total_len desc -> pair IDs asc
    sorted_indices = sorted(
        range(n_pairs),
        key=lambda idx: (
            -pairs[idx]["score"],
            -pairs[idx]["evidence_count"],
            -pairs[idx]["total_len"],
            pairs[idx]["a"],
            pairs[idx]["b"],
        )
    )
    
    unassigned = set(range(n_pairs))
    clusters = []
    pair_assignments = {}  # idx -> (cluster_id, is_rep, sim_to_rep)
    
    for rep_idx in sorted_indices:
        if rep_idx not in unassigned:
            continue
        
        cluster_id = len(clusters) + 1
        unassigned.remove(rep_idx)
        pair_assignments[rep_idx] = (cluster_id, True, 1.0)
        
        members = [rep_idx]
        for neighbor_idx, sim_val in redundant_neighbors[rep_idx]:
            if neighbor_idx in unassigned:
                unassigned.remove(neighbor_idx)
                members.append(neighbor_idx)
                pair_assignments[neighbor_idx] = (cluster_id, False, round(sim_val, 4))
                
        clusters.append({
            "cluster_id": cluster_id,
            "rep_idx": rep_idx,
            "size": len(members),
            "members": members,
        })

    elapsed = time.time() - t0
    
    # 5. Calculate Metrics
    n_rep = len(clusters)
    rep_proteins = set()
    for c in clusters:
        rep_p = pairs[c["rep_idx"]]
        rep_proteins.add(rep_p["a"])
        rep_proteins.add(rep_p["b"])
    m_rep = len(rep_proteins)
    
    redundancy_removed_pct = ((n_pairs - n_rep) / n_pairs) * 100.0 if n_pairs > 0 else 0.0
    cluster_sizes = [c["size"] for c in clusters]
    avg_size = sum(cluster_sizes) / len(cluster_sizes) if cluster_sizes else 0.0
    max_size = max(cluster_sizes) if cluster_sizes else 0
    min_size = min(cluster_sizes) if cluster_sizes else 0
    non_singleton_sizes = [s for s in cluster_sizes if s > 1]
    min_non_singleton = min(non_singleton_sizes) if non_singleton_sizes else 0
    singletons_count = sum(1 for s in cluster_sizes if s == 1)
    multi_pair_clusters_count = len(non_singleton_sizes)

    metrics = {
        "theta": theta,
        "cmin": cmin,
        "total_pairs": n_pairs,
        "representative_pairs": n_rep,
        "unique_proteins_kept": m_rep,
        "redundancy_removed_pct": round(redundancy_removed_pct, 2),
        "total_clusters": len(clusters),
        "avg_cluster_size": round(avg_size, 3),
        "max_cluster_size": max_size,
        "min_cluster_size": min_size,
        "min_non_singleton_size": min_non_singleton,
        "singletons": singletons_count,
        "multi_pair_clusters": multi_pair_clusters_count,
        "redundant_edges_found": redundant_edges_count,
        "candidate_checks": candidate_checks,
        "elapsed_seconds": round(elapsed, 2),
    }

    return metrics, clusters, pair_assignments


def export_cluster_results(output_dir, theta, cmin, pairs, clusters, pair_assignments):
    """Export representative TSV and pair-to-cluster mapping TSV."""
    os.makedirs(output_dir, exist_ok=True)
    mappings_dir = os.path.join(output_dir, "mappings")
    os.makedirs(mappings_dir, exist_ok=True)
    
    theta_str = f"{theta:.2f}".replace(".", "_")
    cmin_str = f"{cmin:.2f}".replace(".", "_")
    
    rep_tsv = os.path.join(output_dir, f"ara_positive_pairs_nr_id{theta_str}_cov{cmin_str}.tsv")
    mapping_tsv = os.path.join(mappings_dir, f"ara_pair_clusters_mapping_id{theta_str}_cov{cmin_str}.tsv")
    
    # 1. Export Representative Pairs TSV
    with open(rep_tsv, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["participant_a_id", "participant_b_id", "cluster_id", "hippie_score", "evidence_count"])
        for c in clusters:
            p = pairs[c["rep_idx"]]
            writer.writerow([p["a"], p["b"], c["cluster_id"], f"{p['score']:.4f}", p["evidence_count"]])
            
    # 2. Export Mappings TSV
    with open(mapping_tsv, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["pair_id", "participant_a_id", "participant_b_id", "cluster_id", "is_representative", "similarity_to_rep"])
        for idx, p in enumerate(pairs):
            cluster_id, is_rep, sim_val = pair_assignments[idx]
            writer.writerow([p["pair_id"], p["a"], p["b"], cluster_id, 1 if is_rep else 0, sim_val])

    print(f"  -> Saved representatives to: {rep_tsv}")
    print(f"  -> Saved cluster mappings to: {mapping_tsv}")
    return rep_tsv, mapping_tsv


def write_summary_tables(output_dir, sweep_metrics):
    """Write comprehensive summary metrics to TSV and Markdown."""
    summary_tsv = os.path.join(output_dir, "summary_clustering_metrics.tsv")
    summary_md = os.path.join(output_dir, "summary_clustering_metrics.md")
    
    headers = [
        "Identity_Threshold_theta",
        "Coverage_Threshold_cmin",
        "Total_Input_Pairs",
        "Representative_Pairs",
        "Unique_Proteins_Kept",
        "Redundancy_Removed_Pct",
        "Total_Clusters",
        "Avg_Cluster_Size",
        "Max_Cluster_Size",
        "Min_Cluster_Size",
        "Min_Non_Singleton_Size",
        "Singletons_Count",
        "Multi_Pair_Clusters",
        "Elapsed_Seconds"
    ]
    
    with open(summary_tsv, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(headers)
        for m in sweep_metrics:
            writer.writerow([
                f"{m['theta']:.2f}",
                f"{m['cmin']:.2f}",
                m["total_pairs"],
                m["representative_pairs"],
                m["unique_proteins_kept"],
                f"{m['redundancy_removed_pct']:.2f}%",
                m["total_clusters"],
                f"{m['avg_cluster_size']:.3f}",
                m["max_cluster_size"],
                m["min_cluster_size"],
                m["min_non_singleton_size"],
                m["singletons"],
                m["multi_pair_clusters"],
                f"{m['elapsed_seconds']:.2f}s"
            ])
            
    with open(summary_md, "w", encoding="utf-8") as f:
        f.write("# Summary of Arabidopsis Pair Redundancy Clustering Parameter Sweep\n\n")
        f.write("| Identity Threshold ($\\theta$) | Coverage ($c_{\\min}$) | Non-Redundant Pairs | Unique Proteins | Redundancy Removed (%) | Total Clusters | Avg Cluster Size | Max Cluster Size | Min Cluster Size | Multi-Pair Clusters | Singletons | Runtime |\n")
        f.write("|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|\n")
        for m in sweep_metrics:
            f.write(
                f"| **{m['theta']:.2f}** | **{m['cmin']:.2f}** | "
                f"{m['representative_pairs']:,} | {m['unique_proteins_kept']:,} | "
                f"**{m['redundancy_removed_pct']:.2f}%** | {m['total_clusters']:,} | "
                f"{m['avg_cluster_size']:.3f} | {m['max_cluster_size']} | "
                f"{m['min_cluster_size']} | {m['multi_pair_clusters']:,} | "
                f"{m['singletons']:,} | {m['elapsed_seconds']:.1f}s |\n"
            )
        f.write("\n\n*Note: Redundancy removed is defined as (Total Input Pairs - Non-Redundant Pairs) / Total Input Pairs.*\n")

    print(f"\nSaved overall summary TSV: {summary_tsv}")
    print(f"Saved overall summary Markdown: {summary_md}")


def main():
    parser = argparse.ArgumentParser(description="MMseqs2 Pair Redundancy Clustering for Arabidopsis PPIs")
    parser.add_argument("--pairs_tsv", type=str, default="data/dataset_construction/01_positive_pairs/arabidopsis/ara_positive_pairs_high_only_hetero.tsv",
                        help="Path to positive pairs TSV")
    parser.add_argument("--alignments_tsv", type=str, default="data/dataset_construction/03_homology_clustering/arabidopsis/ara_prot_all_vs_all.tsv",
                        help="Path to MMseqs2 convertalis TSV")
    parser.add_argument("--metadata_tsv", type=str, default="data/dataset_construction/01_positive_pairs/arabidopsis/ara_positive_pairs_deduplicated.tsv",
                        help="Path to deduplicated pairs TSV for score/evidence metadata lookup")
    parser.add_argument("--fasta", type=str, default="data/dataset_construction/01_positive_pairs/arabidopsis/ara_positive_pairs_high_only_hetero.fasta",
                        help="Path to FASTA for sequence length lookup")
    parser.add_argument("--output_dir", type=str, default="data/dataset_construction/03_homology_clustering/arabidopsis",
                        help="Output directory")
    parser.add_argument("--theta", type=float, default=0.40, help="Sequence identity threshold (default: 0.40)")
    parser.add_argument("--cmin", type=float, default=0.50, help="Alignment coverage threshold (default: 0.50)")
    parser.add_argument("--sweep", action="store_true", help="Run full parameter grid sweep: theta in [0.40, 0.70] x cmin in [0.40, 0.50, 0.70]")
    
    args = parser.parse_args()

    # Dynamic path resolution if run from repo root or scripts dir
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    
    def resolve_path(p):
        if not os.path.isabs(p) and not os.path.exists(p):
            alt = os.path.join(repo_root, p)
            if os.path.exists(alt):
                return alt
        return p

    pairs_tsv = resolve_path(args.pairs_tsv)
    alignments_tsv = resolve_path(args.alignments_tsv)
    metadata_tsv = resolve_path(args.metadata_tsv)
    fasta_path = resolve_path(args.fasta)
    output_dir = resolve_path(args.output_dir)

    print("=====================================================================")
    print("      AraPPLM: MMseqs2 Accelerated Pair Redundancy Clustering        ")
    print("=====================================================================")
    print(f"Input Pairs TSV:   {pairs_tsv}")
    print(f"Alignments TSV:    {alignments_tsv}")
    print(f"Metadata TSV:      {metadata_tsv}")
    print(f"FASTA Sequence:    {fasta_path}")
    print(f"Output Directory:  {output_dir}")
    print("=====================================================================")

    if not os.path.isfile(pairs_tsv):
        print(f"Error: Pairs file not found: {pairs_tsv}", file=sys.stderr)
        sys.exit(1)
    if not os.path.isfile(alignments_tsv):
        print(f"Error: Alignments file not found: {alignments_tsv}", file=sys.stderr)
        print("Please run MMseqs2 all-vs-all search first using run_pair_clustering_pipeline.sh", file=sys.stderr)
        sys.exit(1)

    # 1. Load supporting data
    print("Loading protein sequence lengths...")
    fasta_lengths = load_fasta_lengths(fasta_path)
    print(f"Loaded lengths for {len(fasta_lengths):,} proteins")

    print("Loading metadata scores (HIPPIE confidence & evidence)...")
    metadata_map = load_metadata_scores(metadata_tsv)
    print(f"Loaded metadata for {len(metadata_map):,} pairs")

    print("Loading positive pairs...")
    pairs = load_positive_pairs(pairs_tsv, metadata_map, fasta_lengths)
    print(f"Loaded {len(pairs):,} unique positive heteromeric pairs")

    # 2. Load MMseqs2 alignments
    raw_alignments = load_raw_alignments(alignments_tsv)

    # 3. Determine run grid
    if args.sweep:
        grid = [
            (0.40, 0.40),
            (0.40, 0.50),
            (0.40, 0.70),
            (0.70, 0.40),
            (0.70, 0.50),
            (0.70, 0.70),
        ]
        print(f"\n--- Starting Parameter Sweep Across {len(grid)} Grid Points ---")
    else:
        grid = [(args.theta, args.cmin)]

    sweep_metrics = []
    for theta_val, cmin_val in grid:
        print(f"\nEvaluating: theta = {theta_val:.2f}, cmin = {cmin_val:.2f} ...")
        metrics, clusters, pair_assignments = run_pair_clustering(
            pairs, raw_alignments, theta=theta_val, cmin=cmin_val
        )
        sweep_metrics.append(metrics)
        
        print(f"  Results for theta={theta_val:.2f}, cmin={cmin_val:.2f}:")
        print(f"    Representative Pairs Kept: {metrics['representative_pairs']:,} / {metrics['total_pairs']:,}")
        print(f"    Unique Proteins Kept:      {metrics['unique_proteins_kept']:,}")
        print(f"    Redundancy Pruned:         {metrics['redundancy_removed_pct']:.2f}%")
        print(f"    Total Clusters:            {metrics['total_clusters']:,}")
        print(f"    Avg Cluster Size:          {metrics['avg_cluster_size']:.3f}")
        print(f"    Max Cluster Size:          {metrics['max_cluster_size']}")
        print(f"    Multi-Pair Clusters:       {metrics['multi_pair_clusters']:,}")
        print(f"    Singletons:                {metrics['singletons']:,}")
        print(f"    Elapsed Time:              {metrics['elapsed_seconds']:.2f}s")
        
        export_cluster_results(output_dir, theta_val, cmin_val, pairs, clusters, pair_assignments)

    # 4. Save overall summary if sweep
    if args.sweep or len(grid) > 1:
        write_summary_tables(output_dir, sweep_metrics)

    print("\n=====================================================================")
    print("Clustering completed successfully!")
    print("=====================================================================")


if __name__ == "__main__":
    main()
