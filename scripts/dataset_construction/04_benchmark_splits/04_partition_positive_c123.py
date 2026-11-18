#!/usr/bin/env python3
"""
04_partition_positive_c123.py
=============================
Partitions non-redundant positive PPI pairs into C1/C2/C3/C4 benchmark splits:
- C1 (Train): 85% of pairs (~19,462 pairs)
- C2 (One Unseen Test): 5% of pairs (~1,145 pairs)
- C3 (Both Unseen Test): 5% of pairs (~1,145 pairs)
- C4 (Both Seen Test): 5% of pairs (~1,145 pairs)

Dataset:
--------
Input: data/dataset_construction/03_homology_clustering/arabidopsis/ara_positive_pairs_nr_id0_50_cov0_80.tsv
Total non-redundant positive pairs: 22,897

Benchmark Split Definitions (Park & Marcotte 2012 / Standard ML PPI Benchmarks):
--------------------------------------------------------------------------------
Let G = (V, E) be the positive interaction graph.
The vertex set is partitioned into disjoint sets V = V_train ∪ V_test (V_train ∩ V_test = ∅).

1. C1 (Train Set, 85%):
   Interactions where both interactors are in V_train: u ∈ V(C1) ⊆ V_train, v ∈ V(C1) ⊆ V_train.
   Used for training the PPI models.

2. C2 (One Unseen Test Set, 5%):
   Interactions crossing the partition boundary where exactly one interactor is in V_train
   (previously appeared in C1) and one interactor is completely unseen (in V_test).
   Evaluates model generalization when predicting partners for partially known proteins.

3. C3 (Both Unseen Test Set, 5%):
   Interactions where BOTH interactors are in V_test: u ∈ V_test, v ∈ V_test.
   Zero protein overlap with C1 train set: V(C3) ∩ V(C1) = ∅.
   Evaluates model generalization to completely novel, uncharacterized proteins.

4. C4 (Both Seen Test Set, 5%):
   Held-out test interactions where BOTH interactors have previously appeared in the C1 train set:
   u ∈ V(C1) and v ∈ V(C1), but the interaction pair (u, v) itself was NOT in C1 (edge-disjoint).
   Evaluates model capability to predict novel combinations among familiar interactors.

Strict Partition Verification Criteria:
----------------------------------------
1. Edge Partition Exhaustiveness: N_C1 + N_C2 + N_C3 + N_C4 == N_total (22,897)
2. Pairwise Edge Disjointness: C1 ∩ C2 = ∅, C1 ∩ C3 = ∅, C1 ∩ C4 = ∅, C2 ∩ C3 = ∅, C2 ∩ C4 = ∅, C3 ∩ C4 = ∅
3. Disjoint Vertex Sets: V_train ∩ V_test = ∅ and V_train ∪ V_test = V
4. C1 Vertex Containment: V(C1) ⊆ V_train
5. C4 Both Seen Requirement: For all (u, v) in C4: u ∈ V(C1) and v ∈ V(C1)
6. C3 Vertex Containment: V(C3) ⊆ V_test
7. Zero Protein Overlap: V(C1) ∩ V(C3) = ∅ and V(C4) ∩ V(C3) = ∅
8. C2 Boundary Cut: For all (u, v) in C2: exactly one node in V_train and one node in V_test
"""

import os
import sys
import csv
import time
import random
import argparse
from collections import defaultdict


def load_pairs(tsv_path):
    """Load non-redundant pairs with metadata."""
    pairs = []
    with open(tsv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            a = row["participant_a_id"].strip()
            b = row["participant_b_id"].strip()
            pairs.append({
                "participant_a_id": a,
                "participant_b_id": b,
                "cluster_id": row.get("cluster_id", ""),
                "hippie_score": row.get("hippie_score", ""),
                "evidence_count": row.get("evidence_count", "")
            })
    return pairs


def partition_graph_c1234(pairs, c1_ratio=0.85, c2_ratio=0.05, c3_ratio=0.05, c4_ratio=0.05,
                           seed=42, max_restarts=35, max_refine_steps=600):
    """
    Two-stage community-aware graph partitioning:
    Stage 1: Bisect V into V_train and V_test targeting |C2| = 5% and |C3| = 5%.
    Stage 2: Select 5% of edges from E[V_train] for C4 such that both endpoints retain degree >= 1 in C1 (85%).
    """
    nodes = list(set([p["participant_a_id"] for p in pairs] + [p["participant_b_id"] for p in pairs]))
    adj = defaultdict(set)
    for p in pairs:
        a = p["participant_a_id"]
        b = p["participant_b_id"]
        adj[a].add(b)
        adj[b].add(a)

    total_edges = len(pairs)
    target_c2 = int(round(c2_ratio * total_edges))
    target_c3 = int(round(c3_ratio * total_edges))
    target_c4 = int(round(c4_ratio * total_edges))
    target_c1 = total_edges - target_c2 - target_c3 - target_c4

    print(f"Target Edge Counts: C1={target_c1} ({c1_ratio*100:.1f}%), "
          f"C2={target_c2} ({c2_ratio*100:.1f}%), "
          f"C3={target_c3} ({c3_ratio*100:.1f}%), "
          f"C4={target_c4} ({c4_ratio*100:.1f}%)")

    def evaluate_loss(c2, c3):
        return (c2 - target_c2) ** 2 + (c3 - target_c3) ** 2

    rng = random.Random(seed)
    candidates_seed = [n for n in nodes if len(adj[n]) >= 6]
    if not candidates_seed:
        candidates_seed = nodes

    best_vtest = None
    best_loss = float("inf")
    best_counts = None

    # Stage 1: Bisection into V_train and V_test
    print("Stage 1: Performing community-aware graph bisection into V_train and V_test...")
    for trial in range(max_restarts):
        seed_node = rng.choice(candidates_seed)
        v_test = {seed_node}
        unassigned = set(nodes) - {seed_node}
        vtest_nbrs = set(adj[seed_node])

        c3 = 0
        c2 = len(adj[seed_node])

        while c3 < target_c3 and vtest_nbrs:
            best_u = None
            best_score = -999999
            cands = list(vtest_nbrs)
            if len(cands) > 50:
                cands = rng.sample(cands, 50)
            for u in cands:
                internal = len(adj[u].intersection(v_test))
                external = len(adj[u]) - internal
                score = internal * 2.8 - external
                if score > best_score:
                    best_score = score
                    best_u = u
            if not best_u:
                break

            v_test.add(best_u)
            unassigned.remove(best_u)
            vtest_nbrs.remove(best_u)
            for nbr in adj[best_u]:
                if nbr in unassigned:
                    vtest_nbrs.add(nbr)

            internal = len(adj[best_u].intersection(v_test - {best_u}))
            external = len(adj[best_u]) - internal
            c3 += internal
            c2 += (external - internal)

        # Boundary refinement
        for step in range(max_refine_steps):
            best_move = None
            curr_loss = evaluate_loss(c2, c3)
            iter_best_loss = curr_loss

            # Move from test to train
            b_test = [u for u in v_test if len(adj[u].intersection(unassigned)) > 0]
            for u in rng.sample(b_test, min(len(b_test), 35)):
                internal = len(adj[u].intersection(v_test))
                external = len(adj[u]) - internal
                new_c3 = c3 - internal
                new_c2 = c2 + internal - external
                loss = evaluate_loss(new_c2, new_c3)
                if loss < iter_best_loss:
                    iter_best_loss = loss
                    best_move = (u, "to_train", new_c2, new_c3)

            # Move from train to test
            b_train = [u for u in unassigned if len(adj[u].intersection(v_test)) > 0]
            for u in rng.sample(b_train, min(len(b_train), 35)):
                internal = len(adj[u].intersection(v_test))
                external = len(adj[u]) - internal
                new_c3 = c3 + internal
                new_c2 = c2 - internal + external
                loss = evaluate_loss(new_c2, new_c3)
                if loss < iter_best_loss:
                    iter_best_loss = loss
                    best_move = (u, "to_test", new_c2, new_c3)

            if best_move and iter_best_loss < curr_loss:
                u, direction, c2, c3 = best_move
                if direction == "to_train":
                    v_test.remove(u)
                    unassigned.add(u)
                else:
                    v_test.add(u)
                    unassigned.remove(u)
            else:
                break

        final_loss = evaluate_loss(c2, c3)
        if final_loss < best_loss:
            best_loss = final_loss
            best_vtest = set(v_test)
            best_counts = (c2, c3)
            if best_loss == 0:
                print(f"  -> Found optimal bisection at trial {trial}: C2={c2}, C3={c3} (loss=0)")
                break

    v_test = best_vtest
    v_train = set(nodes) - v_test
    print(f"Bisection complete: |V_train| = {len(v_train):,}, |V_test| = {len(v_test):,}")

    c2_pairs = []
    c3_pairs = []
    train_edges = []

    for p in pairs:
        a = p["participant_a_id"]
        b = p["participant_b_id"]
        if a in v_train and b in v_train:
            train_edges.append(p)
        elif a in v_test and b in v_test:
            c3_pairs.append(p)
        else:
            c2_pairs.append(p)

    # Stage 2: Partition train_edges into C1 (85%) and C4 (5%)
    print(f"Stage 2: Partitioning {len(train_edges):,} train-induced edges into C1 (85%) and C4 (5%)...")
    deg_in_train = defaultdict(int)
    for p in train_edges:
        deg_in_train[p["participant_a_id"]] += 1
        deg_in_train[p["participant_b_id"]] += 1

    rng_c4 = random.Random(seed)
    shuffled_train = list(train_edges)
    rng_c4.shuffle(shuffled_train)

    c4_pairs = []
    c1_deg = dict(deg_in_train)

    for p in shuffled_train:
        if len(c4_pairs) >= target_c4:
            break
        u = p["participant_a_id"]
        v = p["participant_b_id"]
        # Enforce both endpoints retain degree >= 1 in C1
        if c1_deg[u] > 1 and c1_deg[v] > 1:
            c4_pairs.append(p)
            c1_deg[u] -= 1
            c1_deg[v] -= 1

    c4_id_set = set(id(p) for p in c4_pairs)
    c1_pairs = [p for p in train_edges if id(p) not in c4_id_set]

    return c1_pairs, c2_pairs, c3_pairs, c4_pairs, v_train, v_test


def verify_partition_c1234(pairs, c1_pairs, c2_pairs, c3_pairs, c4_pairs, v_train, v_test):
    """
    Comprehensive verification of all C1/C2/C3/C4 criteria:
    """
    total_edges = len(pairs)
    nodes = set(v_train).union(set(v_test))

    # 1. Edge Exhaustiveness
    n_sum = len(c1_pairs) + len(c2_pairs) + len(c3_pairs) + len(c4_pairs)
    assert n_sum == total_edges, f"Edge count mismatch: sum={n_sum} vs total={total_edges}"

    # 2. Canonical Edge Disjointness
    set_c1 = set((min(p["participant_a_id"], p["participant_b_id"]), max(p["participant_a_id"], p["participant_b_id"])) for p in c1_pairs)
    set_c2 = set((min(p["participant_a_id"], p["participant_b_id"]), max(p["participant_a_id"], p["participant_b_id"])) for p in c2_pairs)
    set_c3 = set((min(p["participant_a_id"], p["participant_b_id"]), max(p["participant_a_id"], p["participant_b_id"])) for p in c3_pairs)
    set_c4 = set((min(p["participant_a_id"], p["participant_b_id"]), max(p["participant_a_id"], p["participant_b_id"])) for p in c4_pairs)

    assert len(set_c1.intersection(set_c2)) == 0, "C1 and C2 share pairs!"
    assert len(set_c1.intersection(set_c3)) == 0, "C1 and C3 share pairs!"
    assert len(set_c1.intersection(set_c4)) == 0, "C1 and C4 share pairs!"
    assert len(set_c2.intersection(set_c3)) == 0, "C2 and C3 share pairs!"
    assert len(set_c2.intersection(set_c4)) == 0, "C2 and C4 share pairs!"
    assert len(set_c3.intersection(set_c4)) == 0, "C3 and C4 share pairs!"

    # 3. Disjoint Vertex Sets
    assert len(v_train.intersection(v_test)) == 0, "V_train and V_test are not disjoint!"

    # Extract unique proteins in each partition
    c1_prots = set([p["participant_a_id"] for p in c1_pairs] + [p["participant_b_id"] for p in c1_pairs])
    c2_prots = set([p["participant_a_id"] for p in c2_pairs] + [p["participant_b_id"] for p in c2_pairs])
    c3_prots = set([p["participant_a_id"] for p in c3_pairs] + [p["participant_b_id"] for p in c3_pairs])
    c4_prots = set([p["participant_a_id"] for p in c4_pairs] + [p["participant_b_id"] for p in c4_pairs])

    # 4. C1 Vertex Membership
    assert c1_prots.issubset(v_train), "C1 contains proteins outside V_train!"

    # 5. C4 Both Seen Requirement: every protein in C4 must have appeared in C1 train set
    assert c4_prots.issubset(c1_prots), f"Found {len(c4_prots - c1_prots)} proteins in C4 that did not appear in C1!"

    # 6. C3 Vertex Membership
    assert c3_prots.issubset(v_test), "C3 contains proteins outside V_test!"

    # 7. Zero Protein Overlap with C3
    c1_c3_overlap = c1_prots.intersection(c3_prots)
    assert len(c1_c3_overlap) == 0, f"Found {len(c1_c3_overlap)} overlapping proteins between C1 and C3!"
    c4_c3_overlap = c4_prots.intersection(c3_prots)
    assert len(c4_c3_overlap) == 0, f"Found {len(c4_c3_overlap)} overlapping proteins between C4 and C3!"

    # 8. C2 Cut Membership
    for p in c2_pairs:
        a, b = p["participant_a_id"], p["participant_b_id"]
        assert (a in v_train and b in v_test) or (a in v_test and b in v_train), f"Invalid C2 pair: ({a}, {b})"

    return {
        "status": "PASS",
        "c1_pairs": len(c1_pairs),
        "c2_pairs": len(c2_pairs),
        "c3_pairs": len(c3_pairs),
        "c4_pairs": len(c4_pairs),
        "total_pairs": total_edges,
        "c1_pct": len(c1_pairs) / total_edges * 100.0,
        "c2_pct": len(c2_pairs) / total_edges * 100.0,
        "c3_pct": len(c3_pairs) / total_edges * 100.0,
        "c4_pct": len(c4_pairs) / total_edges * 100.0,
        "c1_prots": len(c1_prots),
        "c2_prots": len(c2_prots),
        "c3_prots": len(c3_prots),
        "c4_prots": len(c4_prots),
        "total_prots": len(nodes),
        "v_train_size": len(v_train),
        "v_test_size": len(v_test),
        "c2_train_prots": len(c2_prots.intersection(v_train)),
        "c2_test_prots": len(c2_prots.intersection(v_test)),
    }


def write_partition_tsv(pairs, out_path):
    """Write partitioned pairs to TSV with schema."""
    fieldnames = ["participant_a_id", "participant_b_id", "cluster_id", "hippie_score", "evidence_count"]
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for p in pairs:
            writer.writerow({k: p.get(k, "") for k in fieldnames})
    print(f"Saved: {out_path} ({len(pairs):,} pairs)")


def write_summary_reports(output_dir, metrics, theta_str, cmin_str):
    """Write comprehensive summary reports in Markdown and TSV."""
    summary_tsv = os.path.join(output_dir, "summary_positive_c123_partitions.tsv")
    summary_tsv_c1234 = os.path.join(output_dir, "summary_positive_c1234_partitions.tsv")
    summary_md = os.path.join(output_dir, "summary_positive_c123_partitions.md")
    summary_md_c1234 = os.path.join(output_dir, "summary_positive_c1234_partitions.md")

    # TSV Header and Row
    headers = [
        "Identity_Theta", "Coverage_Cmin", "Total_Pairs",
        "C1_Train_Pairs", "C1_Train_Pct", "C1_Train_Prots",
        "C2_Test_Pairs", "C2_Test_Pct", "C2_Test_Prots", "C2_Train_Side_Prots", "C2_Test_Side_Prots",
        "C3_Test_Pairs", "C3_Test_Pct", "C3_Test_Prots",
        "C4_Test_Pairs", "C4_Test_Pct", "C4_Test_Prots",
        "Total_Unique_Prots", "V_Train_Prots", "V_Test_Prots", "Verification_Status"
    ]
    row = [
        theta_str.replace("_", "."), cmin_str.replace("_", "."), metrics["total_pairs"],
        metrics["c1_pairs"], f"{metrics['c1_pct']:.2f}%", metrics["c1_prots"],
        metrics["c2_pairs"], f"{metrics['c2_pct']:.2f}%", metrics["c2_prots"], metrics["c2_train_prots"], metrics["c2_test_prots"],
        metrics["c3_pairs"], f"{metrics['c3_pct']:.2f}%", metrics["c3_prots"],
        metrics["c4_pairs"], f"{metrics['c4_pct']:.2f}%", metrics["c4_prots"],
        metrics["total_prots"], metrics["v_train_size"], metrics["v_test_size"], metrics["status"]
    ]

    for p in [summary_tsv, summary_tsv_c1234]:
        with open(p, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(headers)
            w.writerow(row)

    # Markdown Report
    theta_disp = theta_str.replace("_", ".")
    cmin_disp = cmin_str.replace("_", ".")
    vtrain_pct = metrics['v_train_size'] / metrics['total_prots'] * 100.0
    vtest_pct = metrics['v_test_size'] / metrics['total_prots'] * 100.0

    md_content = """# Arabidopsis Positive PPI Benchmark C1/C2/C3/C4 Partition Summary

This report documents the rigorous graph-theoretic benchmark partitioning of the high-confidence, non-redundant *Arabidopsis thaliana* positive protein-protein interactions (at sequence identity $\\theta = {theta_disp}$ and alignment coverage $c_{{\\min}} = {cmin_disp}$) into **C1 (Train, 85%)**, **C2 (One Unseen Test, 5%)**, **C3 (Both Unseen Test, 5%)**, and **C4 (Both Seen Test, 5%)**.

---

## 1. Partition Statistics Table

| Partition | Role | Number of Pairs | % of Total Pairs | Unique Proteins | Vertex Membership Property |
|:---:|:---:|:---:|:---:|:---:|:---|
| **C1** | **Training Set** | **{c1_pairs:,}** | **{c1_pct:.2f}%** | **{c1_prots:,}** | $u, v \\in V_{{\\text{{train}}}}$ (Both interactors seen in train) |
| **C2** | **One Unseen Test** | **{c2_pairs:,}** | **{c2_pct:.2f}%** | **{c2_prots:,}** | $u \\in V_{{\\text{{train}}}}, v \\in V_{{\\text{{test}}}}$ (One seen, one unseen) |
| **C3** | **Both Unseen Test** | **{c3_pairs:,}** | **{c3_pct:.2f}%** | **{c3_prots:,}** | $u, v \\in V_{{\\text{{test}}}}$ (Both interactors completely unseen; $V(C_3) \\cap V(C_1) = \\emptyset$) |
| **C4** | **Both Seen Test** | **{c4_pairs:,}** | **{c4_pct:.2f}%** | **{c4_prots:,}** | $u, v \\in V(C_1)$ (Both interactors seen in C1; novel test interaction pair) |
| **Total** | **Full NR Dataset** | **{total_pairs:,}** | **100.00%** | **{total_prots:,}** | Disjoint graph union $V = V_{{\\text{{train}}}} \\cup V_{{\\text{{test}}}}$ |

---

## 2. Vertex & Graph Partition Properties

- **Total Unique Proteins**: **{total_prots:,}**
- **Disjoint Train Vertex Set ($|V_{{\\text{{train}}}}|$)**: **{v_train_size:,}** proteins ({vtrain_pct:.1f}%)
- **Disjoint Test Vertex Set ($|V_{{\\text{{test}}}}|$)**: **{v_test_size:,}** proteins ({vtest_pct:.1f}%)
- **C2 Cross-Boundary Node Distribution**: {c2_train_prots:,} proteins from $V_{{\\text{{train}}}}$ and {c2_test_prots:,} proteins from $V_{{\\text{{test}}}}$
- **C4 Seen Node Coverage**: **100.0%** of proteins in C4 ({c4_prots:,} / {c4_prots:,}) strictly appear in the C1 training set.
- **C3 Unseen Node Isolation**: **0** proteins in C3 overlap with C1 or C4 ($V(C_3) \\cap V(C_1) = \\emptyset$).

---

## 3. Strict Verification & Integrity Audit

| Verification Check | Target Requirement | Evaluated Result | Status |
|:---|:---|:---:|:---:|
| **1. Edge Exhaustiveness** | $|C_1| + |C_2| + |C_3| + |C_4| == N_{{\\text{{total}}}}$ | {c1_pairs:,} + {c2_pairs:,} + {c3_pairs:,} + {c4_pairs:,} = {total_pairs:,} | **PASS** |
| **2. Pairwise Edge Disjointness** | Zero shared interaction pairs across all 4 sets | $\\max_{{i \\neq j}} |C_i \\cap C_j| = 0$ | **PASS** |
| **3. Vertex Disjointness** | $V_{{\\text{{train}}}} \\cap V_{{\\text{{test}}}} = \\emptyset$ | $|V_{{\\text{{train}}}} \\cap V_{{\\text{{test}}}}| = 0$ | **PASS** |
| **4. C1 Vertex Containment** | $V(C_1) \\subseteq V_{{\\text{{train}}}}$ | All {c1_prots:,} nodes $\\in V_{{\\text{{train}}}}$ | **PASS** |
| **5. C4 Both Seen Requirement** | $\\forall (u, v) \\in C_4: u, v \\in V(C_1)$ | All {c4_prots:,} nodes $\\in V(C_1)$ (100%) | **PASS** |
| **6. C3 Vertex Containment** | $V(C_3) \\subseteq V_{{\\text{{test}}}}$ | All {c3_prots:,} nodes $\\in V_{{\\text{{test}}}}$ | **PASS** |
| **7. C1 / C3 Zero Protein Overlap** | $V(C_1) \\cap V(C_3) = \\emptyset$ | $|V(C_1) \\cap V(C_3)| = 0$ | **PASS** |
| **8. C4 / C3 Zero Protein Overlap** | $V(C_4) \\cap V(C_3) = \\emptyset$ | $|V(C_4) \\cap V(C_3)| = 0$ | **PASS** |
| **9. C2 Cut Boundary Property** | Exactly one protein $\\in V_{{\\text{{train}}}}$, one $\\in V_{{\\text{{test}}}}$ | 100% of C2 pairs cross cut boundary | **PASS** |

---

## 4. Output Partition Files

All partition files have been written to `data/dataset_construction/04_benchmark_splits/arabidopsis/`:
- `ara_positive_c1_train_id{theta_str}_cov{cmin_str}.tsv` ({c1_pairs:,} pairs)
- `ara_positive_c2_test_id{theta_str}_cov{cmin_str}.tsv` ({c2_pairs:,} pairs)
- `ara_positive_c3_test_id{theta_str}_cov{cmin_str}.tsv` ({c3_pairs:,} pairs)
- `ara_positive_c4_test_id{theta_str}_cov{cmin_str}.tsv` ({c4_pairs:,} pairs)
""".format(
        theta_disp=theta_disp,
        cmin_disp=cmin_disp,
        theta_str=theta_str,
        cmin_str=cmin_str,
        vtrain_pct=vtrain_pct,
        vtest_pct=vtest_pct,
        **metrics
    )

    for p in [summary_md, summary_md_c1234]:
        with open(p, "w", encoding="utf-8") as f:
            f.write(md_content)

    print(f"Saved summary TSV: {summary_tsv}")
    print(f"Saved summary MD:  {summary_md}")


def main():
    parser = argparse.ArgumentParser(description="Partition positive PPIs into C1/C2/C3/C4 splits.")
    parser.add_argument("--input_tsv", type=str,
                        default="data/dataset_construction/03_homology_clustering/arabidopsis/ara_positive_pairs_nr_id0_50_cov0_80.tsv",
                        help="Input non-redundant positive TSV file.")
    parser.add_argument("--output_dir", type=str,
                        default="data/dataset_construction/04_benchmark_splits/arabidopsis",
                        help="Output directory for partitioned TSVs and summary reports.")
    parser.add_argument("--c1_pct", type=float, default=0.85, help="Target C1 fraction (default: 0.85)")
    parser.add_argument("--c2_pct", type=float, default=0.05, help="Target C2 fraction (default: 0.05)")
    parser.add_argument("--c3_pct", type=float, default=0.05, help="Target C3 fraction (default: 0.05)")
    parser.add_argument("--c4_pct", type=float, default=0.05, help="Target C4 fraction (default: 0.05)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility.")

    args = parser.parse_args()

    # Dynamic path resolution
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    def resolve_path(p):
        if not os.path.isabs(p) and not os.path.exists(p):
            alt = os.path.join(repo_root, p)
            if os.path.exists(alt):
                return alt
        return p

    input_tsv = resolve_path(args.input_tsv)
    output_dir = resolve_path(args.output_dir)

    print("=" * 80)
    print("    AraPPLM: Positive PPI C1/C2/C3/C4 Benchmark Partitioning")
    print("=" * 80)
    print(f"Input TSV:        {input_tsv}")
    print(f"Output Directory: {output_dir}")
    print(f"Target Splits:    85% C1 (Train), 5% C2 (One Unseen), 5% C3 (Both Unseen), 5% C4 (Both Seen)")
    print("=" * 80)

    if not os.path.isfile(input_tsv):
        print(f"ERROR: Input TSV file not found: {input_tsv}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)

    # Infer theta and cmin from filename
    fname = os.path.basename(input_tsv)
    theta_str = "0_50"
    cmin_str = "0_80"
    if "id" in fname and "cov" in fname:
        parts = fname.replace(".tsv", "").split("_")
        for i, part in enumerate(parts):
            if part.startswith("id"):
                theta_str = parts[i][2:] + "_" + parts[i+1]
            elif part.startswith("cov"):
                cmin_str = parts[i][3:] + "_" + parts[i+1]

    # 1. Load pairs
    print(f"Loading pairs from: {input_tsv} ...")
    pairs = load_pairs(input_tsv)
    print(f"Loaded {len(pairs):,} non-redundant pairs.")

    # 2. Partition graph
    c1, c2, c3, c4, v_train, v_test = partition_graph_c1234(
        pairs,
        c1_ratio=args.c1_pct,
        c2_ratio=args.c2_pct,
        c3_ratio=args.c3_pct,
        c4_ratio=args.c4_pct,
        seed=args.seed
    )

    # 3. Verify partitions
    print("\nAuditing and verifying partition integrity...")
    metrics = verify_partition_c1234(pairs, c1, c2, c3, c4, v_train, v_test)
    print(f"Verification Result: {metrics['status']}")
    print(f"  - C1 (Train):            {metrics['c1_pairs']:,} pairs ({metrics['c1_pct']:.2f}%), {metrics['c1_prots']:,} unique proteins")
    print(f"  - C2 (One Unseen Test):  {metrics['c2_pairs']:,} pairs ({metrics['c2_pct']:.2f}%), {metrics['c2_prots']:,} unique proteins")
    print(f"  - C3 (Both Unseen Test): {metrics['c3_pairs']:,} pairs ({metrics['c3_pct']:.2f}%), {metrics['c3_prots']:,} unique proteins")
    print(f"  - C4 (Both Seen Test):   {metrics['c4_pairs']:,} pairs ({metrics['c4_pct']:.2f}%), {metrics['c4_prots']:,} unique proteins")
    print(f"  - Total Evaluated:       {metrics['total_pairs']:,} pairs (100.00%), {metrics['total_prots']:,} unique proteins")

    # 4. Write partition TSV files
    print("\nExporting partition TSVs...")
    write_partition_tsv(c1, os.path.join(output_dir, f"ara_positive_c1_train_id{theta_str}_cov{cmin_str}.tsv"))
    write_partition_tsv(c2, os.path.join(output_dir, f"ara_positive_c2_test_id{theta_str}_cov{cmin_str}.tsv"))
    write_partition_tsv(c3, os.path.join(output_dir, f"ara_positive_c3_test_id{theta_str}_cov{cmin_str}.tsv"))
    write_partition_tsv(c4, os.path.join(output_dir, f"ara_positive_c4_test_id{theta_str}_cov{cmin_str}.tsv"))

    # 5. Write summary reports
    print("\nGenerating summary reports...")
    write_summary_reports(output_dir, metrics, theta_str, cmin_str)

    print("\n" + "=" * 80)
    print("C1/C2/C3/C4 Partitioning Pipeline Completed Successfully!")
    print("=" * 80)


if __name__ == "__main__":
    main()
