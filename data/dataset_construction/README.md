# Dataset Construction Data Warehouse

This directory houses all progressive artifacts of the AraPPLM plant PPI dataset construction pipeline. The pipeline is structured into four sequential stages, starting from audited raw database interactions and culminating in strict graph-theoretic benchmark splits.

---

## Pipeline Stages

```
data/dataset_construction/
├── README.md                      # Pipeline data architecture & status (this file)
│
├── 01_positive_pairs/             # [STAGE 1: COMPLETE] Master verified physical PPI pairs
│   ├── arabidopsis/               # 120,873 unique pairs (115,012 heteromers; 30,376 high-confidence hetero)
│   ├── maize/                     # 50,984 unique pairs (50,983 pure PPIs)
│   ├── rice/                      # 998 unique pairs (667 pure PPI pairs)
│   ├── tomato/                    # 123 unique pairs (123 pure PPI pairs)
│   └── soybean/                   # 48 unique pairs (48 pure PPI pairs)
│
├── 02_negative_sampling/          # [STAGE 2: COMPLETE] Candidate negative pool generation & filtering
│   ├── ara_all_databases_positive_blacklist.parquet # ~3.8M positive blacklist
│   ├── candidates/                # 101.88M candidate pairs across 10 TAIR compartments (102 chunks)
│   ├── filtered_candidates/       # 101.85M clean pairs free of positive interologs (102 chunks)
│   ├── final_candidates_nr/       # 50.31M non-redundant negative candidate pool (51 chunks)
│   ├── summary_candidate_negatives_filtered.md # Stage 2.5 interolog audit report
│   └── summary_candidate_negatives_nr.md       # Stage 2.6 redundancy reduction audit report
│
├── 03_homology_clustering/        # [STAGE 3: COMPLETE] Sequence & pair redundancy reduction
│   └── arabidopsis/
│       ├── ara_positive_pairs_nr_id0_50_cov0_80.tsv # 22,897 non-redundant positive pairs
│       ├── summary_clustering_metrics.md           # Parameter sweep report across theta in [0.40, 0.80]
│       └── top20_cluster_density_report_0.5.md     # Cluster density analysis for theta=0.50
│
└── 04_benchmark_splits/           # [STAGE 4: COMPLETE] Strict C1/C2/C3/C4 partitions (1:10 ratio)
    ├── arabidopsis/
    │   ├── positive_data/         # Positive partition TSVs & summaries (C1: 19,462; C2: 1,145; C3: 1,145; C4: 1,145)
    │   ├── negative_data/         # Negative partition TSVs, audit JSON & summary (C1: 194,620; C2: 11,450; C3: 10,965; C4: 11,450)
    │   ├── benchmark_data/        # Combined full benchmark TSVs (Pos=1, Neg=0; Total: 251,382 pairs)
    │   └── README.md              # Organization guide & file catalog
    └── figures/
        ├── c1234_network_id0_50_cov0_80.svg     # SVG network diagram
        └── c1234_interactive_network_viewer.html # Cytoscape interactive graph visualization
```

---

## 1. Stage 1: Verified Positive Pairs (`01_positive_pairs/`)

Foundational positive interaction ground truth extracted from BioGRID v5.0.261, IntAct r252, STRING v12.0, and PRIDE PhoX XL-MS (2026).
- Arabidopsis Ground Truth: **120,873 unique pairs** (115,012 heteromers, 5,861 homomers, 20,836 proteins).
- High-Confidence Heteromeric Positive Dataset: **30,376 pairs** with continuous HIPPIE score $\ge 0.72$ across **8,970 unique proteins**.

---

## 2. Stage 2: Negative Candidate Pool (`02_negative_sampling/`)

A 5-phase screening pipeline generating non-interacting controls:
1. **Enumeration**: 409,222,436 non-self pairs across 28,609 canonical loci.
2. **Positive Blacklist Purge**: Filtered against ~3.8M known interactions across 6 databases.
3. **Subcellular Co-Localization Filter**: 10 distinct TAIR compartments; pairs sharing any compartment discarded $\to$ **101,884,330 candidate pairs** in 102 chunks.
4. **Interolog Homology Filter (Stage 2.5)**: 32,788 positive interologs ($S \ge 0.50$ against 22,897 confirmed positive interactions) purged $\to$ **101,851,542 clean pairs** (99.97% retention).
5. **Pair-Level Redundancy Reduction (Stage 2.6)**: Component decoupling and greedy CD-HIT representative selection $\to$ **50,308,730 non-redundant candidate negatives** across 51 Parquet chunks.

---

## 3. Stage 3: Homology Clustering (`03_homology_clustering/`)

Pair-level sequence redundancy reduction evaluated via MMseqs2 ($S(P_i, P_j) \ge 0.50, c_{\min} \ge 0.80$):
- Starting Input: 30,376 high-confidence heteromeric pairs.
- Final Selected Threshold: $\theta = 0.50$, $c_{\min} = 0.80$.
- Retained Representatives: **22,897 non-redundant positive pairs** (24.62% redundancy pruned) across **8,426 unique proteins**.

---

## 4. Stage 4: Benchmark Partitions (`04_benchmark_splits/`)

Rigorous graph-theoretic partitioning into **C1 (Train, 85%)**, **C2 (One Unseen Test, 5%)**, **C3 (Both Unseen Test, 5%)**, and **C4 (Both Seen Test, 5%)** at a strict **1:10 positive-to-negative ratio**:

| Partition | Role | Positive Pairs | Negative Pairs (10x) | Combined Pairs | Pos:Neg Ratio | Mathematical Property |
| :---: | :---: | :---: | :---: | :---: | :---: | :---|
| **C1** | **Training Set** | **19,462** | **194,620** | **214,082** | **1 : 10.00** | $u, v \in V(C_1) \subseteq V_{\text{train}}$ |
| **C2** | **One Unseen Test** | **1,145** | **11,450** | **12,595** | **1 : 10.00** | $u \in V_{\text{train}}, v \in V_{\text{test}}$ (Cut boundary) |
| **C3** | **Both Unseen Test** | **1,145** | **10,965** | **12,110** | **1 : 9.58** | $u, v \in V_{\text{test}}$ ($V(C_3) \cap V(C_1) = \emptyset$, zero leakage) |
| **C4** | **Both Seen Test** | **1,145** | **11,450** | **12,595** | **1 : 10.00** | $u, v \in V(C_1)$ (Novel test pairs among seen proteins) |
| **Total** | **Full Benchmark** | **22,897** | **228,485** | **251,382** | **1 : 9.98** | Disjoint graph union ($V_{\text{train}} \cap V_{\text{test}} = \emptyset$) |
