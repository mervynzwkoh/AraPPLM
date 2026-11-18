# Dataset Construction Documentation

This directory documents the methodology, verification standards, and operational guidelines for constructing the AraPPLM plant PPI datasets.

The dataset construction pipeline begins with audited, verified physical interaction pairs and progresses through negative sampling, pair-level homology clustering, and strict graph-theoretic benchmark partitioning.

---

## Roadmap & Pipeline Structure

```
docs/dataset_construction/
├── README.md                              # Master overview & pipeline roadmap (this file)
│
├── 01_positive_pairs/                     # [STAGE 1: COMPLETE] Verified positive interactomes
│   ├── 01_arabidopsis_positive_pairs.md   # Arabidopsis thaliana (120,873 deduplicated pairs; 30,376 high-confidence hetero)
│   ├── 02_maize_positive_pairs.md         # Zea mays (50,984 deduplicated pairs; 50,983 pure PPIs)
│   └── 03_rice_tomato_soybean_positive_pairs.md # Rice (998 / 667 PPI), Tomato (123), Soybean (48)
│
├── 02_negative_sampling/                  # [STAGE 2: COMPLETE] Negative candidate pool generation & filtering
│   └── 01_arabidopsis_negative_dataset_construction.md # 6-phase negative construction, CD-HIT clustering, & 1:10 partitioning
│
├── 03_homology_clustering/                # [STAGE 3: COMPLETE] Pair redundancy reduction (MMseqs2 & CD-HIT)
│   ├── 01_pair_redundancy_reduction.md    # Pair similarity formula & MMseqs2 inverted index methodology
│   ├── summary_clustering_metrics.md      # Final identity sweep report across theta in [0.40, 0.80] at cmin=0.80
│   ├── summary_clustering_metrics_04_07_09.md # Full 8-condition grid sweep across theta in {0.40, 0.70} x cmin in [0.50, 0.90]
│   └── top20_cluster_density_report_0.5.md # Subgraph density sweep across top 20 gene family clusters
│
└── 04_benchmark_splits/                   # [STAGE 4: COMPLETE] Finalized C1/C2/C3/C4 partitions (1:10 ratio)
    ├── 01_c1234_benchmark_partitions.md   # Authoritative master benchmark report (positives, negatives, & combined TSVs)
    ├── summary_positive_c1234_partitions.md # Positive partition audit (85% C1, 5% C2, 5% C3, 5% C4; 22,897 pairs)
    ├── summary_negative_c1234_partitions.md # Negative partition audit (10x ratio; 228,485 pairs; zero leakage)
    └── summary_positive_c123_partitions.md # Legacy parameter sweep partition report (8 permutations)
```

---

## 1. Stage 1: Verified Positive Pairs Ground Truth (`01_positive_pairs/`)

Foundational positive interactions extracted directly from raw primary database releases (BioGRID v5.0.261, IntAct r252, STRING v12.0, PRIDE PhoX XL-MS):
- Arabidopsis Ground Truth: **120,873 unique pairs** (115,012 heteromers, 5,861 homomers, 20,836 proteins).
- High-Confidence Heteromeric Positive Dataset: **30,376 pairs** ($S_{\text{HIPPIE}} \ge 0.72$) across **8,970 unique proteins**.

| Document | Target Species | Taxonomy ID | Total Deduplicated Pairs | Pure PPI Pairs | Key Data Sources |
| :--- | :--- | :--- | :--- | :--- | :--- |
| [`01_arabidopsis_positive_pairs.md`](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/AraPPLM/docs/dataset_construction/01_positive_pairs/01_arabidopsis_positive_pairs.md) | *Arabidopsis thaliana* | `3702` | **120,873** | **120,873** | BioGRID, IntAct, STRING Full, PhoX XL-MS |
| [`02_maize_positive_pairs.md`](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/AraPPLM/docs/dataset_construction/01_positive_pairs/02_maize_positive_pairs.md) | *Zea mays* | `4577` | **50,984** | **50,983** | IntAct, BioGRID |
| [`03_rice_tomato_soybean_positive_pairs.md`](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/AraPPLM/docs/dataset_construction/01_positive_pairs/03_rice_tomato_soybean_positive_pairs.md) | *Oryza sativa*<br>*Solanum lycopersicum*<br>*Glycine max* | `39947` / `4530`<br>`4081`<br>`3847` | **998** (Rice)<br>**123** (Tomato)<br>**48** (Soybean) | **667** (Rice)<br>**123** (Tomato)<br>**48** (Soybean) | BioGRID, IntAct, STRING Full |

---

## 2. Stage 2: Negative Candidate Pool Generation (`02_negative_sampling/`)

Full 6-phase screening pipeline documented in [`01_arabidopsis_negative_dataset_construction.md`](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/AraPPLM/docs/dataset_construction/02_negative_sampling/01_arabidopsis_negative_dataset_construction.md):
1. **Phase 1 (Pair Enumeration)**: 409,222,436 non-self pairs across 28,609 canonical loci.
2. **Phase 2 (Positive Blacklist Purge)**: Filtered against ~3.8M interactions across 6 databases.
3. **Phase 3 (Subcellular Compartment Filter)**: Disjoint subcellular compartments across 10 TAIR compartments $\to$ **101,884,330 candidate pairs** in 102 Parquet chunks.
4. **Phase 4 (Positive Interolog Filter, Stage 2.5)**: 32,788 positive interologs ($S \ge 0.50$ against 22,897 positive pairs) purged $\to$ **101,851,542 clean pairs** (99.97% retention).
5. **Phase 5 (Paralog Redundancy Reduction, Stage 2.6)**: Component decoupling & greedy CD-HIT selection $\to$ **50,308,730 non-redundant representative negatives** across 51 Parquet chunks.
6. **Phase 6 (Benchmark Partitioning, Stage 2.7)**: Option A quality sampling yielding **228,485 negatives** partitioned across C1 (194,620), C2 (11,450), C3 (10,965), and C4 (11,450).

---

## 3. Stage 3: Homology Clustering (`03_homology_clustering/`)

Eliminates pair-level sequence redundancy from the 30,376 high-confidence positive pairs:
- Methodology: Inverted homology indexing and greedy CD-HIT selection documented in [`01_pair_redundancy_reduction.md`](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/AraPPLM/docs/dataset_construction/03_homology_clustering/01_pair_redundancy_reduction.md).
- Final Benchmark Condition: $\theta = 0.50$, $c_{\min} = 0.80$, retaining **22,897 non-redundant positive pairs** across **8,426 unique proteins** (24.62% redundancy pruned).
- Comparative Parameter Sweeps: Documented in [`summary_clustering_metrics.md`](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/AraPPLM/docs/dataset_construction/03_homology_clustering/summary_clustering_metrics.md) and [`top20_cluster_density_report_0.5.md`](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/AraPPLM/docs/dataset_construction/03_homology_clustering/top20_cluster_density_report_0.5.md).

---

## 4. Stage 4: Benchmark Splits (`04_benchmark_splits/`)

Rigorous graph-theoretic partitioning into **C1 (Train, 85%)**, **C2 (One Unseen Test, 5%)**, **C3 (Both Unseen Test, 5%)**, and **C4 (Both Seen Test, 5%)** at a strict **1:10 positive-to-negative ratio**:
- Master Report: [`01_c1234_benchmark_partitions.md`](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/AraPPLM/docs/dataset_construction/04_benchmark_splits/01_c1234_benchmark_partitions.md).
- Total Combined Benchmark: **251,382 pairs** (22,897 Positives + 228,485 Negatives).
- Mathematical Guarantees: Strict vertex bisection ($V_{\text{train}} \cap V_{\text{test}} = \emptyset$), zero protein leakage ($V(C_3) \cap V(C_1) = \emptyset$), 100% C4 seen containment, and pairwise edge disjointness.
- Visualizations: Cytoscape interactive viewer (`c1234_interactive_network_viewer.html`) and vector network graphs (`c1234_network_id0_50_cov0_80.svg`).
