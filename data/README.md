# AraPPLM Data Directory Architecture

This directory houses all data assets for the AraPPLM project, organized to maintain a strict separation between raw downloads, verified positive interaction pairs, negative sampling pools, non-redundant clusters, benchmark partitions, and legacy evaluation datasets.

---

## Directory Structure Overview

```
data/
├── dataset_construction/          # [DATASET CONSTRUCTION] Complete 4-stage pipeline artifacts
│   ├── README.md                  # Dataset construction pipeline data documentation
│   ├── 01_positive_pairs/         # [Stage 1] Master verified physical positive interactions
│   │   ├── arabidopsis/           # 120,873 deduplicated pairs (100% pure PPI; 30,376 high-confidence hetero)
│   │   ├── maize/                 # 50,984 deduplicated pairs (50,983 pure PPIs)
│   │   ├── rice/                  # 998 unique pairs (667 pure PPI pairs)
│   │   ├── soybean/               # 48 deduplicated pairs
│   │   └── tomato/                # 123 deduplicated pairs
│   │
│   ├── 02_negative_sampling/      # [Stage 2] Negative sampling pipeline outputs
│   │   ├── ara_all_databases_positive_blacklist.parquet # Positive blacklist (~3.8M pairs)
│   │   ├── candidates/            # Subcellular disjoint candidate pairs (101.88M pairs, 102 chunks)
│   │   ├── filtered_candidates/   # Positive interolog-filtered pairs (101.85M pairs, 102 chunks)
│   │   └── final_candidates_nr/   # Non-redundant negative candidate pool (50.31M pairs, 51 chunks)
│   │
│   ├── 03_homology_clustering/    # [Stage 3] Pair-level CD-HIT redundancy clustering (theta=0.50, cmin=0.80)
│   │   └── arabidopsis/           # 22,897 non-redundant positive pairs & parameter sweep reports
│   │
│   └── 04_benchmark_splits/       # [Stage 4] Finalized C1/C2/C3/C4 benchmark partitions (1:10 ratio)
│       ├── arabidopsis/           # Positive, negative, & combined TSVs + audit reports
│       └── figures/               # Cytoscape interactive HTML and SVG network graphs
│
├── legacy_benchmarks/             # [LEGACY BENCHMARKS] Historical datasets & missing ID rescue assets
│   ├── DeepAraPPI/                # DeepAraPPI C1 (train), C2 (test), C3 (test), Rice test
│   ├── ESMAraPPI/                 # ESMAraPPI C1 (train), C2 (pred), C3 (pred)
│   ├── arabidopsis_missing/       # Rescued UniParc sequences & UniProt Arabidopsis reference stream
│   └── rice_missing/              # Rescued UniParc sequences for rice accessions
│
└── raw/                           # [IMMUTABLE] Raw, unmodified primary database downloads
    ├── arabidopsis/               # BioGRID v5.0.261, IntAct r252, STRING v12.0, PRIDE PhoX XL-MS, TAIR, SUBA5
    ├── rice/                      # BioGRID, IntAct, STRING, PoPPIN (PMC XML)
    ├── external_holdout/          # Maize, Tomato, and Soybean raw downloads
    └── legacy_benchmarks/         # Raw source benchmark archives from prior literature
```

---

## 1. Verified Positive Pairs (`data/dataset_construction/01_positive_pairs/`)

This directory is the **authoritative starting point** for the dataset construction pipeline. All records have been independently verified against raw primary database dumps with PCA strictly retained across all sources.

### Master Ground-Truth Summary

| Species | Taxonomy ID | Total Deduplicated Pairs | Pure PPI Pairs | Heteromers | Homomers | Unique Proteins | Sources Included |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Arabidopsis** | `3702` | **120,873** | **120,873** | 115,012 | 5,861 | 20,836 | BioGRID, IntAct, STRING Full, PhoX XL-MS |
| **Maize** | `4577` | **50,984** | **50,983** | 50,409 | 575 | 14,888 | IntAct, BioGRID |
| **Rice** | `39947` / `4530` | **998** | **667** | 651 | 16 | 609 *(892 all)* | BioGRID, IntAct, STRING Full |
| **Tomato** | `4081` | **123** | **123** | 108 | 15 | 61 | BioGRID, IntAct |
| **Soybean** | `3847` | **48** | **48** | 47 | 1 | 51 | BioGRID, IntAct |

---

## 2. Negative Candidate Pool (`data/dataset_construction/02_negative_sampling/`)

Houses the multi-stage negative sampling artifacts for *Arabidopsis thaliana*:
- **Positive Blacklist**: `ara_all_databases_positive_blacklist.parquet` (~3.8M pairs aggregated across all databases).
- **Candidates**: 101,884,330 candidate pairs across 10 distinct subcellular compartments.
- **Filtered Candidates**: 101,851,542 clean candidate pairs with zero interolog homology ($S < 0.50$) to confirmed positive interactions.
- **Final Non-Redundant Candidates**: 50,308,730 representative pairs across 51 Parquet chunks satisfying zero pair-level redundancy.

---

## 3. Homology-Clustered Positives (`data/dataset_construction/03_homology_clustering/`)

Contains non-redundant positive interaction pairs and parameter sweep reports:
- Final representative dataset at $\theta = 0.50, c_{\min} = 0.80$: **22,897 non-redundant positive pairs** across **8,426 unique proteins**.
- Detailed density and parameter sweep audits across identity ($\theta \in [0.40, 0.80]$) and coverage ($c_{\min} \in [0.50, 0.90]$).

---

## 4. Benchmark Partitions (`data/dataset_construction/04_benchmark_splits/`)

Contains the finalized C1/C2/C3/C4 benchmark splits organized into three dedicated datasets under `arabidopsis/`:
- **Positive Data (`positive_data/`)**: C1 (19,462), C2 (1,145), C3 (1,145), C4 (1,145) positive interactome partitions with HIPPIE confidence scores.
- **Negative Data (`negative_data/`)**: C1 (194,620), C2 (11,450), C3 (10,965, capped at theoretical non-redundant capacity of $V_{\text{test}}$), C4 (11,450) with subcellular compartment annotations and audit logs.
- **Benchmark Data (`benchmark_data/`)**: Full combined benchmark datasets containing both positive (`label=1`) and negative (`label=0`) samples (C1: 214,082, C2: 12,595, C3: 12,110, C4: 12,595; Total: 251,382 pairs).
- **Network Visualizations (`figures/`)**: SVG diagrams and Cytoscape interactive HTML viewer.

---

## 5. Legacy Benchmarks & Raw Data (`data/legacy_benchmarks/` & `data/raw/`)

- **`data/legacy_benchmarks/`**: Preserved for historical comparisons (DeepAraPPI, ESMAraPPI) and rescued UniParc sequence assets.
- **`data/raw/`**: Immutable primary database downloads with cryptographic checksum verification.
