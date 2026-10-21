# Dataset Construction Documentation

This directory documents the methodology, verification standards, and operational guidelines for constructing the new AraPPLM plant PPI datasets.

The dataset construction process starts from audited, verified physical interaction pairs and follows a modular multi-stage roadmap.

---

## Roadmap & Pipeline Structure

```
docs/dataset_construction/
├── README.md                              # Master overview & pipeline roadmap (this file)
│
├── 01_positive_pairs/                     # [STAGE 1: STARTING POINT] Verified positive interactomes
│   ├── 01_arabidopsis_positive_pairs.md   # Arabidopsis thaliana (122,683 deduplicated pairs)
│   ├── 02_maize_positive_pairs.md         # Zea mays (50,984 deduplicated pairs)
│   └── 03_rice_tomato_soybean_positive_pairs.md # Rice (998 / 667 PPI), Tomato (123), Soybean (48)
│
├── 02_negative_sampling/                  # [STAGE 2: UPCOMING] Decoy generation & non-interacting controls
│
├── 03_homology_clustering/                # [STAGE 3: UPCOMING] Redundancy reduction (CD-HIT / MMseqs2)
│
└── 04_benchmark_splits/                   # [STAGE 4: UPCOMING] C1/C2/C3 splits & cross-species benchmarks
```

---

## Stage 1: Verified Positive Pairs Ground Truth (`01_positive_pairs/`)

All downstream steps (negative generation, sequence clustering, and train/val/test partitioning) build directly upon the verified positive pairs documented in this section.

### Core Principles
1. **Raw Database Provenance**: Every positive interaction is traced back to a specific raw file in `data/raw/` (BioGRID v5.0.261, IntAct r252, STRING v12.0, or PRIDE PhoX XL-MS).
2. **Strict Physical Semantic Filtering**:
   - IntAct: Core physical interactions (`MI:0915`, `MI:0914`, `MI:0407` direct interaction, `MI:0217` physical/phosphorylation) with PCA retained (`MI:0112`, `MI:0809`, `MI:1204`, `MI:1037`, `MI:0090`). Strictly exclude proximity (`MI:2364`, `MI:0403`, `MI:1314` BioID).
   - BioGRID: `sys_type == 'physical'` with PCA retained (`exp_sys == 'PCA'`). Strictly exclude proximity labeling (`Proximity Label-MS`), co-localization, and genetic interactions.
   - STRING: Direct experimental channel only (`experiments > 0`). Strictly exclude transferred interologs (`experiments == 0 and experiments_transferred > 0`) and computational co-expression.
   - PRIDE PhoX XL-MS: High-confidence intra- and inter-protein cross-links under 1% PSM FDR (<35 Å Euclidean constraint).
3. **Pure Protein Interactions**: Strictly exclude non-protein interactors (DNA `MI:0319`, RNA `MI:0320`, small molecules `MI:0328`, gene loci).
4. **Canonical Deduplication**: Interactions are normalized as undirected pairs (`participant_a <= participant_b`) with provenance tracked.

### Species Overview

| Document | Target Species | Taxonomy ID | Total Deduplicated Pairs | Pure PPI Pairs | Key Data Sources |
| :--- | :--- | :--- | :--- | :--- | :--- |
| [`01_arabidopsis_positive_pairs.md`](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/AraPPLM/docs/dataset_construction/01_positive_pairs/01_arabidopsis_positive_pairs.md) | *Arabidopsis thaliana* | `3702` | **122,683** | **122,683** | BioGRID, IntAct, STRING Full, PhoX XL-MS |
| [`02_maize_positive_pairs.md`](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/AraPPLM/docs/dataset_construction/01_positive_pairs/02_maize_positive_pairs.md) | *Zea mays* | `4577` | **50,984** | **50,983** | IntAct, BioGRID |
| [`03_rice_tomato_soybean_positive_pairs.md`](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/AraPPLM/docs/dataset_construction/01_positive_pairs/03_rice_tomato_soybean_positive_pairs.md) | *Oryza sativa*<br>*Solanum lycopersicum*<br>*Glycine max* | `39947` / `4530`<br>`4081`<br>`3847` | **998** (Rice)<br>**123** (Tomato)<br>**48** (Soybean) | **667** (Rice)<br>**123** (Tomato)<br>**48** (Soybean) | BioGRID, IntAct, STRING Full |
