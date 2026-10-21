# Dataset Construction Data Warehouse

This directory houses all progressive artifacts of the fresh AraPPLM dataset construction pipeline. The dataset construction process is organized sequentially into numbered stages, starting from audited, verified positive interaction pairs.

---

## Pipeline Stages

```
data/dataset_construction/
├── README.md                      # Pipeline data architecture & status (this file)
│
├── 01_positive_pairs/             # [STAGE 1: STARTING POINT] Verified physical PPI pairs
│   ├── arabidopsis/               # 122,683 unique pairs (122,683 pure PPIs)
│   ├── maize/                     # 50,984 unique pairs (50,983 pure PPIs)
│   ├── rice/                      # 998 unique pairs (667 pure PPI pairs)
│   ├── tomato/                    # 123 unique pairs (123 pure PPI pairs)
│   └── soybean/                   # 48 unique pairs (48 pure PPI pairs)
│
├── 02_negative_sampling/          # [STAGE 2: UPCOMING] Verified negative pairs & decoy selection
│
├── 03_homology_clustering/        # [STAGE 3: UPCOMING] Sequence redundancy reduction & clustering
│
└── 04_benchmark_splits/           # [STAGE 4: UPCOMING] Strict C1/C2/C3 partitions & holdouts
```

---

## Stage 1: Verified Positive Pairs (`01_positive_pairs/`)

The foundational positive interaction sets were extracted directly from raw primary database releases (BioGRID v5.0.261, IntAct r252, STRING v12.0, PRIDE PhoX XL-MS) and independently audited.

Each species subfolder contains:
- `<species>_positive_evidence.parquet` & `.tsv`: Comprehensive provenance-tracked evidence records.
- `<species>_positive_pairs_deduplicated.parquet` & `.tsv`: Canonical deduplicated pair entries (`participant_a <= participant_b`) with recalculated continuous HIPPIE scores and confidence tiers.
- `<species>_positive_pairs_high_only.parquet` & `.tsv`: Streamlined 2-column minimal dataset (`participant_a_id`, `participant_b_id`) for High Quality pairs ($S \ge 0.72$).
- `<species>_positive_pairs_high_only_hetero.parquet` & `.tsv`: High Quality pairs excluding homomers (heteromers only).
- `<species>_positive_pairs_high_medium.parquet` & `.tsv`: Streamlined 2-column minimal dataset (`participant_a_id`, `participant_b_id`) for High + Medium Quality pairs ($S \ge 0.629$).
- `<species>_positive_pairs_high_medium_hetero.parquet` & `.tsv`: High + Medium Quality pairs excluding homomers (heteromers only).
- `<species>_positive_audit.json`: Machine-readable audit summary including filter counts, interaction topology, and HIPPIE scoring statistics.

### Summary Metrics

| Species | Taxonomy ID | Total Deduplicated Pairs | Pure PPI Pairs | Heteromers | Homomers | Unique Interactors | Primary Sources |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Arabidopsis** | `3702` | **120,873** | **120,873** | 115,012 | 5,861 | 20,836 | BioGRID, IntAct, STRING Full, PhoX XL-MS |
| **Maize** | `4577` | **50,984** | **50,983** | 50,409 | 575 | 14,888 | IntAct, BioGRID |
| **Rice** | `39947` / `4530` | **998** | **667** | 651 | 16 | 609 *(892 all)* | BioGRID, IntAct, STRING Full |
| **Tomato** | `4081` | **123** | **123** | 108 | 15 | 61 | BioGRID, IntAct |
| **Soybean** | `3847` | **48** | **48** | 47 | 1 | 51 | BioGRID, IntAct |
