# Dataset Construction Scripts Pipeline

This directory contains the computational scripts driving the construction of plant PPI benchmark datasets for AraPPLM.

---

## Pipeline Overview

```
scripts/dataset_construction/
├── README.md                      # Pipeline script architecture (this file)
│
├── 01_verify_positive_pairs.py    # [STAGE 1: STARTING POINT] Raw database parser & verification engine
│
├── 02_sample_negatives.py         # [STAGE 2: UPCOMING] Decoy sampling & negative control generation
│
├── 03_cluster_homology.py         # [STAGE 3: UPCOMING] Sequence redundancy clustering (CD-HIT / MMseqs2)
│
└── 04_generate_splits.py          # [STAGE 4: UPCOMING] Partitioning into C1 (train), C2 (semi-unseen), C3 (unseen)
```

---

## Stage 1: Positive Pair Verification (`01_verify_positive_pairs.py`)

Parses raw database downloads from `data/raw/` across all target plant species (*Arabidopsis*, *Maize*, *Rice*, *Tomato*, *Soybean*):
- Applies explicit PSI-MI physical interaction filters.
- Filters out non-protein participants (DNA, RNA, small molecules).
- Computes canonical undirected deduplicated pairs.
- Audits and matches the ground-truth records stored in `data/dataset_construction/01_positive_pairs/`.

### Usage
```bash
# Verify positive pair parsers
py -3.11 scripts/dataset_construction/01_verify_positive_pairs.py
```
