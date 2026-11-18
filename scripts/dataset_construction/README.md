# Dataset Construction Scripts Pipeline

This directory contains the computational scripts driving the end-to-end construction of plant PPI benchmark datasets for AraPPLM from raw database intake to final benchmark splits.

---

## Pipeline Overview

```
scripts/dataset_construction/
├── README.md                                      # Master script architecture (this file)
│
├── raw_acquisition/                               # [STAGE 0] Raw database downloaders & checksum verifiers
│   ├── phase1_download_sources.py                 # Automated downloader for BioGRID, IntAct, STRING, PhoX XL-MS, TAIR
│   ├── verify_phase1_integrity.py                 # SHA-256 checksum & integrity validator for downloaded archives
│   └── verify_raw_positive_pairs.py               # Convenience wrapper for positive pair verification
│
├── parsing_and_harmonization/                     # [STAGE 1 INTAKE] Standardized parsers & PSI-MI assay ontologies
│   ├── parse_biogrid.py                           # BioGRID TAB 3.0 parser
│   ├── parse_intact.py                            # IntAct PSI-MI TAB 2.7 parser
│   ├── parse_string.py                            # STRING v12.0 physical links parser
│   ├── parse_xlms.py                              # PRIDE PhoX XL-MS CSV parser (Trinh et al., 2026)
│   ├── parse_bipseq.py                            # BIP-seq / PoPPIN XML protocol auditor
│   ├── parse_legacy.py                            # DeepAraPPI / ESMAraPPI benchmark parser
│   ├── harmonize_assays.py                        # Maps experimental assay names to PSI-MI ontologies
│   ├── build_assay_mapping.py                     # Controlled vocabulary table builder
│   ├── combine_evidence.py                        # Normalized evidence combiner
│   ├── generate_manifest_and_checksums.py         # Checksum manifest compiler
│   └── validate_phase2.py                         # QC gate validator
│
├── 01_positive_pairs/                             # [STAGE 1] Positive interactome verification & scoring
│   ├── 01_verify_positive_pairs.py                # Positive interactome verification engine across species
│   ├── apply_hippie_scoring.py                    # Continuous confidence scoring based on HIPPIE formula
│   ├── update_arabidopsis_verified_pairs.py       # Master positive deduplicated interactome builder (120,873 pairs)
│   └── build_high_only_hetero_fasta.py            # Extracts FASTA for high-confidence positive proteins
│
├── 02_negative_sampling/                          # [STAGE 2] 6-phase candidate negative generation & filtering
│   ├── 02_build_positive_blacklist_hpc.py         # Exhaustive positive blacklist compiler (~3.8M pairs)
│   ├── 02_generate_candidate_negatives_hpc.py     # Subcellular localization candidate generator (101.88M)
│   ├── 03_filter_negative_positive_interologs_hpc.py # Positive interolog filter (101.85M clean pairs)
│   ├── 04_cluster_candidate_negatives_mmseqs_hpc.py # Paralog redundancy reduction among negatives (50.31M NR pool)
│   ├── run_filter_negative_interologs_hpc.sh      # Bash wrapper for Stage 2.5 interolog filtering
│   ├── submit_filter_negative_interologs_hpc.pbs  # PBS submission script for Stage 2.5
│   ├── run_cluster_candidate_negatives_hpc.sh     # Bash wrapper for Stage 2.6 CD-HIT redundancy reduction
│   ├── submit_cluster_candidate_negatives_hpc.pbs # PBS submission script for Stage 2.6
│   └── submit_negative_sampling_hpc.pbs           # PBS submission script for Stage 2.3 candidate generator
│
├── 03_homology_clustering/                        # [STAGE 3] Pair-level CD-HIT redundancy clustering for positives
│   ├── cluster_positive_pairs_mmseqs.py           # Inverted homology indexing & greedy CD-HIT clustering (22,897 pairs)
│   ├── evaluate_cluster_density.py                # Parameter sweep density evaluator across top-20 clusters
│   ├── run_pair_clustering_pipeline.sh            # Bash wrapper for Stage 3 clustering
│   ├── submit_pair_clustering_hpc.pbs             # PBS submission script for Stage 3 clustering
│   ├── run_cluster_density_hpc.sh                 # Bash wrapper for Stage 3 density analysis
│   └── submit_cluster_density_hpc.pbs             # PBS submission script for Stage 3 density analysis
│
└── 04_benchmark_splits/                           # [STAGE 4] C1/C2/C3/C4 benchmark graph partitioning (1:10 ratio)
    ├── 04_partition_positive_c123.py              # Graph-theoretic C1/C2/C3/C4 bisection for positive pairs
    ├── 05_partition_negative_c1234_hpc.py         # Negative C1/C2/C3/C4 benchmark partitioner (1:10 ratio)
    ├── run_partition_negative_hpc.sh              # Bash wrapper for Stage 2.7 / Stage 4 negative partitioning
    ├── submit_partition_negative_hpc.pbs          # PBS submission script for Stage 4 negative partitioning
    └── 05_generate_c123_visualizations.py         # SVG network graphs & interactive Cytoscape network visualizer
```

---

## 1. Stage 0 & 1: Raw Acquisition, Parsing, and Positive Verification

### Raw Acquisition (`raw_acquisition/`)
Downloads and validates primary data archives from authoritative public repositories:
- **`phase1_download_sources.py`**: Automated download manager downloading all 32 primary database archives (BioGRID 5.0.261, IntAct r252, STRING v12.0, PRIDE PhoX XL-MS, TAIR GAF, SUBA5, and holdout species).
- **`verify_phase1_integrity.py`**: Verifies SHA-256 checksums of downloaded files against `checksums/PHASE1_SHA256SUMS` and registers metadata in `manifests/phase1_manifest.json`.

### Raw Parsing & Assay Harmonization (`parsing_and_harmonization/`)
Extracts standardized interactions and harmonizes assay ontologies across databases:
- **`parse_biogrid.py`**: Extracts physical interactions, strictly separating proximity labeling.
- **`parse_intact.py`**: Classifies PSI-MI interaction types (physical association, direct interaction, colocalization).
- **`parse_string.py`**: Isolates direct experimental links (`experiments > 0`), stripping transferred interologs.
- **`parse_xlms.py`**: Extracts high-confidence cross-linked peptide pairs under 1% PSM FDR.
- **`harmonize_assays.py`**: Standardizes assay descriptions across databases to PSI-MI controlled vocabularies.
- **`combine_evidence.py`**: Merges records into a unified multi-column evidence table.

### Positive Pair Auditing & Scoring (`01_positive_pairs/`)
- `01_verify_positive_pairs.py`: Audits and matches the ground-truth records stored in `data/dataset_construction/01_positive_pairs/`.
- `apply_hippie_scoring.py`: Computes calibrated continuous confidence scores ($0.0 - 1.0$) using the HIPPIE scoring model.
- `update_arabidopsis_verified_pairs.py`: Builds the master deduplicated positive interactomes (120,873 total pairs; 30,376 high-confidence heteromeric pairs).

---

## 2. Stage 2: Negative Candidate Pool Generation & Filtering (`02_negative_sampling/`)

A 6-phase negative screening pipeline executed on HPC:
1. `02_build_positive_blacklist_hpc.py`: Aggregates ~3.8M known interactions across BioGRID, IntAct, STRING, PhoX XL-MS, DeepAraPPI, and ESMAraPPI into an indexed parquet blacklist.
2. `02_generate_candidate_negatives_hpc.py`: Screens all non-self pairs across 28,609 canonical loci, enforcing disjoint subcellular compartments across 10 TAIR cellular compartments (101.88M candidate pairs in 102 Parquet chunks).
3. `03_filter_negative_positive_interologs_hpc.py`: Screens candidate negatives against 22,897 confirmed positive interactions, purging all pairs with interolog homology $S(P_{\text{neg}}, P_{\text{pos}}) \ge 0.50$ (101.85M clean pairs retained).
4. `04_cluster_candidate_negatives_mmseqs_hpc.py`: Applies the Component Decoupling Theorem and greedy CD-HIT representative selection to prune paralogous redundancy among negatives (50.31M non-redundant pairs retained in 51 chunks).

---

## 3. Stage 3: Homology Clustering & Redundancy Reduction (`03_homology_clustering/`)

Applies the pair similarity metric $S(P_i, P_j) = \max[\min(s(A,C), s(B,D)), \min(s(A,D), s(B,C))]$ across the 30,376 high-confidence heteromeric positive pairs:
- `cluster_positive_pairs_mmseqs.py`: Runs all-vs-all MMseqs2 search (`-s 7.5`, reciprocal coverage mode) and performs greedy CD-HIT clustering prioritized by HIPPIE score $\to$ evidence count $\to$ combined sequence length.
- `evaluate_cluster_density.py`: Audits cluster density across parameter sweeps.
- Final Selected Condition: $\theta = 0.50, c_{\min} = 0.80$ (**22,897 non-redundant positive pairs** across 8,426 unique proteins).

---

## 4. Stage 4: Benchmark Graph Partitioning (C1/C2/C3/C4 at 1:10 Ratio) (`04_benchmark_splits/`)

Partitions the positive and negative interactomes into **C1 (Train, 85%)**, **C2 (One Unseen Test, 5%)**, **C3 (Both Unseen Test, 5%)**, and **C4 (Both Seen Test, 5%)**:
- `04_partition_positive_c123.py`: Two-stage graph bisection separating $V$ into disjoint $V_{\text{train}}$ and $V_{\text{test}}$, isolating C3 (both unseen) and extracting C4 (both seen novel edges).
- `05_partition_negative_c1234_hpc.py`: Partitions negatives at a 1:10 ratio matching positive vertex rules with zero homology leakage ($V(C_3) \cap V(C_1) = \emptyset$).
- `05_generate_c123_visualizations.py`: Exports publication-ready SVG network graphs and standalone Cytoscape.js interactive HTML network viewers.
