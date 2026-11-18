# AraPPLM Scripts Directory Architecture

This directory houses all computational scripts for the AraPPLM project. All scripts are organized into modular, purpose-driven subdirectories categorized strictly by the process carried out.

---

## Directory Overview

```
scripts/
├── README.md                                  # Master scripts architecture (this file)
│
├── dataset_construction/                      # [DATASET CONSTRUCTION] Complete pipeline for AraPPLM benchmarks
│   ├── README.md                              # Dataset construction script guide & execution order
│   ├── raw_acquisition/                       # [STAGE 0] Raw database downloaders & checksum verifiers
│   │   ├── phase1_download_sources.py         # Automated downloader for BioGRID, IntAct, STRING, PhoX XL-MS, TAIR
│   │   ├── verify_phase1_integrity.py         # SHA-256 checksum & integrity validator for downloaded archives
│   │   └── verify_raw_positive_pairs.py       # Convenience wrapper for positive pair verification
│   │
│   ├── parsing_and_harmonization/             # [STAGE 1 INTAKE] Standardized parsers & PSI-MI assay ontologies
│   │   ├── parse_biogrid.py                   # BioGRID TAB 3.0 parser
│   │   ├── parse_intact.py                    # IntAct PSI-MI TAB 2.7 parser
│   │   ├── parse_string.py                    # STRING v12.0 physical links parser
│   │   ├── parse_xlms.py                      # PRIDE PhoX XL-MS CSV parser (Trinh et al., 2026)
│   │   ├── parse_bipseq.py                    # BIP-seq / PoPPIN XML protocol auditor
│   │   ├── parse_legacy.py                    # DeepAraPPI / ESMAraPPI benchmark parser
│   │   ├── harmonize_assays.py                # Maps experimental assay names to PSI-MI ontologies
│   │   ├── build_assay_mapping.py             # Controlled vocabulary table builder
│   │   ├── combine_evidence.py                # Normalized evidence combiner
│   │   ├── generate_manifest_and_checksums.py # Checksum manifest compiler
│   │   └── validate_phase2.py                 # QC gate validator
│   │
│   ├── 01_positive_pairs/                     # [STAGE 1] Positive interactome verification & scoring
│   │   ├── 01_verify_positive_pairs.py        # Positive interactome verification engine across species
│   │   ├── apply_hippie_scoring.py            # Continuous confidence scoring based on HIPPIE formula
│   │   ├── update_arabidopsis_verified_pairs.py # Master positive deduplicated interactome builder (120,873 pairs)
│   │   └── build_high_only_hetero_fasta.py    # Extracts FASTA for high-confidence positive proteins
│   │
│   ├── 02_negative_sampling/                  # [STAGE 2] 6-phase candidate negative generation & filtering
│   │   ├── 02_build_positive_blacklist_hpc.py # Exhaustive positive blacklist compiler (~3.8M pairs)
│   │   ├── 02_generate_candidate_negatives_hpc.py # Subcellular localization candidate generator (101.88M)
│   │   ├── 03_filter_negative_positive_interologs_hpc.py # Positive interolog filter (101.85M clean pairs)
│   │   ├── 04_cluster_candidate_negatives_mmseqs_hpc.py # Paralog redundancy reduction among negatives (50.31M NR pool)
│   │   ├── run_filter_negative_interologs_hpc.sh # Bash wrapper for Stage 2.5 interolog filtering
│   │   ├── submit_filter_negative_interologs_hpc.pbs # PBS submission script for Stage 2.5
│   │   ├── run_cluster_candidate_negatives_hpc.sh # Bash wrapper for Stage 2.6 CD-HIT redundancy reduction
│   │   ├── submit_cluster_candidate_negatives_hpc.pbs # PBS submission script for Stage 2.6
│   │   └── submit_negative_sampling_hpc.pbs   # PBS submission script for Stage 2.3 candidate generator
│   │
│   ├── 03_homology_clustering/                # [STAGE 3] Pair-level CD-HIT redundancy clustering for positives
│   │   ├── cluster_positive_pairs_mmseqs.py   # Inverted homology indexing & greedy CD-HIT clustering (22,897 pairs)
│   │   ├── evaluate_cluster_density.py        # Parameter sweep density evaluator across top-20 clusters
│   │   ├── run_pair_clustering_pipeline.sh    # Bash wrapper for Stage 3 clustering
│   │   ├── submit_pair_clustering_hpc.pbs     # PBS submission script for Stage 3 clustering
│   │   ├── run_cluster_density_hpc.sh         # Bash wrapper for Stage 3 density analysis
│   │   └── submit_cluster_density_hpc.pbs     # PBS submission script for Stage 3 density analysis
│   │
│   └── 04_benchmark_splits/                   # [STAGE 4] C1/C2/C3/C4 benchmark graph partitioning (1:10 ratio)
│       ├── 04_partition_positive_c123.py      # Graph-theoretic C1/C2/C3/C4 bisection for positive pairs
│       ├── 05_partition_negative_c1234_hpc.py # Negative C1/C2/C3/C4 benchmark partitioner (1:10 ratio)
│       ├── run_partition_negative_hpc.sh      # Bash wrapper for Stage 2.7 / Stage 4 negative partitioning
│       ├── submit_partition_negative_hpc.pbs  # PBS submission script for Stage 4 negative partitioning
│       └── 05_generate_c123_visualizations.py # SVG network graphs & interactive Cytoscape network visualizer
│
├── training/                                  # [MODEL TRAINING] PPI classification head training & evaluation
│   ├── ppi_model.py                           # PyTorch 5-layer MLP head + author 9-metric evaluate()
│   ├── dataset.py                             # Dataset_PPI_Interact (15 pooled features, order_mode) & CV splits
│   ├── train_ppi_head.py                      # Per-fold 10-fold CV training (AdamW lr=5e-5, 15 epochs, BCE)
│   ├── test_ppi_head.py                       # Batched multi-model multi-pooling evaluation (order_mode=both)
│   ├── extract_features.py                    # Feature extractor from paired sequence representations
│   ├── select_top_models.py                   # Checkpoint ranker (column-14 AUC_PR) and ensemble packager
│   └── hpc/                                   # HPC batch scripts for model training & feature extraction
│       ├── extract_ara_pplm_embeddings.py     # PPLM paired embedding extractor for positive pairs
│       ├── run_extract_ara_embeddings_nscc.pbs # PBS job script for PPLM embedding extraction
│       ├── run_extract_araplm_c1_shard_nscc.pbs # Sharded C1 feature extraction for curated AraPPLM
│       ├── run_benchmark_araplm_nscc.pbs      # Zero-shot + author-aligned retraining benchmark on AraPPLM
│       ├── run_train_deeparappi_nscc.pbs      # PBS job script for training DeepAraPPI models
│       ├── run_train_esmarappi_nscc.pbs       # PBS job script for training ESMAraPPI models
│       └── run_resume_deeparappi_nscc.pbs     # PBS job script to resume DeepAraPPI training
│
└── benchmarking/                              # [EVALUATION] Large-scale inference & benchmark evaluation
    ├── batch_predict.py                       # HPC-optimized batched inference engine
    ├── evaluate_pplm.py                       # Model scoring against legacy test suites
    ├── data_preparation/                     # Proteome FASTA parsing & legacy sequence rescue lookups
    │   ├── build_sequence_db.py               # Compiles proteome FASTAs into indexed pickle lookups
    │   ├── merge_rice_uniparc.py              # Merges rescued UniParc sequences for rice accessions
    │   ├── merge_esmarappi_sequences.py       # Merges rescued secondary UniProt accessions
    │   └── verify_coverage.py                 # Validates 100.00% sequence coverage
    └── hpc/                                   # HPC batch scripts for benchmark evaluation & inference
        ├── run_batch_predict_nscc.pbs         # PBS job script for batched prediction
        ├── run_all_benchmarks_nscc.pbs        # PBS job script evaluating all benchmarks
        ├── run_esmarappi_benchmarks_nscc.pbs  # PBS job script evaluating ESMAraPPI benchmarks
        ├── run_remaining_benchmarks_nscc.pbs  # PBS job script evaluating remaining test suites
        └── run_cluster_and_calibrate_nscc.pbs # PBS job script for embedding clustering & calibration
```

---

## 1. Dataset Construction Pipeline (`scripts/dataset_construction/`)

Drives the 4-stage dataset construction:
- **`raw_acquisition/`**: Downloads 32 external database archives and validates SHA-256 integrity.
- **`parsing_and_harmonization/`**: Normalizes database formats, filters physical interactions, and harmonizes assays into PSI-MI ontologies.
- **`01_positive_pairs/`**: Verifies positive pairs across plant species, applies continuous HIPPIE scoring, and compiles FASTA sequences.
- **`02_negative_sampling/`**: Generates negative candidate pairs, enforces subcellular disjointness, and purges interologs and paralog redundancy.
- **`03_homology_clustering/`**: Performs pair-level CD-HIT redundancy clustering on positive pairs via MMseqs2.
- **`04_benchmark_splits/`**: Partitions positive and negative interactomes into C1/C2/C3/C4 benchmark splits at 1:10 ratio with zero leakage.

---

## 2. Model Training Pipeline (`scripts/training/`)

Handles per-fold model training, top-K selection, and batched multi-pooling testing matching the original PPLM author's workflow:
- **`train_ppi_head.py`**: Trains the top 5-layer MLP classification head on frozen PPLM features (10-fold CV, fixed 15 epochs, AdamW `lr=5e-5`, `--seed_round 1`).
- **`test_ppi_head.py`**: Batched test-set inference across all 6 pooling combinations (`mean`, `max`, `min`, `combine_mean_max`, `combine_mean_max_min`, `combine_mean_min`) with symmetric averaging (`order_mode=both`).
- **`extract_features.py`**: Feature extractor generating 15 pooled attention/embedding tensors per pair for both train and test splits.
- **`hpc/`**: Production PBS job scripts for feature extraction and multi-epoch model training on NSCC ASPIRE 2A.

---

## 3. Benchmarking Suite (`scripts/benchmarking/`)

Zero-shot evaluation and legacy comparative benchmarks:
- **`batch_predict.py`**: Streamlined inference engine over large benchmark test sets.
- **`data_preparation/`**: Scripts restoring 100% sequence coverage for legacy benchmarks (DeepAraPPI, ESMAraPPI).
- **`hpc/`**: Production PBS job scripts for large-scale inference and automated benchmark scoring.
