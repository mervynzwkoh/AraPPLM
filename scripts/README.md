# AraPPLM Scripts Directory Architecture

This directory houses all computational scripts for the AraPPLM project. Following codebase sanitization, outdated and non-reproducible dataset construction scripts have been purged. The retained scripts are organized into modular, purpose-driven subdirectories.

---

## Directory Overview

```
scripts/
├── README.md                      # Master scripts guide (this file)
│
├── dataset_construction/          # [DATASET CONSTRUCTION] Progressive dataset construction pipeline
│   ├── README.md                  # Dataset construction pipeline script guide
│   └── 01_verify_positive_pairs.py # [STAGE 1: STARTING POINT] Raw DB parser & verification engine
│
├── verify_raw_positive_pairs.py   # Convenience wrapper for 01_verify_positive_pairs.py
├── phase1_download_sources.py     # Download manager for primary databases (BioGRID, IntAct, STRING, etc.)
├── verify_phase1_integrity.py     # Checksum & integrity validator for downloaded raw archives
│
├── data_preparation/              # Proteome FASTA parsing, UniParc legacy rescue, & sequence DB generation
│   ├── build_sequence_db.py       # Compiles proteome FASTAs into indexed pickle lookups
│   ├── merge_rice_uniparc.py      # Merges rescued UniParc sequences for rice accessions
│   ├── merge_esmarappi_sequences.py # Merges rescued secondary UniProt accessions for ESMAraPPI
│   └── verify_coverage.py         # Validates 100.00% sequence coverage against benchmark datasets
│
├── phase2/                        # Robust database parsers & evidence harmonizers
│   ├── parse_biogrid.py           # BioGRID TAB 3.0 parser
│   ├── parse_intact.py            # IntAct PSI-MI TAB 2.7 parser
│   ├── parse_string.py            # STRING v12.0 physical links parser
│   ├── parse_xlms.py              # PRIDE PhoX XL-MS CSV parser (Trinh et al., 2026)
│   ├── parse_bipseq.py            # BIP-seq / PoPPIN XML protocol auditor
│   ├── parse_legacy.py            # DeepAraPPI / ESMAraPPI benchmark parser
│   ├── harmonize_assays.py        # Maps experimental assay names to PSI-MI ontologies
│   ├── build_assay_mapping.py     # Controlled vocabulary table builder
│   ├── combine_evidence.py        # Normalized evidence combiner
│   ├── generate_manifest_and_checksums.py # Checksum manifest compiler
│   └── validate_phase2.py         # QC gate validator
│
├── training/                      # PPI classification head training & evaluation
│   ├── ppi_model.py               # PyTorch classification head (MLP over pooled embeddings)
│   ├── dataset.py                 # PyTorch Dataset implementation for pair representations
│   ├── train_ppi_head.py          # Training loop with early stopping & CosineAnnealing LR
│   ├── test_ppi_head.py           # Comprehensive metric evaluation (AUC, AUPRC, F1, MCC)
│   ├── extract_features.py        # Feature extractor from paired sequence representations
│   └── select_top_models.py       # Checkpoint ranker and selector
│
├── benchmarking/                  # Large-scale inference & benchmark evaluation
│   ├── batch_predict.py           # HPC-optimized batched inference engine
│   └── evaluate_pplm.py           # Model scoring against legacy test suites
│
└── hpc/                           # NSCC PBS batch scripts for GPU execution
    ├── run_batch_predict_nscc.pbs
    ├── run_train_deeparappi_nscc.pbs
    ├── run_train_esmarappi_nscc.pbs
    ├── run_esmarappi_benchmarks_nscc.pbs
    ├── run_all_benchmarks_nscc.pbs
    └── extract_ara_pplm_embeddings.py
```

---

## 1. Dataset Construction Pipeline (`scripts/dataset_construction/`)

All scripts driving the progressive plant PPI dataset construction pipeline reside in `scripts/dataset_construction/`:
- **`01_verify_positive_pairs.py`**: Directly parses raw database files in `data/raw/` to compute and verify the authoritative positive interaction numbers for all five plant species (*Arabidopsis*, *Maize*, *Rice*, *Tomato*, *Soybean*). Applies exact physical interaction semantic filters, filters out non-protein participants, and computes canonical deduplicated pair metrics.
- Subsequent stages (`02_sample_negatives.py`, `03_cluster_homology.py`, `04_generate_splits.py`) will be sequentially integrated here.


Directly parses raw database files in `data/raw/` to compute and verify the authoritative positive interaction numbers for all five plant species:
- *Arabidopsis thaliana* (`3702`)
- *Zea mays* (`4577`)
- *Oryza sativa* (`39947` / `4530`)
- *Solanum lycopersicum* (`4081`)
- *Glycine max* (`3847`)

The parser evaluates raw record volume, applies exact physical interaction semantic filters, filters out non-protein participants (DNA/RNA/small molecules), and computes canonical deduplicated pair metrics.

---

## 2. Sequence Preparation (`scripts/data_preparation/`)

Responsible for generating the indexed sequence databases required by PPLM:
- **`build_sequence_db.py`**: Compiles proteome FASTAs into indexed pickle lookup maps (`data/arabidopsis/uniprot_final.pkl` and `data/rice/uniprot_rice_final.pkl`).
- **`merge_rice_uniparc.py` & `merge_esmarappi_sequences.py`**: Incorporates legacy or secondary accession sequences fetched from UniParc/UniProtKB.
- **`verify_coverage.py`**: Ensures 100% accession coverage across DeepAraPPI (Tasks 1, 2, 3, Rice) and ESMAraPPI (Tasks C2, C3).

---

## 3. Database Parsers (`scripts/phase2/`)

Provides individual parsers that extract standardized interactions from each source repository:
- **`parse_biogrid.py`**: Extracts physical interactions and separates proximity labeling.
- **`parse_intact.py`**: Classifies PSI-MI interaction types (physical association, direct interaction, colocalization).
- **`parse_string.py`**: Isolates direct experimental links (`experiments > 0`), stripping transferred interologs.
- **`parse_xlms.py`**: Extracts high-confidence cross-linked peptide pairs under 1% PSM FDR.
- **`harmonize_assays.py`**: Standardizes assay descriptions across databases.

---

## 4. PPI Head Retraining & Benchmarking (`scripts/training/` & `scripts/benchmarking/`)

Scripts for training and benchmarking the 2-layer MLP classification head on frozen PPLM representations:
- **`train_ppi_head.py`**: Trains on balanced plant datasets with cross-entropy loss, dropout, and learning rate scheduling.
- **`test_ppi_head.py`**: Evaluates test sets across Accuracy, Precision, Recall, Specificity, F1-score, MCC, AUROC, and AUPRC.
- **`batch_predict.py`**: High-throughput batched inference on GPU with expandable memory segment allocation and sequence length budgeting.

---

## 5. HPC Cluster Execution (`scripts/hpc/`)

Contains ready-to-submit PBS job submission scripts configured for NSCC Singapore (National Supercomputing Centre) GPU nodes (NVIDIA V100/A100).
