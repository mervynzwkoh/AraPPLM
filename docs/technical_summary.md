# Technical Summary: PPLM Benchmarking for Plant Protein-Protein Interaction

## Project Overview

**Objective:** Benchmark the pre-trained PPLM (Paired Protein Language Model) on *Arabidopsis thaliana* and *Oryza sativa* (Rice) protein-protein interaction (PPI) datasets to evaluate zero-shot cross-species generalization from animal/human to plant proteins.

**Key Challenge:** PPLM was trained on human, yeast, E. coli, C. elegans, D. melanogaster, and mouse PPI data. The goal is to evaluate its performance against the DeepAraPPI plant benchmarks without fine-tuning, and lay the groundwork for subsequent plant-specific fine-tuning.

---

## Project Architecture

```
AraPPLM/
├── PPLM/                               # [SUBMODULE] Original PPLM codebase (junliu621/PPLM)
│   ├── pplm/                           # Core PPLM model (33-layer Transformer)
│   ├── pplm_ppi/                       # PPI classifier (5-layer MLP)
│   ├── pplm_affinity/                  # Affinity prediction module
│   ├── pplm_contact/                   # Contact prediction module
│   └── weights/                        # pplm_t33_650M.pt & ppi_models.pkl (HPC)
│
├── scripts/                            # Custom execution & data tools
│   ├── batch_predict.py                # Core batch inference script
│   ├── evaluate_pplm.py                # Evaluation & DeepAraPPI benchmark comparison
│   ├── run_all_benchmarks_nscc.pbs     # PBS batch script for all 4 benchmark tasks
│   ├── run_batch_predict_nscc.pbs      # PBS script for single-task execution
│   └── data_preparation/
│       ├── build_sequence_db.py        # Unified FASTA -> PKL builder with regex
│       ├── verify_coverage.py          # Sequence database coverage verifier
│       ├── fetch_dataset_sequences.py  # Automated UniProt REST sequence fetcher
│       └── merge_rice_uniparc.py       # Rice UniParc merger script
│
├── data/
│   ├── DeepAraPPI/                     # Benchmark interaction datasets
│   │   ├── c1_ppi_sample_DeepAraPPI.txt # Task 1 (C1: ~1k pairs)
│   │   ├── c2_ppi_sample_DeepAraPPI.txt # Task 2 (C2: ~3k pairs, one unseen)
│   │   ├── c3_ppi_sample_DeepAraPPI.txt # Task 3 (C3: ~1.5k pairs, both unseen)
│   │   ├── total_positive_negative_samples_DeepAraPPI.txt # Full 130k pairs
│   │   └── all_rice_positive_negative_DeepAraPPI.txt      # Task 4 (Rice: 6.7k pairs)
│   │
│   ├── arabidopsis/                    # Arabidopsis sequence database
│   │   ├── missing_ids.txt             # 348 rescued IDs
│   │   ├── uniparc_348.fasta           # Rescued UniParc sequences
│   │   ├── uniparc_id_mapping.json     # ID mapping JSON
│   │   └── uniprot_final.pkl           # 100% complete database (~109k entries)
│   │
│   └── rice/                           # Rice sequence database
│       ├── missing_rice_ids.txt        # 163 rescued Rice IDs
│       ├── missing_rice.fasta          # Rescued UniParc sequences
│       └── uniprot_rice_final.pkl      # 100% complete database (~100k entries)
│
├── results/                            # Benchmark predictions & summaries
│
├── docs/                               # Comprehensive project documentation
│   ├── technical_summary.md            # Architecture & progress tracker
│   ├── data_processing_pipeline.md     # Detailed data preparation documentation
│   └── plant-pplm-methodology.md       # Plant-PPLM methodology design
│
└── PPLM_NSCC_A100.yml                  # A100-optimized Conda environment specification
```

---

## Benchmarking Progress

### Phase 1: Environment Setup ✓
- [x] Verified PPLM installation on NSCC ASPIRE 2A HPC (RHEL 8 + PBS Pro).
- [x] Configured `PPLM_NSCC_A100.yml` with PyTorch 2.x and CUDA 12.4.
- [x] Verified NVIDIA A100-SXM4-40GB GPU accessibility.
- [x] Model weights linked in `PPLM/weights/`.

### Phase 2: Sequence Database Preparation ✓
- [x] **Arabidopsis thaliana:**
  - Resolved Swiss-Prot/TrEMBL header formatting using regex parsing.
  - Rescued 348 legacy UniParc sequences.
  - **100.00% coverage verified** (130,478 / 130,478 pairs).
- [x] **Oryza sativa (Rice):**
  - Resolved *japonica* vs *indica* subspecies mismatch.
  - Rescued 163 missing accessions via UniParc integration (`merge_rice_uniparc.py`).
  - **100.00% coverage verified** (6,721 / 6,721 pairs).

### Phase 3: Batch Prediction Pipeline ✓
- [x] Path resolution configured in `scripts/batch_predict.py` for repository root.
- [x] GPU feature extraction (inter-chain attention, intra-chain attention, embeddings).
- [x] 10-weight classifier ensemble inference.
- [x] Configured PBS batch job runners for single-task and multi-task execution on A100 GPU.

### Phase 4: Evaluation Framework ✓
- [x] Built `scripts/evaluate_pplm.py` supporting:
  - **AUPRC** (Primary plant PPI metric).
  - **AUROC**, F1-Score, MCC, Accuracy, Sensitivity, Specificity.
  - Direct comparison tables against DeepAraPPI baselines (DeepAraPPI, GO2vec, Domain2vec, RCNN).
  - Multi-task automated consolidated reports (`results/benchmark_summary.csv`).

---

## Benchmark Suite Targets (DeepAraPPI)

| Task | Test Dataset | Difficulty | Primary Metric (AUPRC Baseline) |
| :--- | :--- | :--- | :--- |
| **Task 1** | `c1_ppi_sample_DeepAraPPI.txt` | Low (Seen domain) | DeepAraPPI: **0.965** |
| **Task 2** | `c2_ppi_sample_DeepAraPPI.txt` | Medium (One unseen) | DeepAraPPI: **0.897** |
| **Task 3** | `c3_ppi_sample_DeepAraPPI.txt` | High (Both unseen) | DeepAraPPI: **0.825** |
| **Task 4** | `all_rice_positive_negative_DeepAraPPI.txt` | Cross-Species (Rice) | DeepAraPPI: **0.305** |