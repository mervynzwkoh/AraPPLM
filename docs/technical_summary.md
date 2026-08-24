# Technical Summary: PPLM Benchmarking for Plant Protein-Protein Interaction

## Project Overview

**Objective:** Benchmark the pre-trained PPLM (Paired Protein Language Model) on *Arabidopsis thaliana* and *Oryza sativa* (Rice) protein-protein interaction (PPI) datasets to evaluate zero-shot cross-species generalization from animal/human to plant proteins, and establish a foundation for plant-specific fine-tuning.

**Key Challenge:** PPLM was trained exclusively on human, yeast, E. coli, C. elegans, D. melanogaster, and mouse PPI data. The goal is to evaluate its zero-shot performance against plant benchmarks (DeepAraPPI, ESMAraPPI, AraCoFusion) without fine-tuning, and subsequently fine-tune Plant-PPLM.

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
│   ├── batch_predict.py                # Core batch inference script (with combined pair cropping)
│   ├── evaluate_pplm.py                # Evaluation & benchmark comparison engine
│   ├── run_all_benchmarks_nscc.pbs     # PBS batch script for complete 4-task execution
│   ├── run_remaining_benchmarks_nscc.pbs # PBS continuation script (Tasks 3, 4 & eval)
│   ├── run_batch_predict_nscc.pbs      # PBS script for single-task execution
│   └── data_preparation/
│       ├── build_sequence_db.py        # Unified FASTA -> PKL builder with regex
│       ├── verify_coverage.py          # Sequence database coverage verifier
│       ├── fetch_dataset_sequences.py  # Automated UniProt REST sequence fetcher
│       └── merge_rice_uniparc.py       # Rice UniParc merger script
│
├── data/
│   ├── DeepAraPPI/                     # DeepAraPPI benchmark interaction datasets
│   │   ├── c1_ppi_sample_DeepAraPPI.txt # Task 1 (C1: 31,284 pairs)
│   │   ├── c2_ppi_sample_DeepAraPPI.txt # Task 2 (C2: 66,055 pairs, one unseen)
│   │   ├── c3_ppi_sample_DeepAraPPI.txt # Task 3 (C3: 33,099 pairs, both unseen)
│   │   ├── total_positive_negative_samples_DeepAraPPI.txt # Full 130k pairs
│   │   └── all_rice_positive_negative_DeepAraPPI.txt      # Task 4 (Rice: 6,721 pairs)
│   │
│   ├── ESMAraPPI/                      # [UPCOMING] ESMAraPPI benchmark datasets
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
├── results/
│   ├── DeepAraPPI/                     # DeepAraPPI benchmark predictions & evaluations
│   │   ├── deepara_c1_scores.csv
│   │   ├── deepara_c2_scores.csv
│   │   ├── deepara_c3_scores.csv
│   │   ├── deepara_rice_scores.csv
│   │   ├── benchmark_summary.csv
│   │   └── task*_metrics.txt
│   │
│   └── ESMAraPPI/                      # [UPCOMING] ESMAraPPI benchmark results
│
├── docs/                               # Comprehensive project documentation
│   ├── technical_summary.md            # Architecture & progress tracker
│   ├── data_processing_pipeline.md     # Detailed data preparation documentation
│   ├── plant-pplm-methodology.md       # Plant-PPLM methodology design
│   ├── benchmark_analysis_deeparappi_vs_pplm.md # DeepAraPPI technical report
│   └── lit_review/                     # Paper information extraction templates
│       ├── LitReview_DeepAraPPI.md
│       ├── LitReview_ESMAraPPI.md
│       ├── LitReview_AraCoFusion.md
│       └── Lit_Review_PPLM.md
│
└── PPLM_NSCC_A100.yml                  # A100-optimized Conda environment specification
```

---

## Benchmarking Progress & Empirical Results

### Phase 1: Environment & Tooling Setup ✓
- [x] Verified PPLM on NSCC ASPIRE 2A HPC (A100 GPU).
- [x] Implemented combined paired length cropping ($L_A + L_B \le 1020$) matching PPLM paper pretraining.
- [x] Implemented 10-fold classifier ensemble with symmetric pair averaging.

### Phase 2: Sequence Database Preparation ✓
- [x] **Arabidopsis:** 100.00% sequence coverage (109,640 entries).
- [x] **Rice:** 100.00% sequence coverage (100,297 entries).

### Phase 3: DeepAraPPI Benchmark Completed ✓

Across all 137,159 test pairs in the DeepAraPPI suite:

| Task | Test Dataset | Difficulty | DeepAraPPI Baseline (AUPRC) | PPLM Zero-Shot (AUPRC) | PPLM AUROC | PPLM Specificity |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Task 1** | `c1_ppi_sample_DeepAraPPI.txt` | Low (Seen domain) | **0.9650** | **0.6095** | **0.8991** | 99.55% |
| **Task 2** | `c2_ppi_sample_DeepAraPPI.txt` | Medium (One unseen) | **0.8970** | **0.5738** | **0.8828** | 99.57% |
| **Task 3** | `c3_ppi_sample_DeepAraPPI.txt` | High (Both unseen) | **0.8250** (Seq-only: 0.4810) | **0.5525** | **0.8710** | 99.50% |
| **Task 4** | `all_rice_positive_negative_DeepAraPPI.txt` | Cross-Species (Rice) | **0.3050** | **0.4297** | **0.7561** | 98.35% |

### Key Findings:
1. **Cross-Species Transfer (Rice):** PPLM beats DeepAraPPI by **+40.9%** (0.4297 vs 0.3050) and beats RCNN sequence baseline (0.2480) by **+73.3%**.
2. **Hard Unseen Pairs (Task 3):** PPLM (**0.5525**) outperforms DeepAraPPI's sequence-only RCNN model (**0.4810**) and Random Forest (**0.4340**).
3. **Arabidopsis In-Domain:** High AUROC (0.87–0.90) and specificity (>99.5%) provide a solid baseline for Phase 4 fine-tuning.

---

## Next Steps: ESMAraPPI Benchmark & Fine-Tuning
1. Download / prepare ESMAraPPI datasets into `data/ESMAraPPI/`.
2. Verify sequence ID coverage against existing sequence databases.
3. Run batch prediction and evaluate performance in `results/ESMAraPPI/`.
4. Proceed to Plant-PPLM LoRA Fine-Tuning.