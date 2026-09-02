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
├── scripts/
│   ├── batch_predict.py                # Core batch inference script (with combined pair cropping)
│   ├── evaluate_pplm.py                # Evaluation & benchmark comparison engine
│   ├── run_all_benchmarks_nscc.pbs     # PBS: 4-task zero-shot benchmarking
│   ├── run_esmarappi_benchmarks_nscc.pbs # PBS: ESMAraPPI zero-shot benchmarking
│   ├── run_train_deeparappi_nscc.pbs   # PBS: End-to-end DeepAraPPI C1 training
│   ├── run_train_esmarappi_nscc.pbs    # PBS: End-to-end ESMAraPPI C1 training
│   ├── run_resume_deeparappi_nscc.pbs  # PBS: Resume DeepAraPPI from Stage 3
│   ├── benchmarking/                   # Benchmark script copies (for training pipeline)
│   │   ├── batch_predict.py
│   │   └── evaluate_pplm.py
│   ├── training/                       # PPI head retraining pipeline
│   │   ├── extract_features.py         # Stage 1: PPLM backbone → pooled feature .pkl files
│   │   ├── train_ppi_head.py           # Stage 2: 10-fold stratified CV training
│   │   ├── select_top_models.py        # Stage 3: Rank by AUPRC, package top-5 ensemble
│   │   ├── test_ppi_head.py            # Stage 4: Ensemble inference on test sets
│   │   ├── ppi_model.py                # PPI MLP classifier definition
│   │   └── dataset.py                  # PyTorch Dataset for pre-extracted features
│   └── data_preparation/
│       ├── build_sequence_db.py        # Unified FASTA -> PKL builder
│       ├── verify_coverage.py          # Sequence database coverage verifier
│       └── merge_rice_uniparc.py       # Rice UniParc merger script
│
├── data/
│   ├── DeepAraPPI/                     # DeepAraPPI benchmark datasets
│   │   ├── c1_ppi_sample_DeepAraPPI.txt # C1 training set (31,284 pairs)
│   │   ├── c2_ppi_sample_DeepAraPPI.txt # C2 test (66,055 pairs, one unseen)
│   │   ├── c3_ppi_sample_DeepAraPPI.txt # C3 test (33,099 pairs, both unseen)
│   │   └── all_rice_positive_negative_DeepAraPPI.txt # Rice (6,721 pairs)
│   ├── ESMAraPPI/                      # ESMAraPPI benchmark datasets
│   │   ├── c1_Train.txt                # C1 training set (38,709 pairs)
│   │   ├── c2Pred.txt                  # C2 test (37,444 pairs, one unseen)
│   │   └── c3Pred.txt                  # C3 test (8,866 pairs, both unseen)
│   ├── arabidopsis/                    # Arabidopsis sequence database
│   │   └── uniprot_final.pkl           # 100% complete (109,994 entries)
│   └── rice/                           # Rice sequence database
│       └── uniprot_rice_final.pkl      # 100% complete (100,297 entries)
│
├── features/                           # Pre-extracted PPLM backbone features
│   ├── DeepAraPPI_C1/                  # Feature .pkl files (31,284 pairs)
│   └── ESMAraPPI_C1/                   # Feature .pkl files (38,709 pairs)
│
├── models/                             # Trained PPI head checkpoints
│   ├── DeepAraPPI/                     # 10-fold CV checkpoints + ensemble weights
│   └── ESMAraPPI/                      # 10-fold CV checkpoints + ensemble weights
│
├── results/
│   ├── DeepAraPPI/                     # Zero-shot + plant-trained predictions
│   └── ESMAraPPI/                      # Zero-shot + plant-trained predictions
│
├── logs/                               # HPC PBS job logs
│   ├── benchmarking/                   # Zero-shot benchmark run logs
│   └── training/                       # PPI head training run logs
│
├── docs/                               # Project documentation
│   ├── technical_summary.md            # Architecture & progress tracker
│   ├── data_processing_pipeline.md     # Data preparation documentation
│   ├── plant-pplm-methodology.md       # Plant-PPLM methodology design
│   ├── lit_review/                     # Literature review templates
│   ├── ppi_head_retraining/            # Head retraining documentation & results
│   │   ├── ppi_head_retraining_methodology.md
│   │   ├── ppi_head_retraining_walkthrough.md
│   │   ├── benchmark_analysis_pretrained_pplm_deeparappi.md
│   │   └── benchmark_analysis_pretrained_pplm_esmarappi.md
│   └── zero_shot_benchmarking/         # Zero-shot benchmarking reports
│       ├── benchmark_analysis_deeparappi_vs_pplm.md
│       └── benchmark_analysis_esmarappi_vs_pplm.md
│
└── PPLM_NSCC_A100.yml                  # A100-optimized Conda environment
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

Across all **105,875 held-out test pairs** in the DeepAraPPI suite:

| Task | Test Dataset | Difficulty | DeepAraPPI Baseline (AUPRC) | PPLM Zero-Shot (AUPRC) | PPLM AUROC | PPLM Specificity |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Task 2** | `c2_ppi_sample_DeepAraPPI.txt` (66,055 pairs) | Medium (One unseen) | **0.8970** | **0.5738** | **0.8828** | 99.57% |
| **Task 3** | `c3_ppi_sample_DeepAraPPI.txt` (33,099 pairs) | High (Both unseen) | **0.8250** (Seq-only: 0.4810) | **0.5525** | **0.8710** | 99.50% |
| **Task 4** | `all_rice_positive_negative_DeepAraPPI.txt` (6,721 pairs) | Cross-Species (Rice) | **0.3050** | **0.4297** | **0.7561** | 98.35% |

---

### Phase 4: ESMAraPPI Benchmark Completed ✓

Across all **46,310 held-out test pairs** in the ESMAraPPI suite (with 40% sequence redundancy filtering):

| Task | Test Dataset | Difficulty | ESMAraPPI Baseline (AUPRC) | TAGPPI AF2 (AUPRC) | PPLM Zero-Shot (AUPRC) | PPLM AUROC | PPLM Specificity |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Task C2** | `c2Pred.txt` (37,444 pairs) | Medium (One unseen) | **0.8340** | **0.7000** | **0.5092** | **0.8563** | 99.55% |
| **Task C3** | `c3Pred.txt` (8,866 pairs) | High (Both unseen) | **0.8100** | **0.5540** | **0.5610** | **0.8657** | 99.45% |

### Key Benchmark Takeaways:
1. **New SOTA on Rice Cross-Species Transfer (Task 4):** PPLM (**0.4297**) outperforms ARACoFusion (**0.3519**), DeepAraPPI (**0.3050**), and ESMAraPPI on Rice (**0.2938**).
2. **Beats AlphaFold2-Based TAGPPI on Hard Unseen Proteins (ESMAraPPI C3):** Zero-shot PPLM (**0.5610**) surpasses TAGPPI (**0.5540**), DeepAraPPI RCNN (**0.3310**), PIPR (**0.3870**), and RAPPPID (**0.3710**).
3. **Exceptional Precision & Specificity Across All Datasets:** Specificity $>99.4\%$ and precision $>82.7\%$ across all **152,185 held-out test pairs** (105,875 DeepAraPPI + 46,310 ESMAraPPI).

---

### Phase 5: PPI Head Retraining (Completed on Both Suites) ✓

Retrained the PPLM-PPI MLP classifier head from random initialisation on plant C1 training sets with a frozen 650M PPLM backbone. See [`docs/ppi_head_retraining/ppi_head_retraining_methodology.md`](ppi_head_retraining/ppi_head_retraining_methodology.md) for full technical methodology.

- [x] **DeepAraPPI C1 Training (31,284 pairs):** 10-fold CV (mean AUPRC 0.8965 / max 0.8954)
  * C2 Held-Out (66,055 pairs): **0.8738 AUPRC** (vs. Zero-Shot 0.5738, +52.3%; DeepAraPPI: 0.8970)
  * C3 Held-Out (33,099 pairs): **0.8118 AUPRC** (vs. Zero-Shot 0.5525, +46.9%; DeepAraPPI: 0.8250)
  * Rice Transfer (6,721 pairs): **0.3555 AUPRC** (vs. Zero-Shot 0.4297; DeepAraPPI: 0.3050)
  * Full Report: [`docs/ppi_head_retraining/benchmark_analysis_pretrained_pplm_deeparappi.md`](ppi_head_retraining/benchmark_analysis_pretrained_pplm_deeparappi.md)

- [x] **ESMAraPPI C1 Training (38,709 pairs, 40% redundancy filtered):** 10-fold CV (mean AUPRC 0.8776 / max 0.8768)
  * Task C2 (37,444 pairs): **0.8408 AUPRC** (vs. Zero-Shot 0.5092, +65.1%; **Beats ESMAraPPI: 0.8340**)
  * Task C3 (8,866 pairs): **0.8103 AUPRC** (★ **New SOTA**; vs. Zero-Shot 0.5610, +44.4%; **Beats ESMAraPPI: 0.8100**, **Beats ARACoFusion: 0.8066**, **Beats DeepAraPPI: 0.7850**)
  * Full Report: [`docs/ppi_head_retraining/benchmark_analysis_pretrained_pplm_esmarappi.md`](ppi_head_retraining/benchmark_analysis_pretrained_pplm_esmarappi.md)

---

## Next Steps

1. **Focal loss ablation** — Implement focal loss ($\gamma=2$, $\alpha=0.25$) as an alternative to BCE to address 1:10 class imbalance.
2. **Fine-tuning from human weights** — Test transfer learning by initialising from the human-trained checkpoint (`ppi_models.pkl`) rather than random weights.
3. **LoRA backbone fine-tuning** — Apply Low-Rank Adaptation ($r=8$, $\alpha=16$) to PPLM cross-attention layers.
4. **Joint plant model training** — Merge DeepAraPPI C1 + ESMAraPPI C1 (69,993 pairs) for a unified plant foundation classifier.
5. **Auxiliary feature fusion** — GO terms and domain interaction priors (per `plant-pplm-methodology.md` §2.3).