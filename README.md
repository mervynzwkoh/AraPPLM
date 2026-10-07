# AIS5281: Plant-Adapted PPLM for Protein-Protein Interaction Prediction

## Overview

Benchmarking and fine-tuning [PPLM](https://github.com/junliu621/PPLM)
(Paired Protein Language Model, 650M params) for plant protein-protein
interaction (PPI) prediction.

**Key research question:** Can a paired protein language model trained on
human/model-organism PPI data generalize to plant proteins, and how much
does plant-specific domain adaptation improve performance?

## Project Status

- [x] Literature review (18 plant PPI papers)
- [x] Sequence database preparation (100% coverage on all datasets)
- [x] Batch prediction pipeline (combined cropping, ensemble inference)
- [x] Zero-shot benchmarking on DeepAraPPI C1/C2/C3 + Rice
- [x] Zero-shot benchmarking on ESMAraPPI C2/C3
- [x] Evaluation framework (AUPRC, AUROC, F1, MCC, Specificity)
- [x] PPI head retraining on DeepAraPPI C1 & test-set evaluation (C2, C3, Rice)
- [x] PPI head retraining on ESMAraPPI C1 & test-set evaluation (C2, C3)
- [ ] Domain-adaptive fine-tuning (LoRA on backbone)
- [ ] Focal loss / class-weighted loss ablation
- [ ] Auxiliary feature fusion (GO terms, domain interactions)

## Key Results

### Zero-Shot Benchmarking (No Plant-Specific Training)

| Task | Difficulty | PPLM Zero-Shot (AUPRC) | Best Supervised Baseline |
|------|-----------|------------------------|--------------------------|
| DeepAraPPI C2 | Medium (one unseen) | 0.5738 | DeepAraPPI: 0.8970 |
| DeepAraPPI C3 | High (both unseen) | 0.5525 | DeepAraPPI: 0.8250 |
| DeepAraPPI Rice | Cross-species | **0.4297** ★ | DeepAraPPI: 0.3050 |
| ESMAraPPI C3 | High (both unseen) | **0.5610** | TAGPPI: 0.5540 |

**Highlights:** Zero-shot PPLM set a new SOTA on Rice cross-species transfer (0.4297) and beat AlphaFold2-based TAGPPI on hard unseen proteins without plant training.

### Plant-Pretrained PPI Head (Frozen Backbone + Plant C1 MLP Retraining)

| Benchmark Suite | Task / Partition | Difficulty | Zero-Shot AUPRC | Pretrained AUPRC | Best Published Baseline |
|:---|:---|:---|:---:|:---:|:---:|
| **DeepAraPPI** | Task 2 (C2) | Medium (1 Unseen) | 0.5738 | **0.8738** | DeepAraPPI: 0.8970 |
| **DeepAraPPI** | Task 3 (C3) | High (Both Unseen) | 0.5525 | **0.8118** | DeepAraPPI: 0.8250 |
| **DeepAraPPI** | Task 4 (Rice) | Cross-Species | **0.4297** | 0.3555 | ARACoFusion: 0.3519 |
| **ESMAraPPI** | Task C2 | Medium (1 Unseen) | 0.5092 | **0.8408** | ESMAraPPI: 0.8340 |
| **ESMAraPPI** | Task C3 | High (Both Unseen) | 0.5610 | **0.8103** ★ | ESMAraPPI: 0.8100 |

**Highlights:**
* **New SOTA on ESMAraPPI Task C3 (0.8103 AUPRC):** Surpasses the native ESMAraPPI model (0.8100), ARACoFusion (0.8066), and DeepAraPPI (0.7850).
* **Near-Parity with DeepAraPPI on Sequence-Only Inputs:** 0.8738 on C2 and 0.8118 on C3 without requiring Gene Ontology or domain graph annotations.

## Quick Start

### Prerequisites
- Python 3.10, PyTorch ≥2.0.0 with CUDA, BioPython, scikit-learn

### Setup
```bash
git clone --recurse-submodules https://github.com/mervynzwkoh/AraPPLM.git
cd AraPPLM

# Create conda environment
conda env create -f PPLM_NSCC_A100.yml
conda activate PPLM_NSCC_A100

# Download PPLM model weights (~650MB)
cd PPLM/weights && bash download_weights.sh && cd ../..

# Build sequence database (see docs/data_processing_pipeline.md)
python scripts/data_preparation/build_sequence_db.py
```

### Run Benchmarking (Zero-Shot)
```bash
# Submit PBS jobs on NSCC
qsub scripts/hpc/run_all_benchmarks_nscc.pbs
qsub scripts/hpc/run_esmarappi_benchmarks_nscc.pbs
```

### Run PPI Head Retraining
```bash
# End-to-end retraining on DeepAraPPI C1
qsub scripts/hpc/run_train_deeparappi_nscc.pbs

# End-to-end retraining on ESMAraPPI C1
qsub scripts/hpc/run_train_esmarappi_nscc.pbs
```

## Project Structure
```
├── PPLM/                  # Original PPLM (submodule)
├── checksums/             # SHA-256 provenance checksums (Phase 1 & Phase 2)
├── manifests/             # Machine-readable manifests and progress tracking
├── metadata/              # Source registry, licenses, and BibTeX citations
├── scripts/
│   ├── hpc/               # PBS batch job scripts for NSCC ASPIRE 2A
│   ├── phase2/            # Streaming ETL evidence normalization pipeline
│   ├── training/          # PPI head retraining pipeline (6 scripts)
│   ├── benchmarking/      # Zero-shot benchmark scripts
│   └── data_preparation/  # Sequence DB builders & verifiers
├── data/                  # Benchmark datasets, interim tables, & sequence data
├── features/              # Pre-extracted PPLM backbone features
├── models/                # Trained model checkpoints
├── results/               # Benchmark & evaluation outputs
├── logs/                  # HPC PBS job logs
├── tests/                 # Unit test suite (pytest)
└── docs/                  # Technical documentation & scientific audit reports
```

## Documentation

| Document | Description |
|----------|-------------|
| [technical_summary.md](docs/technical_summary.md) | Architecture, benchmarking results & progress tracker |
| [project_overview.md](docs/project_overview.md) | High-level synthesis of benchmark results and methodology |
| [data_processing_pipeline.md](docs/data_processing_pipeline.md) | Sequence retrieval & database construction |
| [plant-pplm-methodology.md](docs/plant-pplm-methodology.md) | Proposed fine-tuning methodology (LoRA, focal loss) |
| **Dataset Construction (Phase 1 & Phase 2)** | |
| [P1_source_freeze_report.md](docs/dataset_construction/P1_source_freeze_report.md) | Phase 1 raw data freeze & provenance registry |
| [P2_evidence_normalization_report.md](docs/dataset_construction/P2_evidence_normalization_report.md) | Phase 2 normalized canonical Parquet warehouse audit |
| [ara_data_breakdown.md](docs/dataset_construction/ara_data_breakdown.md) | Arabidopsis source-by-source interaction breakdown |
| [rice_data_breakdown.md](docs/dataset_construction/rice_data_breakdown.md) | Rice source-by-source interaction breakdown |
| **Head Retraining** | |
| [ppi_head_retraining_methodology.md](docs/ppi_head_retraining/ppi_head_retraining_methodology.md) | Retraining technical methodology, design decisions & training protocol |
| [ppi_head_retraining_walkthrough.md](docs/ppi_head_retraining/ppi_head_retraining_walkthrough.md) | Step-by-step run instructions for HPC training |
| [benchmark_analysis_pretrained_pplm_deeparappi.md](docs/ppi_head_retraining/benchmark_analysis_pretrained_pplm_deeparappi.md) | DeepAraPPI pretrained results report (C2, C3, Rice) |
| [benchmark_analysis_pretrained_pplm_esmarappi.md](docs/ppi_head_retraining/benchmark_analysis_pretrained_pplm_esmarappi.md) | ESMAraPPI pretrained results report (C2, C3) |
| **Zero-Shot Benchmarking** | |
| [benchmark_analysis_deeparappi_vs_pplm.md](docs/zero_shot_benchmarking/benchmark_analysis_deeparappi_vs_pplm.md) | DeepAraPPI zero-shot report |
| [benchmark_analysis_esmarappi_vs_pplm.md](docs/zero_shot_benchmarking/benchmark_analysis_esmarappi_vs_pplm.md) | ESMAraPPI zero-shot report |

## References
- PPLM: Liu, Chen & Zhang, Nature Communications 2026
- DeepAraPPI: Zheng et al., The Plant Journal 2023
- ESMAraPPI: Zhou et al., Plant Methods 2023
- ARACoFusion: Sarkar & Sarkar, bioRxiv 2026
