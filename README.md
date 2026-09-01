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
- [x] PPI head retraining pipeline (Stages 1–2 complete on DeepAraPPI C1)
- [ ] PPI head test-set evaluation (Stages 3–5 pending)
- [ ] ESMAraPPI C1 training
- [ ] Domain-adaptive fine-tuning (LoRA)
- [ ] Focal loss / class-weighted loss ablation
- [ ] Auxiliary feature fusion (GO terms, domain interactions)

## Key Results

### Zero-Shot Benchmarking (No Plant-Specific Training)

| Task | Difficulty | PPLM (AUPRC) | Best Supervised Baseline |
|------|-----------|--------------|--------------------------|
| DeepAraPPI C2 | Medium (one unseen) | 0.5738 | DeepAraPPI: 0.8970 |
| DeepAraPPI C3 | High (both unseen) | 0.5525 | DeepAraPPI: 0.8250 |
| DeepAraPPI Rice | Cross-species | **0.4297** | DeepAraPPI: 0.3050 |
| ESMAraPPI C3 | High (both unseen) | **0.5610** | TAGPPI: 0.5540 |

**Highlights:** New SOTA on Rice cross-species transfer; beats AlphaFold2-based TAGPPI on hard unseen proteins.

### PPI Head Retraining (DeepAraPPI C1, 10-fold CV)

| Pooling Mode | Avg Validation AUPRC |
|-------------|---------------------|
| Mean | **0.8965** |
| Max | **0.8954** |

Test-set evaluation (C2/C3/Rice) pending — see `docs/ppi_head_retraining_methodology.md`.

## Quick Start

### Prerequisites
- Python 3.10, PyTorch ≥2.0.0 with CUDA, BioPython, scikit-learn

### Setup
```bash
git clone --recurse-submodules https://github.com/YOUR_USERNAME/AIS5281-Plant-PPLM.git
cd AIS5281-Plant-PPLM

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
# Submit PBS job on NSCC
qsub scripts/run_all_benchmarks_nscc.pbs
```

### Run PPI Head Retraining
```bash
# Full pipeline: extract features → train → select → test → evaluate
qsub scripts/run_train_deeparappi_nscc.pbs

# Or resume from model selection if Stages 1–2 done
qsub scripts/run_resume_deeparappi_nscc.pbs
```

## Project Structure
```
├── PPLM/                  # Original PPLM (submodule)
├── scripts/
│   ├── training/          # PPI head retraining pipeline (6 scripts)
│   ├── benchmarking/      # Zero-shot benchmark scripts
│   └── data_preparation/  # Sequence DB builders
├── data/                  # Benchmark datasets & sequence data
├── features/              # Pre-extracted PPLM backbone features
├── models/                # Trained model checkpoints
├── results/               # Benchmark & evaluation outputs
├── logs/                  # HPC PBS job logs
│   ├── benchmarking/
│   └── training/
└── docs/                  # Technical documentation
```

## Documentation

| Document | Description |
|----------|-------------|
| [technical_summary.md](docs/technical_summary.md) | Architecture, benchmarking results & progress tracker |
| [ppi_head_retraining_methodology.md](docs/ppi_head_retraining_methodology.md) | Retraining technical methodology, design decisions & CV results |
| [PPLM-PPI_pretraining_walkthrough.md](docs/PPLM-PPI_pretraining_walkthrough.md) | Step-by-step run instructions for the training pipeline |
| [data_processing_pipeline.md](docs/data_processing_pipeline.md) | Sequence retrieval & database construction |
| [plant-pplm-methodology.md](docs/plant-pplm-methodology.md) | Proposed fine-tuning methodology (LoRA, focal loss) |
| [benchmark_analysis_deeparappi_vs_pplm.md](docs/benchmark_analysis_deeparappi_vs_pplm.md) | DeepAraPPI benchmark analysis |
| [benchmark_analysis_esmarappi_vs_pplm.md](docs/benchmark_analysis_esmarappi_vs_pplm.md) | ESMAraPPI benchmark analysis |

## References
- PPLM: Liu, Chen & Zhang, Nature Communications 2026
- DeepAraPPI: Zheng et al., The Plant Journal 2023
- ESMAraPPI: Zhou et al., Plant Methods 2023
- ARACoFusion: Sarkar & Sarkar, bioRxiv 2026
