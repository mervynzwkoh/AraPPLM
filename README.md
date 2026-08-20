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
- [x] Sequence database preparation (100% coverage on DeepAraPPI)
- [x] Batch prediction pipeline
- [ ] Zero-shot benchmarking on DeepAraPPI C1/C2/C3
- [ ] Evaluation framework
- [ ] Cross-species benchmarking (Rice, Maize)
- [ ] Domain-adaptive fine-tuning
- [ ] Auxiliary feature fusion (GO terms, domain interactions)
- [ ] Full ablation suite with statistical testing

## Quick Start

### Prerequisites
- Python 3.8+, PyTorch with CUDA, BioPython

### Setup
```bash
git clone --recurse-submodules https://github.com/YOUR_USERNAME/AIS5281-Plant-PPLM.git
cd AIS5281-Plant-PPLM

# Download PPLM model weights (~650MB)
cd PPLM/weights && bash download_weights.sh && cd ../..

# Build sequence database
# 1. Download Arabidopsis proteome from UniProt (uniprot_arabidopsis.fasta)
# 2. Place in data/
cd scripts/data_preparation
python uniprot_index.py
python merge_uniparc.py
```

## Project Structure
├── PPLM/                  # Original PPLM (submodule)
├── scripts/               # Custom prediction & data prep scripts
├── data/                  # Benchmark datasets & sequence data
├── results/               # Benchmark outputs
└── docs/                  # Technical documentation

## References
- PPLM: Liu, Chen & Zhang, Nature Communications 2026
- DeepAraPPI: Zheng et al. 2023
