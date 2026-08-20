# Technical Summary: PPLM Benchmarking for Plant Protein-Protein Interaction

## Project Overview

**Objective:** Benchmark the pre-trained PPLM (Paired Sequence Language Model) on Arabidopsis thaliana protein-protein interaction (PPI) datasets to evaluate cross-species generalization from animal/human to plant proteins.

**Key Challenge:** PPLM was trained on human, yeast, E. coli, C. elegans, D. melanogaster, and mouse PPI data. The goal is to test its performance on plant (Arabidopsis) PPI without fine-tuning.

---

## Project Architecture
AIS5281/
├── PPLM/ # Original PPLM codebase
│ ├── pplm/ # Core PPLM model
│ │ ├── pplm.py # Main PPLM architecture (33-layer Transformer)
│ │ ├── modules.py # Attention layers, LM heads
│ │ ├── data.py # Alphabet, batch converter, FASTA parsing
│ │ └── ...
│ ├── pplm_ppi/ # PPI classifier
│ │ └── model.py # PPLM_PPI: 5-layer MLP on extracted features
│ ├── weights/
│ │ ├── pplm_t33_650M.pt # Pre-trained PPLM backbone (650M params)
│ │ └── ppi_models.pkl # 10 PPI classifier weight sets
│ ├── run_pplm-ppi.py # Original single-pair prediction script
│ ├── batch_predict.py # [NEW] Batch prediction script (user-created)
│ └── evaluate_pplm.py # [NEW] Evaluation metrics script (user-created)
│
├── data/ # Data and sequence databases
│ ├── DeepAraPPI/ # Protein interaction datasets
│ │ ├── total_positive_negative_samples_DeepAraPPI.txt # 130K pairs
│ │ ├── c1_ppi_sample_DeepAraPPI.txt # ~1K pairs
│ │ ├── c2_ppi_sample_DeepAraPPI.txt # ~3K pairs
│ │ ├── c3_ppi_sample_DeepAraPPI.txt # ~1.5K pairs
│ │ └── all_rice_positive_negative_DeepAraPPI.txt # Rice PPI
│ └── uniprot_final.pkl # [FINAL] Complete sequence database (~110K IDs)
│
├── results/ # Benchmark outputs
│ ├── c1_test_scores.csv # Test predictions
│ ├── c1_test_metrics.txt # Test metrics
│ └── ...
│
└── technical_summary.md # Comprehensive documentation

---

## Data Format

**Input Data (DeepAraPPI):**
Protein1 Protein2 label
Q9LUI9 F4IDF1 0
F4IHD3 Q9M8S6 0
...


- **Format:** Tab-separated, 3 columns
- **ID types:** Mixed - Swiss-Prot (Q9LUI9) and TrEMBL (A0A178VKB9)
- **Label:** 0 = non-interacting, 1 = interacting

**Sequence Database:**
- **Original:** 54,646 UniProt sequences (Swiss-Prot + TrEMBL)
- **Final:** ~110,000 ID mappings (109,292 from UniProt + 348 from UniParc)
- **Coverage:** 97.6% → 100% after UniParc merge

---

## Benchmarking Progress

### Phase 1: Environment Setup ✓
- [x] Verified PPLM installation on HPC (CentOS7 + PBS)
- [x] Confirmed GPU availability (CUDA enabled)
- [x] Downloaded model weights (pplm_t33_650M.pt, ppi_models.pkl)
- [x] Created conda environment with required dependencies

### Phase 2: Sequence Database Preparation ✓
- [x] Downloaded Arabidopsis thaliana protein sequences from UniProt
- [x] Identified ID format mismatch problem:
  - Data uses simple IDs: `Q9LUI9`, `F4IDF1`
  - FASTA has full IDs: `sp|Q9LUI9|PIAL1_ARATH`
- [x] Created regex-based ID extraction in PKL generation
- [x] Discovered 2.4% missing sequences (TrEMBL/UniParc entries)
- [x] Downloaded missing 348 sequences from UniParc
- [x] Created ID mapping JSON file
- [x] Merged all sequences into final PKL database
- [x] **Verified 100% coverage** on DeepAraPPI dataset

### Phase 3: Batch Prediction Script ✓
- [x] Created `batch_predict.py` with proper path handling:
  ```python
  SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
  PPLM_DIR = os.path.dirname(SCRIPT_DIR)  # Goes up to AraPPLM/PPLM
  sys.path.insert(0, PPLM_DIR)
 - [x] Implemented PPLM feature extraction (mean/max attention + embeddings)
 - [x] Implemented PPI classification (ensemble of 10 weight sets)
 - [x] Added batch processing for GPU efficiency
 - [x] Tested on small dataset (c1_ppi_sample, ~1K pairs)
 
### Phase 4: Evaluation Framework ✓
- [x] Created evaluate_pplm.py with comprehensive metrics:
  - AUROC: Area under ROC curve
  - AUPRC: Area under precision-recall curve
  - FPR@95TPR: False positive rate at 95% true positive rate
  - Threshold-based metrics: Accuracy, sensitivity, specificity, F1
- [x] Output saved to CSV and TXT formats

## Key Technical Details
### PPLM Architecture
Backbone: 33-layer Transformer with cross-chain attention
Input: Paired protein sequences with inter-chain mask
Features extracted:
Mean/Max inter-protein attention (660 dims)
Mean/Max intra-protein attention (660 × 2 dims)
Mean/Max embeddings (1280 × 2 dims)
Classifier: 5-layer MLP (5 × 660 → 1024 → 512 → 256 → 128 → 1)
Ensemble: 10 weight sets (5 mean + 5 max), symmetric averaging