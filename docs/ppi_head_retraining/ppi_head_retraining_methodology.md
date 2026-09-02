# PPLM-PPI Prediction Head Retraining: Technical Methodology

## 1. Motivation & Objective

The PPLM-PPI prediction head was originally trained on human/model-organism PPI data from the D-SCRIPT dataset (Chen et al., 2019). Our zero-shot benchmarking showed PPLM achieves strong precision (>82.7%) and specificity (>99.4%) on plant PPI data, but AUPRC lags behind supervised plant-specific models (e.g., 0.5738 vs. 0.8970 for DeepAraPPI Task 2). 

**Objective:** Retrain the PPLM-PPI MLP classifier head on plant *Arabidopsis thaliana* C1 training data while keeping the PPLM backbone frozen, following the original PPLM training protocol to enable a fair comparison.

**Scientific rationale:** The PPLM backbone (ESM2-initialized, 33-layer Transformer with paired cross-attention) encodes generalizable biophysical protein-pair representations. By retraining only the MLP prediction head (~3.9M parameters out of 650M total), we test whether the backbone's learned representations are sufficient for plant PPI classification when given plant-specific supervision.

---

## 2. Training Data

### 2.1 DeepAraPPI C1 Dataset

| Property | Value |
|----------|-------|
| **Source** | Zheng et al., "DeepAraPPI" (The Plant Journal, 2023) |
| **File** | `data/DeepAraPPI/c1_ppi_sample_DeepAraPPI.txt` |
| **Total pairs** | 31,284 |
| **Positive pairs** | 2,844 (high-quality PPIs from IntAct, MIscore ≥ 0.45) |
| **Negative pairs** | 28,440 (random non-interacting pairs) |
| **Pos:Neg ratio** | 1:10 |
| **Species** | *Arabidopsis thaliana* |
| **Format** | TSV: `Protein1\tProtein2\tlabel` (UniProt accession IDs) |

### 2.2 ESMAraPPI C1 Dataset

| Property | Value |
|----------|-------|
| **Source** | Zhou et al., "ESMAraPPI" (Plant Methods, 2023) |
| **File** | `data/ESMAraPPI/c1_Train.txt` |
| **Total pairs** | 38,709 |
| **Positive pairs** | 3,519 (IntAct, "direct interaction" or "physical association", MIscore ≥ 0.45) |
| **Negative pairs** | 35,190 (40% sequence identity filtering on negative candidates) |
| **Pos:Neg ratio** | 1:10 |
| **Species** | *Arabidopsis thaliana* |

### 2.3 Data Splitting Strategy

Both datasets follow Park & Marcotte's (2012) pair-input evaluation scheme:

- **C1** = Training set (used for model training via 10-fold cross-validation)
- **C2** = Test set where exactly *one* protein in each pair appeared in C1 (medium difficulty)
- **C3** = Test set where *both* proteins are completely unseen in C1 (high difficulty, zero-shot)
- **Rice** = Cross-species transfer test (DeepAraPPI only)

Each dataset is trained **independently** to follow the same evaluation protocol used by the published baseline models (DeepAraPPI, ESMAraPPI, ARACoFusion), ensuring directly comparable results.

---

## 3. Model Architecture

### 3.1 PPLM Backbone (Frozen)

The backbone is a 33-layer Transformer (650M parameters) initialised from ESM-2 (`esm2_t33_650M_UR50D`) with paired cross-attention:

- **Intra-chain attention**: Standard self-attention with Rotary Position Embeddings (RoPE)
- **Inter-chain attention**: Cross-attention between protein A and B (no positional encoding, treating inter-chain contacts as position-independent)
- **Input**: Two protein sequences concatenated with separator tokens, combined cropping at $L_A + L_B \le 1020$

The backbone produces:
1. **Attention tensor**: `[33 layers × 20 heads, L, L]` — capturing pairwise residue-residue interaction patterns
2. **Per-residue embeddings**: `[L, 1280]` from layer 33 — encoding residue-level structural context

### 3.2 Feature Extraction (Stage 1)

For each protein pair, the backbone's raw outputs are **pooled** into fixed-size feature vectors:

| Feature | Dimension | Description |
|---------|-----------|-------------|
| `inter_attn` | `[660]` | Cross-attention between A and B (symmetrized: `(A→B + B→A) / 2`), pooled across spatial dims |
| `attn_AA` | `[660]` | Self-attention within protein A |
| `attn_BB` | `[660]` | Self-attention within protein B |
| `embed_A` | `[1280]` | Per-residue embedding for A, pooled across residues |
| `embed_B` | `[1280]` | Per-residue embedding for B, pooled across residues |

Dimension 660 = 33 layers × 20 heads. Three pooling strategies are computed: **mean**, **max**, and **min**, yielding 15 feature vectors per pair.

These are saved as individual `.pkl` files in `features/{Dataset}_C1/{ProtA}@{ProtB}.pkl`.

### 3.3 PPI Classifier MLP (Trainable — Stage 2)

The classifier is a 5-layer MLP operating on features from **one pooling mode at a time**:

```
  inter_attn [660] ─→ Linear(660→660)  ─┐
  attn_AA    [660] ─→ Linear(660→660)  ─┤  (shared weights for A & B)
  attn_BB    [660] ─→ Linear(660→660)  ─┤
  embed_A    [1280]─→ Linear(1280→660) ─┤
  embed_B    [1280]─→ Linear(1280→660) ─┘
                                         │
  Concatenate ──────────────────── [3300] │
                                         │
  FC 3300→1024 + LayerNorm + ReLU       │
  FC 1024→512  + LayerNorm + ReLU       │
  FC 512→256   + LayerNorm + ReLU       │
  FC 256→128   + LayerNorm + ReLU       │
  FC 128→1     + Sigmoid ─────── P(interaction)
```

**Total trainable parameters**: ~3.9M per model instance.

> **Design Decision**: The architecture is kept **identical** to the original PPLM-PPI model (`pplm_ppi/model.py`) to ensure that any performance difference is attributable purely to the training data (plant vs. human), not to architectural changes. This enables a controlled comparison.

---

## 4. Training Protocol

### 4.1 Hyperparameters

All hyperparameters match the original PPLM-PPI paper:

| Hyperparameter | Value | Rationale |
|----------------|-------|-----------|
| **Optimizer** | AdamW | Following PPLM paper |
| **Learning rate** | 5 × 10⁻⁵ | Following PPLM paper |
| **Weight decay** | 0.01 | Following PPLM paper |
| **Batch size** | 512 | Following PPLM paper |
| **Epochs** | 15 per fold | Following PPLM paper |
| **Loss function** | Binary Cross-Entropy (BCE) | See §4.3 |
| **Cross-validation** | 10-fold stratified | Maintains 1:10 class ratio per fold |
| **Pooling modes** | mean, max | Following PPLM paper (min dropped for efficiency) |
| **Data augmentation** | 50% random A↔B swapping | Enforces interaction symmetry |
| **Random seed** | 26240761 | Deterministic reproducibility |

### 4.2 Cross-Validation Procedure

1. All C1 training pairs are split into 10 stratified folds using `sklearn.StratifiedKFold`, preserving the 1:10 positive-to-negative ratio in each fold
2. For each fold: 9 folds for training, 1 for validation
3. A fresh model is initialised for each fold (no weight transfer between folds)
4. Per-epoch checkpoints are saved with validation metrics
5. The process is repeated for each pooling mode (mean, max)

Total training runs per dataset: **10 folds × 2 pooling modes = 20 training runs**

### 4.3 Loss Function Decision

> **Design Decision — BCE vs. Focal Loss**: The original PPLM paper uses standard Binary Cross-Entropy (BCE) loss. Our plant-pplm-methodology document (§3.2) proposes focal loss ($\gamma=2$) with label smoothing to handle the 1:10 class imbalance more aggressively. 
>
> We chose to start with **plain BCE** for the initial retraining to establish a clean controlled baseline — any performance gain over zero-shot PPLM can then be attributed purely to plant-specific training data, not to loss function engineering.
>
> **Future work**: Focal loss ($\gamma=2$, $\alpha=0.25$) and class-weighted BCE ($w_{pos}=10.0$) will be implemented as ablation experiments to quantify their effect on recall of the minority positive class.

### 4.4 Model Selection & Ensemble

After all training runs complete:

1. **Ranking**: All 150 checkpoints per pooling mode (10 folds × 15 epochs) are ranked by **validation AUPRC**
2. **Selection**: Top-5 checkpoints per pooling mode are selected
3. **Packaging**: The 10 model state_dicts (5 mean + 5 max) are saved as `ppi_plant_models.pkl`

### 4.5 Inference Protocol

At test time, for each protein pair:
1. Run PPLM backbone to extract features (same as Stage 1)
2. For each of the 10 ensemble classifiers:
   - Run forward pass: A→B prediction
   - Run swapped pass: B→A prediction
   - Average the two (symmetric pair averaging)
3. Final score = mean of all 10 averaged predictions

This is identical to the original PPLM-PPI inference protocol, ensuring a fair comparison.

---

## 5. Training Infrastructure

### 5.1 Hardware

| Resource | Specification |
|----------|---------------|
| **Platform** | NSCC ASPIRE 2A (National Supercomputing Centre, Singapore) |
| **GPU** | NVIDIA A100-SXM4-40GB |
| **CPU** | 8 cores per job |
| **RAM** | 64 GB |
| **CUDA** | 12.8, Driver 570.124.06 |

### 5.2 Software Environment

| Package | Version |
|---------|---------|
| **Python** | 3.10 |
| **PyTorch** | ≥2.0.0 (CUDA 12.4) |
| **fair-esm** | ≥2.0.0 |
| **scikit-learn** | (latest) |
| **Conda environment** | `PPLM_NSCC_A100` (spec: `PPLM_NSCC_A100.yml`) |

### 5.3 Execution Times (DeepAraPPI C1, 31,284 pairs)

| Stage | Description | Wall Time |
|-------|-------------|-----------|
| **Stage 1** | Feature extraction (frozen backbone inference) | ~4–6 hours |
| **Stage 2** | MLP training (10-fold CV × 2 pooling modes) | ~1–2 hours |
| **Stage 3** | Model selection | <1 minute |
| **Stage 4** | Test-set evaluation (C2: 66k + C3: 33k + Rice: 6.7k) | ~4–8 hours |
| **Total** | End-to-end | ~10–16 hours |

---

## 6. Results

Training results and test-set evaluation are documented separately:

- **DeepAraPPI pretrained results:** [`benchmark_analysis_pretrained_pplm_deeparappi.md`](benchmark_analysis_pretrained_pplm_deeparappi.md)
- **Zero-shot baseline results:** [`benchmark_analysis_deeparappi_vs_pplm.md`](benchmark_analysis_deeparappi_vs_pplm.md)

---

## 7. Pipeline Scripts

### 7.1 Script Reference

| Script | Stage | Description |
|--------|-------|-------------|
| [`extract_features.py`](../scripts/training/extract_features.py) | 1 | Frozen PPLM backbone → pooled feature `.pkl` files |
| [`dataset.py`](../scripts/training/dataset.py) | 2 | PyTorch Dataset for loading `.pkl` features with A↔B augmentation |
| [`ppi_model.py`](../scripts/training/ppi_model.py) | 2 | MLP classifier architecture + evaluation metrics |
| [`train_ppi_head.py`](../scripts/training/train_ppi_head.py) | 2 | K-fold stratified CV training loop |
| [`select_top_models.py`](../scripts/training/select_top_models.py) | 3 | Rank by AUPRC, package top-5 into ensemble |
| [`test_ppi_head.py`](../scripts/training/test_ppi_head.py) | 4 | Ensemble inference on held-out test sets |
| [`evaluate_pplm.py`](../scripts/benchmarking/evaluate_pplm.py) | 5 | Metrics computation + baseline comparison |

### 7.2 PBS Job Scripts

| Script | Purpose |
|--------|---------|
| [`run_train_deeparappi_nscc.pbs`](../scripts/run_train_deeparappi_nscc.pbs) | Full pipeline: Stages 1–5 on DeepAraPPI C1 |
| [`run_train_esmarappi_nscc.pbs`](../scripts/run_train_esmarappi_nscc.pbs) | Full pipeline: Stages 1–5 on ESMAraPPI C1 |
| [`run_resume_deeparappi_nscc.pbs`](../scripts/run_resume_deeparappi_nscc.pbs) | Resume from Stage 3 (model selection → evaluation) |

---

## 8. Future Directions

1. **Focal loss ablation** ($\gamma=2$, $\alpha=0.25$): Test whether focal loss improves recall for the minority positive class without sacrificing the exceptional specificity
2. **LoRA backbone fine-tuning** ($r=8$, $\alpha=16$): Apply Low-Rank Adaptation to the PPLM cross-attention layers to adapt the backbone representations to plant protein biophysics
3. **Joint training**: Train on merged DeepAraPPI C1 + ESMAraPPI C1 to produce a unified plant model
4. **Extended cross-species**: Evaluate rice, maize, and soybean transfer with the plant-trained model
5. **Auxiliary feature fusion**: Integrate GO term and domain interaction features as proposed in `plant-pplm-methodology.md` §2.3

---

## References

1. Liu, J., Chen, X. & Zhang, Y. PPLM: Paired Protein Language Model for Protein-Protein Interaction Prediction. *Nature Communications* (2026).
2. Zheng, J. et al. DeepAraPPI: A Benchmark Dataset and Comprehensive Study for Arabidopsis Protein-Protein Interaction Prediction. *The Plant Journal* (2023).
3. Zhou, K. et al. Pre-trained protein language model sheds new light on the prediction of Arabidopsis protein-protein interactions. *Plant Methods* 19:141 (2023).
4. Park, Y. & Marcotte, E.M. Flaws in evaluation schemes for pair-input computational predictions. *Nature Methods* 9:1134–1136 (2012).
