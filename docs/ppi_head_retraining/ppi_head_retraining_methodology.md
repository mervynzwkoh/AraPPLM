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
| **File** | `data/legacy_benchmarks/DeepAraPPI/c1_ppi_sample_DeepAraPPI.txt` |
| **Total pairs** | 31,284 |
| **Positive pairs** | 2,844 (high-quality PPIs from BioGRID, DIP, IntAct, MINT, TAIR; HIPPIE score ≥ 0.72) |
| **Negative pairs** | 28,440 (random non-interacting pairs) |
| **Pos:Neg ratio** | 1:10 |
| **Species** | *Arabidopsis thaliana* |
| **Format** | TSV: `Protein1\tProtein2\tlabel` (UniProt accession IDs) |

### 2.2 ESMAraPPI C1 Dataset

| Property | Value |
|----------|-------|
| **Source** | Zhou et al., "ESMAraPPI" (Plant Methods, 2023) |
| **File** | `data/legacy_benchmarks/ESMAraPPI/c1Train.txt` |
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

Following the original PPLM training and evaluation setup (`ppi_training_code/`), `.pkl` features are pre-extracted for **both** C1 training pairs and held-out test partitions (`C2`, `C3`, and `Rice`) into `features/{Dataset}_{Split}/{ProtA}@{ProtB}.pkl`. This allows both training and held-out testing to run via fast batched `DataLoader` (`batch_size=512`).

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

All hyperparameters and determinism flags match the original PPLM-PPI training codebase (`ppi_training_code/training.PPI.single_gpu.cross_validation.py` and `run_training.sh`):

| Hyperparameter | Value | Rationale |
|----------------|-------|-----------|
| **Optimizer** | AdamW | Following PPLM paper |
| **Learning rate** | 5 × 10⁻⁵ | Following PPLM paper |
| **Weight decay** | 0.01 | Following PPLM paper |
| **Batch size** | 512 | Following PPLM paper (`num_workers=8`, `prefetch_factor=4`, `pin_memory=True`) |
| **Epochs** | 15 per fold | Following PPLM paper |
| **Loss function** | Binary Cross-Entropy (BCE) | See §4.3 |
| **Cross-validation** | 10-fold stratified pre-split files | Persists 1:10 class ratio per fold for independent `--cross_round 0..9` execution |
| **Pooling modes** | `mean`, `max`, `min` | All 3 pooling modes trained matching `run_training.sh` |
| **Data augmentation** | 50% random A↔B swapping (`order_mode="random"`) | Enforces interaction symmetry |
| **Random seed** | `73591024` (`SEEDS[1]`, `--seed_round 1`) | Author's `run_training.sh` seed + `torch.backends.cudnn.enabled = False` |

### 4.2 Cross-Validation Procedure

1. All C1 training pairs are partitioned once into 10 stratified folds using `sklearn.StratifiedKFold(n_splits=10, shuffle=True, random_state=73591024)` and persisted on disk as `{prefix}.cross_validation_{0..9}.trainset` and `.testset` files, preserving the 1:10 positive-to-negative ratio in every fold
2. Each fold (`--cross_round 0..9`) runs as an independent, resumable process (9 folds train, 1 fold validation; checkpoint resumption via `model.pkl`)
3. A fresh model is initialised for each fold (no weight transfer between folds)
4. Per-epoch checkpoints (`model_1.pkl`–`model_15.pkl`) and exact author-format `.recording` logs (`AUC_PR` at column 14) are saved
5. The process is repeated for all three pooling modes (`mean`, `max`, `min`)

Total training runs per dataset: **10 folds × 3 pooling modes = 30 training runs** (450 candidate checkpoints per dataset)

### 4.3 Loss Function Decision

> **Design Decision — BCE vs. Focal Loss**: The original PPLM paper uses standard Binary Cross-Entropy (BCE) loss. Our plant-pplm-methodology document (§3.2) proposes focal loss ($\gamma=2$) with label smoothing to handle the 1:10 class imbalance more aggressively. 
>
> We chose to start with **plain BCE** for the initial retraining to establish a clean controlled baseline — any performance gain over zero-shot PPLM can then be attributed purely to plant-specific training data, not to loss function engineering.
>
> **Future work**: Focal loss ($\gamma=2$, $\alpha=0.25$) and class-weighted BCE ($w_{pos}=10.0$) will be implemented as ablation experiments to quantify their effect on recall of the minority positive class.

### 4.4 Model Selection & Ensemble

After all 30 training runs complete:

1. **Ranking**: All 150 checkpoints per pooling mode (10 folds × 15 epochs) are ranked globally by **validation `AUC_PR`** (column 14 of `========== Accuarcy Epoch:` lines)
2. **Selection**: Top-5 checkpoints per pooling mode are selected and written to author-compatible path lists (`PPLM-PPI.inter_intra_attn_embed.{mode}.round_1.top5_list`)
3. **Packaging**: The 15 selected model `state_dict`s (`5 mean + 5 max + 5 min`) are also packaged into `ppi_plant_models.pkl`

### 4.5 Inference Protocol

Following `ppi_training_code/testing.PPI.single_gpu.multi-model.multi_pooling.py` and `run_testing.sh`, held-out evaluation runs over **pre-extracted `.pkl` test features** via `Dataset_PPI_Interact` and `DataLoader(batch_size=512)`:
1. For each test pair in the batch and each selected classifier in the ensemble:
   - Run forward pass: $A \to B$ prediction
   - Run swapped pass: $B \to A$ prediction (`order_mode="both"`)
   - Average the two (symmetric pair averaging)
2. Final score = mean across all models in the evaluated pooling combination:
   - **Single modes (5 models each):** `mean`, `max`, `min`
   - **Combined ensembles:** `combine_mean_max` (10 models), `combine_mean_max_min` (15 models), `combine_mean_min` (10 models)

An on-the-fly backbone fallback (`--seq_db`) remains supported if pre-extracted test features are unavailable.

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

### 5.3 Execution Times (DeepAraPPI Suite, Train + Test Pre-Extraction)

| Stage | Description | Wall Time |
|-------|-------------|-----------|
| **Stage 1** | Feature extraction (frozen backbone on C1: 31.3k + C2: 66.1k + C3: 33.1k + Rice: 6.7k) | ~10–14 hours |
| **Stage 2** | MLP training (10 folds × 3 pooling modes = 30 runs) | ~1.5–2.5 hours |
| **Stage 3** | Model selection | <1 minute |
| **Stage 4** | Batched test-set evaluation (all 6 pooling combos via `DataLoader` `batch_size=512`) | <2 minutes |
| **Total** | End-to-end | ~12–17 hours |

---

## 6. Results

Training results and test-set evaluation are documented separately:

- **DeepAraPPI pretrained results:** [`benchmark_analysis_pretrained_pplm_deeparappi.md`](benchmark_analysis_pretrained_pplm_deeparappi.md)
- **ESMAraPPI pretrained results:** [`benchmark_analysis_pretrained_pplm_esmarappi.md`](benchmark_analysis_pretrained_pplm_esmarappi.md)
- **Zero-shot baseline results (DeepAraPPI):** [`../zero_shot_benchmarking/benchmark_analysis_deeparappi_vs_pplm.md`](../zero_shot_benchmarking/benchmark_analysis_deeparappi_vs_pplm.md)
- **Zero-shot baseline results (ESMAraPPI):** [`../zero_shot_benchmarking/benchmark_analysis_esmarappi_vs_pplm.md`](../zero_shot_benchmarking/benchmark_analysis_esmarappi_vs_pplm.md)

---

## 7. Pipeline Scripts

### 7.1 Script Reference

| Script | Stage | Description |
|--------|-------|-------------|
| [`extract_features.py`](../../scripts/training/extract_features.py) | 1 | Frozen PPLM backbone → pooled feature `.pkl` files (train & test splits) |
| [`dataset.py`](../../scripts/training/dataset.py) | 2, 4 | `Dataset_PPI_Interact` (all 15 feature keys, `order_mode`) + stratified CV split file generator |
| [`ppi_model.py`](../../scripts/training/ppi_model.py) | 2, 4 | MLP classifier architecture + author 9-metric `evaluate()` & `evaluate_extended()` |
| [`train_ppi_head.py`](../../scripts/training/train_ppi_head.py) | 2 | Per-fold CV training (`--cross_round 0..9`, `--seed_round 1`, `cudnn.enabled=False`) |
| [`select_top_models.py`](../../scripts/training/select_top_models.py) | 3 | Rank by column-14 `AUC_PR`, generate `.top5_list` path files & `ppi_plant_models.pkl` |
| [`test_ppi_head.py`](../../scripts/training/test_ppi_head.py) | 4 | Fast batched multi-model multi-pooling evaluation across all 6 pooling configs |
| [`evaluate_pplm.py`](../../scripts/benchmarking/evaluate_pplm.py) | 5 | Metrics computation, baseline comparison, and multi-pooling ablation summary |

### 7.2 PBS Job Scripts

| Script | Purpose |
|--------|---------|
| [`run_train_deeparappi_nscc.pbs`](../../scripts/training/hpc/run_train_deeparappi_nscc.pbs) | Full aligned pipeline: Stages 1–5 on DeepAraPPI (3 modes × 10 folds, 6 test combos) |
| [`run_train_esmarappi_nscc.pbs`](../../scripts/training/hpc/run_train_esmarappi_nscc.pbs) | Full aligned pipeline: Stages 1–5 on ESMAraPPI (3 modes × 10 folds, 6 test combos) |
| [`run_resume_deeparappi_nscc.pbs`](../../scripts/training/hpc/run_resume_deeparappi_nscc.pbs) | Resume from Stage 3 (model selection → batched multi-pooling testing → evaluation) |

---

## 8. Future Directions

1. **Focal loss ablation** ($\gamma=2$, $\alpha=0.25$): Test whether focal loss improves recall for the minority positive class without sacrificing the exceptional specificity
2. **LoRA backbone fine-tuning** ($r=8$, $\alpha=16$): Apply Low-Rank Adaptation to the PPLM cross-attention layers to adapt the backbone representations to plant protein biophysics
3. **Joint training**: Train on merged DeepAraPPI C1 + ESMAraPPI C1 to produce a unified plant model
4. **Extended cross-species**: Evaluate rice, maize, and soybean transfer with the plant-trained model
5. **Auxiliary feature fusion**: Integrate GO term and domain interaction features as proposed in `plant-pplm-methodology.md` §2.3

---

## References

1. Liu, J., Chen, H. & Zhang, Y. A paired sequence language model for protein-protein interaction modeling. *Nature Communications* 17:3733 (2026).
2. Zheng, J., Yang, X., Huang, Y., Yang, S., Wuchty, S. & Zhang, Z. Deep learning-assisted prediction of protein–protein interactions in Arabidopsis thaliana. *The Plant Journal* 114:984–994 (2023).
3. Zhou, K. et al. Pre-trained protein language model sheds new light on the prediction of Arabidopsis protein-protein interactions. *Plant Methods* 19:141 (2023).
4. Park, Y. & Marcotte, E.M. Flaws in evaluation schemes for pair-input computational predictions. *Nature Methods* 9:1134–1136 (2012).
