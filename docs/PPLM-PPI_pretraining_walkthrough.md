# Plant-PPLM PPI Head Retraining: Technical Report & Run Instructions

## 1. Summary of Changes

Built a complete pipeline to retrain the PPLM-PPI classifier MLP head on plant *Arabidopsis thaliana* C1 training datasets, replicating the original PPLM-PPI training protocol but substituting plant data for human data. The pipeline is designed to run end-to-end on NSCC ASPIRE 2A (A100 GPU).

### What was built

| # | File | Purpose |
|---|------|---------|
| 1 | [extract_features.py](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/AraPPLM/scripts/training/extract_features.py) | Stage 1: Extract PPLM backbone features for C1 training pairs → `.pkl` files |
| 2 | [ppi_model.py](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/AraPPLM/scripts/training/ppi_model.py) | PPI classifier architecture (identical to original PPLM) + extended metrics |
| 3 | [dataset.py](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/AraPPLM/scripts/training/dataset.py) | PyTorch Dataset for loading pre-extracted `.pkl` features |
| 4 | [train_ppi_head.py](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/AraPPLM/scripts/training/train_ppi_head.py) | Stage 2: K-fold stratified CV training of the MLP head |
| 5 | [select_top_models.py](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/AraPPLM/scripts/training/select_top_models.py) | Stage 3: Rank checkpoints by validation AUPRC, package top-5 into ensemble |
| 6 | [test_ppi_head.py](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/AraPPLM/scripts/training/test_ppi_head.py) | Stage 4: Full inference (PPLM backbone + plant-trained head) on test sets |
| 7 | [run_train_deeparappi_nscc.pbs](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/AraPPLM/scripts/run_train_deeparappi_nscc.pbs) | PBS: End-to-end DeepAraPPI C1 training → C2/C3/Rice evaluation |
| 8 | [run_train_esmarappi_nscc.pbs](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/AraPPLM/scripts/run_train_esmarappi_nscc.pbs) | PBS: End-to-end ESMAraPPI C1 training → C2/C3 evaluation |

### Repository reorganisation

```diff
 AraPPLM/
 ├── PPLM/                               # [UNCHANGED] Original PPLM submodule
 ├── data/                               # [UNCHANGED] Benchmark datasets
 ├── docs/                               # [UNCHANGED] Documentation
 ├── results/                            # [UNCHANGED] Benchmark results
+├── features/                           # [NEW] Pre-extracted PPLM features
+│   ├── DeepAraPPI_C1/                  #   Feature .pkl files for DeepAraPPI C1
+│   └── ESMAraPPI_C1/                   #   Feature .pkl files for ESMAraPPI C1
+├── models/                             # [NEW] Trained model checkpoints
+│   ├── DeepAraPPI/                     #   Checkpoints from DeepAraPPI C1 training
+│   └── ESMAraPPI/                      #   Checkpoints from ESMAraPPI C1 training
 ├── scripts/
 │   ├── batch_predict.py                # [KEPT] Original zero-shot benchmark (backward compat)
 │   ├── evaluate_pplm.py                # [KEPT] Original evaluation engine (backward compat)
+│   ├── benchmarking/                   # [NEW] Copies of benchmark scripts
+│   │   ├── batch_predict.py
+│   │   └── evaluate_pplm.py
+│   ├── training/                       # [NEW] Complete PPI head retraining pipeline
+│   │   ├── extract_features.py
+│   │   ├── train_ppi_head.py
+│   │   ├── select_top_models.py
+│   │   ├── test_ppi_head.py
+│   │   ├── ppi_model.py
+│   │   └── dataset.py
 │   ├── data_preparation/               # [UNCHANGED]
 │   ├── run_all_benchmarks_nscc.pbs     # [UNCHANGED] Existing zero-shot PBS scripts
+│   ├── run_train_deeparappi_nscc.pbs   # [NEW] PBS: Train on DeepAraPPI C1
+│   └── run_train_esmarappi_nscc.pbs    # [NEW] PBS: Train on ESMAraPPI C1
 └── .gitignore                          # [UPDATED] Added features/, models/, *.log
```

---

## 2. How the Original PPLM-PPI Training Works

The original PPLM-PPI training pipeline operates in **two completely decoupled stages**, matching the codebase in [`ppi_training_code/`](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/ppi_training_code):

### Stage 1: Offline Feature Extraction (Frozen PPLM Backbone)

The 650M-parameter PPLM backbone (33 Transformer layers, ESM2 initialization, paired cross-attention with RoPE for intra-chain and non-positional for inter-chain) is run **once** in inference mode on every protein pair. For each pair, it produces:

1. **Full attention tensor**: `[33 layers, 20 heads, L, L]` where `L = len(seqA) + len(seqB) + 4` (special tokens)
2. **Per-residue embeddings**: `[L, 1280]` from the final transformer layer

These raw tensors are immediately **pooled** into compact fixed-size vectors:

| Feature | Dimension | Description |
|---------|-----------|-------------|
| `mean_inter_attn` | `[660]` | Inter-protein cross-attention (mean over spatial dims) |
| `mean_attn_AA` | `[660]` | Intra-protein self-attention for protein A |
| `mean_attn_BB` | `[660]` | Intra-protein self-attention for protein B |
| `mean_embed_A` | `[1280]` | Mean-pooled per-residue embedding for A |
| `mean_embed_B` | `[1280]` | Mean-pooled per-residue embedding for B |

The same is repeated with max-pooling and min-pooling, yielding 15 feature vectors per pair. The dimension 660 = 33 layers × 20 heads. All features are saved as `.pkl` files.

### Stage 2: MLP Head Training

The classifier is a **5-layer MLP** (`PPI_inter_intra_attn_embed_single_pooling`) trained on the pre-extracted features. Key architectural details:

```
Input (one pooling mode at a time):
  inter_attn [660] ─→ linear_inter(660→660) ─┐
  attn_AA    [660] ─→ linear_intra(660→660) ─┤  (shared weights for A/B)
  attn_BB    [660] ─→ linear_intra(660→660) ─┤
  embed_A    [1280]─→ linear_embed(1280→660)─┤
  embed_B    [1280]─→ linear_embed(1280→660)─┘
                                              │
  Concatenate ─────────────────────────── [3300]
                                              │
  FC 3300→1024 + LayerNorm + ReLU            │
  FC 1024→512  + LayerNorm + ReLU            │
  FC 512→256   + LayerNorm + ReLU            │
  FC 256→128   + LayerNorm + ReLU            │
  FC 128→1     + Sigmoid ─────────────── [1]  (probability)
```

Training protocol:
- **Cross-validation**: 10-fold on the human D-SCRIPT training set
- **Loss**: Binary Cross-Entropy (BCE)
- **Optimizer**: AdamW (lr=5e-5, weight_decay=0.01)
- **Epochs**: 15 per fold
- **Batch size**: 512
- **Data augmentation**: 50% random protein-order swapping (A↔B) to enforce PPI symmetry
- **Seeds**: Deterministic from a pool of 10 pre-defined seeds

### Model Selection & Ensemble Inference

After training 10 folds × 3 pooling modes (mean, max, min):
1. All 150 checkpoints (10 folds × 15 epochs × 1 mode) are ranked by validation AUPRC
2. Top-5 per pooling mode are selected
3. At inference, all 10 classifiers (5 mean + 5 max) run on each pair
4. Each classifier also runs with A↔B swapped, then averages (symmetric pair averaging)
5. Final score = mean of all 10 × 2 = 20 forward passes

---

## 3. Plant Retraining Pipeline Architecture

Our pipeline follows the exact same 2-stage architecture but substitutes plant C1 training data for human D-SCRIPT data:

```
┌─────────────────────────────────────────────────────────────────┐
│ Stage 1: Feature Extraction (extract_features.py)               │
│   Frozen PPLM backbone → pooled features → .pkl per pair        │
│   Input: data/DeepAraPPI/c1_ppi_sample_DeepAraPPI.txt           │
│   Output: features/DeepAraPPI_C1/{ProtA}@{ProtB}.pkl            │
│   Resources: 1× A100 GPU, ~3-6 hours for 31,284 pairs          │
└──────────────────────────┬──────────────────────────────────────┘
                           │
┌──────────────────────────▼──────────────────────────────────────┐
│ Stage 2: MLP Head Training (train_ppi_head.py)                  │
│   10-fold stratified CV × 2 pooling modes (mean, max)           │
│   BCE loss, AdamW lr=5e-5, 15 epochs, batch 512                 │
│   50% random A↔B swapping                                       │
│   Output: models/DeepAraPPI/plant_ppi.{mode}.cv_{0-9}/          │
│   Resources: any GPU, ~1-2 hours total                          │
└──────────────────────────┬──────────────────────────────────────┘
                           │
┌──────────────────────────▼──────────────────────────────────────┐
│ Stage 3: Model Selection (select_top_models.py)                 │
│   Parse recording files → rank by validation AUPRC              │
│   Package top-5 mean + top-5 max → ppi_plant_models.pkl         │
│   Output: models/DeepAraPPI/ppi_plant_models.pkl                │
│   Resources: CPU only, seconds                                  │
└──────────────────────────┬──────────────────────────────────────┘
                           │
┌──────────────────────────▼──────────────────────────────────────┐
│ Stage 4: Ensemble Testing (test_ppi_head.py)                    │
│   Full PPLM backbone + plant-trained head on test sets           │
│   10-model ensemble with symmetric pair averaging               │
│   Output: results/DeepAraPPI/deepara_c{2,3}_plant_scores.csv    │
│   Resources: 1× A100 GPU, ~2-4 hours per test set              │
└──────────────────────────┬──────────────────────────────────────┘
                           │
┌──────────────────────────▼──────────────────────────────────────┐
│ Stage 5: Evaluation (evaluate_pplm.py)                          │
│   AUPRC, AUROC, Precision, Recall, F1, MCC, Specificity         │
│   Comparison against published baselines                        │
│   Output: results/DeepAraPPI/benchmark_summary.csv              │
└─────────────────────────────────────────────────────────────────┘
```

---

## 4. Per-File Technical Details

### 4.1 `extract_features.py`

Runs the frozen PPLM backbone on each protein pair and saves 15 pooled feature vectors + label as a `.pkl` file.

**Key design decisions:**
- Combined cropping `LA + LB ≤ 1020` matching PPLM paper pretraining
- OOM resilience: two-stage fallback (1020 → 800 token budget)
- `--resume` flag for crash recovery: skips pairs with existing `.pkl` files
- Periodic VRAM cleanup every 200 pairs

**Feature file format** (each `.pkl`):
```python
{
    "mean_inter_attn": tensor([660]),   # Mean-pooled inter-protein cross-attention
    "mean_attn_AA":    tensor([660]),   # Mean-pooled intra-protein attention A
    "mean_attn_BB":    tensor([660]),   # Mean-pooled intra-protein attention B
    "mean_embed_A":    tensor([1280]),  # Mean-pooled embedding A
    "mean_embed_B":    tensor([1280]),  # Mean-pooled embedding B
    "max_inter_attn":  tensor([660]),   # Max-pooled (same structure)
    # ... (same for max and min)
    "label": int,                       # Ground truth (0 or 1)
    "lens": (lenA, lenB),              # Sequence lengths after cropping
}
```

### 4.2 `dataset.py`

**`PlantPPIDataset`**: PyTorch Dataset that loads pre-extracted `.pkl` files. Key features:
- Takes a list of `(pair_id, label)` tuples + feature directory path
- Selects correct pooling mode (mean/max/min) via constructor arg
- 50% random A↔B swapping during training (matching PPLM protocol)
- Resilient to corrupted files (re-samples on error, up to 10 attempts)

**`load_pair_list()`**: Reads a TSV dataset file and returns `(pair_id, label)` tuples. Handles both DeepAraPPI format (3 cols with header) and ESMAraPPI format (4 cols, no header).

### 4.3 `train_ppi_head.py`

Main training script with stratified K-fold CV. Key details:
- Uses `StratifiedKFold` from scikit-learn to maintain 1:10 class ratio in each fold
- Fresh model initialization per fold (no transfer between folds)
- Logs per-epoch train/valid loss and validation metrics to `.recording` files
- Saves per-epoch checkpoints (`model_1.pkl` through `model_15.pkl`) and latest (`model.pkl`)
- Debug mode: `--debug` limits to 5000 samples for rapid iteration

### 4.4 `select_top_models.py`

Post-training utility that:
1. Parses all `.recording` files using regex to extract per-epoch AUPRC
2. Ranks all fold×epoch checkpoints by validation AUPRC descending
3. Loads top-5 `model_state_dict` per pooling mode
4. Saves them as `ppi_plant_models.pkl` in format: `{'mean': [state_dict, ...], 'max': [state_dict, ...]}`

This format is **directly compatible** with the existing benchmarking pipeline (`batch_predict.py`), so you could drop the plant weights into the original PPLM weights directory as a replacement.

### 4.5 `test_ppi_head.py`

End-to-end test-time inference:
1. Loads frozen PPLM backbone + plant-trained ensemble weights
2. For each test pair: runs PPLM backbone → extracts features → runs all 10 classifiers
3. Symmetric pair averaging (forward + swapped A↔B, divided by 2) per classifier
4. Final score = mean of all ensemble predictions
5. Outputs CSV compatible with `evaluate_pplm.py`

### 4.6 `ppi_model.py`

Architecture-identical to the original PPLM `PPLM_PPI` class with an extended `evaluate()` function that adds specificity and MCC to the metrics.

---

## 5. How to Run

### 5.1 Prerequisites

1. **PPLM weights** must be downloaded to `PPLM/weights/`:
   ```bash
   cd PPLM/weights && bash download_weights.sh && cd ../..
   ```
2. **Sequence databases** must be built (see [data_processing_pipeline.md](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/AraPPLM/docs/data_processing_pipeline.md)):
   - `data/arabidopsis/uniprot_final.pkl` (109,994 entries)
   - `data/rice/uniprot_rice_final.pkl` (100,297 entries)
3. **Conda environment** must be set up:
   ```bash
   conda env create -f PPLM_NSCC_A100.yml
   conda activate PPLM_NSCC_A100
   ```

### 5.2 Option A: One-Command PBS Submission (Recommended)

Submit the end-to-end PBS job to run all 5 stages automatically:

```bash
# Train on DeepAraPPI C1 → evaluate on C2/C3/Rice
cd AraPPLM
qsub scripts/run_train_deeparappi_nscc.pbs

# Monitor progress
tail -f pplm_train_deeparappi.log
```

For ESMAraPPI (after placing the C1 training file):
```bash
# Requires: data/ESMAraPPI/c1_ppi_sample_ESMAraPPI.txt
qsub scripts/run_train_esmarappi_nscc.pbs
tail -f pplm_train_esmarappi.log
```

> [!IMPORTANT]
> The DeepAraPPI PBS script requests 24 hours walltime. Feature extraction for 31,284 pairs takes ~3-6 hours on A100. MLP training is fast (~1-2 hours). Testing on C2 (66k pairs) + C3 (33k pairs) + Rice (6.7k pairs) adds ~4-6 hours. Total: ~8-14 hours.

### 5.3 Option B: Manual Stage-by-Stage Execution

If you prefer to run each stage separately (useful for debugging or if you need to split across multiple PBS jobs):

#### Stage 1: Extract Features
```bash
python -u scripts/training/extract_features.py \
    --input data/DeepAraPPI/c1_ppi_sample_DeepAraPPI.txt \
    --output_dir features/DeepAraPPI_C1/ \
    --seq_db data/arabidopsis/uniprot_final.pkl \
    --gpu_id 0 --resume
```

#### Stage 2: Train MLP Head
```bash
# Mean pooling (10-fold CV)
python -u scripts/training/train_ppi_head.py \
    --feat_dir features/DeepAraPPI_C1/ \
    --output_dir models/DeepAraPPI/ \
    --mode mean --n_folds 10 --epochs 15 --batch_size 512 --lr 5e-5

# Max pooling (10-fold CV)
python -u scripts/training/train_ppi_head.py \
    --feat_dir features/DeepAraPPI_C1/ \
    --output_dir models/DeepAraPPI/ \
    --mode max --n_folds 10 --epochs 15 --batch_size 512 --lr 5e-5
```

#### Stage 3: Select Top-5 Models
```bash
python -u scripts/training/select_top_models.py \
    --model_dir models/DeepAraPPI/ \
    --output models/DeepAraPPI/ppi_plant_models.pkl \
    --top_k 5 --modes mean max
```

#### Stage 4: Evaluate on Test Sets
```bash
# C2 (One Unseen Protein)
python -u scripts/training/test_ppi_head.py \
    --input data/DeepAraPPI/c2_ppi_sample_DeepAraPPI.txt \
    --output results/DeepAraPPI/deepara_c2_plant_scores.csv \
    --model_weights models/DeepAraPPI/ppi_plant_models.pkl \
    --seq_db data/arabidopsis/uniprot_final.pkl

# C3 (Both Unseen Proteins)
python -u scripts/training/test_ppi_head.py \
    --input data/DeepAraPPI/c3_ppi_sample_DeepAraPPI.txt \
    --output results/DeepAraPPI/deepara_c3_plant_scores.csv \
    --model_weights models/DeepAraPPI/ppi_plant_models.pkl \
    --seq_db data/arabidopsis/uniprot_final.pkl

# Rice Cross-Species
python -u scripts/training/test_ppi_head.py \
    --input data/DeepAraPPI/all_rice_positive_negative_DeepAraPPI.txt \
    --output results/DeepAraPPI/deepara_rice_plant_scores.csv \
    --model_weights models/DeepAraPPI/ppi_plant_models.pkl \
    --seq_db data/rice/uniprot_rice_final.pkl
```

#### Stage 5: Evaluation Metrics
```bash
python -u scripts/benchmarking/evaluate_pplm.py \
    --c2 results/DeepAraPPI/deepara_c2_plant_scores.csv \
    --c3 results/DeepAraPPI/deepara_c3_plant_scores.csv \
    --rice results/DeepAraPPI/deepara_rice_plant_scores.csv \
    --output_dir results/DeepAraPPI/
```

### 5.4 Debug / Quick Test
```bash
# Quick test with 2 folds, 2 epochs, limited samples
python scripts/training/train_ppi_head.py \
    --feat_dir features/DeepAraPPI_C1/ \
    --output_dir models/DeepAraPPI/ \
    --mode mean --n_folds 2 --epochs 2 --debug
```

---

## 6. Expected Outputs

After training completes:

```
models/DeepAraPPI/
├── plant_ppi.mean.cv_0/           # Fold 0 checkpoints (model_1.pkl ... model_15.pkl)
├── plant_ppi.mean.cv_0.recording  # Training log with per-epoch metrics
├── ...
├── plant_ppi.mean.cv_9/
├── plant_ppi.mean.cv_9.recording
├── plant_ppi.max.cv_0/
├── ...
├── plant_ppi.max.cv_9.recording
├── plant_ppi.mean.top5_list       # Top-5 model paths with AUPRC values
├── plant_ppi.max.top5_list
└── ppi_plant_models.pkl           # Ensemble weights (5 mean + 5 max state_dicts)

results/DeepAraPPI/
├── deepara_c2_plant_scores.csv    # Plant-trained predictions on C2
├── deepara_c3_plant_scores.csv    # Plant-trained predictions on C3
├── deepara_rice_plant_scores.csv  # Plant-trained predictions on Rice
└── benchmark_summary.csv          # Comparison vs. published baselines
```

---

## 7. Comparison to Zero-Shot Baselines

The purpose of this retraining is to compare plant-trained PPLM-PPI against the zero-shot baselines you've already established:

| Task | Zero-Shot PPLM (AUPRC) | Plant-Trained PPLM (AUPRC) | DeepAraPPI Baseline |
|------|----------------------|---------------------------|---------------------|
| C2 (One Unseen) | 0.5738 | *TBD — run training* | 0.8970 |
| C3 (Both Unseen) | 0.5525 | *TBD — run training* | 0.8250 |
| Rice (Cross-Species) | 0.4297 | *TBD — run training* | 0.3050 |

> [!NOTE]
> The key hypothesis: retraining the MLP head on Arabidopsis C1 data (while keeping the PPLM backbone frozen) should significantly improve C2 performance (where one protein was seen during training), while C3 performance (both proteins unseen) depends primarily on the quality of the backbone's learned representations. The Rice cross-species transfer is the most informative signal about whether plant-specific training helps generalization.

---

## 8. Reminder: ESMAraPPI C1 Data

> [!WARNING]
> The ESMAraPPI C1 training file (`data/ESMAraPPI/c1_ppi_sample_ESMAraPPI.txt`) is **not yet in the repository**. You need to obtain it from the [ESMAraPPI paper supplementary data](https://doi.org/10.1186/s13007-023-01119-6) or its GitHub repository, and place it at that path in TSV format (`Protein1\tProtein2\tlabel`). The ESMAraPPI PBS script will verify this file exists before proceeding.
