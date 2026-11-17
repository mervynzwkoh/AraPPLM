# Plant-PPLM PPI Head Retraining: Technical Report & Run Instructions

## 1. Summary of Changes

Built and aligned the complete pipeline to retrain the PPLM-PPI classifier MLP head on plant *Arabidopsis thaliana* C1 training datasets, following the original PPLM author's training and evaluation codebase ([`ppi_training_code/`](../../ppi_training_code/)) 1:1 while substituting plant data for human data. The pipeline runs end-to-end on NSCC ASPIRE 2A (A100 GPU).

### What was built & aligned

| # | File | Purpose |
|---|------|---------|
| 1 | [extract_features.py](../../scripts/training/extract_features.py) | Stage 1: Extract PPLM backbone features for C1 train AND C2/C3/Rice test pairs → `.pkl` files |
| 2 | [ppi_model.py](../../scripts/training/ppi_model.py) | PPI classifier architecture (identical to original PPLM) + author 9-metric `evaluate()` & `evaluate_extended()` |
| 3 | [dataset.py](../../scripts/training/dataset.py) | `Dataset_PPI_Interact` (all 15 feature keys, `"sample"`, `order_mode`) + stratified `.trainset`/`.testset` split generator |
| 4 | [train_ppi_head.py](../../scripts/training/train_ppi_head.py) | Stage 2: Per-fold CV training (`--cross_round 0..9`, `--seed_round 1`, `cudnn.enabled=False`) across `mean`, `max`, `min` |
| 5 | [select_top_models.py](../../scripts/training/select_top_models.py) | Stage 3: Rank checkpoints by column-14 `AUC_PR`, write `.top5_list` files and `ppi_plant_models.pkl` |
| 6 | [test_ppi_head.py](../../scripts/training/test_ppi_head.py) | Stage 4: Fast batched evaluation (`batch_size=512`) over pre-extracted test `.pkl` features across all 6 pooling combos |
| 7 | [run_train_deeparappi_nscc.pbs](../../scripts/training/hpc/run_train_deeparappi_nscc.pbs) | PBS: End-to-end DeepAraPPI training & 6-configuration evaluation |
| 8 | [run_train_esmarappi_nscc.pbs](../../scripts/training/hpc/run_train_esmarappi_nscc.pbs) | PBS: End-to-end ESMAraPPI training & 6-configuration evaluation |

---

## 2. How the Original PPLM-PPI Training Works

The original PPLM-PPI training pipeline ([`ppi_training_code/`](../../ppi_training_code/)) operates in **two completely decoupled stages**:

### Stage 1: Offline Feature Extraction (Frozen PPLM Backbone)

The 650M-parameter PPLM backbone (33 Transformer layers, ESM2 initialization, paired cross-attention with RoPE for intra-chain and non-positional for inter-chain) is run **once** in inference mode on every protein pair (both training and held-out test sets). For each pair, it produces:

1. **Full attention tensor**: `[33 layers, 20 heads, L, L]` where `L = len(seqA) + len(seqB) + 4` (special tokens)
2. **Per-residue embeddings**: `[L, 1280]` from layer 33

These raw tensors are immediately **pooled** into compact fixed-size vectors:

| Feature | Dimension | Description |
|---------|-----------|-------------|
| `inter_attn` | `[660]` | Inter-protein cross-attention (symmetrized `(AB + BA^T)/2`, pooled over spatial dims) |
| `attn_AA` | `[660]` | Intra-protein self-attention for protein A |
| `attn_BB` | `[660]` | Intra-protein self-attention for protein B |
| `embed_A` | `[1280]` | Pooled per-residue embedding for A |
| `embed_B` | `[1280]` | Pooled per-residue embedding for B |

Three pooling operators (**mean**, **max**, and **min**) are saved in each `.pkl` file (15 feature tensors + `label` + `lens`).

### Stage 2: Per-Fold MLP Head Training

The classifier is a **5-layer MLP** (`PPI_inter_intra_attn_embed_single_pooling`) trained on the pre-extracted features:

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

Training protocol (`ppi_training_code/training.PPI.single_gpu.cross_validation.py` & `run_training.sh`):
- **Cross-validation**: 10 pre-generated fold splits (`--cross_round 0..9`)
- **Loss**: Binary Cross-Entropy (BCE)
- **Optimizer**: AdamW (`lr=5e-5`, `weight_decay=0.01`)
- **Epochs**: 15 per fold
- **Batch size**: 512 (`num_workers=8`, `prefetch_factor=4`, `pin_memory=True`)
- **Data augmentation**: 50% random protein-order swapping (A↔B) to enforce PPI symmetry
- **Seeds & Determinism**: `--seed_round 1` (`seeds[1] = 73591024`) with `torch.backends.cudnn.enabled = False`

### Model Selection & Multi-Pooling Testing (`run_testing.sh`)

After training **10 folds × 3 pooling modes (mean, max, min) = 30 runs**:
1. All 150 checkpoints per pooling mode (10 folds × 15 epochs) are ranked globally by validation `AUC_PR` (column 14 of `.recording` files).
2. Top-5 checkpoints per pooling mode are listed into `.top5_list` files (and packaged into `ppi_plant_models.pkl`).
3. Held-out test evaluation runs via batched `DataLoader(batch_size=512)` over pre-extracted test `.pkl` features with `order_mode="both"` (symmetric pair averaging $(f(A, B) + f(B, A)) / 2$) across all **6 configurations**:
   - `mean` (5 models), `max` (5 models), `min` (5 models)
   - `combine_mean_max` (10 models), `combine_mean_max_min` (15 models), `combine_mean_min` (10 models)

---

## 3. Aligned Plant Retraining Pipeline Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│ Stage 1: Feature Extraction (extract_features.py)               │
│   Frozen PPLM backbone → .pkl per pair (mean, max, min)         │
│   Run on C1 Train AND C2, C3, Rice Test Sets                    │
│   Output: features/{Dataset}_{C1,C2,C3,Rice}/{ProtA}@{ProtB}.pkl│
└──────────────────────────┬──────────────────────────────────────┘
                           │
┌──────────────────────────▼──────────────────────────────────────┐
│ Stage 2: Per-Fold MLP Head Training (train_ppi_head.py)         │
│   Stratified 10-fold split files × 3 modes (mean, max, min)     │
│   --cross_round 0..9, --seed_round 1 (73591024), cudnn=False    │
│   Output: models/{Dataset}/PPLM-PPI.inter_intra_attn_embed.*    │
└──────────────────────────┬──────────────────────────────────────┘
                           │
┌──────────────────────────▼──────────────────────────────────────┐
│ Stage 3: Model Selection (select_top_models.py)                 │
│   Parse `Accuarcy Epoch:` (col 14 AUC_PR) → top-5 per mode      │
│   Output: .top5_list files + ppi_plant_models.pkl               │
└──────────────────────────┬──────────────────────────────────────┘
                           │
┌──────────────────────────▼──────────────────────────────────────┐
│ Stage 4: Fast Batched Ensemble Testing (test_ppi_head.py)       │
│   DataLoader(batch_size=512) over test .pkl features            │
│   Evaluates all 6 pooling configs with order_mode="both"        │
│   Output: results/{Dataset}/pretrained/*_scores.csv (<2 mins)   │
└──────────────────────────┬──────────────────────────────────────┘
                           │
┌──────────────────────────▼──────────────────────────────────────┐
│ Stage 5: Evaluation & Ablation (evaluate_pplm.py)               │
│   AUPRC, AUC_PR, AUROC, Precision, Recall, F1, MCC, Top-K       │
│   Output: benchmark_summary.csv + pooling_ablation_summary.csv  │
└─────────────────────────────────────────────────────────────────┘
```

---

## 4. Per-File Technical Details

### 4.1 `extract_features.py`
Runs the frozen PPLM backbone on each protein pair and saves 15 pooled feature vectors + label as `{ProtA}@{ProtB}.pkl`.
- Combined cropping `LA + LB <= 1020` matching PPLM pretraining.
- Two-stage OOM fallback (1020 → 800 token budget), `--resume` skip logic, and periodic VRAM cleanup every 200 pairs.

### 4.2 `dataset.py`
- **`Dataset_PPI_Interact`**: Matches `ppi_training_code/dataset.py`. Reads split files (`.trainset`/`.testset`) or raw TSVs, loads all 15 feature keys (`mean_*`, `max_*`, `min_*`), returns `"sample"` ID + `"label"`, and supports `order_mode` (`"random"`, `"both"`, `"normal"`, `"reverse"`).
- **`generate_cv_split_files()`**: Persists stratified 10-fold `.trainset` and `.testset` files on disk (`StratifiedKFold`, `seed=73591024`) preserving the 1:10 ratio in every fold.

### 4.3 `train_ppi_head.py`
- Implements the author's `SEEDS` pool defaulting to `--seed_round 1` (`73591024`) and `torch.backends.cudnn.enabled = False`.
- Supports per-fold execution (`--cross_round 0..9`) matching `run_training.sh` (or `--all_folds` for single-command execution) with `model.pkl` checkpoint resumption.
- Writes exact author-format `.recording` logs (`========== Accuarcy Epoch:` with `AUC_PR` at column 14).

### 4.4 `select_top_models.py`
- Parses `.recording` logs across `mean`, `max`, and `min`, ranks globally by column-14 `AUC_PR`, writes `PPLM-PPI.inter_intra_attn_embed.{mode}.round_1.top5_list` files (one checkpoint path per line), and packages `ppi_plant_models.pkl`.

### 4.5 `test_ppi_head.py`
- Loads `--mean_path`, `--max_path`, `--min_path` (or `--model_weights`) and runs batched `DataLoader(batch_size=512)` over pre-extracted test `.pkl` features (`--feat_dir`) with `--order_mode both` symmetric pair averaging.
- Supports fallback on-the-fly backbone extraction (`--seq_db`) when `--feat_dir` is omitted.
- Outputs both `pair_id,pred_score,true_label` CSVs and optional author `.preds.csv` / `.preds.statistic` files.

### 4.6 `ppi_model.py`
- Identical `PPI_inter_intra_attn_embed_single_pooling` MLP head.
- `evaluate()` returns the author's exact 9-tuple `(precision, recall, accuracy, F1_score, Top10, Top20, Top50, AUC_ROC, AUC_PR)`; `evaluate_extended()` additionally returns `specificity` and `MCC`.

---

## 5. How to Run

### 5.1 Prerequisites
1. **PPLM weights** downloaded to `PPLM/weights/` (`pplm_t33_650M.pt`, `ppi_models.pkl`).
2. **Sequence databases**:
   - `data/arabidopsis/uniprot_final.pkl`
   - `data/rice/uniprot_rice_final.pkl`
3. **Conda environment**: `conda activate PPLM_NSCC_A100`

### 5.2 Option A: One-Command PBS Submission (Recommended)

```bash
# Train & evaluate on DeepAraPPI (all 3 modes, all 6 combinations)
qsub scripts/training/hpc/run_train_deeparappi_nscc.pbs

# Train & evaluate on ESMAraPPI (all 3 modes, all 6 combinations)
qsub scripts/training/hpc/run_train_esmarappi_nscc.pbs
```

### 5.3 Option B: Manual Stage-by-Stage Execution

#### Stage 1: Extract Features (Train & Test Sets)
```bash
for split in c1 c2 c3; do
    python -u scripts/training/extract_features.py \
        --input data/legacy_benchmarks/DeepAraPPI/${split}_ppi_sample_DeepAraPPI.txt \
        --output_dir features/DeepAraPPI_${split^^}/ \
        --seq_db data/arabidopsis/uniprot_final.pkl \
        --gpu_id 0 --resume
done
```

#### Stage 2: Train MLP Head (10 folds × 3 modes, `--seed_round 1`)
```bash
for mode in mean max min; do
    for round in 0 1 2 3 4 5 6 7 8 9; do
        python -u scripts/training/train_ppi_head.py \
            models/DeepAraPPI/PPLM-PPI.inter_intra_attn_embed.${mode} \
            --feat_dir features/DeepAraPPI_C1/ \
            --split_dir models/DeepAraPPI/splits \
            --pair_list data/legacy_benchmarks/DeepAraPPI/c1_ppi_sample_DeepAraPPI.txt \
            --dataset_prefix deeparappi_c1 \
            --mode PPI_inter_intra_attn_embed_${mode} \
            --cross_round ${round} --seed_round 1 \
            --epochs 15 --batch_size 512 --learning_rate 5e-5 --gpu_id 0
    done
done
```

#### Stage 3: Select Top-5 Models
```bash
python -u scripts/training/select_top_models.py \
    --model_dir models/DeepAraPPI/ \
    --output models/DeepAraPPI/ppi_plant_models.pkl \
    --top_k 5 --modes mean max min --seed_round 1
```

#### Stage 4: Fast Batched Test Evaluation
```bash
python -u scripts/training/test_ppi_head.py \
    --input data/legacy_benchmarks/DeepAraPPI/c2_ppi_sample_DeepAraPPI.txt \
    --feat_dir features/DeepAraPPI_C2/ \
    --mean_path models/DeepAraPPI/PPLM-PPI.inter_intra_attn_embed.mean.round_1.top5_list \
    --max_path  models/DeepAraPPI/PPLM-PPI.inter_intra_attn_embed.max.round_1.top5_list \
    --min_path  models/DeepAraPPI/PPLM-PPI.inter_intra_attn_embed.min.round_1.top5_list \
    --output results/DeepAraPPI/pretrained/deepara_c2_combine_mean_max_min_scores.csv \
    --order_mode both --batch_size 512 --gpu_id 0
```

#### Stage 5: Evaluation & Multi-Pooling Ablation Summary
```bash
python -u scripts/benchmarking/evaluate_pplm.py \
    --c2 results/DeepAraPPI/pretrained/deepara_c2_combine_mean_max_min_scores.csv \
    --c3 results/DeepAraPPI/pretrained/deepara_c3_combine_mean_max_min_scores.csv \
    --rice results/DeepAraPPI/pretrained/deepara_rice_combine_mean_max_min_scores.csv \
    --pooling_ablation_dir results/DeepAraPPI/pretrained/ \
    --suite DeepAraPPI \
    --output_dir results/DeepAraPPI/pretrained/
```

---

## 6. Expected Outputs

```
models/DeepAraPPI/
├── splits/
│   ├── deeparappi_c1.cross_validation_{0..9}.trainset
│   └── deeparappi_c1.cross_validation_{0..9}.testset
├── PPLM-PPI.inter_intra_attn_embed.{mean,max,min}.round_1.cv_{0..9}/
├── PPLM-PPI.inter_intra_attn_embed.{mean,max,min}.round_1.cv_{0..9}.recording
├── PPLM-PPI.inter_intra_attn_embed.{mean,max,min}.round_1.top5_list
└── ppi_plant_models.pkl

results/DeepAraPPI/pretrained/
├── deepara_{c2,c3,rice}_{mean,max,min,combine_mean_max,combine_mean_max_min,combine_mean_min}_scores.csv
├── benchmark_summary.csv
└── pooling_ablation_summary.csv
```

---

## 7. Prior 2-Mode Run Reference & Pending Aligned Rerun

> [!NOTE]
> The benchmark numbers summarized below are from the **initial 2-mode retraining run** (`mean + max` pooling only, `seed = 26240761`, 10-model ensemble). After aligning `scripts/training/` and `scripts/benchmarking/` with the original author's full 3-mode (`mean + max + min`, `seed_round 1 = 73591024`) protocol and 6-combination testing setup, the NSCC ASPIRE 2A rerun is scheduled to populate updated tables in [`benchmark_analysis_pretrained_pplm_deeparappi.md`](benchmark_analysis_pretrained_pplm_deeparappi.md) and [`benchmark_analysis_pretrained_pplm_esmarappi.md`](benchmark_analysis_pretrained_pplm_esmarappi.md).

### DeepAraPPI Reference (Initial 2-Mode Run)

| Task | Zero-Shot PPLM (AUPRC) | Prior 2-Mode Plant PPLM (AUPRC) | DeepAraPPI Baseline |
|------|------------------------|----------------------------------|---------------------|
| C2 (One Unseen) | 0.5738 | **0.8738** | 0.8970 |
| C3 (Both Unseen) | 0.5525 | **0.8118** | 0.8250 |
| Rice (Cross-Species) | 0.4297 | **0.3555** | 0.3050 |

### ESMAraPPI Reference (Initial 2-Mode Run)

| Evaluation Set | Zero-Shot PPLM (AUPRC) | Prior 2-Mode Plant PPLM (AUPRC) | ESMAraPPI Baseline |
| :--- | :--- | :--- | :--- |
| **Test 1 (`c2Pred.txt` - Homology/Random)** | 0.5092 | **0.8408** | 0.8340 |
| **Test 2 (`c3Pred.txt` - Strict $\le 40\%$ Non-redundant)** | 0.5610 | **0.8103** | 0.8100 |
