# Benchmark Report: PPLM-PPI (Zero-Shot & Plant-Retrained) on the Curated AraPPLM Dataset

**Author / Project:** AIS5281 - Plant Protein Language Modeling (`AraPPLM`)  
**Dataset:** Curated *Arabidopsis thaliana* Benchmark Splits ($\theta = 0.50$, $c_{\min} = 0.80$)  
**Dataset Location:** [`data/dataset_construction/04_benchmark_splits/arabidopsis/benchmark_data/`](../../data/dataset_construction/04_benchmark_splits/arabidopsis/benchmark_data/)  
**Training Methodology:** Aligned with [`ppi_training_code/`](../../ppi_training_code/) (see [`ppi_head_retraining_methodology.md`](ppi_head_retraining_methodology.md))

> [!IMPORTANT]
> **Pending NSCC Execution:** All extraction, zero-shot evaluation, and training scripts for this benchmark are prepared and verified. Submit the sharded C1 feature extraction jobs followed by the benchmark job on NSCC ASPIRE 2A (`run_extract_araplm_c1_shard_nscc.pbs` → `run_benchmark_araplm_nscc.pbs`) to populate the metric tables below.

---

## 1. Dataset Partition Summary

| Split | File Name | Positives | Negatives | Total Pairs | Pos:Neg Ratio | Partition Definition |
| :--- | :--- | :---: | :---: | :---: | :---: | :--- |
| **C1 (Train)** | `ara_benchmark_c1_train_id0_50_cov0_80.tsv` | 19,462 | 194,620 | 214,082 | 1 : 10.00 | Both interactors in $V(C_1)$ |
| **C2 (Test)** | `ara_benchmark_c2_test_id0_50_cov0_80.tsv` | 1,145 | 11,450 | 12,595 | 1 : 10.00 | Cut boundary (1 seen, 1 unseen) |
| **C3 (Test)** | `ara_benchmark_c3_test_id0_50_cov0_80.tsv` | 1,145 | 10,965 | 12,110 | 1 : 9.58 | Both unseen ($u, v \in V_{\text{test}}$, strict zero-shot) |
| **C4 (Test)** | `ara_benchmark_c4_test_id0_50_cov0_80.tsv` | 1,145 | 11,450 | 12,595 | 1 : 10.00 | Both seen ($u, v \in V(C_1)$, novel test edge) |
| **Total** | — | **22,897** | **228,485** | **251,382** | **1 : 9.98** | **100.00% sequence coverage** |

---

## 2. Benchmark Configurations

### 2.1 PPLM-PPI Zero-Shot
- **Backbone:** Frozen PPLM (`pplm_t33_650M.pt`, 650M params)
- **Head:** Pre-trained human D-SCRIPT ensemble (`PPLM/weights/ppi_models.pkl`, 5 `mean` + 5 `max` classifiers)
- **Execution:** Fast batched `DataLoader` (`batch_size=512`, `order_mode=both`) over pre-extracted `features/AraPPLM_{C2,C3,C4}/`
- **Outputs:** `results/AraPPLM/zero-shot/benchmark_summary.csv`

### 2.2 PPLM-PPI Plant-Retrained (Author Protocol)
- **Backbone:** Frozen PPLM (`pplm_t33_650M.pt`)
- **Head:** 5-layer MLP (`PPI_inter_intra_attn_embed_single_pooling`, ~3.9M params) trained from random initialization on `AraPPLM C1` (214,082 pairs)
- **Protocol:** 10-fold stratified CV (`--seed_round 1`, `seed = 73591024`, `cudnn.enabled = False`), 3 pooling modes (`mean`, `max`, `min` = 30 runs), AdamW (`lr=5e-5`, `weight_decay=0.01`, `epochs=15`, `batch_size=512`), plain BCE loss
- **Primary Summary (`benchmark_summary.csv`):** Top-5 `mean` + Top-5 `max` (`combine_mean_max`, 10 models) with symmetric pair averaging (`order_mode=both`)
- **Outputs:** `results/AraPPLM/pretrained/benchmark_summary.csv` & `pooling_ablation_summary.csv`

---

## 3. Zero-Shot vs. Plant-Retrained Comparison (`combine_mean_max`)

| Metric | C2 Zero-Shot | C2 Retrained | Δ | C3 Zero-Shot | C3 Retrained | Δ | C4 Zero-Shot | C4 Retrained | Δ |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **AUPRC (AP)** | *Pending* | *Pending* | — | *Pending* | *Pending* | — | *Pending* | *Pending* | — |
| **AUC_PR (Trapz)** | *Pending* | *Pending* | — | *Pending* | *Pending* | — | *Pending* | *Pending* | — |
| **AUROC** | *Pending* | *Pending* | — | *Pending* | *Pending* | — | *Pending* | *Pending* | — |
| **Precision ($\tau=0.5$)** | *Pending* | *Pending* | — | *Pending* | *Pending* | — | *Pending* | *Pending* | — |
| **Recall ($\tau=0.5$)** | *Pending* | *Pending* | — | *Pending* | *Pending* | — | *Pending* | *Pending* | — |
| **Specificity ($\tau=0.5$)** | *Pending* | *Pending* | — | *Pending* | *Pending* | — | *Pending* | *Pending* | — |
| **F1 Score ($\tau=0.5$)** | *Pending* | *Pending* | — | *Pending* | *Pending* | — | *Pending* | *Pending* | — |
| **MCC** | *Pending* | *Pending* | — | *Pending* | *Pending* | — | *Pending* | *Pending* | — |

---

## 4. Multi-Pooling Ablation Summary (Plant-Retrained)

Populated from `results/AraPPLM/pretrained/pooling_ablation_summary.csv`:

| Pooling Configuration | Models in Ensemble | C2 AUPRC | C3 AUPRC | C4 AUPRC |
| :--- | :---: | :---: | :---: | :---: |
| `mean` (top-5) | 5 (`mean`) | *Pending* | *Pending* | *Pending* |
| `max` (top-5) | 5 (`max`) | *Pending* | *Pending* | *Pending* |
| `min` (top-5) | 5 (`min`) | *Pending* | *Pending* | *Pending* |
| `combine_mean_max` | 10 (`mean` + `max`) | *Pending* | *Pending* | *Pending* |
| `combine_mean_max_min` | 15 (`mean` + `max` + `min`) | *Pending* | *Pending* | *Pending* |
| `combine_mean_min` | 10 (`mean` + `min`) | *Pending* | *Pending* | *Pending* |

---

## 5. Execution Commands (NSCC ASPIRE 2A)

```bash
# 1. Run C1 feature extraction across 2 parallel shards (214,082 pairs total)
qsub -v NUM_SHARDS=2,SHARD_ID=0 -o pplm_extract_araplm_c1_s0.log scripts/training/hpc/run_extract_araplm_c1_shard_nscc.pbs
qsub -v NUM_SHARDS=2,SHARD_ID=1 -o pplm_extract_araplm_c1_s1.log scripts/training/hpc/run_extract_araplm_c1_shard_nscc.pbs

# 2. Once both C1 extraction shards finish, run Zero-Shot + Retraining + Evaluation
qsub scripts/training/hpc/run_benchmark_araplm_nscc.pbs
```
