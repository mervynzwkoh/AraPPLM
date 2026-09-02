# Technical Report: Plant-Pretrained PPLM-PPI on the DeepAraPPI Benchmark Suite

**Author / Project:** AIS5281 - Plant Protein Language Modeling (`AraPPLM`)  
**Date:** September 2026  
**Model Configuration:** PPLM backbone (frozen, 650M params) + MLP head retrained from random initialisation on DeepAraPPI C1  
**Datasets:** DeepAraPPI *Arabidopsis thaliana* Partitions (C2, C3) & *Oryza sativa* (Rice Monocot Benchmark)  
**Evaluation Standard:** Park & Marcotte (2012) Pair-Input Partition Scheme at 1:10 Positive-to-Negative Ratio  
**Training Methodology:** See [`ppi_head_retraining_methodology.md`](ppi_head_retraining_methodology.md)  

---

## 1. Data Provenance & Study Lineage

```
=============================================================================================================================
                                          DATA PROVENANCE & BENCHMARK METHODOLOGY MATRIX
=============================================================================================================================
Model / Benchmark                 Architecture Authors         Evaluation Performed By               Training Domain
-----------------------------------------------------------------------------------------------------------------------------
DeepAraPPI (C2/C3/Rice)           Zheng et al. (2023)          Zheng et al. (The Plant Journal)      Arabidopsis C1
ARACoFusion (Rice)                Sarkar & Sarkar (2026)       Sarkar & Sarkar (bioRxiv 2026)        Arabidopsis C1
ESMAraPPI on Rice*                Zhou et al. (2023)           Sarkar & Sarkar (bioRxiv 2026)*       Arabidopsis C1
-----------------------------------------------------------------------------------------------------------------------------
PPLM (Zero-Shot)                  Liu et al. (2026)            This Work (AIS5281, Aug 2026)         Non-Plant Pretraining
PPLM (Plant-Pretrained)           Liu et al. (2026)            This Work (AIS5281, Sep 2026)         Arabidopsis C1 (head only)
=============================================================================================================================
*ESMAraPPI Rice: re-implemented by Sarkar & Sarkar (2026) on the Zheng et al. (2023) Rice dataset.
```

**Key distinction — "Plant-Pretrained" vs. supervised baselines:** All supervised baseline models (DeepAraPPI, ESMAraPPI, ARACoFusion) trained their *entire* architecture on plant C1 data. In the PPLM plant-pretrained configuration, only the 5-layer MLP prediction head (~3.9M params) was trained on C1, with the 650M-parameter backbone frozen. The MLP weights were **randomly initialised** (not fine-tuned from the human-trained checkpoint). See [`ppi_head_retraining_methodology.md`](ppi_head_retraining_methodology.md) §4 for full training protocol details.

---

## 2. Executive Summary

This report presents the results of **plant-pretrained PPLM-PPI** — where the PPLM backbone remains frozen and a new MLP prediction head is trained from random initialisation on the DeepAraPPI C1 training set (31,284 pairs). Results are compared against the zero-shot PPLM baseline and all published supervised models.

```
=============================================================================================================================
                                AUPRC PERFORMANCE ON DEEPARAPPI HELD-OUT TEST SUITE
=============================================================================================================================
Task                               DeepAraPPI (2023)   ARACoFusion (2026)   PPLM Zero-Shot   PPLM Plant-Pretrained   Δ vs Zero-Shot
-----------------------------------------------------------------------------------------------------------------------------
Task 2: C2 (One Unseen Protein)         0.8970               —                 0.5738              0.8738              +0.3000 (+52.3%)
Task 3: C3 (Both Unseen Proteins)       0.8250               —                 0.5525              0.8118              +0.2593 (+46.9%)
Task 4: Rice (Cross-Species Transfer)   0.3050             0.3519              0.4297              0.3555              -0.0742 (-17.3%)
=============================================================================================================================
```

### Key Findings

1. **Massive improvement on Arabidopsis tasks:** Plant-pretraining the MLP head delivers +52.3% relative AUPRC improvement on C2 (0.5738 → 0.8738) and +46.9% on C3 (0.5525 → 0.8118), nearly closing the gap with fully supervised DeepAraPPI.

2. **Near-parity with DeepAraPPI on C2 and C3:** The plant-pretrained PPLM achieves **0.8738 AUPRC** on C2 (vs. DeepAraPPI's 0.8970, a gap of only 0.0232) and **0.8118 AUPRC** on C3 (vs. DeepAraPPI's 0.8250, a gap of only 0.0132). This is remarkable given that PPLM uses only primary sequence features, while DeepAraPPI additionally leverages Gene Ontology (`GO2vec`: 0.8710) and protein domain (`Domain2vec`: 0.7800) graph embeddings.

3. **Rice cross-species performance decreases:** The plant-pretrained model scores **0.3555 AUPRC** on Rice, down from the zero-shot **0.4297**. Training on Arabidopsis C1 data introduces dicot-specific bias that reduces cross-species transfer, consistent with findings reported for all other supervised models.

---

## 3. Detailed Empirical Results

### 3.1 Plant-Pretrained PPLM Performance Summary

Across all **105,875 held-out test pairs** in the DeepAraPPI benchmark suite:

| Metric | Task 2: C2 (One Unseen) | Task 3: C3 (Both Unseen) | Task 4: Rice (Cross-Species) |
| :--- | :--- | :--- | :--- |
| **Total Test Pairs** | **66,055** | **33,099** | **6,721** |
| **Positives / Negatives (1:10)** | 6,005 / 60,050 | 3,009 / 30,090 | 611 / 6,110 |
| **AUPRC (Primary Metric)** | **0.8738** | **0.8118** | **0.3555** |
| **AUROC** | **0.9678** | **0.9561** | **0.7122** |
| **Accuracy ($\tau = 0.5$)** | 96.58% | 95.52% | 86.58% |
| **Precision ($\tau = 0.5$)** | **90.98%** | **88.03%** | **30.63%** |
| **Specificity ($\tau = 0.5$)** | **99.31%** | **99.20%** | **91.47%** |
| **Sensitivity / Recall ($\tau = 0.5$)** | 69.24% | 58.66% | 37.64% |
| **F1 Score ($\tau = 0.5$)** | 0.7864 | 0.7040 | 0.3377 |
| **Matthews Correlation (MCC)** | 0.7766 | 0.6973 | 0.2657 |
| **Optimal Threshold ($\tau^*$)** | **0.3556** | **0.3108** | **0.8227** |
| **Optimal F1 Score ($F_1^*$)** | **0.8038** | **0.7318** | **0.3642** |

### 3.2 Head-to-Head: Zero-Shot vs. Plant-Pretrained PPLM

| Metric | C2 Zero-Shot | C2 Pretrained | Δ | C3 Zero-Shot | C3 Pretrained | Δ |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **AUPRC** | 0.5738 | **0.8738** | **+0.3000** | 0.5525 | **0.8118** | **+0.2593** |
| **AUROC** | 0.8828 | **0.9678** | +0.0850 | 0.8710 | **0.9561** | +0.0851 |
| **Precision** | 85.14% | **90.98%** | +5.84pp | 82.72% | **88.03%** | +5.31pp |
| **Recall** | 24.53% | **69.24%** | **+44.71pp** | 24.03% | **58.66%** | **+34.63pp** |
| **F1** | 0.3809 | **0.7864** | **+0.4055** | 0.3724 | **0.7040** | **+0.3316** |
| **MCC** | 0.4339 | **0.7766** | +0.3427 | 0.4218 | **0.6973** | +0.2755 |
| **Optimal F1** | 0.5452 | **0.8038** | +0.2586 | 0.5305 | **0.7318** | +0.2013 |
| **Specificity** | 99.57% | 99.31% | -0.26pp | 99.50% | 99.20% | -0.30pp |

**The most significant gain is in recall:** Plant-pretraining increased sensitivity from ~24% to ~69% (C2) and ~59% (C3) at the default threshold, meaning the model now correctly identifies a far greater proportion of true interactions while maintaining high precision (>88%).

### 3.3 Comprehensive Multi-Model Benchmark Comparison (AUPRC)

```
=============================================================================================================================
                                COMPREHENSIVE BENCHMARK COMPARISON TABLE (METRIC: AUPRC)
=============================================================================================================================
Model / Method                     Evaluated By                    Task 2 (C2)   Task 3 (C3)   Task 4 (Rice)
-----------------------------------------------------------------------------------------------------------------------------
DeepAraPPI (Integrated LR)         Zheng et al. (2023)               0.8970         0.8250         0.3050
GO2vec (GO Graph MLP)              Zheng et al. (2023)               0.8710         0.8030         0.2650
PPLM (Plant-Pretrained, Ours)      This Study (AIS5281, Sep 2026)    0.8738         0.8118         0.3555
Domain2vec (Domain Graph MLP)      Zheng et al. (2023)               0.7800         0.6810         0.2790
RCNN (Sequence Word2Vec+GRU)       Zheng et al. (2023)               0.7460         0.4810         0.2480
ARACoFusion (ESM-1b + CrossAttn)   Sarkar & Sarkar (2026)              —               —            0.3519
ESMAraPPI on Rice                  Sarkar & Sarkar (2026)*              —               —            0.2938
-----------------------------------------------------------------------------------------------------------------------------
PPLM (Zero-Shot, Ours)             This Study (AIS5281, Aug 2026)    0.5738         0.5525         0.4297 ★
-----------------------------------------------------------------------------------------------------------------------------
                                                                                 ★ = Best score for that task
=============================================================================================================================
```

---

## 4. In-Depth Analysis & Discussion

### 4.1 Arabidopsis Task 2 (C2): One Unseen Protein

Plant-pretrained PPLM achieves **0.8738 AUPRC** on C2, within 0.0232 of DeepAraPPI's 0.8970. This is a significant result because:

- **Sequence-only features:** PPLM uses only primary protein sequences, while DeepAraPPI fuses `GO2vec` (0.8710), `Domain2vec` (0.7800), and `RCNN` (0.7460) — three separate feature modalities — via logistic regression.
- **PPLM outperforms all individual DeepAraPPI components except GO2vec:** The plant-pretrained PPLM (0.8738) exceeds `Domain2vec` (0.7800) and `RCNN` (0.7460), and matches `GO2vec` (0.8710) using only sequence information.
- **Recall is the primary driver:** The zero-shot model's recall was only 24.5% (the model was too conservative). After plant-pretraining, recall jumped to 69.2% while precision *increased* from 85.1% to 91.0%.

### 4.2 Arabidopsis Task 3 (C3): Both Proteins Unseen

On the hardest Arabidopsis task where both proteins are completely absent from C1 training, plant-pretrained PPLM scores **0.8118 AUPRC** — within 0.0132 of DeepAraPPI's 0.8250. Key observations:

- **Outperforms all sequence-only methods by a wide margin:** PPLM (0.8118) exceeds RCNN (0.4810) by +68.8% relative and surpasses Domain2vec (0.6810) by +19.2% relative.
- **Gap to DeepAraPPI Integrated is only 1.6%:** DeepAraPPI's advantage on C3 comes from `GO2vec` (0.8030) — GO annotations provide interaction evidence even when the proteins themselves are unseen. PPLM nearly matches this with sequence alone.
- **Comparison to zero-shot:** The +0.2593 AUPRC gain from C1 plant-pretraining on this fully unseen task indicates that the PPLM backbone's learned representations benefit from even indirect supervision — seeing other Arabidopsis proteins during MLP training helps calibrate the decision boundary.

### 4.3 Rice Cross-Species Transfer (Task 4): Performance Trade-off

Plant-pretrained PPLM scores **0.3555 AUPRC** on Rice, compared to zero-shot **0.4297** — a 17.3% relative decrease. This is an important finding:

- **Dicot-specific training bias:** By training the MLP head exclusively on *Arabidopsis* (a dicot), the decision boundary becomes tuned to dicot-specific interaction signatures. This reduces generalization to *Oryza sativa* (a monocot), consistent with findings for DeepAraPPI (0.3050), ARACoFusion (0.3519), and ESMAraPPI on Rice (0.2938) — all of which were trained on Arabidopsis and performed poorly on Rice.
- **Zero-shot PPLM remains the best Rice model (0.4297):** The zero-shot configuration, which learned general biophysical interaction rules from human/yeast/E. coli data, still achieves the highest Rice AUPRC of any published model.
- **Threshold shift explains part of the gap:** The optimal threshold for the plant-pretrained model on Rice shifts to 0.8227 (from 0.4625 zero-shot), indicating the model has become more confident overall but less calibrated for the Rice protein distribution.

```
=============================================================================================================================
                              RICE CROSS-SPECIES: ALL MODELS RANKED BY AUPRC
=============================================================================================================================
Model / Method                        Training Domain         AUPRC     AUROC    Precision   Recall    Specificity   F1
-----------------------------------------------------------------------------------------------------------------------------
PPLM Zero-Shot (Ours, Aug 2026)       Non-Plant Pretrain      0.4297    0.7561    0.6367     0.2897     0.9835      0.3982
PPLM Plant-Pretrained (Ours, Sep)     Arabidopsis C1 (head)   0.3555    0.7122    0.3063     0.3764     0.9147      0.3377
ARACoFusion                           Arabidopsis C1          0.3519    0.6864    0.3176     0.3748     0.9195      0.3438
DeepAraPPI Integrated                 Arabidopsis C1          0.3050      —         —          —          —           —
ESMAraPPI on Rice                     Arabidopsis C1          0.2938    0.7034    0.3638     0.2995     0.9476      0.3285
=============================================================================================================================
```

**Implication:** For maximum cross-species generalization, the zero-shot PPLM configuration should be used. For Arabidopsis-specific tasks (C2, C3), plant-pretraining provides substantial gains.

### 4.4 Calibration & Threshold Analysis

| Setting | C2 $\tau^*$ | C2 Best F1 | C3 $\tau^*$ | C3 Best F1 | Rice $\tau^*$ | Rice Best F1 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Zero-Shot** | 0.0924 | 0.5452 | 0.0963 | 0.5305 | 0.4625 | 0.4048 |
| **Plant-Pretrained** | 0.3556 | **0.8038** | 0.3108 | **0.7318** | 0.8227 | **0.3642** |

The optimal threshold shifted from ~0.09 (zero-shot) to ~0.31–0.36 (plant-pretrained) on Arabidopsis tasks, indicating that plant-pretraining has substantially improved the model's calibration. The model now assigns higher scores to true interactions, making the 0.5 default threshold more practical (the zero-shot model required threshold tuning to ~0.09 for reasonable F1).

---

## 5. Summary & Recommendations

### 5.1 When to Use Each Configuration

| Configuration | Best For | AUPRC Range |
|--------------|----------|-------------|
| **PPLM Plant-Pretrained** | Arabidopsis PPI prediction (C2, C3 tasks) | 0.81–0.87 |
| **PPLM Zero-Shot** | Cross-species PPI prediction (Rice, other non-Arabidopsis) | 0.43 (Rice) |
| **DeepAraPPI Integrated** | Arabidopsis PPI when GO + domain annotations available | 0.83–0.90 |

### 5.2 Next Steps

1. **ESMAraPPI C1 training:** Run the same pipeline on ESMAraPPI C1 (38,709 pairs with 40% sequence redundancy filtering) to evaluate whether the more stringent negative sampling improves generalization.
2. **Transfer learning from human weights:** The current plant-pretrained MLP was initialised randomly. Initialising from the human-trained `ppi_models.pkl` checkpoint and fine-tuning on plant data may provide additional gains.
3. **Focal loss / class-weighted BCE ablation:** The current plain BCE loss may under-weight the minority positive class. Focal loss ($\gamma=2$) or class-weighted BCE ($w_{pos}=10$) could improve recall without sacrificing specificity.
4. **LoRA backbone adaptation:** Apply Low-Rank Adaptation to the PPLM cross-attention layers to learn plant-specific interaction patterns within the backbone itself.
5. **Joint/mixed strategy for Rice:** Train on a combination of Arabidopsis C1 + synthetic rice augmentation, or use an ensemble of zero-shot + plant-pretrained models for Rice tasks.

---

## 6. Artifact & File References

* **Benchmark Summary (Pretrained):** [`results/DeepAraPPI/pretrained/benchmark_summary.csv`](../results/DeepAraPPI/pretrained/benchmark_summary.csv)
* **Task 2 Predictions (Pretrained):** [`results/DeepAraPPI/pretrained/deepara_c2_plant_scores.csv`](../results/DeepAraPPI/pretrained/deepara_c2_plant_scores.csv)
* **Task 3 Predictions (Pretrained):** [`results/DeepAraPPI/pretrained/deepara_c3_plant_scores.csv`](../results/DeepAraPPI/pretrained/deepara_c3_plant_scores.csv)
* **Task 4 Predictions (Pretrained):** [`results/DeepAraPPI/pretrained/deepara_rice_plant_scores.csv`](../results/DeepAraPPI/pretrained/deepara_rice_plant_scores.csv)
* **Zero-Shot Results:** [`results/DeepAraPPI/zero-shot/`](../results/DeepAraPPI/zero-shot/)
* **Training Methodology:** [`docs/ppi_head_retraining_methodology.md`](ppi_head_retraining_methodology.md)
* **Zero-Shot Benchmark Report:** [`docs/benchmark_analysis_deeparappi_vs_pplm.md`](benchmark_analysis_deeparappi_vs_pplm.md)
* **Training Logs:** [`logs/training/`](../logs/training/)
