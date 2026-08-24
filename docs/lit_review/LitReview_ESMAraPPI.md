# Paper Information Extraction Template

## 1. Paper Metadata

- **Title:** Pre-trained protein language model sheds new light on the prediction of Arabidopsis protein–protein interactions
- **Authors:** Kewei Zhou, Chenping Lei, Jingyan Zheng, Yan Huang, Ziding Zhang
- **Publication Year:** 2023 (received 29 August 2023; accepted 28 November 2023; published online 7 December 2023)
- **Journal/Conference:** Plant Methods, 19:141 (doi: 10.1186/s13007-023-01119-6)
- **Model Name:** ESMAraPPI

---

## 2. Dataset Information

### 2.1 Data Sources

- **Positive PPI data source(s):** IntAct database (https://www.ebi.ac.uk/intact/home); only PPIs of type "direct interaction" or "physical association" retained; PPIs with MIscore <0.45 removed (i.e., MIscore ≥0.45 retained)
- **Negative sample generation method:** Proteins in positive samples first removed from the complete Arabidopsis protein list; remaining Arabidopsis proteins sharing ≥40% sequence identity with positive-sample proteins were filtered out; redundant proteins further removed via a 40% sequence identity cutoff, retaining 8,382 proteins; a protein list was then formed by mixing these 8,382 proteins with proteins in positive samples, and negative pairs were constructed via random pairing, controlling the positive:negative ratio at 1:10, yielding 77,290 random protein pairs not experimentally identified as PPIs
- **Proteome source:** Complete Arabidopsis protein list (specific database not explicitly named beyond IntAct-derived protein set)

### 2.2 Dataset Construction

- **Total number of positive PPI pairs:** 7,729
- **Total number of negative pairs:** 77,290
- **Positive-to-negative ratio:** 1:10
- **Species covered:** Arabidopsis thaliana
- **Sequence similarity threshold applied:** 40% sequence identity cutoff applied twice — once to filter proteins similar to positive-sample proteins, and again to remove redundant proteins among remaining candidates

### 2.3 Train/Test Split

- **Training dataset size and composition:** C1: 3,519 positive samples (1,415 proteins involved), 35,190 negative samples (7,068 proteins involved)
- **Test/validation dataset size and composition:** C2: 3,404 positive samples (1,781 proteins involved), 34,040 negative samples (10,586 proteins involved); C3: 806 positive samples (551 proteins involved), 8,060 negative samples (3,534 proteins involved)
- **Split strategy:** Followed Park and Marcotte's advice (2012) for pair-input evaluation scheme; C1 = training dataset, C2 and C3 = two independent test datasets
- **Specific split criteria:** Only one protein in each pair from C2 is allowed to appear in C1 (partially unseen); both proteins in each pair from C3 are unseen in C1 (fully unseen — most stringent)
- **Any temporal split:** Not specified

### 2.4 Cross-Species Datasets (if applicable)

- **Other species tested on:** Not specified — this paper does not report cross-species (e.g., rice) evaluation; it focuses solely on Arabidopsis thaliana. (Note: cross-species evaluation on rice appears in the related DeepAraPPI and ARACoFusion papers, but is not part of ESMAraPPI's own reported experiments.)
- **Source of cross-species data:** Not applicable
- **Number of samples per cross-species dataset:** Not applicable

---

## 3. Model Architecture

### 3.1 Protein Language Model (PLM)

- **PLM used:** ESM-1b (esm1b_t33_650M_UR50S), selected as the best-performing of nine ESM pLM variants tested; also compared against ProtTrans (ProtT5-XL-U50), UniRep, and TAPE
- **PLM source/training:** ESM models obtained from https://github.com/facebookresearch/esm/tree/v1.0.2; ESM-1b was pretrained on UniRef50, learning multi-scale representations including biochemical properties, remote homology, and within-family alignment, via masked language modeling with a deep Transformer architecture (Rives et al. 2021)
- **Embedding dimension:** 1280-dimensional (all nine ESM models tested produce 1280-dim vectors; note ProtT5 = 1024-dim, TAPE = 768-dim, UniRep = 1900-dim, per Methods — though the Results text describes TAPE embedding size as 768 in one place and the Methods section states TAPE=1900 and UniRep=768; the Methods section values are taken as authoritative: TAPE → 1900, UniRep → 768)
- **How embeddings are extracted:** Final layer's hidden parameters extracted, then the matrix was averaged along the first dimension (i.e., mean pooling across residues) to generate a fixed-length (1280-dim for ESM) feature vector per sequence
- **Which protein representation is used:** Mean-pooled (averaged) per-residue embeddings from the final transformer layer

### 3.2 Downstream Architecture

- **Architecture type:** Multilayer perceptron (MLP) — also benchmarked against Random Forest (RF) and Support Vector Machine (SVM) as alternative downstream classifiers, with MLP performing best across all pLMs tested
- **Detailed layer structure:** 4-layer MLP with 1024, 512, 128, and 16 nodes (implemented in PyTorch)
- **Activation functions:** Not explicitly specified for hidden layers; sigmoid function applied at the final output layer
- **How protein pairs are combined:** Hadamard (element-wise) product of the two protein feature vectors — chosen specifically to avoid order bias from protein pair input, rather than concatenation
- **Output layer:** Sigmoid function producing a prediction score between 0 and 1; a prediction score ≥0.5 corresponds to a positive interaction
- **Loss function:** Binary cross entropy (BCE) loss
- **Optimizer and training details:** Optimizer not explicitly named; training epochs = 40 (per Table 1, "Training epoch: 40"); total training time = 56 seconds (on the specified hardware); total prediction time on C3 test set = 0.1 seconds. Batch size and learning rate not specified.

### 3.3 Any Additional Components

- **Attention mechanisms:** None in the downstream model itself (attention is internal to the pretrained ESM-1b Transformer, not part of the custom downstream architecture)
- **Feature fusion methods:** Hadamard product used to fuse/combine the two protein embeddings into a single pair-representation vector before MLP input
- **Regularization techniques:** Not specified for the MLP (e.g., no dropout or weight decay mentioned for the ESMAraPPI model itself); for the SVM baseline, regularization parameter C was set to 1 with 'rbf' kernel and 'scale' kernel coefficient (optimized via grid search); for RF, n_estimators=100, max_depth=None (optimized via grid search)

---

## 4. Results

### 4.1 Performance on Primary Test Set (same species, stringent split)

- **Metric values (ESMAraPPI, best ESM-1b+MLP combination):**
    - C2: AUPR = 0.834, AUROC = 0.966, Accuracy = 0.957, Specificity = 0.994, MCC = 0.708, Recall = 0.589, Precision = 0.901
    - C3: AUPR = 0.810 (0.824 reported in Table 3 for a separate run — see note below), AUROC = 0.960, Accuracy = 0.954, Specificity = 0.994, MCC = 0.688, Recall = 0.557, Precision = 0.902
    - Note: There is a minor inconsistency in the paper — the abstract and Fig. 3/4 report AUPR = 0.834 (C2) / 0.810 (C3), while Table 3 reports AUPR = 0.824 (C2) / 0.810 (C3) for ESMAraPPI. Both are reported here as they appear in the source.
- **Comparison to baseline methods:**
    - **Other pLMs (with MLP), AUPR (C2/C3):** esm1_t34_670M_UR50D: 0.794/0.772; esm1_t34_670M_UR50S: 0.813/0.762; esm1_t34_670M_UR100: 0.807/0.769; esm1b_t33_650M_UR50S (ESM-1b, best): 0.834/0.810; esm1v_t33_650M_UR90S_1: 0.823/0.799; _2: 0.823/0.796; _3: 0.820/0.798; _4: 0.826/0.809; _5: 0.828/0.796
    - **Other pLM families (with MLP), AUPR/AUROC:** ProtTrans (ProtT5-XL-U50): AUPR 0.813 (C2)/0.756 (C3), AUROC 0.958 (C2)/0.942 (C3); UniRep: AUPR 0.761/0.687, AUROC 0.934/0.910; TAPE: AUPR 0.720/0.658, AUROC 0.920/0.893
    - **Baseline sequence encodings:** AAC+SVM (optimal for AAC): AUPR 0.519 (C2)/0.481 (C3), AUROC 0.852/0.824; DPC+RF (optimal for DPC): AUPR 0.646 (C2)/0.564 (C3), AUROC 0.884/0.845 — both much inferior to ESM-1b
    - **Machine learning algorithm comparison:** MLP > RF > SVM in combination with all pLMs tested, judged by AUPR/AUROC

### 4.2 Performance on Cross-Species Test Sets (if applicable)

- Not specified — this paper does not report any cross-species (e.g., rice) evaluation experiments. All experiments are within Arabidopsis thaliana (C1/C2/C3 splits).

### 4.3 Ablation Studies (if performed)

- **What was ablated:** Comparison across 9 different ESM pLM variants combined with 3 different ML algorithms (MLP, RF, SVM); comparison of ESM-1b against 3 other pLM families (ProtTrans, UniRep, TAPE); comparison of ESM-1b against 2 baseline encoding schemes (AAC, DPC) each combined with SVM/RF
- **Results of ablation:**
    - Across all 9×3 pLM/algorithm combinations, MLP consistently yielded the highest AUPR, followed by RF then SVM
    - ESM-1b (esm1b_t33_650M_UR50S) + MLP achieved the best overall performance among all 27 pLM×algorithm combinations (AUPR = 0.834 on C2, 0.810 on C3)
    - ESM-1b outperformed ProtTrans, UniRep, and TAPE when each was combined with MLP
    - ESM-1b (pLM-based) vastly outperformed traditional AAC/DPC encoding schemes regardless of downstream ML algorithm

### 4.4 Comparison with Existing Methods

- **Methods compared against:**
    - Generic sequence/structure-based PPI predictors: D-SCRIPT, RAPPPID, PIPR (all sequence-based), and TAGPPI (structure-based, using AlphaFold2 contact maps)
    - Arabidopsis-specific PPI predictors: AraPPINet (structure/function-based, interolog mapping-derived), DeepAraPPI (RCNN, Domain2vec, Go2vec individual predictors + integrated LR model)
- **Performance comparison:**
    - **vs. generic predictors (AUPR/AUROC, C2/C3):** ESMAraPPI: 0.834/0.810 (AUPR), 0.966/0.960 (AUROC); TAGPPI: 0.700/0.554 (AUPR), 0.925/0.873 (AUROC); PIPR: 0.588/0.387 (AUPR), 0.871/0.780 (AUROC); RAPPPID: 0.516/0.371 (AUPR), 0.852/0.800 (AUROC); D-SCRIPT: 0.292/0.291 (AUPR), 0.781/0.739 (AUROC). ESMAraPPI considerably outperformed all four existing generic predictors on both C2 and C3.
    - **Computational efficiency (Table 1):** ESM-1b+MLP: 40 training epochs, 56s total training time, 0.1s total prediction time (on C3); TAGPPI: 10 epochs, 9.29h training, 583s prediction; RAPPPID: 20 epochs, 1.12h training, 18s prediction; PIPR: 20 epochs, 700s training, 5s prediction; D-SCRIPT: 10 epochs, 7.22h training, 82s prediction. ESMAraPPI was markedly faster than all four in both training and prediction.
    - **vs. AraPPINet (Table 2, C2/C3):** ESMAraPPI: Accuracy 0.957/0.954, Specificity 0.994/0.994, MCC 0.708/0.688, Recall 0.589/0.557, Precision 0.901/0.902. AraPPINet: Accuracy 0.939/0.937, Specificity 0.999/0.999, MCC 0.551/0.534, Recall 0.337/0.318, Precision 0.966/0.966. ESMAraPPI outperformed AraPPINet on MCC (the authors' preferred comprehensive metric) in both datasets, despite AraPPINet's higher precision/specificity.
    - **vs. DeepAraPPI (Table 3, AUPR/AUROC, C2/C3):** DeepAraPPI_RCNN: 0.541/0.331 (AUPR), 0.852/0.778 (AUROC); DeepAraPPI_Domain2vec: 0.706/0.639, 0.884/0.845; DeepAraPPI_Go2vec: 0.771/0.709, 0.942/0.917; DeepAraPPI (integrated): 0.871/0.785, 0.978/0.944; ESMAraPPI: 0.824/0.810, 0.966/0.960. ESMAraPPI outperformed all three individual DeepAraPPI predictors on both test sets. On C2, the integrated DeepAraPPI model outperformed ESMAraPPI (0.871 vs 0.824 AUPR; 0.978 vs 0.966 AUROC), but on C3 (the more stringent, fully-unseen-protein test set), ESMAraPPI surpassed the integrated DeepAraPPI model (0.810 vs 0.785 AUPR; 0.960 vs 0.944 AUROC), indicating stronger generalization/extrapolation ability for ESMAraPPI on unseen proteins.

---

## 5. Evaluation Metrics Used

- **List all metrics reported in the paper:** Accuracy, Specificity, Precision, Recall (TPR), MCC (Matthews Correlation Coefficient), AUPR/AUPRC (area under precision-recall curve), AUROC (area under ROC curve)
- **Primary metric (the one emphasized by authors):** AUPR (area under the precision-recall curve) — emphasized due to the highly imbalanced positive:negative sample ratio (1:10); MCC also emphasized as "a more comprehensive measurement" when comparing against AraPPINet (since AUPR/AUROC could not be computed for that comparison)
- **Threshold used for binary classification:** 0.5 (a prediction score ≥0.5 corresponds to a positive interaction); this same default threshold (0.5) was also used for AraPPINet in the head-to-head comparison, as reported by AraPPINet's own web server

---

## 6. Key Findings & Claims

- **Main conclusion about model performance:** Sequence representations directly generated by large-scale pretrained pLMs (specifically ESM-1b), without any further feature engineering, can be successfully combined with a simple MLP to build a highly accurate Arabidopsis PPI predictor (ESMAraPPI). This dramatically outperforms models built on baseline sequence encoding schemes (AAC, DPC) and also outperforms several state-of-the-art generic (D-SCRIPT, RAPPPID, PIPR, TAGPPI) and plant-specific (AraPPINet, DeepAraPPI) PPI predictors, while also being computationally more efficient (faster training and prediction).
- **Generalization claims:** ESMAraPPI yielded an AUPR of 0.810 on the C3 independent test set, where both proteins in each pair are unseen in the training dataset, which the authors state "suggest[s] its strong generalization and extrapolating ability." Notably, on this most stringent C3 test set, ESMAraPPI outperformed the integrated DeepAraPPI model (0.810 vs 0.785 AUPR), which the authors interpret as evidence that ESMAraPPI is "more competitive and will be more reliable in practical applications," particularly for previously uncharacterized proteins. No cross-species (e.g., rice) generalization was tested in this paper.
- **Limitations mentioned by authors:** Not extensively discussed as explicit limitations in the paper; however, the authors implicitly note that on C2 (less stringent test set), the integrated DeepAraPPI model outperforms ESMAraPPI, suggesting ESMAraPPI's relative advantage is most pronounced specifically in the hardest (unseen-protein) generalization setting rather than uniformly across all test conditions. No explicit "Limitations" section is present.

---

## 7. Implementation Details (if available)

- **Code availability:** Yes — GitHub repository: https://github.com/keiwo/ESMAraPPI (code and datasets freely available)
- **Hardware used for training:** High-performance computer with 20 cores CPU, 256GB RAM, and a Tesla V100 GPU
- **Training time:** 56 seconds total training time (40 epochs) for ESM-1b+MLP; total prediction time on the C3 test dataset = 0.1 seconds (compared to 9.29h/7.22h/1.12h/700s training and 583s/82s/18s/5s prediction times for TAGPPI/D-SCRIPT/RAPPPID/PIPR respectively, per Table 1)