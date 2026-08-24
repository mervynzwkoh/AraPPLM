# Paper Information Extraction Template

## 1. Paper Metadata

- **Title:** ARACoFusion: Uncertainty-aware calibrated deep learning for protein-protein interaction network prediction in Arabidopsis thaliana
- **Authors:** Dipayan Sarkar, Chiranjib Sarkar
- **Publication Year:** 2026 (preprint posted May 26, 2026)
- **Journal/Conference:** bioRxiv preprint (doi: 10.64898/2026.05.22.727120) — not peer-reviewed
- **Model Name:** ARACoFusion

---

## 2. Dataset Information

### 2.1 Data Sources

- **Positive PPI data source(s):** IntAct database (obtained from Zhou et al. [21], i.e., the ESMAraPPI paper), filtered for physical association interactions with MI score < 0.45
- **Negative sample generation method:** Three-step process — (1) positive samples removed from the complete protein pool, and remaining sequences with >40% sequence identity to positive samples removed; (2) proteins with similar sequence identity removed using a 40% cutoff, yielding 8,382 proteins; (3) negative pairs constructed by random pairing among those 8,382 proteins and proteins from positive samples
- **Proteome source:** Not explicitly named (proteome pool derived from the IntAct-sourced protein set obtained via Zhou et al.); no specific proteome database (e.g., UniProt) named for this step

### 2.2 Dataset Construction

- **Total number of positive PPI pairs:** 7,729
- **Total number of negative pairs:** 77,290
- **Positive-to-negative ratio:** 1:10
- **Species covered:** Arabidopsis thaliana (primary); Oryza sativa (rice) for cross-species testing
- **Sequence similarity threshold applied:** 40% sequence identity cutoff used when constructing the negative sample pool (to exclude near-duplicates of positive-sample proteins)

### 2.3 Train/Test Split

- **Training dataset size and composition:** C1 = 38,709 total samples (3,519 positive, 35,190 negative; 1:10 ratio; 90.91% negative class imbalance)
- **Test/validation dataset size and composition:** C2 = 37,444 total samples (3,404 positive, 34,040 negative, 90.90% negative); C3 = 8,866 total samples (806 positive, 8,060 negative, 90.91% negative)
- **Split strategy:** Park and Marcotte's method (2012) for pair-input data partitioning; also stratified 5-fold cross-validation (on merged C1+C2+C3) for robustness assessment
- **Specific split criteria:** C1 used for training; C2 and C3 used for validation (consistent with the Park & Marcotte C1/C2/C3 scheme used in DeepAraPPI — C2 typically shares one protein with C1 "medium difficulty," C3 shares no proteins with C1 "high difficulty," though the paper does not restate this distinction explicitly beyond citing the method)
- **Any temporal split:** Not specified

### 2.4 Cross-Species Datasets (if applicable)

- **Other species tested on:** Oryza sativa (rice)
- **Source of cross-species data:** Rice PPI dataset downloaded from four public databases (DIP, MINT, BioGRID, IntAct), obtained from Zheng et al. [20] (i.e., the DeepAraPPI paper); self, non-physical, and redundant interactions removed to yield 611 positive interactions; negatives randomly selected at 1:10 ratio
- **Number of samples per cross-species dataset:** Rice dataset total = 6,721 samples (611 positive, 6,110 negative; 90.90% negative class imbalance)

Additionally, a network-level validation dataset was used: derived from STRING database (Arabidopsis thaliana, physical links only, confidence ≥0.900, AB unidirectional interactions retained), yielding a high-confidence network of 2,477 proteins and 9,549 physical interactions; a connected subnetwork of 43 nodes and 103 edges was used for case-study evaluation.

---

## 3. Model Architecture

### 3.1 Protein Language Model (PLM)

- **PLM used:** ESM-1b (ESM 1b-t33-650M parameter model); also benchmarked against ESM2 (T6, T12, T30, T33 variants) and ProtT5 (T5-T33) for comparison, with ESM-1b selected as optimal
- **PLM source/training:** Developed by Meta AI; pretrained with a masked language model objective on the UniRef50 dataset, comprising 40 million proteins
- **Embedding dimension:** 1280-dimensional per-residue embeddings from the final layer; averaged to a 1280-dimensional per-protein embedding
- **How embeddings are extracted:** Mean pooling — per-residue contextual embeddings (dimension d=1280) from the final layer of ESM-1b-650M are averaged across sequence length L to produce a fixed-size per-protein embedding: z = (1/L)Σxᵢ
- **Which protein representation is used:** Mean-pooled (average) of all per-residue token embeddings from the final transformer layer (33 layers total, 20 attention heads per the figure, feed-forward intermediate dim 5120, hidden dim 1280)

### 3.2 Downstream Architecture

- **Architecture type:** Custom deep neural network combining a reciprocal cross-attention encoder (Siamese-like, bidirectional cross-attention), a latent interaction projector (nonlinear feature adaptation), and a multi-stage fusion-based classifier head (feedforward MLP)
- **Detailed layer structure:**
    - Input: two 1280-dim protein embeddings (p1, p2) from ESM-1b
    - Reciprocal cross-attention encoder: multi-head cross-attention where each protein's query attends to the other's key-value pairs (Q, K, V generated via learnable weight matrices Wq, Wkv), producing contextual embeddings c1, c2 ∈ ℝ^1280; optimal number of attention heads = 8 (Optuna-tuned; default was 4)
    - Latent interaction projector: two-layer fully connected network per protein embedding with GELU activation, projecting to 512-dim latent representations a1, a2 ∈ ℝ^512 (Latent(pᵢ) = W2·GELU(W1·pᵢ + b1) + b2)
    - Explicit pairwise relational features: element-wise product P_prod = p1⊙p2 (∈ℝ^1280) and element-wise absolute difference P_diff = |p1−p2| (∈ℝ^1280)
    - Final concatenated representation: Z = [p1; p2; c1; c2; P_prod; P_diff; a1; a2]
    - Classifier head: h1 = Dropout(GELU(Linear(BatchNorm(Z)))); h2 = Dropout(GELU(Linear(h1))); ŷ = σ(Linear(h2)) ∈ (0,1)
    - Exact hidden-layer dimensions of the classifier head (sizes of h1, h2) are not explicitly specified numerically
- **Activation functions:** GELU (Gaussian Error Linear Unit) used throughout feature adaptation and classifier head; sigmoid (σ) at output layer
- **How protein pairs are combined:** Multiple fusion strategies combined via concatenation — reciprocal cross-attention contextual embeddings, element-wise product, element-wise absolute difference, latent projections, and original embeddings all concatenated into one interaction representation Z; skip-connections preserve original embeddings and element-wise features
- **Output layer:** Sigmoid function producing an interaction probability score in (0,1)
- **Loss function:** Focal loss with label smoothing (L_Focal-LS), combined with a variance-based uncertainty regularization term: Loss_total = L_Focal-LS + λ·V[ŷ1,ŷ2,ŷ3] (variance computed across 3 Monte Carlo forward passes)
- **Optimizer and training details:** AdamW optimizer; learning rate tuned via Optuna (best value = 0.000681383, vs. default 0.001); other Optuna-tuned hyperparameters: attention heads = 8 (default 4), focal loss gamma = 3.892573 (default 2.0), dropout = 0.118479 (default 0.3), label smoothing = 0.0514500 (default 0.1), uncertainty weight λ = 0.579057 (default 0.1). Batch size and number of epochs not specified.

### 3.3 Any Additional Components

- **Attention mechanisms:** Reciprocal multi-head cross-attention encoder — each protein embedding attends to the key-value pairs derived from the other protein's embedding (bidirectional/reciprocal cross-attention), computed with H heads and scaled dot-product attention (Attention(Q,K,V) formula with √(d/H) scaling)
- **Feature fusion methods:** Multi-source fusion combining raw embeddings, cross-attention contextual embeddings, element-wise product, element-wise absolute difference, and nonlinear latent projections, all concatenated (Z) before the classifier head; skip-connections included
- **Regularization techniques:** Dropout (tuned via Optuna, optimal 0.118479); uncertainty-aware variance regularization (Monte Carlo sampling, 3 forward passes, variance penalty term λ); focal loss with label smoothing to address class imbalance; temperature scaling (post-hoc calibration, not training regularization but a calibration technique) — optimal T = 0.58, reducing Expected Calibration Error (ECE) from 0.034 (raw) to 0.020 (scaled); BatchNorm applied before the classifier head

---

## 4. Results

### 4.1 Performance on Primary Test Set (same species, stringent split)

- **Metric values (C2, with Optuna-optimized hyperparameters):** ACC = 0.9657, SPEC = 0.9875, PREC = 0.8568, SN (Sensitivity/Recall) = 0.7471, F1 = 0.7982, MCC = 0.7817, NPV = 0.975, AP = 0.8546, AUROC = 0.9548, BACC = 0.8673, AUPRC = 0.8546
- **Metric values (C3, with Optuna-optimized hyperparameters):** ACC = 0.9594, SPEC = 0.9897, PREC = 0.8644, SN = 0.6563, F1 = 0.7461, MCC = 0.7326, NPV = 0.9664, AP = 0.8067, AUROC = 0.9308, BACC = 0.823, AUPRC = 0.8066
- **Metric values (5-fold cross-validation on merged C1+C2+C3):** ACC = 0.9978, SPEC = 0.9987, PREC = 0.9867, SN = 0.9893, F1 = 0.988, MCC = 0.9868, NPV = 0.9989, AP = 0.9967, AUROC = 0.9993, BACC = 0.994, AUPRC = 0.9967
- **Comparison to baseline methods:**
    - vs. ESMAraPPI (C2): ACC 0.9518, SPEC 0.9966, PREC 0.9361, SN 0.5038, F1 0.6551, MCC 0.6669, NPV 0.9526, AP 0.8359, AUROC 0.9665, BACC 0.7502, AUPRC 0.8358
    - vs. ESMAraPPI (C3): ACC 0.9501, SPEC 0.9967, PREC 0.9354, SN 0.4851, F1 0.6389, MCC 0.6534, NPV 0.9509, AP 0.8087, AUROC 0.9593, BACC 0.7409, AUPRC 0.8085
    - vs. AraPPINet (C2): ACC 0.939, SPEC 0.999, PREC 0.966, SN 0.337, MCC 0.551 (other metrics N/A — taken from original publication)
    - vs. AraPPINet (C3): ACC 0.937, SPEC 0.999, PREC 0.966, SN 0.318, MCC 0.534 (others N/A)
    - vs. DeepAraPPI variants (AUPRC/AUROC): RCNN 0.541/0.852 (C2), 0.331/0.778 (C3); Domain2vec 0.706/0.884 (C2), 0.639/0.845 (C3); GO2vec 0.771/0.942 (C2), 0.709/0.917 (C3)
    - vs. general PPI methods (AUPRC/AUROC, C2/C3): D-SCRIPT 0.292/0.781 (C2), 0.291/0.739 (C3); RAPPPID 0.516/0.852 (C2), 0.371/0.800 (C3); PIPR 0.588/0.871 (C2), 0.387/0.780 (C3); TAGPPI 0.700/0.925 (C2), 0.554/0.873 (C3)
    - ARACoFusion outperformed all compared methods on AUPRC and generally on sensitivity/F1/MCC/BACC across both C2 and C3

### 4.2 Performance on Cross-Species Test Sets (if applicable)

- **Species: Oryza sativa (rice)**
    - **AUROC:** 0.6864
    - **AUPRC:** 0.3519
    - **Other metrics:** ACC = 0.87, SPEC = 0.9195, PREC = 0.3176, SN (Recall) = 0.3748, F1 = 0.3438, MCC = 0.2734, NPV = 0.9363, AP = 0.3525, BACC = 0.6471
    - Comparison — ESMAraPPI on rice: ACC 0.8887, SPEC 0.9476, PREC 0.3638, SN 0.2995, F1 0.3285, MCC 0.27, NPV 0.9312, AP 0.2949, AUROC 0.7034, BACC 0.6236, AUPRC 0.2938
    - Comparison — other cross-species/general methods (AUPRC only): RCNN 0.248, Domain2vec 0.279, GO2vec 0.265, RF+DPC 0.171, LR (DeepAraPPI integrated) 0.305, ESMAraPPI 0.2938; ARACoFusion 0.3519 (highest)
    - Confusion matrix (rice, ARACoFusion): TP=220, FP=408, FN=391, TN=5702 (vs. ESMAraPPI: TP=183, FP=320, FN=428, TN=5790)

### 4.3 Ablation Studies (if performed)

- **What was ablated:** Four components removed one at a time from the full architecture: (a) reciprocal cross-attention encoder; (b) element-wise product and absolute difference features; (c) temperature scaling module; (d) uncertainty-aware training (variance regularization). All ablation variants trained on C1, evaluated on C2 and C3.
- **Results of ablation:**
    - **C2:** Full model (original): SN 0.7471, F1 0.7982, MCC 0.7817, BACC 0.8673, AUPRC 0.8546, AUROC 0.9548. Removing cross-attention: SN drops to 0.7051, F1 to 0.787 (largest recall drop, though AUPRC slightly higher at 0.8515). Removing element-wise product/difference: MCC drops to 0.7686, AUPRC 0.8553. Removing temperature scaling: AUROC drops marginally to 0.9559→0.9548 comparison note (AUPRC actually increases slightly to 0.8555, calibration ability reduced though). Removing uncertainty-aware training: minor negative impact across most metrics (AUPRC 0.8635).
    - **C3:** Full model: SN 0.6563, F1 0.7461, MCC 0.7326, AUPRC 0.8066, AUROC 0.9308. Removing cross-attention gives lowest SN (0.6079) and F1 (0.7275), confirming its centrality. Removing temperature scaling and uncertainty regularization slightly increased AUPRC (0.8138 and 0.8207, respectively) versus the full model's 0.8066, though these components' primary benefit is calibration/reliability rather than raw discrimination.
    - Overall conclusion: the reciprocal cross-attention encoder contributes the most to sensitivity/recall and overall discrimination; element-wise features encode explicit pairwise compatibility; uncertainty training and temperature scaling primarily improve calibration/stability rather than raw classification metrics.

### 4.4 Comparison with Existing Methods

- **Methods compared against:** AraPPINet, DeepAraPPI (RCNN, Domain2vec, GO2vec, LR variants), ESMAraPPI, D-SCRIPT, RAPPPID, PIPR, TAGPPI (general-purpose sequence-based PPI predictors)
- **Performance comparison:**
    - ARACoFusion significantly outperforms ESMAraPPI and AraPPINet on C2/C3 splits (higher sensitivity/F1/AUPRC/BACC), though ESMAraPPI and AraPPINet show higher precision/specificity (more conservative predictions)
    - ARACoFusion outperforms all DeepAraPPI variants on AUPRC and AUROC on both C2 (0.8546 vs. best DeepAraPPI GO2vec 0.771) and C3 (0.8066 vs. 0.709)
    - ARACoFusion vastly outperforms general-purpose PPI predictors (D-SCRIPT, RAPPPID, PIPR, TAGPPI) on AUPRC/AUROC for both C2 and C3, with TAGPPI as the strongest general baseline (still substantially below ARACoFusion)
    - On 5-fold CV, ARACoFusion (AUROC 0.9993, AUPRC 0.9967) outperforms ESMAraPPI (AUROC 0.9938, AUPRC 0.9571)
    - On a STRING-derived Arabidopsis subnetwork (43 nodes, 103 edges), ARACoFusion achieved 70% accuracy vs. ESMAraPPI's 66%

---

## 5. Evaluation Metrics Used

- **List all metrics reported in the paper:** Accuracy (ACC), Sensitivity/Recall (SN), Specificity (SP), Precision/Positive Predictive Value (PREC/PPV), Average Precision (AP), Negative Predictive Value (NPV), F1-score, Matthews Correlation Coefficient (MCC), Balanced Accuracy (BACC), AUROC, AUPRC, Expected Calibration Error (ECE, for calibration assessment)
- **Primary metric (the one emphasized by authors):** AUPRC (explicitly used as the Optuna objective function to maximize; emphasized as most informative under class imbalance)
- **Threshold used for binary classification:** Not explicitly specified as a fixed value (e.g., no stated 0.5 cutoff); calibrated probabilities derived via temperature scaling (T=0.58) but the specific decision threshold for converting probability to binary label is not stated

---

## 6. Key Findings & Claims

- **Main conclusion about model performance:** ARACoFusion, combining ESM-1b embeddings with a reciprocal cross-attention encoder, latent interaction projection, and multi-source feature fusion, substantially outperforms existing Arabidopsis-specific predictors (AraPPINet, DeepAraPPI, ESMAraPPI) and general-purpose sequence-based PPI predictors (D-SCRIPT, RAPPPID, PIPR, TAGPPI) across AUPRC, BACC, MCC, and other imbalance-sensitive metrics, while also providing calibrated, uncertainty-aware probability outputs
- **Generalization claims:** The model shows "robust cross-species generalization," outperforming ESMAraPPI and other cross-species baselines on rice (AUPRC 0.3519 vs. ESMAraPPI's 0.2938), though the authors note the cross-species performance is still much lower than same-species (Arabidopsis) performance, reflecting a persistent generalization gap. t-SNE visualizations show ARACoFusion restructures the latent space to improve class separability even on unseen rice sequences.
- **Limitations mentioned by authors:** Some false predictions persist in network-level evaluation, "likely attributable to the negative-class bias in the training data (C1), which was generated with a large proportion of presumed non-interacting pairs," introducing conservatism in classifying positive interactions especially in densely connected modules. The authors also note that "domain adaptation remains a challenge" for cross-species prediction, and cross-species performance, while improved over baselines, remains far below within-species performance.

---

## 7. Implementation Details (if available)

- **Code availability:** Web server at https://ARAcofusion.compbiosysnbu.in/; documentation and usage instructions stated to be "publicly available on GitHub" (specific repository URL not given in the text)
- **Hardware used for training:** Workstation with Intel Core i5-12400F CPU and NVIDIA RTX 3060 GPU (12 GB VRAM)
- **Training time:** Not specified