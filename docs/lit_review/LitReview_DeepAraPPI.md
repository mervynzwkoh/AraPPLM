# Paper Information Extraction Template

## 1. Paper Metadata

- **Title:** Deep learning-assisted prediction of protein–protein interactions in _Arabidopsis thaliana_
- **Authors:** Jingyan Zheng, Xiaodi Yang, Yan Huang, Shiping Yang, Stefan Wuchty, Ziding Zhang
- **Publication Year:** 2023
- **Journal/Conference:** The Plant Journal, 114, 984–994 (doi: 10.1111/tpj.16188)
- **Model Name:** DeepAraPPI

---

## 2. Dataset Information

### 2.1 Data Sources

- **Positive PPI data source(s):** BioGRID, DIP, IntAct, MINT, and TAIR (no version numbers given, only URLs)
- **Negative sample generation method:** Randomly selected protein pairs from the reference proteome pool, constrained so sampled proteins neither share the same subcellular localization nor belong to the pool of known interactions
- **Proteome source:** Arabidopsis reference proteome sequences downloaded from UniProt (28,361 sequences after removing proteins with fewer than 40 amino acids or non-standard amino acids)

### 2.2 Dataset Construction

- **Total number of positive PPI pairs:** 49,398 experimentally verified PPIs between 10,330 Arabidopsis proteins (before quality filtering); high-quality subset = 11,858 PPIs (score ≥0.72 via HIPPIE scoring scheme); low-quality subset = 37,540 PPIs (score <0.72)
- **Total number of negative pairs:** Varies by task (e.g., Task1: 118,580 negative pairs)
- **Positive-to-negative ratio:** 1:10
- **Species covered:** Primarily Arabidopsis thaliana; cross-species testing on Oryza sativa (rice)
- **Sequence similarity threshold applied:** Not specified (no explicit sequence identity/redundancy threshold mentioned)

### 2.3 Train/Test Split

- **Training dataset size and composition:** Task1: random 80% of 11,858 high-quality PPIs as positives, plus corresponding negatives (1:10 ratio)
- **Test/validation dataset size and composition:** Task1: remaining 20% of PPIs as independent test set
- **Split strategy:** Three difficulty-level tasks (Task1, Task2, Task3) based on dataset partition method proposed by Park & Marcotte (2012)
- **Specific split criteria:** High-quality PPI set segmented into three subsets C1 (2,844 PPIs), C2 (6,005 PPIs), C3 (3,009 PPIs). C1 used as positive training set for Task2 and Task3. C2 = positive test samples for Task2 (each PPI in C2 shares only one protein with C1 — "medium difficulty," partially unseen proteins). C3 = positive test samples for Task3 (each PPI in C3 shares no protein with C1 — "high difficulty," fully unseen proteins). Task1 uses a simple random 80/20 split of all 11,858 high-quality PPIs (low difficulty). Negative samples in C1/C2/C3 follow the same 1:10 sampling strategy as Task1.
- **Any temporal split:** Yes — for comparison with three pre-2018 existing methods (AraPPINet, AtPIN, AtPID), the 11,858 high-quality PPIs were repartitioned: 8,997 PPIs published before 2018 as training positives, 2,861 PPIs published after 2018 as independent test positives.

### 2.4 Cross-Species Datasets (if applicable)

- **Other species tested on:** Oryza sativa (rice)
- **Source of cross-species data:** Four public databases (DIP, MINT, BioGRID, IntAct) plus literature (Wierbowski et al., 2020)
- **Number of samples per cross-species dataset:** 611 rice PPIs between 555 proteins as positive samples (after removing self-interactions, non-physical interactions, redundant interactions); negatives sampled at 1:10 ratio. A further split: 80% training / remainder as independent test set was used when building a hybrid Arabidopsis+rice training set.

---

## 3. Model Architecture

### 3.1 Protein Language Model (PLM)

- **PLM used:** DeepAraPPI itself does **not** use a pretrained protein language model (e.g., ESM/ProtT5). It uses a custom-trained **word2vec** embedding for amino acids (RCNN branch), **Domain2vec** (node2vec-based embedding on a domain-interaction network), and **GO2vec** (node2vec-based embedding on a GO term graph). PLMs are only discussed as related work/future direction (e.g., D-SCRIPT, which uses a pretrained protein language model, is used as a comparison baseline).
- **PLM source/training:** word2vec was trained on protein sequences from the Uniref50 database using a continuous bag-of-words (CBOW) architecture, implemented via the Python Gensim library
- **Embedding dimension:** 32-dimensional per amino acid (word2vec, window size = 3, chosen via fivefold cross-validation on Task1)
- **How embeddings are extracted:** Sequences truncated/zero-padded to fixed length L = 2000, representing each protein as an L × 32 array
- **Which protein representation is used:** Not a single pooled vector — full per-residue L×32 array fed into CNN-GRU layers

### 3.2 Downstream Architecture

- **Architecture type:** Ensemble of three separate predictors combined via logistic regression: (i) Siamese RCNN (CNN + bidirectional GRU) on sequence/word2vec input; (ii) MLP on Domain2vec embeddings; (iii) MLP on GO2vec embeddings
- **Detailed layer structure:**
    - RCNN: 1D convolution (kernel size 3) → max pooling → bidirectional GRU (50 channels) → element-wise multiplication layer → MLP with three fully connected layers → softmax output
    - Domain2vec: node2vec-derived protein embeddings concatenated into a pair vector → MLP with three fully connected layers → interaction probability
    - GO2vec: node2vec-derived GO-graph embeddings concatenated into a pair vector → MLP with three fully connected layers → interaction probability
    - Final: individual scores (S_RCNN, S_Domain2vec, S_GO2vec) combined into a vector → Logistic Regression model → final prediction score
- **Activation functions:** Not specified in detail (softmax used at RCNN output layer; specific activations for MLP hidden layers not stated)
- **How protein pairs are combined:** Element-wise multiplication of the pair of protein embedding vectors (RCNN); concatenation of embedding vectors into a pair vector (Domain2vec and GO2vec)
- **Output layer:** Softmax function (RCNN); MLP outputs interaction probability (Domain2vec, GO2vec); final LR model outputs comprehensive prediction score
- **Loss function:** Not specified
- **Optimizer and training details:** Not specified (learning rate, batch size, epochs not given). LR model used L2 penalty and a linear solver (scikit-learn), with hyperparameters tuned via GridSearchCV with fivefold cross-validation

### 3.3 Any Additional Components

- **Attention mechanisms:** None described (GRU-based sequential modeling only)
- **Feature fusion methods:** Late fusion — individual predictor scores (S_RCNN, S_Domain2vec, S_GO2vec) combined into a vector for logistic regression
- **Regularization techniques:** L2 penalty (in the logistic regression model only); no dropout/weight decay specified for RCNN or MLP components

---

## 4. Results

### 4.1 Performance on Primary Test Set (same species, stringent split)

- **Metric values (AUPRC, Table 1):**

||Task1 (low)|Task2 (medium)|Task3 (high)|
|---|---|---|---|
|RCNN|0.925|0.746|0.481|
|Domain2vec|0.868|0.780|0.681|
|GO2vec|0.939|0.871|0.803|
|Logistic Regression (DeepAraPPI)|0.965|0.897|0.825|

- Other metrics reported: TPR, FPR, precision (via PR curves); AUC-ROC not explicitly reported (AUPRC is the primary metric)
- **Comparison to baseline methods:** RF models with AC/CT/DPC encodings for Task1–3 — RF+AC: 0.875/0.657/0.276; RF+CT: 0.892/0.712/0.392; RF+DPC: 0.903/0.720/0.434 (vs. DeepAraPPI: 0.965/0.897/0.825). DeepAraPPI outperformed all RF-based methods across all three tasks (Figure 3a–c).

### 4.2 Performance on Cross-Species Test Sets (if applicable)

- **Species: Oryza sativa (rice)**
    - **AUPRC:** RCNN = 0.248, Domain2vec = 0.279, GO2vec = 0.265, Logistic Regression (integrated) = 0.305, RF+DPC baseline = 0.171
    - **AUC-ROC:** Not specified
    - **Other metrics:** When GO2vec was retrained on a hybrid Arabidopsis+rice training set, AUPRC improved to 0.561 on the independent rice test set (vs. 0.285 for the original Arabidopsis-only GO2vec model on the same test set)

### 4.3 Ablation Studies (if performed)

- **What was ablated:** Comparison of the three individual baseline predictors (RCNN, Domain2vec, GO2vec) versus the integrated LR model, across all three difficulty tasks and cross-species prediction
- **Results of ablation:** The integrated LR model consistently outperformed any single predictor in every task (Task1: 0.965 vs best individual 0.939; Task2: 0.897 vs 0.871; Task3: 0.825 vs 0.803; rice: 0.305 vs 0.279), demonstrating the integrative strategy's value. GO2vec was uniformly the best single predictor; RCNN was most sensitive to increasing prediction difficulty (largest performance drop from Task1 to Task3); Domain2vec was the most robust individual model across difficulty levels.

### 4.4 Comparison with Existing Methods

- **Methods compared against:** RF+AC, RF+CT, RF+DPC (traditional ML baselines); AraPPINet, AtPIN, AtPID (existing Arabidopsis-specific PPI predictors); D-SCRIPT (a human PPI predictor using a pretrained protein language model)
- **Performance comparison:**
    - Against AraPPINet: at fixed FPR = 0.063%, DeepAraPPI achieved TPR = 24.1% vs AraPPINet's TPR = 9.1%
    - Against AtPIN and AtPID: DeepAraPPI achieved much higher TPR at matched FPR thresholds (0.031% and 0.025% respectively) — exact TPR values shown only in Figure 3d, not stated numerically in text
    - Against D-SCRIPT: DeepAraPPI AUPRC = 0.828 vs D-SCRIPT AUPRC = 0.708 (on a repartitioned dataset used for the three-method comparison)

---

## 5. Evaluation Metrics Used

- **List all metrics reported in the paper:** TPR (recall), FPR, precision, Precision-Recall (PR) curves, Area Under the PR Curve (AUPRC)
- **Primary metric (the one emphasized by authors):** AUPRC
- **Threshold used for binary classification:** Not a fixed single threshold — the paper instead tunes/reports performance across multiple FPR thresholds (e.g., 0.063%, 0.05%, 0.031%, 0.025%) for comparison purposes and provides FPR threshold options (e.g., 0.05%) in the webserver

---

## 6. Key Findings & Claims

- **Main conclusion about model performance:** DeepAraPPI, by integrating sequence (RCNN/word2vec), domain (Domain2vec), and GO (GO2vec) information via logistic regression, outperforms individual baseline models and existing state-of-the-art Arabidopsis PPI prediction methods (AraPPINet, AtPIN, AtPID) as well as RF-based and non-plant deep learning methods (D-SCRIPT) across difficulty levels
- **Generalization claims:** DeepAraPPI shows better cross-species predictive ability in rice than traditional ML methods (RF), but cross-species performance is "dramatically inferior" to same-species (Arabidopsis) performance. Authors attribute this partly to insufficient PPI data in the rice test set and to models being highly biased toward Arabidopsis proteins, limiting generalizability. Using a hybrid Arabidopsis+rice training set substantially improved rice prediction (AUPRC 0.285 → 0.561), suggesting training-set similarity to the test set matters for cross-species performance. Authors state generalizability of the model for cross-species application "remains an open issue."
- **Limitations mentioned by authors:** (1) Cross-species prediction performance is far behind that of Arabidopsis; (2) Domain2vec, GO2vec, and the integrated LR model can only predict proteins present in their pre-trained corpus — proteins not in the corpus cannot be predicted by those models (only RCNN, sequence-only, can handle novel proteins); (3) authors note future improvements could come from incorporating AlphaFold2-derived structural information and pretrained protein language models (e.g., Rives et al. 2021 ESM-type models) to improve cross-species prediction, which they had not yet incorporated in this work.

---

## 7. Implementation Details (if available)

- **Code availability:** Yes — GitHub: https://github.com/zjy1125/DeepAraPPI; also downloadable from online platform http://zzdlab.com/deeparappi/ (webserver implemented with CentOS 7.4 and Apache 2.4.6)
- **Hardware used for training:** Not specified
- **Training time:** Not specified