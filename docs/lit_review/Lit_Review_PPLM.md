# Paper Information Extraction Template

#### 1. The Pretraining Distribution Match

During PPLM pretraining, the cross-chain inter-protein attention mechanism was trained on paired sequences with **full combined lengths $L_A + L_B \le 1024$**.

- When evaluating an uncropped pair where $L_A + L_B \gg 1024$ (e.g., $3,000$ residues), two problems occur:
    1. **Out-of-Distribution Attention Scaling:** The model's inter-chain attention weights are evaluated on token distances far beyond what the pre-trained attention heads saw during pretraining.
    2. **Memory Blowup:** The $33 \times 20 \times L_{total}^2$ attention tensor exceeds 40 GB VRAM.

#### 2. Combined Length Cropping Strategy ($L_A + L_B \le 1020$)

To strictly honor the paper's methodology and maintain full compatibility with the pre-trained weights:

- **Maximum Combined Tokens:** $\mathbf{1020}$ residues (leaving 4 tokens for the `<cls>`, `<sep>`, and boundary markers to make the total input tensor exactly $\le 1024$).
- **Proportional Cropping (Preserving Sequence Representation):** If $\text{len}(\text{seqA}) + \text{len}(\text{seqB}) > 1020$: $$\text{Budget}_A = \max\left(50, \left\lfloor 1020 \times \frac{\text{len}(\text{seqA})}{\text{len}(\text{seqA}) + \text{len}(\text{seqB})} \right\rfloor\right)$$ $$\text{Budget}_B = 1020 - \text{Budget}_A$$ $$\text{seqA} = \text{seqA}[:\text{Budget}_A], \quad \text{seqB} = \text{seqB}[:\text{Budget}_B]$$

This ensures:

1. **Biological balance:** Neither protein is completely clipped if one is substantially longer than the other (each gets at least a proportional share of the context).
2. **Total token budget $\le 1024$**: Perfectly adheres to the PPLM paper pretraining specification.
3. **Guaranteed Peak VRAM $\le 6 \text{ GB}$**: Attention matrices for $L \le 1024$ require only $\sim 2.6 \text{ GB}$ of VRAM, running with zero OOM on the 40 GB A100 GPU.

## 1. Paper Metadata

- **Title:** A paired sequence language model for protein-protein interaction modeling
- **Authors:** Jun Liu, Hungyu Chen, Yang Zhang
- **Publication Year:** 2026 (received 18 July 2025; accepted 17 February 2026; published online 10 March 2026)
- **Journal/Conference:** Nature Communications, 17:3733 (doi: 10.1038/s41467-026-70457-5)
- **Model Name:** PPLM (Protein Pair Language Model), with three task-specific downstream variants: PPLM-PPI (binary interaction prediction), PPLM-Affinity (binding affinity prediction), PPLM-Contact / PPLM-Contact2 (inter-protein contact and interface residue prediction)

_Note: This paper is not plant-specific — it covers general (non-plant) species: H. sapiens, M. musculus, D. melanogaster, C. elegans, S. cerevisiae, and E. coli for PPI prediction, plus PDB/STRING complexes broadly for pretraining, and antibody–antigen/TCR–pMHC complexes for affinity prediction. Included here for architectural/methodological comparison purposes within the broader PPI-PLM literature review._

---

## 2. Dataset Information

### 2.1 Data Sources

**For PPLM pretraining (language model):**

- **Positive PPI data source(s):** Protein structure complexes from the Protein Data Bank (PDB) released before January 1, 2024, and interaction sequences from the STRING database (2023 version referenced)
- **Negative sample generation method:** Not applicable — this is unsupervised masked language model pretraining on interacting pairs only, not a binary classification dataset
- **Proteome source:** Not species-specific; broad multi-species protein pairs from PDB and STRING

**For PPLM-PPI (binary interaction prediction):**

- **Positive PPI data source(s):** Benchmark datasets adopted from D-SCRIPT, spanning six species (H. sapiens, M. musculus, D. melanogaster, C. elegans, S. cerevisiae, E. coli)
- **Negative sample generation method:** Random negative sampling (inherited from D-SCRIPT benchmark construction); duplicate, erroneous, and invalid samples arising from random negative sampling were identified and removed (per Supplementary Method S8)
- **Proteome source:** Not separately specified beyond the D-SCRIPT benchmark species

**For PPLM-Affinity:**

- **Positive PPI data source(s):** PPB-Affinity dataset (Liu et al. 2024) — a large curated resource compiled from multiple public affinity databases
- **Negative sample generation method:** Not applicable (regression task on binding affinity values, not classification)
- **Proteome source:** PDB entries (3,027 distinct entries)

**For PPLM-Contact / PPLM-Contact2:**

- **Positive PPI data source(s):** Same datasets as used in DeepInter (Lin et al. 2023) — non-redundant homodimers and heterodimers curated from PDB
- **Negative sample generation method:** Not applicable (contact prediction, not binary classification of interacting/non-interacting pairs)
- **Proteome source:** PDB-derived homodimer/heterodimer complexes; additional CASP13–CASP16 targets for independent test sets

### 2.2 Dataset Construction

**PPLM pretraining dataset:**

- **Total number of positive PPI pairs:** Over 3.3 million protein sequence pairs; after clustering/redundancy removal via pair-level clustering (MMseqs2 + custom clustering): 25,245 heteromeric clusters + 23,082 homomeric clusters from PDB, plus 629,045 clusters from STRING = 672,372 total clusters used for training; 4,678 single-pair clusters held out for validation
- **Total number of negative pairs:** Not applicable (unsupervised pretraining)
- **Positive-to-negative ratio:** Not applicable
- **Species covered:** Broad, multi-species (not specified exhaustively; PDB + STRING spans many organisms)
- **Sequence similarity threshold applied:** Non-redundant clustering via MMseqs2 plus a custom clustering procedure; validation set (4,678 pairs) was non-redundant to any training pairs

**PPLM-PPI dataset (six species):**

- **Total number of positive PPI pairs / negative pairs:** Not explicitly stated in the main text (detailed in Table S19, not provided); ratio of negative to positive maintained at 10:1
- **Positive-to-negative ratio:** 10:1 (negative:positive)
- **Species covered:** H. sapiens (train/validation), M. musculus, D. melanogaster, C. elegans, S. cerevisiae, E. coli (test species)
- **Sequence similarity threshold applied:** For the sequence-similarity-stratified generalization evaluation, test interactions were grouped by maximum single-sequence identity of either protein to any human protein in the training set (identity intervals analyzed, not a hard exclusion cutoff)

**PPB-Affinity dataset (for PPLM-Affinity):**

- **Total number of positive PPI pairs:** 12,052 interaction samples from 3,027 distinct PDB entries
- **Sequence similarity threshold applied:** All PDB entries clustered using US-align with a complex-level TM-score cutoff of 0.8 to prevent homology-based leakage; complexes with different numbers of chains clustered separately

**Contact prediction dataset (for PPLM-Contact):**

- **Training set:** 3,504 homodimers + 1,881 heterodimers
- **Validation set:** 296 homodimers + 96 heterodimers
- **Test sets:** 300 homodimers (Homodimer300), 99 heterodimers (Heterodimer99), plus independent CASP-derived sets: 43 homodimer targets (CASP_Homodimer43) and 20 heterodimer targets (CASP_Heterodimer20) from CASP13–CASP16

### 2.3 Train/Test Split

- **Training dataset size and composition:** Varies by task (see above for each of the four datasets: PPLM pretraining, PPLM-PPI, PPLM-Affinity, PPLM-Contact)
- **Test/validation dataset size and composition:** Varies by task (see above)
- **Split strategy:**
    - PPLM pretraining: cluster-based non-redundant hold-out (4,678 validation pairs non-redundant to training)
    - PPLM-PPI: trained/validated via 10-fold cross-validation on H. sapiens dataset only; tested cross-species on the other five species (train on one species, test on five others — this is itself a form of cross-species generalization test built into the primary evaluation)
    - PPLM-Affinity: five-fold cross-validation, with samples sharing the same PDB ID grouped into the same fold (to prevent leakage); folds also regrouped by structural similarity (TM-score cutoff 0.8) for the final PPB-Affinity comparison
    - PPLM-Contact: fixed train/validation/test split by dimer type (homodimer-specific and heterodimer-specific training sets), plus independent CASP-derived test sets for robustness checking
- **Specific split criteria:** Not a C1/C2/C3 Park & Marcotte-style scheme (unlike the plant PPI papers); instead uses (a) cross-species holdout for PPLM-PPI, (b) structural-similarity-based fold grouping for PPLM-Affinity, and (c) PDB-curated non-redundant dimer sets plus temporally/community-defined CASP target sets for PPLM-Contact
- **Any temporal split:** Yes — PPLM pretraining used only PDB complexes released before January 1, 2024; CASP_Homodimer43 and CASP_Heterodimer20 test sets are drawn from CASP13 to CASP16, representing chronologically later/independent targets for real-world robustness assessment

### 2.4 Cross-Species Datasets (if applicable)

- **Other species tested on:** For PPLM-PPI: model trained/validated on H. sapiens only, then tested on M. musculus, D. melanogaster, C. elegans, S. cerevisiae, and E. coli (five cross-species test sets) — this is the paper's primary evaluation paradigm, explicitly described as testing "generalizes robustly to unseen proteins under a strict cross-species inductive setting"
- **Source of cross-species data:** D-SCRIPT benchmark datasets (Sledzieski et al. 2021)
- **Number of samples per cross-species dataset:** Not explicitly stated in main text (detailed in Supplementary Table S19, not provided in extracted content); ratio of 10:1 negative:positive maintained across all species datasets

---

## 3. Model Architecture

### 3.1 Protein Language Model (PLM)

- **PLM used:** PPLM — a novel Protein Pair Language Model, initialized from and built upon ESM2 (650M-parameter version)
- **PLM source/training:** Initialized from the 650M-parameter ESM2 model; further pretrained by the authors on a composite dataset of >3.3 million protein sequence pairs from PDB (physical interface complexes, pre-2024) and STRING (interaction sequences), using masked language modeling with pairs sampled from PDB:STRING clusters at a 1:2 ratio
- **Embedding dimension:** Not explicitly stated as a single number in the main text (inherited from ESM2-650M's architecture; dimensionality details for downstream tasks are in Supplementary Tables S20–22, not provided in extracted content)
- **How embeddings are extracted:** Both sequences in a pair are independently tokenized (with BOS/EOS tokens marking chain boundaries), then concatenated and passed through 33 stacked transformer blocks; final-layer embeddings and both intra- and inter-protein attention matrices are extracted as features for downstream tasks
- **Which protein representation is used:** Full per-residue sequence embeddings AND both intra-protein and inter-protein attention matrices (not just a single pooled vector) — for PPLM-PPI specifically, both max-pooling and mean-pooling are applied to these features along the sequence dimension to form two parallel representation branches

### 3.2 Downstream Architecture

**PPLM (core language model) — architecture type:** Transformer-based paired-sequence language model with a hybrid intra-/inter-protein attention mechanism

- **Detailed layer structure:** 33 serially connected transformer blocks, each with a tailored multi-head attention module and a feed-forward network; within each attention module, two separate attention matrices are computed — one using rotary positional embeddings (RoPE) for intra-chain residue pairs, and one without positional encoding (ROW) for inter-protein residue pairs — combined via a learnable weight (W_lh, per layer l and head h) and an inter-protein binary mask (M_ij): A_lhij = W_lh · [M_ij · A_row_lhij + (1−M_ij) · A_rope_lhij] + M_ij
- **Activation functions:** Sigmoid used in the attention-gating formula (mask/weight combination); ReLU used in downstream MLPs (see PPLM-PPI/Affinity below)
- **How protein pairs are combined:** Joint co-representation — both sequences are concatenated (with BOS/EOS boundary tokens) and processed together through shared transformer blocks, allowing explicit modeling of inter-protein residue-residue relationships via the dedicated inter-protein attention pathway (rather than combining two independently-computed single-protein embeddings post hoc)
- **Loss function (pretraining):** Masked language modeling (MLM) loss — weighted average of cross-entropy losses computed separately over masked residues in each of the two sequences: L_MLM = −(1/|M1|)Σlog p(x_i^(1)|x\M) − (1/|M2|)Σlog p(x_i^(2)|x\M)
- **Optimizer and training details (pretraining):** AdamW optimizer, β1=0.9, β2=0.98; learning rate linearly warmed up to 1×10⁻⁶ over first 2,000 steps, then linearly decayed to 5×10⁻⁷ over remaining training; trained for 50,000 total steps on 4 NVIDIA A100 GPUs with gradient accumulation step of 32; masking strategy differs by source — PDB pairs: 30% of interface residues masked vs. 15% of non-interface residues; STRING pairs: uniform 15% random masking; sequences >1024 residues cropped (biased toward interface residues for PDB; proportional continuous fragments for STRING)

**PPLM-PPI downstream architecture:**

- **Architecture type:** Dual-branch (max-pooling + mean-pooling) MLP ensemble on top of frozen/extracted PPLM features
- **Detailed layer structure:** From PPLM, three feature types are extracted (sequence embeddings, intra-protein attention matrices ×2, inter-protein attention matrix); each is pooled via max-pooling and mean-pooling independently (two parallel branches); each pooled representation is projected via a dedicated linear layer, then passed into its own MLP; each MLP consists of 5 linear layers, with the first 4 followed by layer normalization and ReLU activation, ending in sigmoid activation for interaction probability; final score = average of outputs from the two pooling branches (dimensionality details in Supplementary Tables S20–21, not fully provided)
- **Activation functions:** ReLU (in the first four of five MLP layers), sigmoid (output layer)
- **How protein pairs are combined:** Combined jointly within PPLM itself (via inter-protein attention); downstream MLP operates on pooled features already reflecting pair-level information (not a simple post hoc combination like Hadamard product)
- **Output layer:** Sigmoid, producing interaction probability; final prediction = average of max-pooling-branch and mean-pooling-branch outputs
- **Loss function:** Binary cross-entropy: L_PPI = −(1/N)Σ[y_i log p_i + (1−y_i) log(1−p_i)]
- **Optimizer and training details:** AdamW, learning rate 5×10⁻⁵, batch size 32, trained on single NVIDIA A100 GPU; 10-fold cross-validation on H. sapiens dataset, 12 epochs per fold; best-checkpoint (highest validation AUPRC) selected per fold; top 5 best-performing models across all folds ensembled for final inference

**PPLM-Affinity downstream architecture:**

- **Architecture type:** Fine-tuned final transformer block + pooling + fully connected regression head
- **Detailed layer structure:** Final (last) transformer block of PPLM fine-tuned on affinity data; full-length embeddings aggregated via max pooling; passed through two fully connected layers with an intermediate ReLU activation to predict a single scalar binding affinity value
- **Activation functions:** ReLU (intermediate layer)
- **How protein pairs are combined:** Receptor and ligand sequences (concatenated into single sequences if multi-chain) processed jointly through PPLM's paired-sequence architecture; embeddings aggregated via max pooling before the regression head
- **Output layer:** Linear (regression) output — a continuous binding affinity (ΔG) value
- **Loss function:** Mean squared error (MSE): L_Affinity = (1/N)Σ(ŷ_i − y_i)²
- **Optimizer and training details:** AdamW, learning rate 1×10⁻⁴; five-fold cross-validation (samples sharing PDB ID grouped in same fold); 15 epochs per fold; single NVIDIA A100 GPU

**PPLM-Contact downstream architecture:**

- **Architecture type:** Inter-protein transformer with parallel triangle multiplication, cross-attention, and self-attention modules (AlphaFold-inspired geometric deep learning architecture)
- **Detailed layer structure:** Two ResNet modules independently encode intra-protein and inter-protein features (intra-protein ResNet shares tied parameters across both chains); followed by 12 inter-protein transformer blocks, each containing: 2 parallel triangle multiplication modules, 2 parallel cross-attention modules, 2 parallel self-attention modules, and a transition layer (2 linear transformations)
- **Activation functions:** Sigmoid (gating in triangle multiplication/cross-attention/self-attention formulas), softmax (attention weight computation), LayerNorm used throughout (φ function)
- **How protein pairs are combined:** Explicit inter-protein pair representation z_ij updated iteratively via triangle multiplication (integrating intra-chain interactions r_in and cross-chain interactions z_nj under geometric triangle constraints), cross-attention (inter-protein pair as query, intra-protein pairs as key/value), and self-attention (attending over other inter-protein pairs sharing residue i or j, modulated by Gaussian-kernel-transformed intra-protein distances)
- **Output layer:** Predicted inter-protein contact probability p_ij per residue pair
- **Loss function:** Focal Loss (to address contact/non-contact class imbalance): L_Contact = −(1/(L1·L2))ΣΣ α(1−p_ij)^γ log(p_ij), with α=0.25, γ=1.5
- **Optimizer and training details:** AdamW, initial learning rate 1×10⁻³, decayed by factor 0.98 per epoch; 100 epochs total; single NVIDIA A100 GPU; monomers >256 residues cropped to a 256-residue fragment maximizing interface residue coverage; contacts defined as residue pairs with minimum heavy-atom distance <8Å; top 5 models by validation top-L precision ensembled for inference
- **PPLM-Contact2 (enhanced variant):** Integrates inter-protein distance maps extracted from predicted complex structures (AlphaFold2.3, AlphaFold3, DMFold) as an additional feature; predictions from the three structure sources are ensembled as final output

### 3.3 Any Additional Components

- **Attention mechanisms:** (1) Core PPLM: hybrid intra-/inter-protein attention with RoPE (intra-chain) vs. no positional encoding (inter-protein), combined via learnable per-layer/per-head weights and an explicit inter-protein binary mask; (2) PPLM-Contact: triangle multiplication (geometric constraint-based), cross-attention (inter-protein pair attends to intra-protein pairs), and self-attention (inter-protein pair attends to other inter-protein pairs, modulated by Gaussian-kernel intra-protein distances) — all applied in parallel per protein and combined
- **Feature fusion methods:** PPLM-PPI fuses embedding + intra-protein attention (×2 proteins) + inter-protein attention via parallel max/mean pooling branches, later averaged; PPLM-Contact fuses PPLM-derived inter-protein attention with MSA-derived features (PSSMs, DCA row scores, ESM-MSA-1b embeddings/attention) and monomer distance maps (experimental or AlphaFold2-predicted); PPLM-Contact2 further fuses in AlphaFold2.3/AlphaFold3/DMFold-predicted complex-structure distance maps, ensembled across the three structure sources
- **Regularization techniques:** Layer normalization (in PPLM-PPI MLP branches and in PPLM-Contact's φ transformations); Focal Loss (α=0.25, γ=1.5) specifically to counter severe class imbalance in contact prediction; cropping strategies to manage sequence length/GPU memory (implicit regularization); model ensembling (top-5 models per task) as a form of variance reduction; gradient accumulation (step=32) during PPLM pretraining

---

## 4. Results

### 4.1 Performance on Primary Test Set (same species, stringent split)

**PPLM language modeling quality (perplexity, vs. ESM2):**

- Random masking: PPLM average perplexity 7.30 (homomeric PDB), 5.08 (heteromeric PDB), 4.50 (STRING-derived) vs. ESM2's 8.40, 6.73, 5.70 — reductions of 13.1%, 24.5%, 21.0% respectively (median reductions 15.1%, 28.5%, 19.0%; Wilcoxon p = 2.96×10⁻¹³⁹, 8.07×10⁻¹⁷⁷, <10⁻³⁰⁰)
- Head-to-head: PPLM outperformed ESM2 on 91.9% (4298/4678) of sequence pairs under random masking
- Interface masking (Dual mode): PPLM 6.79 avg (median 5.78) vs. ESM2 8.50 (7.75); p=9.90×10⁻²²³; PPLM better in 80.7% of cases
- Interface masking (Single mode): PPLM 7.81 avg (median 5.62) vs. ESM2 10.36 (8.89); p=1.17×10⁻²⁰⁸; PPLM better in 78.9% of cases
- Case study (1Y9B homodimer): top 20 attention-ranked residue pairs — 90% were true heavy-atom contacts (80% Cβ–Cβ contacts); 52/62 (83.9%) experimentally determined interface residues recoverable from attention rankings

**PPLM-PPI (binary interaction prediction, five test species):**

- AUPRC: M. musculus 0.920, D. melanogaster 0.906, C. elegans 0.883, S. cerevisiae 0.745, E. coli 0.784
- F1-score improvements over TUnA (second-best): 4.8%, 8.6%, 4.8%, 16.9%, 15.5% (species order as above)
- Overall mean±SD AUPRC across species: PPLM-PPI 0.848±0.078
- **Comparison to baseline methods:** TUnA 0.778±0.109; ESMDNN-PPI 0.762±0.106; D-SCRIPT 0.532±0.072; Topsy-Turvy 0.548±0.118. Relative to TUnA/ESMDNN-PPI: F1 improved 10.1–10.5%, AUPRC improved 9.6–11.9% (Cohen's d ≈1.742–6.586). Relative to D-SCRIPT/Topsy-Turvy: F1 improved 32.5–72.3%, AUPRC improved 58.7–60.9% (Cohen's d ≈2.062–6.183)

**PPLM-Affinity (binding affinity prediction, PPB-Affinity dataset, entire dataset):**

- PCC = 0.643±0.058, SRCC = 0.636±0.082
- vs. ESM2-Affinity: PCC 0.548±0.061, SRCC 0.547±0.053 (PPLM-Affinity improvements: +17.3% PCC, +16.4% SRCC; Cohen's d = 1.274, 0.958)
- vs. structure-based PPB-Affinity model: PCC 0.545±0.072, SRCC 0.540±0.088 (PPLM-Affinity improvements: +18.0% PCC, +17.8% SRCC; Cohen's d = 1.326, 1.064)
- Mean absolute error: PPLM-Affinity 1.68 (σ=1.44) vs. ESM2-Affinity 1.85 (σ=1.52) vs. PPB-Affinity 1.85 (σ=1.63)
- RMSE: PPLM-Affinity 2.312±0.297 vs. ESM2-Affinity 2.476±0.376 vs. PPB-Affinity 2.463±0.394
- Antibody–antigen subgroup: PPLM-Affinity PCC 0.380, SRCC 0.404 (improvements of 117.1%/111.5% over ESM2-Affinity; 60.3%/57.2% over PPB-Affinity)
- TCR–pMHC subgroup: PPLM-Affinity PCC 0.366, SRCC 0.312 (improvements of 144.0%/127.7% over ESM2-Affinity; 39.2%/36.8% over PPB-Affinity)

**PPLM-Contact (inter-protein contact prediction, top-L precision):**

- Homodimer300 (experimental monomers): PPLM-Contact 77.8% vs. DeepInter 68.9%, PLMGraph-Inter 50.6%, CDPred 63.0%, DeepHomo2.0 51.4%, GLINTER 34.5% (improvements 12.8–125.3%, all p=3.63×10⁻²¹ to 1.29×10⁻⁴⁷)
- Homodimer300 (AlphaFold2-predicted monomers): PPLM-Contact 66.6%, +10.4% over DeepInter (p=3.96×10⁻¹³)
- CASP_Homodimer43 (experimental): 65.2% precision, improvements 6.0–157.7% over baselines
- CASP_Homodimer43 (AlphaFold2-predicted): improvements of 12.9%, 39.5%, 46.3%, 72.3%, 140.0% over baselines
- Heterodimer99: PPLM-Contact 48.9% (experimental) / 45.1% (AlphaFold2-predicted); improvements 37.0–182.7% over competitors
- CASP_Heterodimer20: 45.6% (experimental) / 40.5% (AlphaFold2-predicted)
- **Comparison to baseline methods:** PLMGraph-Inter, DeepInter, CDPred, GLINTER, DeepHomo2.0 (all outperformed by PPLM-Contact)

**PPLM-Contact2 (enhanced with predicted complex structures):**

- Homodimers: top-L precision improved from 65.0% (PPLM-Contact) to 85.1% (PPLM-Contact2); Heterodimers: 44.3% → 88.0%
- vs. AlphaFold2.3/AlphaFold3/DMFold: improvements of 4.8%/8.7%/5.6% (homodimers) and 7.6%/5.6%/8.1% (heterodimers)
- Head-to-head win rates vs AlphaFold2.3: 60.3%/12.8% (homodimers, higher/lower), 63.0%/8.4% (heterodimers)
- F1-score: PPLM-Contact2 0.780 vs. AlphaFold2.3 0.748, AlphaFold3 0.739, DMFold 0.743 (improvements 4.3%, 5.5%, 5.0%)

**Interface residue identification precision:**

- PPLM-Contact (general methods): 0.824 (homodimers), 0.695 (heterodimers) — highest among general contact-prediction methods; improvements of 5.0–28.7% (homodimers) and 7.1–39.7% (heterodimers) over DeepInter/PLMGraph-Inter/CDPred/DeepHomo2.0/GLINTER
- PPLM-Contact2: average/median precision 0.897/0.962 (homodimers), 0.904/0.945 (heterodimers) — highest among all methods including structure-based approaches
- Case study (DNMT3A–DNMT3L complex, PDB 4U7P): PPLM-Contact identified 40/43 (93.0%) interface residues on chain A, 40/44 (90.9%) on chain B; overall precision 0.92 vs. DeepInter 0.35, PLMGraph-Inter 0.75, CDPred 0.69, GLINTER 0.09

### 4.2 Performance on Cross-Species Test Sets (if applicable)

- **Species: M. musculus** — AUPRC 0.920; F1-score improvement over TUnA +4.8%
- **Species: D. melanogaster** — AUPRC 0.906; F1-score improvement over TUnA +8.6%
- **Species: C. elegans** — AUPRC 0.883; F1-score improvement over TUnA +4.8%
- **Species: S. cerevisiae** — AUPRC 0.745; F1-score improvement over TUnA +16.9%
- **Species: E. coli** — AUPRC 0.784; F1-score improvement over TUnA +15.5%
- Other metrics: precision, recall, accuracy, F1-score, AUROC for all species reported in Supplementary Table S1 (not extracted in full here); overall mean±SD AUPRC across all five species: 0.848±0.078
- Sequence-similarity-stratified evaluation (grouped by max single-sequence identity to any human training-set protein): PPLM-PPI maintained strong performance across all homology ranges and delivered the largest gains over existing methods at every identity interval, indicating generalization is not driven by shortcut/memorization effects

### 4.3 Ablation Studies (if performed)

**PPLM-PPI ablations (pooling modules and feature removal):**

- **What was ablated:** Pooling strategy (max-only, mean-only, combined mean-max) and selective removal of PPLM-derived feature types (inter-protein attention, intra-protein attention, embedding)
- **Results of ablation:** Among individual pooling strategies, mean pooling achieved the highest accuracy alone; the combined mean-max pooling strategy produced the most robust and consistent performance across species. At the feature level, incorporating all three PPLM-derived feature types (inter-protein attention, intra-protein attention, embedding) yielded the best results, with embedding features contributing the strongest individual predictive signal (multi-seed results with 95% CI in Supplementary Tables S4–5 and Fig. S4)

**PPLM-Contact ablations (component/feature removal):**

- **What was ablated:** Cross-attention module, self-attention module, triangle-multiplication module, PPLM-derived inter-protein attention features, MSA features, monomer distance maps — each removed individually
- **Results of ablation:** Removing any major network component or feature reduces precision; the largest performance drops occur when excluding the triangle-multiplication module, PPLM inter-protein attention, MSA features, or monomer distance maps (Fig. 4D,E; detailed statistics in Supplementary Method S5)
- **Case study (2PMY, EF-hand domain, human RASEF):** PPLM-Contact without PPLM features (MSA + monomer structure only): 16.7% top-N precision (20.8% top-L); PPLM-Contact without MSA (PPLM features + monomer structure only): 78.2% top-N precision (201/257 ground-truth contacts), 100% top-L precision; full PPLM-Contact (both PPLM and MSA features): 81.3% top-N precision while maintaining 100% top-L precision — demonstrating PPLM-derived features alone can substantially outperform MSA-only features, and combining both gives the best result

### 4.4 Comparison with Existing Methods

**PPLM-PPI comparisons:** TUnA, ESMDNN-PPI, D-SCRIPT, Topsy-Turvy (all reimplemented/run locally except ESMDNN-PPI, which was reimplemented from the publication due to lack of public source code). PPLM-PPI achieved highest AUPRC and F1-score on all five test species, with markedly smoother precision-recall curves indicating superior stability.

**PPLM-Affinity comparisons:** ESM2-Affinity (same architecture/training procedure as PPLM-Affinity but initialized from vanilla ESM2, for controlled comparison), and the structure-based PPB-Affinity model (retrained from released source code under identical cross-validation). PPLM-Affinity outperformed both on PCC, SRCC, MAE, and RMSE, on the full dataset and especially on antibody–antigen and TCR–pMHC subgroups.

**PPLM-Contact / PPLM-Contact2 comparisons:** PLMGraph-Inter, DeepInter, CDPred, GLINTER, DeepHomo2.0 (general sequence/MSA-based contact predictors, all installed/run locally); AlphaFold2.3, AlphaFold3, DMFold (complex structure-prediction methods, compared for PPLM-Contact2). PPLM-Contact outperformed all general contact predictors on top-L precision and interface residue identification across Homodimer300, Heterodimer99, CASP_Homodimer43, and CASP_Heterodimer20. PPLM-Contact2 (incorporating predicted complex structures) outperformed AlphaFold2.3/AlphaFold3/DMFold on both top-L contact precision and F1-score, and on interface residue identification precision.

---

## 5. Evaluation Metrics Used

- **List all metrics reported in the paper:**
    - Perplexity (language modeling quality)
    - Precision, Recall, Accuracy, F1-score, AUROC, AUPRC (PPLM-PPI binary classification)
    - Pearson correlation coefficient (PCC), Spearman rank correlation coefficient (SRCC), mean absolute error (MAE), RMSE (PPLM-Affinity regression)
    - Top-N/Top-L/Top-1/Top-10/Top-50/L-10/L-5 contact precision (PPLM-Contact); interface residue identification precision
    - Effect sizes: Cohen's d; statistical significance: Wilcoxon signed-rank test p-values, Benjamini–Hochberg FDR correction
- **Primary metric (the one emphasized by authors):** AUPRC for PPLM-PPI (explicitly stated as "especially informative for imbalanced PPI datasets," used as the primary metric); top-L contact precision for PPLM-Contact/Contact2; PCC/SRCC jointly emphasized for PPLM-Affinity
- **Threshold used for binary classification:** Not explicitly stated as a fixed value (e.g., no stated 0.5 cutoff for PPLM-PPI's sigmoid output); for contact prediction, a residue-pair distance cutoff of <8Å (heavy atoms) defines a "contact" ground-truth label, and top-N/top-L is used as a ranking-based evaluation rather than a fixed probability threshold

---

## 6. Key Findings & Claims

- **Main conclusion about model performance:** PPLM, a paired-sequence protein language model with a hybrid intra-/inter-protein attention mechanism, substantially outperforms single-chain PLMs (ESM2) in modeling protein pairs (lower perplexity, especially at interfaces) and, when combined with task-specific downstream architectures, achieves state-of-the-art performance across three distinct PPI-related tasks: binary interaction prediction (PPLM-PPI), binding affinity estimation (PPLM-Affinity), and inter-protein contact/interface prediction (PPLM-Contact/Contact2) — the latter even outperforming leading structure-prediction methods (AlphaFold2.3, AlphaFold3, DMFold) when combined with their predicted structures.
- **Generalization claims:** PPLM-PPI, trained/validated only on H. sapiens data, generalizes robustly to five other species (mouse, fly, worm, yeast, E. coli) without species-specific retraining, and a sequence-similarity-stratified analysis confirms this generalization is not driven by "shortcut learning" or memorization of near-duplicate training sequences — performance gains over baselines hold across all homology intervals, described as "strict cross-species inductive setting" generalization. The bias audit (Fig. 2H) similarly shows PPLM's advantage over ESM2 holds robustly across identity bins to the training set.
- **Limitations mentioned by authors:** (1) PPLM's training data (composite PDB + STRING sequence-pair dataset), while diverse, may not fully capture the breadth of protein interactions across organisms and cellular states — transient, weak, or condition-dependent interactions remain underrepresented; (2) a representative failure case (Supplementary Fig. S13) involves a nanobody-antigen complex with a very small, loop-mediated interface and weak co-evolutionary signal, where PPLM-derived features only partially compensate for lack of MSA/structural information and remain insufficient to recover these small, flexible, weakly co-evolving interfaces — highlighting an intrinsic limitation of sequence-based contact prediction for such cases.

---

## 7. Implementation Details (if available)

- **Code availability:** Yes — webserver and source code freely available at https://zhanggroup.org/PPLM/; source code also on GitHub at https://github.com/junliu621/PPLM (MIT License); publication release deposited on Zenodo at https://zenodo.org/records/18256392. Datasets available at https://github.com/junliu621/PPLM/tree/main/data/. Training data sourced from PDB (https://www.rcsb.org/) and STRING (https://string-db.org/); PPLM-PPI interaction dataset from D-SCRIPT GitHub; PPLM-Affinity dataset from PPB-Affinity GitHub; Uniref30_2021_03 database used for MSA search.
- **Hardware used for training:**
    - PPLM pretraining: 4× NVIDIA A100 GPUs
    - PPLM-PPI: single NVIDIA A100 GPU
    - PPLM-Affinity: single NVIDIA A100 GPU (per fold)
    - PPLM-Contact: single NVIDIA A100 GPU
- **Training time:** Not given in wall-clock time; specified instead in training steps/epochs: PPLM pretraining = 50,000 steps (gradient accumulation step of 32); PPLM-PPI = 12 epochs per fold (10-fold CV); PPLM-Affinity = 15 epochs per fold (5-fold CV); PPLM-Contact = 100 epochs total