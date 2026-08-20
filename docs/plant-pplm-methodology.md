# Proposed Methodology: A Plant-Adapted Paired Protein Language Model for PPI Prediction

Base architecture: PPLM (Liu, Chen & Zhang, *Nat. Commun.* 2026 / bioRxiv 2025.07.07.663595) — a paired protein language model initialized from ESM2-650M and further pretrained on 3.3M PDB+STRING sequence pairs with a hybrid intra-/inter-protein cross-attention mechanism (RoPE for intra-chain, non-positional for inter-chain, learnable inter-protein attention weights + explicit inter-protein mask). Downstream PPI head (PPLM-PPI): embeddings + intra-attention + inter-attention → mean+max pooling → concatenation → 5-layer MLP → sigmoid.

Everything below is organized around one governing idea: **PPLM already solves the architecture problem the plant-specific literature never solved (true paired cross-attention). Your contribution is closing the domain gap between PPLM's training distribution (PDB structures + STRING's highest-confidence pairs, overwhelmingly human/model-organism/microbial) and the actual plant PPI problem (Arabidopsis-centric, small, HIPPIE/MI-score-filtered, heavily imbalanced, with a known and severe cross-species generalization failure mode).**

---

## 0. What PPLM changes about your starting position

Before this, my recommendation (per the prior lit-review synthesis) was: build per-residue cross-attention yourself, since nothing in the plant literature had it. That's now moot — PPLM has it, at a scale (33 transformer blocks, ESM2-650M init, 3.3M training pairs) no plant-specific group could plausibly match from scratch. Three things follow:

1. **Do not pretrain a paired PLM from scratch.** Start from the public PPLM checkpoint (`github.com/junliu621/PPLM`) and adapt it. Building competing from-scratch pretraining infrastructure would be reinventing a wheel that a much better-resourced team already built and open-sourced.
2. **Your novel contribution shifts from architecture to domain adaptation + plant-specific auxiliary signal fusion.** This is a legitimate, well-motivated, and currently unfilled research gap: PPLM's paper never touches plants, and every plant-specific paper you've reviewed (ESMAraPPI, DeepAraPPI, ARACoFusion, DWPPI, AraPathogen2.0, the DHT paper) uses either frozen single-sequence embeddings or an from-scratch/lightweight architecture — none use a pretrained *paired* PLM at all.
3. **PPLM-PPI, evaluated as-is on Arabidopsis data with zero adaptation, becomes your most important baseline** — it tells you how much of your eventual gain is "generic paired-attention transfers to plants for free" versus "plant-specific adaptation was necessary."

---

## 1. Data sources

### 1.1 Core Arabidopsis dataset (primary benchmark)

Use **DeepAraPPI's construction** as your base (most rigorous of the plant-specific sources you've reviewed):

- Positives: BioGRID + DIP + IntAct + MINT + TAIR, unified via UniProt ID, HIPPIE score ≥ 0.72 → 11,858 high-quality pairs (from Zheng et al. 2023; DHT paper and DWPPI use a similar but less-curated pool of 28,110 pairs from IntAct+TAIR+BioGRID without HIPPIE scoring — use theirs as a **sensitivity-check / larger-but-noisier alternative**, not primary).
- Negatives: random pairing + explicit subcellular co-localization exclusion, 1:10 ratio (consistent with ESMAraPPI/ARACoFusion convention, which is now the field standard — note DWPPI and the DHT paper use 1:1, which is a weaker, easier, and non-standard setup; don't adopt it for your primary results, but it's worth reporting as a secondary robustness check since real interactomes are far more imbalanced than 1:1).
- Split: Park & Marcotte C1 (train) / C2 (one protein unseen) / C3 (both unseen).

### 1.2 Cross-species datasets (generalization testing)

Use **all three** available plant species rather than just rice, since you now have access to two additional independent sources:

| Species | Source | Positive pairs | Notes |
|---|---|---|---|
| Rice (*O. sativa*) | DeepAraPPI's construction (DIP+MINT+BioGRID+IntAct+literature) | 611 | Small but rigorously filtered (self/non-physical/redundant removed); use as primary rice benchmark for direct comparability to ESMAraPPI/ARACoFusion numbers |
| Rice (*O. sativa*), alternative | PRIN database, via DWPPI or the DHT paper | 51,514 (DWPPI) or 4,800 (DHT paper) | Much larger n (DWPPI's) gives statistical power the 611-pair set lacks; use as a secondary, higher-powered rice test |
| Maize (*Z. mays*) | PPIM database, via DWPPI | 81,989 | Largest single non-Arabidopsis plant PPI resource across all papers reviewed; strongly recommend including — none of ESMAraPPI/DeepAraPPI/ARACoFusion test on maize at all, so this is genuinely novel coverage for the paired-PLM approach |

Negative sampling for cross-species sets: match DWPPI's convention (random pairing across different subcellular localizations) for consistency with the source data, but re-derive at a 1:10 ratio rather than DWPPI's native 1:1 to keep the imbalance regime consistent with your Arabidopsis benchmark.

### 1.3 Optional extension: plant–pathogen interactions

AraPathogen2.0's dataset (1,387 Arabidopsis–pathogen-effector PPIs from PPIN-1/PPIN-2/EffectorK, 8,505 Arabidopsis proteins, 872 effectors) is a **different biological question** (interspecies host–pathogen, not intraspecies plant–plant), but it's worth including as a stretch-goal secondary task for two reasons: (1) it comes with the most rigorous partitioning scheme of any paper you've reviewed (see §5.2), and (2) a paired-PLM approach is unusually well-suited to interspecies pairs, since PPLM was never trained assuming both proteins come from the same organism. If time permits, a small "does plant-adapted PPLM also improve host–pathogen PPI prediction" side-experiment would meaningfully broaden your paper's scope beyond a single narrow benchmark re-run.

### 1.4 External validation network

Follow ARACoFusion's precedent: pull a high-confidence (≥0.900), physical-interaction-only STRING subnetwork for Arabidopsis as an independent, non-training-derived validation set for qualitative network-level case-study analysis (Cytoscape visualization of correctly/incorrectly predicted edges).

---

## 2. Data representation & preprocessing

### 2.1 Base encoder: adapt PPLM, don't replace it

- Load the public PPLM checkpoint (33 transformer blocks, ESM2-650M initialization).
- **Domain-adaptive continued pretraining (recommended core step):** Continue PPLM's own masked-language-modeling objective, using its own masking strategy (30% of interface/known-interacting residues, 15% of others — approximate "interface" using your positive pair labels as a proxy since you won't have per-residue interface annotations for most plant pairs) on your combined plant paired-sequence corpus (Arabidopsis + rice + maize positive pairs, ~90K+ pairs across species once pooled). This directly addresses the domain-shift concern: PPLM's original training corpus (PDB structures + STRING top-0.999-confidence pairs) is overwhelmingly non-plant, so its residue-level statistics likely under-represent plant-specific motifs, domain family expansions (NLR/RLK-type receptors), and plant-typical co-evolutionary signal.
- Budget-conscious alternative if full continued pretraining is infeasible: fine-tune only the **final 2-4 transformer blocks** (mirroring PPLM-Affinity's own strategy of fine-tuning just the final block for its lower-data-volume affinity task) rather than all 33, which is both cheaper and less prone to catastrophic forgetting of PPLM's general paired-representation ability.
- **Ablation-worthy decision point:** compare (a) fully frozen PPLM, (b) final-block-only fine-tuned, (c) continued-pretrained-then-frozen, (d) continued-pretrained-then-fine-tuned. This maps directly onto the frozen-vs-fine-tuned question that every prior plant PPI paper (ESMAraPPI, DeepAraPPI, ARACoFusion) left completely unexplored — all of them only ever used frozen embeddings.

### 2.2 Feature extraction (mirror PPLM-PPI's own pipeline)

For each protein pair, extract exactly what PPLM-PPI extracts:
- Per-protein embeddings (from the final transformer layer)
- Intra-protein attention matrices (both proteins)
- Inter-protein attention matrix

Apply PPLM's own validated pooling strategy: **mean + max pooling combined** (their own ablation study found this outperforms either alone, and mean alone beats max/min alone) — don't deviate from this without a specific reason, since it's already been empirically validated at far larger scale than you can replicate.

### 2.3 Auxiliary plant-specific branches (this is where your genuine novelty lives)

Add branches that address the single most consistent finding across the plant-specific literature: **pure sequence signal generalizes worse than functional/domain signal**, especially cross-species.

**(a) GO-similarity branch.** DeepAraPPI's GO2vec was the best-generalizing individual predictor of any single signal across every difficulty level (Task1/2/3), and no PLM-based model (ESMAraPPI, ARACoFusion) has ever combined GO information with a pretrained embedding. Options in increasing complexity:
- Simplest: direct GO-term semantic similarity (GOSemSim-style, as DeepAraPPI itself used for validation) between the two proteins' annotations, as a small scalar/vector feature.
- Matching DeepAraPPI's original approach: node2vec embedding over the GO annotation graph, concatenated pairwise.
- This branch is inherently species-agnostic (GO terms are shared across organisms), which directly targets the cross-species generalization weakness every paper you've reviewed has documented.

**(b) Domain-interaction branch.** DeepAraPPI's Domain2vec (DDI network + node2vec) was the most *robust* individual predictor across difficulty levels even where it wasn't the best absolute performer — worth including as a second auxiliary signal, built from a plant-relevant domain-domain interaction network (3did, as DeepAraPPI used) plus HMMER-derived domain annotation of your protein set.

**(c) Network-embedding branch (optional, flagged for careful use).** DWPPI showed that combining sequence ("attribute") with a DeepWalk network embedding ("behavior") gave a 7-17% accuracy improvement over either alone, the largest single-paper improvement from feature fusion in the entire literature you've reviewed. **Important caveat you must handle explicitly**: a network embedding is inherently *transductive* — it only exists for proteins already present in the training PPI network, so it cannot help (and shouldn't be used) for genuinely novel proteins in your C3 test set or in true prospective discovery. Implement this as a strictly optional, ablatable branch that is zeroed out or omitted for C3-style evaluation, and be explicit in your writeup that any gain from this branch is a transductive-setting result, not a generalization result — DWPPI itself doesn't flag this distinction, and you doing so would be a genuine methodological improvement over that paper.

**Fusion point:** concatenate (a)/(b)/(optionally c) with PPLM's pooled embedding+attention representation, following ARACoFusion's late-concatenation pattern, before the final MLP.

---

## 3. Model architecture modifications

### 3.1 Downstream head

Start from PPLM-PPI's own validated design (5-layer MLP, LayerNorm + ReLU on the first four layers, sigmoid output) as your baseline head — don't redesign this without cause, it's already been ablated by the PPLM authors themselves. Extend the input dimensionality to accommodate the auxiliary branches from §2.3.

### 3.2 Class imbalance handling

PPLM-PPI itself uses plain binary cross-entropy despite training on a 10:1 imbalanced set (per D-SCRIPT's convention). Your plant data will likely be at least as imbalanced. Adopt **ARACoFusion's proven improvement** here: focal loss + label smoothing in place of plain BCE — this is the one architectural choice in ARACoFusion that showed a clear, real (not leakage-inflated) recall/F1 improvement over the plain-BCE ESMAraPPI baseline, and PPLM-Contact independently validates focal loss's value for a different severely-imbalanced task (contact prediction) in the same paper family. This is a well-supported, low-risk choice.

### 3.3 Calibration

Neither PPLM nor any prior plant-specific paper except ARACoFusion performs calibration. Add temperature scaling (reliability diagram + Expected Calibration Error) as a post-hoc step — cheap, useful for any downstream probabilistic-network use case, and gives you a second axis of comparison against ARACoFusion specifically.

### 3.4 Ensembling / uncertainty

Two established options from the literature, worth comparing directly rather than picking one blind:
- **PPLM-PPI's own approach**: 10-fold CV, keep the 5 best checkpoints by validation AUPRC, ensemble at inference. Already validated at scale by the PPLM authors.
- **ARACoFusion's approach**: MC-dropout-style variance regularization during training (3 forward passes, variance penalty added to loss).
These are not mutually exclusive — you could run both and report which contributes more, which would itself be a useful ablation nobody in either lineage has done (PPLM never tries variance regularization; ARACoFusion never tries checkpoint ensembling).

---

## 4. Training protocol

- **Stage 1 (optional, resource-permitting):** continued MLM pretraining of PPLM on pooled plant paired sequences (§2.1), using PPLM's own masking ratios and loss formulation, for domain adaptation.
- **Stage 2:** train the PPI-specific head (§3.1-3.2) on Arabidopsis C1, using stratified k-fold (following PPLM-PPI's 10-fold protocol) with focal loss, AdamW, and the auxiliary-branch fusion from §2.3.
- **Stage 3:** evaluate on C2/C3, then on the rice/maize cross-species sets (§1.2) without any additional fine-tuning (true zero-shot cross-species test — matching the protocol every plant paper you've read has used, so directly comparable), and optionally with a small amount of target-species fine-tuning (matching DeepAraPPI's "hybrid Arabidopsis+rice training set" experiment) to report both zero-shot and few-shot cross-species numbers.
- **Negative sampling sensitivity check:** as flagged in the earlier lit-review synthesis, every paper in this space (yours included, absent action) inherits an unvalidated negative-sampling heuristic. Retrain once with an alternative negative-sampling scheme (e.g., swap subcellular-exclusion for pure sequence-identity exclusion, or vice versa) and report how much your headline AUPRC shifts — nobody in this specific lineage has done this, and it's a cheap, genuine rigor addition.

---

## 5. Evaluation methodology

### 5.1 Primary benchmark (Arabidopsis)

- Metrics: AUPRC (primary), AUROC, MCC, Balanced Accuracy, F1, NPV — matching ARACoFusion's expanded suite, which is the most complete reporting standard in the literature you've reviewed.
- Report C2 and C3 **separately**, never merged/re-split (avoid ARACoFusion's leakage trap from the merged-5-fold-CV result).
- Add **bootstrap confidence intervals** on every headline metric, and **Wilcoxon signed-rank tests** for pairwise model comparisons — borrowed directly from PPLM's own paper, where this is used extensively and rigorously, and is a gap in every plant-specific paper you've reviewed to date.

### 5.2 Cross-species evaluation — upgrade the partitioning scheme

Rather than the simple "train on Arabidopsis, test on rice" binary used by DeepAraPPI/ARACoFusion, adopt **AraPathogen2.0's more granular partitioning logic**, adapted from host-pathogen to species-species framing:
- **Regular test**: held-out Arabidopsis pairs (proteins may have been seen).
- **Novel species**: pairs where both proteins come from a species never seen in training (e.g., rice, maize).
- This is a direct generalization of AraPathogen2.0's "novel host / novel pathogen / novel host & pathogen" four-way split, and gives a more granular, more informative picture of exactly where generalization breaks down than the single-number cross-species AUPRC every other plant paper reports.

### 5.3 Baselines — retrain, don't cite

Following the reproducibility discipline that differentiates ESMAraPPI/DeepAraPPI from ARACoFusion's partial shortcuts (as flagged earlier):
- **Retrain on your exact training split:** ESMAraPPI (Hadamard+MLP on frozen ESM-1b), DeepAraPPI (or at minimum its GO2vec branch, the strongest individual predictor), D-SCRIPT, RAPPPID, PIPR, TAGPPI.
- **Critically, retrain unmodified PPLM-PPI** on the same plant data with zero domain adaptation — this is your most important baseline, as noted in §0, since it isolates exactly how much your plant-specific adaptation contributes versus generic paired-attention transfer.
- If feasible, ARACoFusion (given its GitHub availability) and DWPPI (given code availability per its paper).

### 5.4 Ablations

Systematically ablate, matching the granularity both ARACoFusion and PPLM's own papers use:
- Frozen vs. continued-pretrained vs. fine-tuned PPLM base (§2.1)
- With/without GO branch, domain branch, network branch (§2.3), individually and combined
- With/without focal loss + label smoothing (§3.2)
- With/without temperature scaling (§3.3) — note per the earlier ARACoFusion review that this should show negligible AUPRC/AUROC change (it's a monotonic rescaling), so use it to sanity-check your own ablation pipeline rather than expecting a discriminative-performance gain
- Checkpoint ensembling vs. variance regularization vs. both (§3.4)

### 5.5 Case studies / external validation

- STRING-derived Arabidopsis subnetwork (ARACoFusion's approach), visualized in Cytoscape.
- Literature-verification of top-K highest-confidence novel predictions (DWPPI's and DeepAraPPI's approach) — DWPPI's maize case study (14/20 top predictions literature-confirmed) is a good template for what a credible version of this looks like.

---

## 6. Phased build plan

1. **Baseline reproduction** — retrain ESMAraPPI and unmodified PPLM-PPI on the DeepAraPPI Arabidopsis C1/C2/C3 split. This alone tells you whether PPLM already beats the plant-specific SOTA with zero adaptation, which calibrates how much headroom your subsequent work realistically has.
2. **Auxiliary branch fusion** — add GO + domain branches to the frozen-PPLM pipeline, evaluate incremental AUPRC/MCC gain on C2/C3, individually and combined.
3. **Domain-adaptive pretraining** — if step 1 shows a meaningful domain gap (frozen general-purpose PPLM underperforms plant-specific frozen-ESM-1b baselines on C3 specifically, where domain-general representations should matter most), invest in continued pretraining; if not, this step may not be worth the compute cost — let step 1's results decide.
4. **Imbalance/calibration/ensembling additions** — layer in focal loss, temperature scaling, ensembling.
5. **Cross-species extension** — rice + maize zero-shot and few-shot evaluation, using the AraPathogen2.0-style granular partitioning.
6. **Full ablation suite + statistical testing** — bootstrap CIs, Wilcoxon tests, negative-sampling sensitivity check.
7. **External validation** — STRING subnetwork case study, literature-verified top-K predictions.
8. **(Stretch) Plant-pathogen extension** — apply the finished pipeline to AraPathogen2.0's host-pathogen data as a secondary demonstration of generality.

Step 1 is the highest-priority, lowest-cost, highest-information experiment — run it first, since its outcome should reshape how much effort you invest in the domain-adaptation-heavy steps versus the fusion/evaluation-heavy steps.
