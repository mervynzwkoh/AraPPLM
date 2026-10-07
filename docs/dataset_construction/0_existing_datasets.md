Yes. For Objective 1, these three papers are especially useful because they represent two distinct dataset lineages:

* **DeepAraPPI (Zheng et al., 2023)** constructed its own multi-database Arabidopsis benchmark and introduced the rice cross-species benchmark.
* **ESMAraPPI (Zhou et al., 2023)** constructed a newer, more stringent IntAct-based Arabidopsis benchmark with explicit sequence-redundancy filtering.
* **ARACoFusion (Sarkar & Sarkar, 2026)** did **not construct a fundamentally new Arabidopsis benchmark**. It reused the ESMAraPPI dataset and the DeepAraPPI rice dataset, which makes it useful for seeing how later models benchmark on established datasets. ARACoFusion is currently a **bioRxiv preprint, not peer-reviewed**. 

---

# 1. DeepAraPPI — Zheng et al., 2023

DeepAraPPI is the broadest of the three in terms of **data provenance**. Rather than relying on one interaction database, Zheng et al. aggregated experimentally reported Arabidopsis PPIs from **BioGRID, DIP, IntAct, MINT, and TAIR**. 

## 1.1 Positive dataset construction

Their starting sources were:

| Source      | Role                              |
| ----------- | --------------------------------- |
| **BioGRID** | Experimentally reported PPIs      |
| **DIP**     | Experimentally reported PPIs      |
| **IntAct**  | Experimentally reported PPIs      |
| **MINT**    | Experimentally reported PPIs      |
| **TAIR**    | Arabidopsis-specific interactions |

Because these resources use different identifier systems, the authors first converted protein identifiers to **UniProt IDs**.

They then removed:

* self-interactions;
* redundant interactions;
* non-physical interactions;
* interactions containing a protein shorter than **40 aa**;
* interactions containing proteins with non-standard amino acids.

After this initial filtering:

**49,398 experimentally verified PPIs involving 10,330 Arabidopsis proteins** remained. 

### Quality filtering using HIPPIE

They did not treat all 49,398 interactions as equally reliable.

The authors applied the **HIPPIE scoring scheme**, which assigns interaction confidence based on factors including:

1. experimental method used;
2. number of publications supporting the interaction;
3. species information associated with the evidence.

They selected a threshold of:

$$
\text{HIPPIE score} \geq 0.72
$$

giving:

* **11,858 high-quality PPIs**
* **37,540 lower-quality PPIs**

The 11,858 high-quality interactions formed the principal positive benchmark dataset. 

This is quite different from ESMAraPPI, which uses an IntAct MIscore threshold rather than integrating five databases and rescoring interactions.

---

# 2. DeepAraPPI negative dataset

The authors downloaded the Arabidopsis **reference proteome from UniProt**.

They removed:

* proteins <40 aa;
* proteins containing non-standard amino acids.

This left:

$$
28,361\text{ proteins}
$$

as the sampling pool.

Negatives were then **randomly generated**, but with an important biological constraint:

> The two sampled proteins must **not share the same subcellular localization**, and the pair must not already belong to the known PPI set.

This was intended to reduce the probability that a randomly generated "negative" is actually an undiscovered interaction. 

They used:

$$
\text{positive : negative} = 1:10.
$$

### Important weakness for your Objective 1

DeepAraPPI does **not report an explicit sequence-identity threshold** for reducing homology/redundancy in this benchmark. 

That matters substantially.

Two proteins that are highly homologous can occur across train and test sets. A sequence model can consequently exploit homology rather than learning interaction determinants.

This is one of the major methodological improvements introduced by ESMAraPPI.

---

# 3. DeepAraPPI train/test design

This is one of the most important parts of the paper.

A naive random pair split is problematic for PPI prediction because the same proteins can appear repeatedly across training and testing.

DeepAraPPI therefore followed the **Park & Marcotte pair-input evaluation framework** and created three difficulty levels. 

### Task 1 — random split

All 11,858 positive PPIs were used.

They generated:

$$
118,580
$$

negative pairs.

Then:

$$
80\% \rightarrow training
$$

$$
20\% \rightarrow independent\ test
$$

This is the easiest setting because proteins appearing in the test set can also occur in training.

---

## Task 2 — one unseen protein

The positive dataset was divided into:

| Partition | Positive PPIs |
| --------- | ------------: |
| C1        |         2,844 |
| C2        |         6,005 |
| C3        |         3,009 |

For Task 2:

**C1 = training**

**C2 = testing**

For every C2 interaction:

$$
(A,B)
$$

only **one** protein can occur in C1.

Conceptually:

```text
Training C1:
A -- B
A -- C
D -- E

C2:
A -- X
B -- Y
E -- Z
```

X, Y and Z are unseen, while their partners are known.

This tests whether the model can generalize an interaction rule learned for a known protein to a new partner.

---

# 4. DeepAraPPI Task 3 — both proteins unseen

Again:

**C1 = training**

but:

**C3 = testing**

Neither protein in a C3 pair occurs in C1.

For example:

```text
C1 training proteins:
A B C D E

C3 test:
X -- Y
P -- Q
R -- S
```

with:

$$
\{X,Y,P,Q,R,S\}\cap \{A,B,C,D,E\}=\varnothing
$$

This is substantially more stringent and is the most relevant of the three settings if your eventual claim is that a model predicts PPIs for **previously unseen plant proteins**.

The complete DeepAraPPI benchmark therefore contains:

| Dataset | Positives | Negatives |  Total |
| ------- | --------: | --------: | -----: |
| C1      |     2,844 |    28,440 | 31,284 |
| C2      |     6,005 |    60,050 | 66,055 |
| C3      |     3,009 |    30,090 | 33,099 |

Your existing project processing has reproduced these exact dataset sizes. 

---

# 5. DeepAraPPI evaluation

Because the benchmark is highly imbalanced at 1:10, the authors emphasized the **precision-recall curve and AUPRC** rather than relying on accuracy.

They also reported:

* TPR/recall;
* FPR;
* precision;
* PR curves.

AUPRC was the primary comparison metric. 

The final integrated model obtained:

| Model                     |    Task 1 |    Task 2 |    Task 3 |
| ------------------------- | --------: | --------: | --------: |
| RCNN                      |     0.925 |     0.746 |     0.481 |
| Domain2vec                |     0.868 |     0.780 |     0.681 |
| GO2vec                    |     0.939 |     0.871 |     0.803 |
| **Integrated DeepAraPPI** | **0.965** | **0.897** | **0.825** |

These results are informative because the sequence-only RCNN collapses from **0.925 → 0.481** between Task 1 and Task 3, showing how dramatically performance can change once protein overlap is eliminated. 

That observation should strongly influence how we design your benchmark.

---

# 6. DeepAraPPI rice dataset

This part is particularly relevant to your question.

The rice positives came from:

* **DIP**
* **MINT**
* **BioGRID**
* **IntAct**
* additional literature, including Wierbowski et al. (2020). 

They removed:

* self-interactions;
* non-physical interactions;
* redundant interactions.

This produced:

$$
611\ PPIs
$$

between:

$$
555\ rice\ proteins.
$$

Negatives were random rice protein pairs not present in the known rice PPI set.

Again:

$$
1:10
$$

giving approximately:

$$
611 + 6,110 = 6,721\ pairs.
$$

---

# 7. Did DeepAraPPI retrain on rice?

### Main cross-species experiment: **No.**

The authors explicitly used the models trained on **Arabidopsis** to predict the rice interactions.

So:

```text
Arabidopsis training
        ↓
DeepAraPPI
        ↓
Rice dataset
        ↓
Cross-species evaluation
```

There was **no rice-specific training for the primary cross-species benchmark**.

The AUPRC values were:

| Model         | Rice AUPRC |
| ------------- | ---------: |
| RCNN          |      0.248 |
| Domain2vec    |      0.279 |
| GO2vec        |      0.265 |
| Integrated LR |  **0.305** |
| RF-DPC        |      0.171 |

These values are reported in the project literature extraction. 

### But they then performed a second experiment

This is the nuance that is easy to miss.

They randomly split the rice dataset:

$$
80\% \ rice \rightarrow training
$$

$$
20\% \ rice \rightarrow independent\ test
$$

and combined the rice training portion with the Arabidopsis training data:

```text
Arabidopsis training
        +
80% Rice training
        ↓
Hybrid training dataset
        ↓
Retrain GO2vec
        ↓
20% held-out Rice
```

The retrained GO2vec obtained:

$$
AUPRC=0.561
$$

compared with:

$$
AUPRC=0.285
$$

for the Arabidopsis-only GO2vec evaluated on that same rice holdout. 

So there are actually **two DeepAraPPI rice experiments**:

| Experiment               | Rice used for training? | Meaning                       |
| ------------------------ | ----------------------- | ----------------------------- |
| Main rice benchmark      | **No**                  | Cross-species transfer        |
| Hybrid GO2vec experiment | **Yes, 80%**            | Rice-adapted supervised model |

For your Objective 2 later, these must not be conflated.

---

# 8. ESMAraPPI — Zhou et al., 2023

ESMAraPPI constructed a different benchmark.

Rather than aggregating five databases, it used **IntAct only**. 

Positive interactions had to satisfy:

1. Arabidopsis;
2. interaction type = **"direct interaction" OR "physical association"**;
3. IntAct:

$$
MIscore \geq 0.45
$$

The paper states that interactions with MIscore <0.45 were removed.

This produced:

$$
7,729\ positive\ PPIs.
$$



This is a smaller but more explicitly filtered positive dataset than DeepAraPPI.

---

# 9. ESMAraPPI negative generation

This is where ESMAraPPI becomes especially relevant to your benchmark design.

The procedure was approximately:

```text
Complete Arabidopsis proteins
           │
           ├── remove proteins participating in positive PPIs
           │
           ↓
Remaining proteins
           │
           ├── remove proteins ≥40% identical
           │   to positive proteins
           ↓
Candidate negatives
           │
           ├── remove internal redundancy
           │   at 40% identity
           ↓
8,382 proteins
           │
           ├── combine with positive proteins
           ↓
Random pair generation
           │
           ├── remove known PPIs
           ↓
77,290 negatives
```

Specifically, proteins sharing:

$$
\geq 40\%
$$

sequence identity with positive proteins were filtered.

They then performed another **40% identity redundancy reduction** within the remaining candidate proteins.

This left:

$$
8,382
$$

background proteins.

They mixed these proteins with proteins appearing in positive interactions and randomly generated pairs not experimentally identified as PPIs.

Final dataset:

$$
7,729\ positives
$$

$$
77,290\ negatives
$$

giving:

$$
1:10.
$$



### This is an important methodological difference

DeepAraPPI mainly attempts to make negatives biologically plausible using **subcellular localization**.

ESMAraPPI mainly attempts to control **sequence redundancy/homology**.

Those address different sources of bias:

$$
\boxed{\text{DeepAraPPI: biological negative plausibility}}
$$

versus

$$
\boxed{\text{ESMAraPPI: sequence-homology leakage control}}
$$

For your dataset, I would not view these as competing alternatives. A stronger benchmark could incorporate **both principles**.

---

# 10. ESMAraPPI C1/C2/C3

Again they followed Park & Marcotte.

| Partition | Positive | Negative |      Total |
| --------- | -------: | -------: | ---------: |
| **C1**    |    3,519 |   35,190 | **38,709** |
| **C2**    |    3,404 |   34,040 | **37,444** |
| **C3**    |      806 |    8,060 |  **8,866** |



The interpretation is explicit:

### C1

Training dataset.

### C2

One member of every test pair may have appeared in C1:

$$
(\text{seen},\text{unseen})
$$

### C3

Neither member appears in C1:

$$
(\text{unseen},\text{unseen})
$$

This is why C3 is the strongest protein-level generalization test.

---

# 11. ESMAraPPI evaluation

The primary model was:

$$
ESM\text{-}1b
\rightarrow
1280D\ embedding/protein
$$

then pair combination via:

$$
z_{pair}=z_A\odot z_B
$$

i.e. the **Hadamard product**, followed by a four-layer MLP.

The MLP had:

$$
1024 \rightarrow 512 \rightarrow 128 \rightarrow 16
$$

hidden units, BCE loss and sigmoid output. 

They evaluated:

* accuracy;
* specificity;
* precision;
* recall;
* MCC;
* PR curve/AUPRC;
* ROC/AUROC.

Because of the 1:10 imbalance, AUPRC is particularly informative.

The main ESM-1b results were approximately:

| Dataset |     AUPRC |     AUROC |
| ------- | --------: | --------: |
| C2      | **0.834** | **0.966** |
| C3      | **0.810** | **0.960** |



---

# 12. An important strength of ESMAraPPI's benchmarking

They did not simply copy published performance numbers for the generic predictors.

For:

* PIPR;
* D-SCRIPT;
* RAPPPID;
* TAGPPI;

they downloaded the source code and **retrained the models on ESMAraPPI C1**, then evaluated them on exactly the same C2 and C3 datasets.

That makes these comparisons substantially cleaner than comparing values reported on different datasets.

Their comparison yielded, for example:

| Model         |   C2 AUPRC |  C3 AUPRC |
| ------------- | ---------: | --------: |
| D-SCRIPT      |      0.292 |     0.291 |
| RAPPPID       |      0.516 |     0.371 |
| PIPR          |      0.588 |     0.387 |
| TAGPPI        |      0.700 |     0.554 |
| **ESMAraPPI** | **~0.834** | **0.810** |



They also tested DeepAraPPI on the ESMAraPPI dataset rather than comparing against DeepAraPPI's original easier dataset. That is why the DeepAraPPI numbers you see in the ESMAraPPI paper differ from Zheng et al.'s original numbers.

This is an important principle for your future benchmark:

$$
\boxed{\text{same training data + same test data + same metric}}
$$

rather than comparing numbers across papers.

---

# 13. Did ESMAraPPI test on rice?

**No.**

This is worth making explicit because later papers make the lineage confusing.

The original Zhou et al. ESMAraPPI paper contains **no rice cross-species experiment**. 

So a number such as:

> ESMAraPPI rice AUPRC = 0.2938

is **not from Zhou et al. (2023)**.

It comes from the later ARACoFusion study.

---

# 14. ARACoFusion — Sarkar & Sarkar, 2026

ARACoFusion is fundamentally different in dataset provenance.

It effectively **inherits the ESMAraPPI Arabidopsis benchmark**.

The paper says the experimentally validated Arabidopsis dataset was downloaded from IntAct and **obtained from Zhou et al.**, giving the same:

$$
7,729\ positives
$$

and:

$$
77,290\ negatives.
$$

It applies the same reported 40% sequence-identity filtering and Park & Marcotte C1/C2/C3 partitioning. 

Thus:

|                 | ESMAraPPI |          ARACoFusion |
| --------------- | --------: | -------------------: |
| Positive source |    IntAct | **ESMAraPPI/IntAct** |
| Positives       |     7,729 |                7,729 |
| Negatives       |    77,290 |               77,290 |
| Identity cutoff |       40% |                  40% |
| C1              |    38,709 |               38,709 |
| C2              |    37,444 |               37,444 |
| C3              |     8,866 |                8,866 |

The dataset is essentially inherited, not newly curated.

---

# 15. Important ARACoFusion paper inconsistency

There is a wording error in the ARACoFusion Methods.

It states that interactions with:

> "MI score of < 0.45 was retained"



But ESMAraPPI—the dataset they explicitly say they obtained from Zhou et al.—does the opposite:

$$
MIscore <0.45 \rightarrow removed
$$

therefore:

$$
\boxed{MIscore \geq 0.45\ retained}
$$



Given that ARACoFusion reports exactly the same **7,729 positive interactions** as ESMAraPPI, the "<0.45 retained" statement is almost certainly a manuscript error.

For your thesis/methodology, I would **not reproduce the ARACoFusion wording as though it were a genuine filtering criterion**. Document the discrepancy.

---

# 16. ARACoFusion evaluation design

The primary evaluation again uses:

```text
C1
 ↓
training
 ↓
ARACoFusion
 ├────────→ C2
 └────────→ C3
```

ARACoFusion reports a much larger metric suite than the previous papers:

* accuracy;
* sensitivity;
* specificity;
* precision;
* AP;
* NPV;
* F1;
* MCC;
* balanced accuracy;
* AUROC;
* AUPRC.

The authors explicitly emphasize **F1, MCC and balanced accuracy** because of the 1:10 class imbalance. 

They report approximately:

|        |     C2 |     C3 |
| ------ | -----: | -----: |
| AUPRC  | 0.8546 | 0.8066 |
| AUROC  | 0.9548 | 0.9308 |
| MCC    | 0.7817 | 0.7326 |
| Recall | 0.7471 | 0.6563 |



---

# 17. ARACoFusion also performs 5-fold CV — but this needs caution

They merge:

$$
C1+C2+C3
$$

and perform **stratified five-fold cross-validation**.

The resulting AUPRC is:

$$
0.9967.
$$



I would treat this result cautiously.

Why?

The key advantage of C2/C3 is **protein-level separation**.

If you merge everything and then perform ordinary stratified pair-level CV, proteins can occur in both training and validation folds:

```text
Fold training:
Protein A -- Protein B

Fold validation:
Protein A -- Protein X
```

The validation model has therefore already seen Protein A.

Consequently:

$$
AUPRC_{5fold}=0.9967
$$

is **not comparable in stringency** to:

$$
AUPRC_{C3}=0.8066.
$$

For your Objective 1 dataset, I would strongly recommend **cluster-/protein-aware splitting rather than ordinary stratified pair CV**.

---

# 18. ARACoFusion rice dataset

ARACoFusion did not construct a new rice benchmark either.

It explicitly takes the rice dataset from **Zheng et al./DeepAraPPI**.

The underlying positives therefore come from:

* DIP;
* MINT;
* BioGRID;
* IntAct;

with:

* self-interactions removed;
* non-physical interactions removed;
* redundant interactions removed.

Final:

$$
611\ positive
$$

$$
6,110\ negative
$$

$$
6,721\ total.
$$

The ARACoFusion Methods have another typo saying "6111 negatives" and even call randomly selected non-PPIs "positive samples"; Table 1 gives the internally consistent value of **6,110 negatives**. 

Again, something to flag rather than inherit blindly.

---

# 19. Did ARACoFusion train on rice?

### No.

This is stated explicitly.

Both:

* **ARACoFusion**
* **their reimplementation of ESMAraPPI**

were:

> **trained exclusively on Arabidopsis C1 data**

and then evaluated on the independent rice dataset.

So the experiment is:

```text
                 Arabidopsis C1
                       │
                supervised training
                       │
              ┌────────┴────────┐
              ↓                 ↓
       ARACoFusion        ESMAraPPI
              │                 │
              └────────┬────────┘
                       ↓
                 Rice dataset
                NO rice training
```

This is a genuine **cross-species transfer experiment**.

Their reported rice results are:

| Model                      | Rice AUPRC |
| -------------------------- | ---------: |
| ARACoFusion                | **0.3519** |
| ESMAraPPI reimplementation | **0.2938** |

The key provenance point is that **0.2938 is Sarkar & Sarkar's evaluation of the ESMAraPPI architecture, not a result reported in Zhou et al.'s original ESMAraPPI paper**. Your existing project documentation records this distinction as well. 

---

# 20. Putting all three dataset pipelines side by side

| Feature                     | DeepAraPPI                                  | ESMAraPPI                       | ARACoFusion                   |
| --------------------------- | ------------------------------------------- | ------------------------------- | ----------------------------- |
| Year                        | 2023                                        | 2023                            | 2026                          |
| Arabidopsis source          | **BioGRID + DIP + IntAct + MINT + TAIR**    | **IntAct**                      | **Reuses ESMAraPPI**          |
| Initial positives           | 49,398                                      | 7,729                           | 7,729                         |
| Positive quality filtering  | HIPPIE ≥0.72 → 11,858                       | Direct/physical + MIscore ≥0.45 | Same dataset as ESMAraPPI     |
| Sequence identity filtering | **None reported**                           | **40%**                         | **40%, inherited**            |
| Negative strategy           | Random + different subcellular localization | Homology-filtered random pairs  | Same as ESMAraPPI             |
| Pos:Neg                     | 1:10                                        | 1:10                            | 1:10                          |
| Main split                  | Task1 + C1/C2/C3                            | C1/C2/C3                        | C1/C2/C3                      |
| C2                          | One protein overlaps train                  | One protein overlaps train      | Same                          |
| C3                          | Neither protein overlaps train              | Neither protein overlaps train  | Same                          |
| Main metric                 | AUPRC                                       | AUPRC + AUROC                   | AUPRC + AUROC + MCC/BACC etc. |
| Rice dataset                | **Constructed here**                        | None                            | **Reuses DeepAraPPI**         |
| Main rice model training    | Arabidopsis only                            | N/A                             | Arabidopsis only              |
| Rice adaptation experiment  | **Yes, separate GO2vec experiment**         | No                              | No                            |

---

# 21. The dataset genealogy is therefore surprisingly simple

Rather than three independent plant benchmarks, what the literature really gives us is approximately:

```text
                    ┌─────────────────────────────┐
                    │      DeepAraPPI 2023        │
                    │ BioGRID/DIP/IntAct/MINT/TAIR│
                    │       11,858 positives      │
                    └──────────────┬──────────────┘
                                   │
                       introduces Rice dataset
                                   │
                                   ↓
                         611 rice positives
                                   │
                                   │ reused by
                                   ↓
                         ARACoFusion 2026


IntAct
  │
  │ direct/physical
  │ MIscore ≥ 0.45
  ↓
7,729 positives
  │
  │ 40% identity filtering
  ↓
77,290 negatives
  │
  ↓
┌─────────────────────┐
│  ESMAraPPI 2023     │
│ C1 / C2 / C3        │
└──────────┬──────────┘
           │
           │ dataset reused
           ↓
┌─────────────────────┐
│ ARACoFusion 2026    │
└─────────────────────┘
```

So **ARACoFusion does not provide a third independent dataset**. It provides a third model evaluated primarily on datasets introduced by the other two papers.

---

# 22. What this tells us about your Objective 1

There are several ideas worth retaining, but also several weaknesses that I would avoid reproducing.

The strongest aspects are **DeepAraPPI's multi-database aggregation**, **ESMAraPPI's sequence-redundancy filtering**, and the **C1/C2/C3 protein-overlap framework** used across the papers. DeepAraPPI also introduced **localization-aware negative sampling** in an effort to avoid spurious random pairings; however, as analyzed in [docs/dataset_construction/2_evidence_gathering_report.md](2_evidence_gathering_report.md), compartment-based negatives introduce a well-documented **"organelle shortcut"** (where models learn individual sequence localization signals like transit peptides rather than true biophysical binding interfaces). The cross-species rice test is also valuable because it measures a qualitatively different problem: species transfer rather than merely unseen proteins within Arabidopsis.

However, neither benchmark completely solves the leakage/negative-label problem. ESMAraPPI's 40% filtering is applied primarily in construction of the negative protein pool; it is not equivalent to guaranteeing that **all proteins across train/test partitions are separated into sequence clusters at ≤40% identity**. Likewise, "not experimentally reported as interacting" is not equivalent to experimentally validated non-interaction. In plant interactomes, where coverage is sparse, random negatives inevitably contain unknown positives, while pure localization-filtered negatives introduce shortcut learning. This motivates degree-matched and hard-paralog negative suites.

For the benchmark we construct for your Objective 1, I therefore think the literature points toward a stronger hybrid design:

$$
\boxed{
\text{multi-source experimental positives}
+
\text{evidence-level provenance}
+
\text{biologically informed negatives}
+
\text{sequence-cluster-aware splitting}
+
\text{protein-disjoint test}
+
\text{cross-species holdout}
}
$$

In particular, I would distinguish **three separate things that these papers partly conflate**: (1) *positive evidence quality*, such as assay type/publication support; (2) *negative reliability*, such as localization/tissue-expression incompatibility; and (3) *train-test leakage*, which should be controlled by sequence clustering **before splitting**. That distinction will become central when we decide whether to build on DeepAraPPI/ESMAraPPI or construct a substantially more rigorous benchmark from IntAct, BioGRID and other current plant resources.
