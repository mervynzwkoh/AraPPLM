# Research Plan: Evidence Mapping for a Modern Plant PPI Benchmark

## 1. Research objective

Conduct a systematic, quantitative investigation of currently available **experimentally supported plant protein–protein interaction datasets**, with particular emphasis on *Arabidopsis thaliana* and *Oryza sativa*, to determine the optimal data sources, species composition, evidence criteria, negative-sampling strategy, and train/validation/test splitting strategy for constructing a new benchmark for protein-language-model-based plant PPI prediction.

The research must answer two overarching questions:

> **Q1. What is the highest-quality and most comprehensive experimentally supported plant PPI dataset that can currently be assembled?**

and

> **Q2. Should the training benchmark remain Arabidopsis-specific, as in DeepAraPPI/ESMAraPPI, or should experimentally supported PPIs from multiple plant species be incorporated?**

Do **not** construct the final dataset during this research stage.

The purpose is to gather the evidence required to make all major dataset-design decisions first.

---

# 2. Background the research agent should understand

Existing plant PPI prediction benchmarks have primarily focused on Arabidopsis.

### DeepAraPPI — Zheng et al., 2023

DeepAraPPI aggregated experimentally reported Arabidopsis PPIs from:

* BioGRID;
* DIP;
* IntAct;
* MINT;
* TAIR.

It obtained 49,398 experimentally reported interactions and subsequently retained 11,858 high-confidence PPIs using a HIPPIE score ≥0.72. Negatives were sampled at a 1:10 positive:negative ratio, with proteins from non-overlapping subcellular localizations used to reduce the likelihood of false negatives. It used the Park & Marcotte C1/C2/C3 framework for protein-overlap-controlled evaluation. 

DeepAraPPI also constructed a rice cross-species dataset containing 611 experimentally supported rice PPIs obtained from DIP, MINT, BioGRID, IntAct and literature. 

### ESMAraPPI — Zhou et al., 2023

ESMAraPPI used IntAct Arabidopsis PPIs classified as `"direct interaction"` or `"physical association"` and retained interactions with MIscore ≥0.45.

This produced:

$$
7,729\ positives.
$$

For negative construction, candidate proteins were filtered at 40% sequence identity relative to positive proteins and internally redundancy-reduced at 40%, producing 8,382 background proteins. Random negative pairs were generated at a 1:10 ratio. 

It used:

* C1 = training;
* C2 = one protein seen during training;
* C3 = both proteins unseen.



### ARACoFusion — Sarkar & Sarkar, 2026

ARACoFusion does not introduce an independent Arabidopsis benchmark. It substantially reuses the ESMAraPPI dataset and uses the DeepAraPPI rice dataset for cross-species evaluation. 

Therefore, there is currently no clearly established modern successor to these 2023 benchmark datasets despite continued growth of experimental PPI resources.

---

# 3. Scope of the research

The investigation should cover at minimum:

* *Arabidopsis thaliana*
* *Oryza sativa*
* *Zea mays*
* *Glycine max*
* *Solanum lycopersicum*

Also investigate other plant species if preliminary searches indicate meaningful experimentally validated PPI datasets.

Potential additional species include:

* *Triticum aestivum*
* *Hordeum vulgare*
* *Nicotiana benthamiana*
* *Medicago truncatula*
* *Brachypodium distachyon*
* *Sorghum bicolor*
* *Chlamydomonas reinhardtii*

Do not assume that database counts alone represent the available literature. Search both **structured databases and primary publications**.

---

# 4. Research Question A — Current plant PPI inventory

For each species, determine how many experimentally supported PPIs currently exist.

Investigate at minimum:

* BioGRID;
* IntAct;
* DIP, if still relevant/current;
* MINT or its successor/integrated resources;
* STRING experimental evidence;
* HitPredict;
* plant-specific databases;
* TAIR/BAR for Arabidopsis;
* POPPIN for rice;
* literature-derived interactomes;
* ProteomeXchange/PRIDE datasets where relevant.

For every source report:

| Field                 | Required information          |
| --------------------- | ----------------------------- |
| Database/resource     | Name                          |
| Version               | Exact release/version         |
| Retrieval date        | Date                          |
| Species               | Scientific name + taxonomy ID |
| Raw interactions      | Number                        |
| Unique protein pairs  | Number                        |
| Unique proteins       | Number                        |
| Physical PPIs         | Number                        |
| Direct PPIs           | Number if available           |
| Publications          | Number represented            |
| Experimental evidence | Assay breakdown               |
| Identifier system     | UniProt/TAIR/RAP/etc.         |
| Download availability | Yes/no                        |
| Programmatic access   | API/bulk download             |
| Evidence metadata     | What is retained?             |

Crucially, distinguish:

$$
\text{evidence records}
$$

from:

$$
\text{unique protein pairs}.
$$

Ten experiments supporting A–B should remain **one PPI with ten evidence records**, not ten PPIs.

---

# 5. Research Question B — What experimental assays constitute the available data?

For Arabidopsis and rice in particular, quantify PPIs according to experimental methodology.

At minimum investigate:

* yeast two-hybrid;
* split-ubiquitin;
* affinity purification–mass spectrometry;
* co-immunoprecipitation;
* pull-down;
* protein-fragment complementation;
* BiFC;
* split-luciferase;
* protein microarray;
* cross-linking mass spectrometry;
* proximity labeling such as TurboID/miniTurbo/PUP-IT;
* biochemical/biophysical direct-binding assays.

Determine whether databases distinguish:

$$
\text{direct interaction}
$$

from:

$$
\text{physical association}
$$

from:

$$
\text{proximity/co-complex association}.
$$

Do not treat these categories as biologically equivalent.

Produce an evidence hierarchy proposal, for example:

```text
Tier 1
Direct physical interaction
+ multiple independent experiments/publications

Tier 2
Direct physical interaction
+ one experimental source

Tier 3
Physical association/co-complex evidence

Tier 4
Proximity-based evidence

Tier 5
Computational/interolog prediction
```

The exact hierarchy should emerge from the evidence rather than being assumed in advance.

Explicitly discuss limitations of plant assays. For example, Y2H is heterologous and may not reproduce plant-specific compartmentalization/PTMs, whereas BiFC occurs *in planta* but irreversible complementation can stabilize weak/spurious interactions. AP-MS identifies complex association rather than necessarily direct binding.

---

# 6. Research Question C — Cross-database redundancy

For Arabidopsis and rice, estimate or calculate overlap among major sources.

For example:

$$
BioGRID\cap IntAct
$$

$$
BioGRID\cap HitPredict
$$

$$
IntAct\cap HitPredict
$$

and unique contributions:

$$
BioGRID-(IntAct\cup others).
$$

Where feasible, normalize protein identifiers before comparing.

Produce a table such as:

| Source  | Unique PPIs | Shared with others | Source-unique |
| ------- | ----------: | -----------------: | ------------: |
| BioGRID |           ? |                  ? |             ? |
| IntAct  |           ? |                  ? |             ? |
| XL-MS   |           ? |                  ? |             ? |
| POPPIN  |           ? |                  ? |             ? |

Determine whether combining databases materially increases biological coverage or primarily duplicates evidence.

---

# 7. Research Question D — Publication and experimental-screen bias

Determine whether the PPI network is dominated by a small number of high-throughput publications.

For each major species calculate, where possible:

$$
N_{\text{PPI/publication}}.
$$

Report:

* median PPIs/publication;
* largest contributing publication;
* top 5 publications;
* percentage of all PPIs contributed by top 1/5/10 publications.

Also identify the assays used by those publications.

We want to determine whether:

$$
\text{dataset}\approx\text{one/few high-throughput screens}.
$$

If so, evaluate whether **publication-held-out testing** should form part of the benchmark.

---

# 8. Research Question E — Evidence replication

Determine how many PPIs have:

* one experimental observation;
* ≥2 experiments;
* ≥2 assay types;
* ≥2 independent publications.

Where possible distinguish:

$$
\text{same assay repeated}
$$

from:

$$
\text{orthogonal validation}.
$$

For example:

$$
Y2H + Y2H
$$

is weaker corroboration than:

$$
Y2H + coIP.
$$

Assess whether there are enough multiply supported PPIs to construct a separate **high-confidence gold-standard test set**.

---

# 9. Research Question F — Recent experimental interactomes

Conduct a targeted search for plant interactome publications from approximately **2023–September 2026**.

Prioritize studies producing hundreds/thousands of PPIs rather than individual interaction studies.

Search particularly for:

* proteome-scale Y2H;
* XL-MS;
* AP-MS;
* protein arrays;
* proximity labeling;
* systematic binary interaction screens.

For each study record:

| Field                         | Information  |
| ----------------------------- | ------------ |
| Citation                      | Authors/year |
| Species                       |              |
| Assay                         |              |
| Number proteins               |              |
| Number interactions           |              |
| Direct vs association         |              |
| Raw data available?           |              |
| Database accession            | PRIDE/etc.   |
| Novel interactions            |              |
| Previously known interactions |              |
| Identifier format             |              |
| Potential benchmark role      |              |

Give particular attention to the 2026 Arabidopsis XL-MS resource and the new rice POPPIN resource identified in preliminary research.

---

# 10. Research Question G — Rice data sufficiency

This requires a dedicated investigation.

Determine the total number of **unique experimentally supported rice PPIs**, not simply database evidence records.

Investigate:

* POPPIN;
* BioGRID;
* IntAct;
* DIP;
* STRING experimental evidence;
* literature interactomes;
* systematic rice Y2H screens.

Determine:

$$
N_{\text{rice PPIs}}
$$

after:

1. removing self-interactions if appropriate;
2. removing duplicates;
3. removing computational predictions;
4. requiring physical experimental evidence;
5. normalizing identifiers.

Break these PPIs down by:

* assay;
* publication;
* number of proteins;
* evidence confidence;
* direct versus association evidence.

Then answer:

> **Is the current rice experimental PPI corpus sufficiently large and diverse for supervised model training?**

Do not answer solely from raw interaction count.

A dataset of 10,000 PPIs from one screen is qualitatively different from 10,000 PPIs supported by 500 independent studies.

---

# 11. Research Question H — Other species data sufficiency

Repeat a lighter version of the rice analysis for:

* maize;
* soybean;
* tomato;
* wheat;
* other promising species.

For each species determine:

$$
N_{\text{usable experimental PPIs}}
$$

and:

$$
N_{\text{unique proteins}}.
$$

Classify each species provisionally as:

**Training-capable**

Enough independent experimentally supported PPIs for supervised learning.

**Evaluation-capable**

Too small for substantial training but potentially large enough for a cross-species test.

**Insufficient**

Too few reliable PPIs even for robust evaluation.

Do not establish universal numerical cutoffs in advance; justify classifications from dataset size, diversity and statistical power.

---

# 12. Research Question I — Negative data availability

Search specifically for **experimentally assayed negative interactions**.

This is different from random negative sampling.

For example, systematic Y2H studies may test:

$$
A-B
$$

and obtain no interaction signal.

Determine whether the raw screening datasets preserve:

* tested positive pairs;
* tested negative pairs;
* replicate information;
* assay quality/control information.

Search Arabidopsis and rice high-throughput screens particularly carefully.

Estimate:

$$
N_{\text{experimentally tested negatives}}.
$$

If substantial, these could provide an unusually valuable negative benchmark.

However, explicitly acknowledge that a negative Y2H result does not prove absence of interaction *in planta*, because false negatives can arise from expression, folding, localization, missing cofactors/PTMs, or unsuitable fusion orientation.

---

# 13. Research Question J — Negative sampling alternatives

Review recent PPI-prediction literature—not only plant-specific papers—for modern negative-sampling strategies.

Compare:

* random negatives;
* subcellular localization-based negatives;
* expression-incompatible negatives;
* experimentally tested negatives;
* degree-matched negatives;
* sequence-similarity-matched negatives;
* hard-negative sampling;
* bait/prey screen-aware negatives.

For each method discuss:

$$
\text{false-negative risk}
$$

and

$$
\text{shortcut-learning risk}.
$$

Determine whether recent benchmark papers recommend **multiple negative test sets** rather than one presumed ground-truth negative set.

---

# 14. Research Question K — Historical-negative validation

Investigate whether historical database snapshots are available.

If feasible, identify pairs that would have been classified as "unknown/non-interacting" in an older release but were subsequently experimentally validated.

Conceptually:

$$
N_{2020}
=
\text{pairs sampled as negatives using 2020 knowledge}
$$

then query 2026 knowledge:

$$
N_{2020}\cap PPI_{2026}.
$$

This could empirically quantify contamination of random negatives.

Determine whether such an analysis is feasible using archived BioGRID/IntAct releases.

---

# 15. Research Question L — Sequence redundancy and paralogy

Review recent PPI benchmark methodology concerning sequence leakage.

Determine commonly used thresholds:

$$
30\%,40\%,50\%,70\%,90\%.
$$

Investigate:

* CD-HIT;
* MMseqs2;
* BLAST-based clustering.

Critically evaluate both:

$$
\text{sequence identity}
$$

and:

$$
\text{alignment coverage}.
$$

Recommend how proteins should be clustered before train/test splitting.

Because plants have extensive paralogy arising from whole-genome duplication, explicitly investigate whether conventional 40% identity filtering sufficiently prevents family-level leakage.

---

# 16. Research Question M — Orthology and cross-species leakage

If Arabidopsis and rice are both used for training/evaluation, quantify the potential problem of orthologous proteins.

Example:

```text
Arabidopsis training:
A -- B

Rice test:
A' -- B'

A' = close ortholog of A
B' = close ortholog of B
```

Technically this is cross-species.

But biologically it may amount to interolog transfer.

Therefore research methods for controlling:

* ortholog leakage;
* paralog leakage;
* conserved-domain leakage.

Determine whether cross-species benchmarks should include:

1. unrestricted species transfer;
2. ortholog-filtered species transfer.

Both may be useful because they answer different questions.

---

# 17. Research Question N — Temporal splitting

Determine whether publication dates can be reliably associated with individual PPIs.

Investigate whether we can construct:

$$
D_{\leq2023}
$$

and:

$$
D_{2024-2026}.
$$

For newer PPIs determine whether they were genuinely absent from earlier databases.

Assess the feasibility of:

```text
TRAIN
PPIs known before cutoff

TEST
PPIs experimentally discovered after cutoff
```

Particularly evaluate whether the 2026 Arabidopsis XL-MS dataset could serve as a temporal test.

Quantify overlap between its interactions and older BioGRID/IntAct/STRING releases.

---

# 18. Research Question O — PPLM-specific pretraining contamination

This is essential for this project.

PPLM's pair-language-model pretraining included interacting pairs from PDB and STRING; therefore a test interaction could have been encountered during unsupervised paired pretraining. The project literature extraction records PDB and STRING as the major PPLM paired-pretraining sources. 

Investigate whether the exact STRING/PDB releases used by PPLM can be reconstructed.

For candidate plant benchmark PPIs classify where possible:

$$
E_0=\text{neither sequence seen}
$$

$$
E_1=\text{sequences seen individually}
$$

$$
E_2=\text{pair seen during paired pretraining}
$$

$$
E_3=\text{pair used for supervised PPI training}.
$$

The most important distinction is \(E_1\) versus \(E_2\).

Do not incorrectly describe ordinary sequence exposure during self-supervised pLM pretraining as label leakage.

---

# 19. Research Question P — Identifier and sequence reference policy

Investigate current identifier standards for each species.

For Arabidopsis consider:

* TAIR locus IDs;
* UniProt accessions;
* current TAIR12 annotations;
* canonical versus alternative isoforms.

For rice consider:

* RAP-DB;
* MSU identifiers;
* UniProt;
* cultivar/subspecies differences.

Determine the mapping coverage and ambiguity between these systems.

The final benchmark should ideally contain:

```text
species
taxonomy_id

gene_id
canonical_protein_id
uniprot_accession

original_database_id
original_database

sequence
sequence_version
proteome_release
```

Research how many interactions would be lost under strict unambiguous mapping.

---

# 20. Research Question Q — Homomer policy

Quantify how many experimental interactions are:

$$
A-A
$$

versus:

$$
A-B,\quad A\neq B.
$$

Determine whether current ML PPI benchmarks generally:

* remove homomers;
* retain them;
* evaluate them separately.

Assess whether the new XL-MS data provide enough Arabidopsis homomers to support a separate benchmark.

---

# 21. Research Question R — Host–pathogen interactions

Determine how many database PPIs involve:

$$
Plant-Plant
$$

versus:

$$
Plant-Pathogen
$$

versus:

$$
Plant-Virus.
$$

The primary dataset should not silently mix these tasks.

Investigate whether there is enough host-pathogen data for a future separate benchmark, but do not automatically include it in the main plant–plant PPI training set.

---

# 22. Research Question S — Network and annotation shortcuts

Review recent literature on shortcut learning in PPI prediction.

Investigate whether models can exploit:

* protein degree;
* sequence length;
* publication count;
* GO annotation richness;
* subcellular localization;
* protein-family identity;
* assay origin.

Look for recommended controls.

An especially useful diagnostic would be whether simple classifiers based only on metadata can predict the labels.

This informs future dataset construction even if that experiment is performed later.

---

# 23. Research Question T — Class imbalance

Do not automatically inherit the DeepAraPPI/ESMAraPPI:

$$
1:10
$$

ratio.

Review modern PPI benchmark practice.

Compare:

$$
1:1,\quad1:5,\quad1:10,\quad1:50,\quad1:100
$$

and realistic proteome-scale prevalence.

Determine whether training and evaluation should use different ratios.

For example:

$$
Train=1:10
$$

while testing could include:

$$
1:10
$$

and a much more imbalanced proteome-like benchmark.

Focus particularly on AUPRC sensitivity to prevalence.

---

# 24. Required analysis of multi-species training

After gathering the species inventory, explicitly compare three potential strategies.

### Strategy A — Arabidopsis only

```text
Train:
Arabidopsis

Test:
Arabidopsis
Rice
other plants
```

Advantages include direct comparability to DeepAraPPI/ESMAraPPI and a clean cross-species test.

### Strategy B — Arabidopsis + rice

```text
Train:
Arabidopsis
Rice

Test:
held-out Ara
held-out Rice
third species
```

Investigate whether there is enough rice data and how species imbalance should be handled.

### Strategy C — multi-plant

```text
Train:
Arabidopsis
Rice
Maize
Tomato
Soybean
...
```

Only recommend this if the non-Arabidopsis/non-rice species have enough experimentally supported data to contribute meaningful supervision.

For each strategy evaluate:

* training sample size;
* species diversity;
* assay diversity;
* sequence diversity;
* taxonomic coverage;
* risk of species imbalance;
* ability to construct independent cross-species tests.

---

# 25. Consider species-balanced sampling

If multi-species training is feasible, research alternatives to simply concatenating datasets.

If:

$$
N_{Ara}\gg N_{Rice},
$$

naive pooling gives the model overwhelmingly Arabidopsis supervision.

Consider:

$$
P(s)=\frac{1}{|S|}
$$

followed by sampling within species.

Also investigate:

* species-weighted losses;
* inverse-frequency sampling;
* temperature-based sampling;
* capped per-species sampling.

The research report should recommend appropriate strategies if multi-species training is pursued.

---

# 26. Do not confuse three different generalization problems

The final research report must preserve this distinction.

### Within-species protein generalization

```text
Train: Arabidopsis proteins A...M
Test:  Arabidopsis proteins N...Z
```

### Cross-species generalization

```text
Train: Arabidopsis
Test: Rice
```

### Cross-family/cross-homology generalization

```text
Train clusters ≠ Test clusters
```

A test can satisfy one without satisfying the others.

Recommend benchmark partitions that allow these effects to be measured independently.

---

# 27. Required source hierarchy

Prioritize sources in approximately this order:

1. primary peer-reviewed publications;
2. official database documentation/releases;
3. primary preprints for very recent resources;
4. associated repositories/datasets;
5. reviews;
6. secondary sources only when primary information is unavailable.

For every numerical dataset claim, record:

* source;
* version/date;
* whether the number refers to evidence records, unique pairs or proteins.

Do not rely on search-result snippets for final conclusions.

---

# 28. Important source-age requirement

Because this research specifically concerns whether the situation has changed since DeepAraPPI/ESMAraPPI, prioritize information current to **September 2026**.

For continuously updated databases, report exact release versions.

Do not quote an old BioGRID count as though it represented the current database.

---

# 29. Required treatment of uncertainty

Use explicit categories:

**Confirmed**

Directly supported by primary source/database documentation.

**Likely**

Supported indirectly but exact underlying data could not be retrieved.

**Unknown**

Required information could not be established.

Never convert "not reported" into zero.

For example:

> Number of maize experimental PPIs: **not determined**

is preferable to inferring that BioGRID contains the entire maize literature.

---

# 30. Required final deliverables

The research should produce **six outputs**.

### Deliverable 1 — Plant PPI inventory

A table such as:

| Species     | Sources | Raw evidence | Unique experimental PPIs | Proteins | Direct PPIs | Publications | Main assays |
| ----------- | ------- | -----------: | -----------------------: | -------: | ----------: | -----------: | ----------- |
| Arabidopsis | ...     |          ... |                      ... |      ... |         ... |          ... | ...         |
| Rice        | ...     |          ... |                      ... |      ... |         ... |          ... | ...         |
| Maize       | ...     |          ... |                      ... |      ... |         ... |          ... | ...         |
| Tomato      | ...     |          ... |                      ... |      ... |         ... |          ... | ...         |
| Soybean     | ...     |          ... |                      ... |      ... |         ... |          ... | ...         |

---

### Deliverable 2 — Source comparison

For each major database/resource:

| Resource   | Species | Strength | Weakness | Evidence detail | Recommended use |
| ---------- | ------- | -------- | -------- | --------------- | --------------- |
| BioGRID    |         |          |          |                 |                 |
| IntAct     |         |          |          |                 |                 |
| HitPredict |         |          |          |                 |                 |
| POPPIN     |         |          |          |                 |                 |
| XL-MS      |         |          |          |                 |                 |
| etc.       |         |          |          |                 |                 |

---

### Deliverable 3 — Evidence-quality framework

Recommend how interactions should eventually be classified, for example:

```text
Gold
Strong direct evidence / replicated

Silver
Single direct experimental evidence

Bronze
Physical/co-complex association

Proximity
Proximity-label evidence

Exclude
Predicted/interolog/text-mined only
```

The exact scheme should be evidence-derived.

---

### Deliverable 4 — Negative-data analysis

Compare all viable negative strategies and identify whether experimentally screened negatives exist at sufficient scale.

Provide a recommended hierarchy rather than immediately choosing one.

---

### Deliverable 5 — Multi-species feasibility assessment

For every investigated plant species assign:

```text
TRAINING-CAPABLE
EVALUATION-CAPABLE
INSUFFICIENT
```

with quantitative justification.

Then explicitly address:

> Is Arabidopsis-only training still justified?

> Is Arabidopsis + rice training feasible?

> Is broader multi-plant training feasible?

> Which species should remain completely held out for cross-species testing?

---

### Deliverable 6 — Dataset design decision matrix

Finish with unresolved choices and the evidence supporting each option:

| Decision           | Option A    | Option B             | Evidence | Recommendation confidence |
| ------------------ | ----------- | -------------------- | -------- | ------------------------- |
| Positive sources   | IntAct only | Multi-database       | ...      |                           |
| Positive evidence  | Direct only | Direct + association | ...      |                           |
| Training species   | Ara         | Ara + rice           | ...      |                           |
| Negatives          | Random      | Biological           | ...      |                           |
| Sequence cutoff    | 40%         | Other                | ...      |                           |
| Temporal split     | No          | Yes                  | ...      |                           |
| Homomers           | Exclude     | Separate benchmark   | ...      |                           |
| Cross-species test | Rice        | Other species        | ...      |                           |

Do **not** present arbitrary choices as settled when the evidence remains ambiguous.

---

# 31. Final synthesis the research agent should provide

The report should finish by proposing **2–3 candidate benchmark architectures**, not one prematurely fixed design.

For example:

### Candidate A — Conservative Arabidopsis benchmark

```text
High-confidence Arabidopsis experimental PPIs
              ↓
sequence clustering
              ↓
Ara train / validation / cluster-disjoint test
              +
Rice external test
              +
2026 temporal XL-MS test
```

### Candidate B — Arabidopsis + rice plant model

```text
Arabidopsis PPIs + Rice PPIs
              ↓
evidence harmonization
              ↓
cross-species sequence clustering
              ↓
species-balanced training
              ↓
Ara/Rice held-outs
              +
third-species external test
```

### Candidate C — Multi-plant benchmark

```text
Ara + Rice + other sufficiently supported plants
              ↓
leave-one-species-out evaluation
              ↓
cluster-disjoint evaluation
              ↓
temporal evaluation
```

For each, state what scientific claim it would allow us to make.

For example, Candidate A primarily supports:

> "generalization to unseen Arabidopsis proteins and transfer to another plant species."

Candidate B supports a stronger claim:

> "learning a PPI decision function from multiple evolutionarily distinct plant species."

Candidate C could support:

> "generalization of a plant-wide PPI predictor to previously unseen plant species."

Those claims should **not** be treated as interchangeable.

---

## Most important instruction to the research agent

The objective is **not to maximize the number of PPIs**.

The objective is to determine how to build a dataset where model performance genuinely measures:

$$
\boxed{\text{interaction generalization}}
$$

rather than:

$$
\text{protein memorization}
+
\text{homology leakage}
+
\text{database bias}
+
\text{assay bias}
+
\text{negative-sampling artifacts}.
$$

Accordingly, a smaller, well-controlled dataset may be scientifically preferable to a much larger but weakly curated one.

Once this reconnaissance is complete, we should have enough information to write the actual **dataset-construction specification**: exact databases/releases → inclusion/exclusion criteria → identifier normalization → evidence aggregation → negative generation → sequence clustering → species allocation → train/validation/test assignment → leakage audits → final benchmark files.
