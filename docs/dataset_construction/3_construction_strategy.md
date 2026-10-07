# Proposed Dataset Construction Strategy for a Modern Plant PPI Benchmark

**Project:** AIS5281 — Plant Protein–Protein Interaction Prediction Using Protein Language Models  
**Prepared:** September 2026  
**Primary objective:** Construct a reproducible, leakage-controlled benchmark for training and evaluating plant PPI predictors, while preserving compatibility with historical Arabidopsis benchmarks and enabling rigorous cross-species generalization tests.

---

## Executive recommendation

The benchmark should be built as a **versioned plant PPI evidence warehouse plus multiple derived benchmark views**, rather than as one irreversible binary interaction table. The primary supervised benchmark (**Benchmark A**) focuses on **high-confidence direct physical interactions from Arabidopsis** with strict sequence-cluster splitting (40% identity, $\ge 70\%$ coverage), preserving Rice and PhoX XL-MS strictly as external holdout benchmarks. A dual-species pan-plant model (**Benchmark B: Arabidopsis + Rice**) is maintained as an exploratory research configuration. Maize, tomato, and soybean remain completely external to model development and are used for zero-shot cross-species evaluation.

The design should separate four questions that previous benchmarks partially conflate:

1. **Evidence quality:** is a positive pair supported by direct binding, co-complex association, proximity, or only prediction?
2. **Negative uncertainty:** is a negative an unobserved pair, an assay non-hit, or a biologically challenging hard negative?
3. **Generalization:** are test proteins, homologous sequence clusters, publications, and species independent of training?
4. **Pretraining exposure:** was a test pair already present in PPLM's PDB/STRING paired pretraining corpus?

The adopted primary model-development dataset is **Benchmark A (Arabidopsis direct PPIs, Tiers 1–2)**, while Tier 3 associations, homomers, publication-held-out interactions, recent XL-MS interactions, and alternative negative types become separate auxiliary or stress-test suites.

The evidence gathering and Phase 2 normalization audits establish approximately **18,500 unique experimentally supported physical Arabidopsis PPIs across 7,820 proteins**, of which approximately 11,200 are direct binary interactions. For Rice, Phase 2 curation identified **923 physical pairs** (spanning 732 proteins across 37 non-overlapping publications) and **356 direct binary pairs** in IntAct and BioGRID (POPPIN bulk data was unavailable as flat-files, and STRING v12.0 `exp >= 700` comprises 2,527 computational interolog transfers). Furthermore, retraining showed that Arabidopsis-specific tuning reduced Rice transfer performance ($0.4297 \rightarrow 0.3555$), confirming that holding Rice out as an external test set (Benchmark A) provides the most rigorous, leakage-free benchmark. In contrast, maize (~850 physical pairs), tomato (~740), and soybean (~620) are substantially smaller and serve as secondary crop holdouts.

---

# 1. Design principles

## 1.1 Build an evidence warehouse first; derive benchmark datasets second

Do not merge databases directly into a final `protein_A, protein_B, label` file. Maintain two linked layers:

### Evidence-level table
One row per experimental evidence record, retaining:

- species and taxonomy ID;
- original protein/gene identifiers;
- normalized identifiers;
- protein A and B canonical sequences/version;
- source database and exact release;
- source record identifier;
- PMID/publication year;
- experimental method and PSI-MI term where available;
- interaction type (direct, physical association, proximity, etc.);
- throughput/screen identifier;
- IntAct MIscore or other source confidence score;
- tissue/compartment if supplied;
- whether evidence is species-native or transferred;
- provenance/retrieval date.

### Pair-level table
Aggregate all evidence for an unordered canonical pair `{A,B}` and derive:

- number of independent publications;
- number of experimental observations;
- number of orthogonal assay classes;
- strongest evidence tier;
- all contributing sources;
- direct-binding flag;
- association flag;
- proximity-only flag;
- first publication year;
- homomer/heteromer status;
- PPLM exposure state where known.

**Rationale.** The source inventory contains extensive database overlap and heterogeneous assays. Collapsing evidence before preserving provenance would make it impossible to distinguish a repeatedly validated binary interaction from a single AP-MS co-complex observation. It would also prevent later temporal or publication-held-out analyses.

---

# 2. Species composition

## 2.1 Primary training species: Arabidopsis (Benchmark A - Primary; Benchmark B - Exploratory)

Construct two supervised configurations from the same warehouse.

### Benchmark A — Arabidopsis-only reference (Adopted Primary)
Train on Arabidopsis Tiers 1–2 with strict sequence-cluster splitting (MMseqs2 40% identity, $\ge 70\%$ coverage). This serves as the methodological bridge to DeepAraPPI/ESMAraPPI, isolates the effect of improved curation/splitting from multi-species distribution shifts, and maintains Rice (356 direct binary / 923 physical pairs) as an uncontaminated monocot external holdout benchmark.

### Benchmark B — Arabidopsis + rice pan-plant model (Exploratory)
Train on Tiers 1–2 from both species, with explicit species balancing.

**Evidence & Reality Check.** Arabidopsis has ~18,500 physical PPIs / ~11,200 direct pairs and >3,800 contributing publications. In Rice, Phase 2 canonical normalization revealed 923 physical pairs and 356 direct binary pairs across IntAct and BioGRID (the previously cited ~2,527 pairs were STRING v12.0 `exp >= 700` computational interolog transfers, and POPPIN bulk flat-files were inaccessible). Furthermore, Arabidopsis-retrained PPLM dropped performance on Rice from 0.4297 (zero-shot) to 0.3555 (-17.3%), demonstrating the hazard of cross-species distribution shift. Benchmark B is preserved for exploratory experiments examining whether adding monocot supervision or interolog priors mitigates this drop.

## 2.2 External species: maize, tomato, soybean

Keep **all usable interactions from these species completely out of model selection, hyperparameter tuning, threshold calibration, negative-sampling design, and training**.

Use them only for locked external evaluation after the model is finalized.

- Maize: ~850 physical PPIs / ~310 proteins.
- Tomato: ~740 / ~290.
- Soybean: ~620 / ~240.

Their limited size gives little training benefit but considerable scientific value as unseen crop species spanning distinct lineages. Avoid unfiltered STRING records because the report finds that non-model plant `experimental` channels can be overwhelmingly composed of transferred interolog evidence rather than species-native assays.

Wheat and *Chlamydomonas* should not be primary headline benchmarks because the report estimates only ~210 and ~410 physical pairs, respectively. They may be retained as exploratory supplementary tests.

---

# 3. Positive interaction construction

## 3.1 Source ingestion

### Arabidopsis
Use, at minimum:

1. **IntAct/IMEx Release 252** — core manually curated source with PSI-MI evidence and MIscore.
2. **BioGRID 5.0.261** — broad literature coverage; retain physical evidence only.
3. **Species-native experimental records from STRING v12.0** only after verifying that the evidence is laboratory-derived and not transferred.
4. **2026 PhoX XL-MS (Trinh et al.)** — ingest into the warehouse but reserve temporally novel pairs from ordinary model development.
5. TAIR/HitPredict as identifier/evidence corroboration layers rather than independent duplicate positives unless they contribute unique primary evidence.

### Rice
Use:

1. experimentally flagged POPPIN records;
2. IntAct/BioGRID species-native physical evidence;
3. primary rice Y2H/AP-MS/other screens represented in POPPIN or literature;
4. STRING only after excluding transferred evidence.

Do not accept LLM-text-mined POPPIN clues as positive labels without traceable experimental evidence.

## 3.2 Canonical pair normalization

Treat PPI classification as symmetric. Convert every heteromeric pair to an unordered canonical representation, e.g. lexicographically sorted stable IDs, so `(A,B)` and `(B,A)` cannot survive as duplicate examples.

Retain directionality separately only when it belongs to the assay design (bait/prey), because this can be useful for diagnosing screen bias but should not redefine the biological PPI label.

## 3.3 Inclusion/exclusion rules

Include a pair in the evidence warehouse when:

- both participants are proteins from the intended plant taxon;
- there is traceable experimental evidence;
- both identifiers can be mapped unambiguously to the chosen reference sequence system;
- the evidence is not purely computational/interolog/text-mined/co-expression based.

Exclude from the primary direct-PPI dataset:

- genetic interactions;
- functional associations without physical evidence;
- transferred STRING interologs;
- pure text-mining predictions;
- plant–pathogen and plant–virus pairs (store separately);
- ambiguous multi-gene identifiers;
- obsolete models without a defensible sequence mapping;
- Tier 3/4 evidence unless the same pair also has Tier 1/2 evidence.

---

# 4. Evidence tiers and label semantics

Use the evidence hierarchy from the research report, but do **not** convert all tiers into the same positive label for the primary task.

| Tier | Definition | Primary direct-PPI training? | Main role |
|---|---|---:|---|
| **1 — Gold** | ≥2 independent publications or ≥2 orthogonal experimental methods supporting direct binding | Yes | High-confidence train/validation; enriched test reporting |
| **2 — Silver** | One well-controlled direct binary assay | Yes | Main source of supervised positives |
| **3 — Bronze** | Physical association/co-complex evidence without direct pairwise binding proof | No, by default | Auxiliary/multi-task experiment |
| **4 — Proximity** | TurboID/miniTurbo/PUP-IT or unconfirmed BiFC-type proximity evidence | No | Separate spatial/proximity task |
| **5 — Excluded** | Interolog-only, computational, text-mining-only, co-expression-only | No | Never a positive label |

### Important refinement: XL-MS

Do not automatically classify every inter-protein XL-MS pair as equivalent to an AP-MS association. Inter-protein crosslinks provide substantially stronger spatial evidence than generic co-complex membership. Preserve residue-level crosslink/FDR/distance metadata and create an `xlms_interprotein` evidence subtype. For the main temporal test, report XL-MS performance separately from conventional direct-binary assays rather than claiming that the two labels are identical.

### Gold interactions should not all be reserved for testing

The report proposes ~3,200 Gold interactions. Reserving every Gold pair for test would distort the training distribution and unnecessarily discard the strongest supervision. Instead, split Tiers 1–2 under the same leakage constraints, then **stratify results by evidence tier**. Additionally create a smaller locked Gold-only challenge set from publications/clusters not represented in training.

---

# 5. Identifier and sequence normalization

## 5.1 Arabidopsis

Primary gene coordinate: current TAIR locus ID (`ATxGxxxxx`). Cross-reference to current UniProt accession and record the exact sequence release/version.

## 5.2 Rice

Primary gene coordinate: RAP-DB ID, cross-referenced to MSU and UniProt. Anchor the main benchmark to *Oryza sativa* subsp. *japonica* reference sequences. Do not silently map *indica* or other cultivar-specific proteins onto japonica sequences when sequence equivalence is uncertain.

## 5.3 Sequence QC

Require:

- valid standard amino-acid sequence;
- length ≥40 aa;
- no unresolved identifier ambiguity;
- explicit sequence checksum/version.

### Do not exclude proteins merely because `L > 1020`

The evidence report proposes restricting each protein to `L ≤ 1020` because of PPLM's token budget. That would make the **biological benchmark model-specific** and systematically remove long plant proteins.

Instead:

1. retain all biologically valid proteins in the master benchmark;
2. store full-length sequences;
3. create a `pplm_compatible` flag and length fields;
4. define model-specific truncation/cropping at inference/training time.

PPLM constrains the **combined pair token budget**, not the biological validity of proteins. The benchmark should therefore remain usable by future models with longer context windows.

---

# 6. Leakage-controlled sequence clustering and splitting

This is the central methodological improvement over the legacy benchmarks.

## 6.1 Cluster proteins before splitting

Run MMseqs2 jointly over the protein universe used for model development (Arabidopsis + rice), rather than independently by species.

Primary clustering criterion:

- sequence identity ≥40%;
- alignment coverage ≥70% of both proteins (bidirectional coverage).

Store cluster membership at several additional thresholds (30%, 50%, 70%, 90%) for sensitivity analysis, but use the prespecified 40%/70% rule for the primary benchmark unless pilot statistics demonstrate unacceptable data collapse.

**Rationale.** Protein-disjoint C3 partitions do not prevent a close paralog of a test protein from occurring in training. This is particularly problematic in plant genomes shaped by whole-genome duplication. Joint clustering also detects cross-species ortholog/homolog leakage.

## 6.2 Do not describe C1/C2/C3 as mutually exclusive 70/15/15 cluster partitions

The Park–Marcotte concepts concern **endpoint exposure**, not three ordinary pair partitions. Construct them explicitly relative to the training protein/cluster universe:

- **C1 / seen–seen:** both endpoint clusters represented in training. This is useful for conventional pair generalization but is not the principal test.
- **C2 / seen–unseen:** exactly one endpoint cluster represented in training.
- **C3 / unseen–unseen:** neither endpoint cluster represented in training.

C2 and C3 should be **locked evaluation sets**, not sources of supervised training.

### Recommended primary emphasis

Use C3 as the principal within-species generalization benchmark; C2 is secondary. A random pair test should be retained only for historical comparability and sanity checking.

## 6.3 Split assignment algorithm

1. Build the positive interaction graph over sequence clusters.
2. Assign clusters—not individual pairs—to train/validation/test pools using a deterministic seeded optimizer.
3. Optimize for approximately 70/10/20 positive-pair availability while maintaining sufficient C2 and C3 examples and preserving species/evidence distributions.
4. Ensure no C3 endpoint cluster occurs in train or validation.
5. Generate negatives **after** positive split assignment so that negative construction cannot reconnect train and test protein universes.
6. Freeze all split manifests before model training.

Exact final percentages should be chosen after observing graph connectivity; preserving leakage constraints is more important than forcing exactly 70/15/15.

---

# 7. Use Tier 3 associations without corrupting the direct-binding label

The association corpus is too large to ignore but too heterogeneous to merge naively.

Run a planned ablation:

### Direct-only model
Train on Tiers 1–2 only.

### Direct + association auxiliary model
Use Tier 3 as an auxiliary task or lower-weight positive class, e.g. a three-class target:

`non-interaction / physical association / direct interaction`.

Alternatively pretrain/adapt the classifier on Tier 3 and fine-tune on Tiers 1–2.

Evaluate both configurations on exactly the same locked direct-PPI tests.

This establishes empirically whether broader co-complex supervision helps or harms prediction of direct binding.

---

# 8. Negative construction

There is no single trustworthy biological negative definition. Therefore separate **training negatives** from **evaluation negative suites**.

## 8.1 Training negatives

Use **degree-aware background sampling** from within the same species and within the allowed training cluster universe.

For each positive `(A,B)`, sample candidate non-edges so that the marginal degree distribution of proteins in negatives approximately matches positives. Also match coarse sequence-length distributions to prevent trivial length shortcuts.

Exclude:

- all known positive pairs across every ingested evidence tier/source, including held-out positives;
- pairs with direct evidence in any external database snapshot used for curation;
- self-pairs when constructing the heteromer task.

Use a **1:10 positive:negative training ratio** initially for continuity with DeepAraPPI/ESMAraPPI and manageable training. Treat these labels as `unknown/non-observed`, not experimentally proven non-interactions.

Do **not** use incompatible subcellular localization as the primary negative generator. Although it lowers false-negative probability, the report identifies a severe shortcut: organelle targeting peptides and localization-associated sequence composition can make the task partly a localization classifier.

## 8.2 Locked evaluation suites

### Suite N1 — Degree-matched background
Evaluate at both:

- 1:10 for comparison with prior plant benchmarks;
- 1:100 as a more stringent low-prevalence stress test.

For fair AUPRC comparison, create the 1:10 test as a deterministic subset of the 1:100 negative pool so the positive set is identical and prevalence is the only designed difference.

### Suite N2 — Paralog/homolog hard negatives
For a positive `(A,B)`, seek a homolog `B'` that is biologically plausible but lacks positive evidence. Prefer experimentally assayed non-binding paralogs when available.

Do not call an untested homolog a verified negative. Maintain `negative_confidence` and `negative_source` fields.

### Suite N3 — Assay non-hits
Use systematic Y2H matrix non-hits where the tested pair and assay controls are recoverable. Keep this suite separate because high-throughput Y2H sensitivity is limited; absence of reporter activation is not proof of no *in planta* interaction.

### Optional Suite N4 — Expression incompatibility
Use only as a biological stress test, not the headline benchmark, because it introduces its own tissue/localization shortcuts.

---

# 9. Publication-held-out benchmark

The evidence report finds severe publication concentration: the Arabidopsis Interactome Mapping Consortium 2011 screen contributes ~5,664 interactions (~30.6% of direct binary interactions), the top five publications ~48.2%, and the top ten >57.5%.

Therefore create a **publication-held-out stress test**.

Recommended design:

1. Select one or more major high-throughput screens, including AI-1, before model development.
2. Remove all evidence uniquely attributable to those publications from training.
3. If a pair also has independent evidence in a training publication, mark it as replicated and either remove it from the held-out set or analyze it separately.
4. Apply the same sequence-cluster leakage audit to the held-out publication set.
5. Evaluate using matched negatives from the screened protein universe where possible.

This tests whether performance transfers across experimental campaigns rather than reproducing screen-specific artifacts.

---

# 10. Temporal 2026 XL-MS benchmark

Reserve the 2026 PhoX XL-MS resource as a **locked temporal benchmark** rather than incorporating it into the main training corpus.

The source report records >3,000 physical PPI complexes, 52,944 cross-linked peptide pairs, 37,531 residue contacts and 5,064 proteins. Its 2026 acquisition postdates PPLM's PDB/STRING paired-pretraining cutoff, making it uniquely valuable for testing pair-level temporal generalization.

However, temporal novelty must be established at the **pair level**, not merely from the publication date.

Construct:

- **T0 — all valid 2026 inter-protein XL-MS pairs**;
- **T1 — temporally novel pairs:** remove any pair already present in pre-cutoff BioGRID/IntAct/STRING/PDB evidence;
- **T2 — temporally novel + cluster-unseen:** additionally require both endpoint sequence clusters to be absent from supervised training.

Report T1 and T2 as the headline temporal results.

### PPLM exposure terminology

The report's claim of guaranteed `E1-only` exposure should be softened. A 2026 discovery guarantees that the **2026 experimental record** was unavailable during PPLM pretraining, but the same protein pair could conceivably have appeared in an older STRING/PDB source. Therefore classify `E2` only after reconstructing/searching the relevant pretraining snapshots.

---

# 11. Cross-species evaluation

## 11.1 Arabidopsis-only baseline

For Benchmark A, evaluate direct transfer to rice before any rice-specific training. This preserves continuity with DeepAraPPI's cross-species question.

Create two rice views:

- **standard transfer:** all eligible rice direct PPIs;
- **strict homology-filtered transfer:** remove pairs whose endpoint clusters are represented by homologous Arabidopsis training proteins under the prespecified similarity rule.

The difference quantifies how much apparent cross-species performance comes from conserved interolog-like relationships.

## 11.2 Pan-plant model

For Benchmark B, rice is no longer an external species because it contributes training data. The true zero-shot species tests become maize, tomato and soybean.

For each external species create:

- standard zero-shot set;
- strict homology/orthology-filtered set.

All external-species decisions and filtering rules must be frozen before seeing model predictions.

---

# 12. Species-balanced training

Naively pooling ~11,200 Arabidopsis direct PPIs with ~2,527 rice direct PPIs would heavily favor Arabidopsis.

Use a two-stage sampler:

1. sample species with probability `P(s) ∝ N_s^τ`;
2. sample a positive/negative example within that species.

Primary setting: **τ = 0.5** (square-root sampling), which increases rice representation without forcing every mini-batch to be exactly 50:50. Include `τ = 1` (natural-frequency pooling) and `τ = 0` (equal-species sampling) as prespecified ablations.

Do not conflate sampling with the benchmark itself: the released data should retain natural counts; balancing belongs to the training protocol.

---

# 13. Homomers

Retain homomer evidence in the warehouse but construct a separate homomer benchmark.

Primary heteromer classifier:

`A != B` only.

Homomer suite:

`A == B`, evaluated separately.

This prevents the trivial sequence identity of `(A,A)` from altering the heteromer distribution while preserving biologically important oligomerization information, including hundreds of homomultimeric observations in the recent XL-MS resource.

---

# 14. PPLM pretraining-exposure audit

Attach a pair-level field:

- `E1`: sequences individually exposed to generic pLM pretraining but pair not identified in PPLM paired corpus;
- `E2`: exact pair found in the reconstructed PPLM PDB/STRING pretraining universe;
- `E3`: supervised downstream PPI-label exposure.

For plant pairs, the report indicates no plant supervised PPLM-PPI training (`E3=0`), but exact `E2` status still requires auditing.

Report model performance separately on `E2` and non-`E2` pairs where reconstruction is possible. Do not call single-sequence exposure label leakage.

---

# 15. Class balance and evaluation metrics

The benchmark should not be summarized by accuracy or AUROC alone.

Primary metric:

- **AUPRC**, reported separately for every negative suite/prevalence.

Also report:

- AUROC;
- MCC;
- balanced accuracy;
- precision;
- recall/sensitivity;
- specificity;
- F1;
- calibration (Brier score and/or ECE) when probabilities are interpreted as confidence.

Threshold-dependent metrics must use a threshold selected on validation data only. Never optimize a threshold on the locked test set.

Because AUPRC depends strongly on class prevalence, never directly compare the numeric AUPRC of a 1:10 test with a 1:100 test without stating the prevalence.

Use bootstrap confidence intervals at the **protein/cluster level**, rather than treating correlated interaction edges as fully independent observations.

---

# 16. Dataset manifests and reproducibility

Release immutable manifests rather than only generated pair files.

Suggested structure:

```text
dataset_release/
├── README.md
├── RELEASE.json
├── evidence/
│   ├── evidence_records.tsv
│   └── pair_evidence.tsv
├── sequences/
│   ├── arabidopsis.fasta
│   ├── rice.fasta
│   └── external_species/
├── mappings/
│   ├── arabidopsis_ids.tsv
│   └── rice_ids.tsv
├── clusters/
│   ├── mmseqs_40id_70cov.tsv
│   └── sensitivity_thresholds/
├── splits/
│   ├── arabidopsis_only/
│   ├── panplant/
│   ├── publication_holdout/
│   ├── temporal_xlms/
│   └── external_species/
├── negatives/
│   ├── degree_matched/
│   ├── hard_paralog/
│   └── assayed_nonhits/
└── audits/
    ├── overlap_audit.tsv
    ├── leakage_audit.tsv
    └── pplm_pretraining_exposure.tsv
```

For every released pair include at minimum:

`pair_id, protein_a, protein_b, species, label, label_semantics, evidence_tier, split, sequence_cluster_a, sequence_cluster_b, negative_source, first_publication_year, sources, pplm_exposure_state`.

Record software versions, command lines, random seeds, source database versions, retrieval dates and file checksums.

---

# 17. Quality-control gates before release

The dataset should not be considered complete until all of the following pass.

### QC1 — Identifier integrity
100% of released pairs resolve to two exact sequence records; no ambiguous identifiers.

### QC2 — Pair uniqueness
No duplicate unordered pairs within a benchmark view; no `(A,B)`/`(B,A)` duplicates.

### QC3 — Label conflicts
No pair is both positive and negative anywhere in the warehouse. A new positive supersedes an unobserved-negative designation.

### QC4 — Split leakage
For C3, neither endpoint's 40%/70%-coverage cluster may occur in training/validation. Audit this programmatically.

### QC5 — Cross-species leakage
Strict external tests must satisfy the predefined homolog/ortholog exclusion rule.

### QC6 — Temporal leakage
Temporal-novel pairs must be absent from all designated pre-cutoff interaction sources.

### QC7 — Negative distribution audit
Compare positive versus negative distributions for protein degree, sequence length, species and annotation availability. Large separability should trigger resampling.

### QC8 — Metadata-only baseline
Train simple classifiers using no protein sequence—only degree, length, publication count/source, etc. Strong performance indicates residual dataset shortcuts and should block release until investigated.

### QC9 — Evidence audit
Randomly sample records from every tier/source and manually trace them back to the source publication/database evidence.

---

# 18. Recommended benchmark suite

The final release should contain several deliberately different tasks rather than one headline split.

| Task | Training | Positive test | Negative test | Question answered |
|---|---|---|---|---|
| **A1 Legacy-compatible** | Arabidopsis direct | Random/seen-aware Ara | 1:10 degree-matched | Comparable reference performance |
| **A2 C2** | Arabidopsis direct | One endpoint cluster unseen | Degree-matched | Partner generalization |
| **A3 C3 (primary Ara)** | Arabidopsis direct | Both endpoint clusters unseen | Degree-matched | Novel-protein/family generalization |
| **A4 Ara→Rice** | Arabidopsis direct | Rice direct | Standard + strict homology-filtered | Cross-species transfer |
| **B1 Pan-plant C3** | Ara + Rice direct | Cluster-disjoint Ara/Rice | Degree-matched | Generalization after dual-species training |
| **B2 External crops** | Ara + Rice direct | Maize/Tomato/Soybean | Standard + strict | True unseen-species transfer |
| **S1 Hard negatives** | Same model | Direct positives | Paralog hard negatives | Binding specificity beyond family identity |
| **S2 Assay non-hits** | Same model | Screen positives | Screen non-hits | Within-screen discrimination |
| **S3 Publication holdout** | Exclude chosen screens | Held-out publication positives | Screen-matched | Assay/publication robustness |
| **T1 Temporal XL-MS** | Pre-2026 development only | Novel 2026 XL-MS | Matched background | Prospective temporal generalization |
| **H1 Homomer** | Defined separately | A–A interactions | Matched proteins | Homooligomerization prediction |

This suite deliberately prevents a single favorable negative construction or split from defining whether a model is successful.

---

# 19. Planned ablations

The dataset design should support, from the beginning, the following controlled experiments:

1. **Arabidopsis-only vs Arabidopsis+rice training.** Tests whether monocot supervision improves general plant transfer.
2. **Natural vs square-root vs equal-species sampling.** Tests sensitivity to species imbalance.
3. **Direct-only vs direct+Tier-3 auxiliary supervision.** Tests whether co-complex evidence improves direct PPI prediction.
4. **40% vs alternative sequence clustering thresholds.** Quantifies sensitivity to homology stringency.
5. **Standard vs strict ortholog-filtered external tests.** Quantifies interolog contribution.
6. **1:10 vs 1:100 prevalence.** Tests precision robustness under stronger imbalance.
7. **Degree-matched vs hard-paralog vs assay-negative suites.** Tests dependence on negative definition.
8. **E2-exposed vs pair-unexposed tests for PPLM.** Quantifies paired-pretraining advantage where auditable.

These should be preregistered in the benchmark documentation to reduce post-hoc benchmark selection.

---

# 20. Decisions deliberately *not* made from the current report alone

Several numerical choices should remain provisional until the raw corpus has been ingested and audited:

- exact train/validation/test percentages after sequence clustering;
- exact number of Gold/Silver pairs after cross-database deduplication;
- exact rice direct-PPI count after POPPIN experimental-only filtering;
- final number of usable external-species PPIs;
- final MMseqs2 clustering threshold if 40%/70% causes severe graph fragmentation/data loss;
- size and validity of paralog hard-negative and assayed-negative sets;
- exact `E2` contamination status of legacy pairs.

This is intentional. The evidence report provides strong design guidance, but approximate inventory counts should not be converted into fixed benchmark sizes before reproducing them from raw records.

---

# 21. Construction workflow

## Phase 1 — Freeze sources

Download and checksum exact database releases and primary supplemental datasets. Create a source registry with version, date, URL/accession, taxon and license.

## Phase 2 — Normalize evidence

Parse all records into the evidence-level schema; harmonize PSI-MI assay terms; normalize taxon and identifiers; preserve original fields.

## Phase 3 — Resolve sequences

Map to reference sequences, retain full-length proteins, flag ambiguous mappings, compute sequence checksums and model compatibility fields.

## Phase 4 — Aggregate pairs and assign evidence tiers

Canonicalize unordered pairs, deduplicate across sources, calculate publication/assay support and assign the strongest defensible evidence tier.

## Phase 5 — Freeze special holdouts

Before ordinary splitting, identify and quarantine:

- 2026 temporal XL-MS pairs;
- external maize/tomato/soybean data;
- publication-held-out screens;
- homomer suite.

This prevents accidental contamination during later development.

## Phase 6 — Joint sequence clustering

Cluster Arabidopsis + rice development proteins with MMseqs2 and generate threshold sensitivity statistics.

## Phase 7 — Construct positive partitions

Create Arabidopsis-only and pan-plant train/validation/C2/C3 manifests under cluster constraints.

## Phase 8 — Generate negatives within each split

Generate degree-/length-matched training and test negatives only from the protein universe legally available to that split. Then generate the hard-negative and assay-non-hit suites separately.

## Phase 9 — Run leakage and shortcut audits

Check pair overlap, protein overlap, cluster overlap, temporal overlap, source/publication overlap and metadata-only predictability.

## Phase 10 — Freeze benchmark v1.0

Only after all QC gates pass, hash all manifests and declare the test sets immutable. Model tuning begins **after** this freeze.

---

# 22. Why this strategy is preferable to reproducing DeepAraPPI/ESMAraPPI

The legacy datasets remain useful comparators, but the proposed benchmark addresses several weaknesses exposed by the evidence review:

- **Larger and more current evidence base:** ~18,500 Arabidopsis physical PPIs are estimated after modern multi-source harmonization, versus 11,858 DeepAraPPI high-confidence pairs and 7,729 ESMAraPPI positives.
- **Rigorous monocot cross-species benchmark:** Rice provides 356 direct binary pairs and 923 physical pairs (across IntAct and BioGRID) as a completely un-leaked monocot test bed to evaluate whether dicot-trained models transfer across evolutionary divergence.
- **Sequence-family leakage control:** clustering is performed before splitting and jointly across training species.
- **Evidence semantics are preserved:** direct binding, co-complex association and proximity are not collapsed into one label.
- **Negative uncertainty is explicit:** multiple negative suites replace one arbitrary presumed-negative distribution.
- **Screen bias is tested:** publication-held-out evaluation addresses the finding that one Arabidopsis screen contributes ~30.6% of direct binary interactions.
- **True external species remain untouched:** maize, tomato and soybean support zero-shot evaluation rather than weakly contributing to training.
- **Temporal testing becomes possible:** the 2026 XL-MS resource can test interactions unavailable to contemporary supervised plant benchmarks, with pair-level novelty auditing.
- **PPLM-specific leakage is auditable:** pair exposure during paired pretraining is separated from ordinary single-sequence pLM exposure.

---

# 23. Final proposed dataset architecture

```text
                         RAW EXPERIMENTAL EVIDENCE
          Arabidopsis                                  Rice
   IntAct / BioGRID / verified STRING          POPPIN / IntAct / BioGRID
                 │                                      │
                 └──────────────┬───────────────────────┘
                                ▼
                    EVIDENCE-LEVEL WAREHOUSE
        assay • publication • source • score • sequence version
                                │
                                ▼
                    CANONICAL PAIR AGGREGATION
                  Tier 1 / Tier 2 / Tier 3 / Tier 4
                                │
             ┌──────────────────┼─────────────────────┐
             │                  │                     │
             ▼                  ▼                     ▼
       SPECIAL HOLDOUTS   DIRECT-PPI CORPUS     AUXILIARY CORPUS
       • 2026 XL-MS       Tier 1 + Tier 2        Tier 3 associations
       • publications            │               Tier 4 proximity
       • homomers                ▼
       • external crops   JOINT MMseqs2 CLUSTERS
                          40% identity + 70% cov
                                │
                 ┌──────────────┴──────────────┐
                 ▼                             ▼
          ARABIDOPSIS-ONLY                PAN-PLANT
             TRAIN/VAL                  ARA + RICE TRAIN/VAL
                 │                             │
        C2 / C3 Ara tests             C2 / C3 Ara/Rice tests
                 │                             │
                 ▼                             ▼
          Ara→Rice transfer        Maize / Tomato / Soybean
                 └──────────────┬──────────────┘
                                ▼
                       STRESS-TEST SUITES
             degree-matched • hard paralog • assay non-hit
                       publication • temporal
```

## Bottom line

The benchmark should **not** be a single Arabidopsis table with randomly sampled negatives. It should be a provenance-rich, multi-view resource in which the primary supervised signal is **direct physical interaction**, Arabidopsis and rice provide complementary dicot/monocot training evidence, and all claims of generalization are tied to explicitly controlled axes: protein/sequence-family novelty, species novelty, publication novelty, temporal novelty, and negative-set difficulty.

The first implementation milestone should therefore be **raw evidence ingestion and pair-level harmonization**, not model training. Once the empirical counts after normalization and tier assignment are known, the remaining provisional parameters—especially split proportions and the final size of the rice/direct and hard-negative sets—can be frozen without relying on approximate literature counts.

---

## Evidence basis

This construction strategy is derived from the supplied September 2026 evidence-mapping report, particularly its current plant PPI inventory, database comparison, five-tier evidence framework, negative-sampling analysis, species-feasibility assessment, publication-bias analysis, PPLM pretraining-exposure audit, and benchmark decision matrix. Numerical values above should therefore be treated as **planning estimates until reproduced from the raw database releases during construction**.
