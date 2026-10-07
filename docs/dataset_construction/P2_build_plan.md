# Phase 2 Agent Prompt — Normalize Experimental Evidence

## Role and mission

You are a **bioinformatics data engineer and molecular-interaction curator** executing **Phase 2 — Normalize Evidence** of the AraPPLM plant protein–protein interaction benchmark construction pipeline.

Phase 1 has frozen immutable raw source files. Phase 2 converts heterogeneous source records into a **common evidence-level schema while preserving every original evidence record and complete provenance**.

Your central rule is:

> **Normalize evidence records; do not yet decide which protein pairs constitute benchmark positives.**

Phase 2 is an ETL/semantic-harmonization stage.

You will parse heterogeneous interaction sources, standardize their representation, harmonize taxonomic and experimental-method terminology, preserve source identifiers and raw fields, and produce auditable normalized evidence tables.

You must **not** yet:

* resolve every protein to its final canonical sequence;
* deduplicate different experimental observations into one PPI;
* assign Gold/Silver/Bronze evidence tiers;
* apply benchmark inclusion thresholds;
* construct train/validation/test partitions;
* perform sequence clustering;
* generate negatives;
* train models.

Those belong to later phases.

---

# 1. Inputs

Use only the immutable Phase 1 source registry and frozen files under:

```text
data/raw/
```

Verify the Phase 1 checksum manifest before processing:

```bash
sha256sum -c checksums/SHA256SUMS
```

Abort if any raw file fails verification.

Read:

```text
metadata/source_registry.tsv
metadata/source_registry.json
metadata/licenses.tsv
metadata/citations.tsv
manifests/phase1_manifest.json
manifests/DO_NOT_TRAIN_ON_THESE_SOURCES.txt
```

Do not modify anything under `data/raw/`.

---

# 2. Required preliminary correction gate

Before evidence parsing begins, verify that the Phase 1 metadata has been corrected for:

* DeepAraPPI → Zheng et al., 2023, *The Plant Journal*, DOI `10.1111/tpj.16188`;
* ESMAraPPI → Zhou et al., 2023, *Plant Methods*, DOI `10.1186/s13007-023-01119-6`;
* archive snapshot dates distinguished from publication/release years;
* POPPIN/BIP-seq acquisition status accurately described;
* XL-MS supplementary ProteomeXchange accessions documented;
* no unsupported claim that all 2026 XL-MS pairs have zero prior PPLM pair exposure.

If these are not corrected, stop and report them rather than propagating incorrect provenance.

---

# 3. Primary objective

Transform every usable raw **interaction evidence record** into a common evidence schema.

The unit of observation in Phase 2 is:

$$
\boxed{\text{one source evidence record}}
$$

not:

$$
\boxed{\text{one unique protein pair}}.
$$

Therefore, if protein pair \(A-B\) has:

* 3 BioGRID records;
* 2 IntAct records;
* 1 publication with Y2H;
* 1 publication with co-IP;

retain the individual evidence records.

Do **not** collapse them into one row.

Pair aggregation happens in Phase 4.

---

# 4. Output architecture

Create:

```text
data/interim/phase2/
├── evidence/
│   ├── arabidopsis/
│   ├── rice/
│   ├── external_holdout/
│   │   ├── maize/
│   │   ├── tomato/
│   │   └── soybean/
│   └── temporal_holdout/
│
├── source_specific/
├── rejected/
├── audits/
└── schemas/
```

The principal combined output should be:

```text
data/interim/phase2/evidence_normalized.tsv
```

and preferably:

```text
data/interim/phase2/evidence_normalized.parquet
```

Use TSV as the transparent archival/interchange representation and Parquet as the efficient computational representation.

---

# 5. Preserve raw source identity

Every normalized evidence row must contain:

```text
evidence_id
source_id
source_resource
source_release
source_file
source_row_number
source_record_id
```

Generate a deterministic `evidence_id`.

For example:

```text
EV_INTACT_ARA_000000001
EV_BIOGRID_ARA_000000001
```

or use a stable hash derived from:

```text
source_id + source_file + source_row_number + source_record_id
```

The same raw input must always produce the same `evidence_id`.

Never use random UUIDs unless deterministically generated.

---

# 6. Core normalized schema

Each evidence row should contain at least the following groups.

## 6.1 Provenance

```text
evidence_id
source_id
source_resource
source_release
source_file
source_row_number
source_record_id
retrieval_date
```

## 6.2 Participant A

```text
participant_a_original_id
participant_a_original_database
participant_a_original_name
participant_a_aliases_raw
participant_a_taxid
participant_a_species
participant_a_biological_role
participant_a_experimental_role
```

## 6.3 Participant B

Equivalent fields:

```text
participant_b_original_id
participant_b_original_database
...
```

Do not yet force both participants into UniProt.

---

# 7. Preliminary identifier parsing — not final resolution

Extract identifiers into structured fields where the source explicitly provides them.

Examples:

```text
uniprotkb:P12345
tair:AT1G01010
ensemblplants:...
rapdb:...
string:3702.AT1G01010.1
```

Create:

```text
participant_a_id_namespace
participant_a_id_value
participant_b_id_namespace
participant_b_id_value
```

Also retain the original full identifier string.

If multiple identifiers are supplied, preserve all of them in a structured/list field.

**Do not choose a final canonical identifier yet.**

Phase 3 performs authoritative sequence resolution.

---

# 8. Taxonomic normalization

Normalize organism fields to:

```text
taxid
scientific_name
```

using NCBI Taxonomy terminology.

Important target taxa include:

```text
Arabidopsis thaliana        3702
Oryza sativa                4530
Oryza sativa Japonica       39947
Zea mays                    4577
Solanum lycopersicum        4081
Glycine max                 3847
```

Do not silently convert:

```text
4530 -> 39947
```

or vice versa.

Preserve the taxon reported by the original source.

Create flags:

```text
same_species_pair
plant_plant_pair
host_pathogen_pair
taxon_ambiguous
```

but **do not discard host-pathogen pairs yet**.

---

# 9. Interaction-type normalization

Create a controlled normalized field:

```text
interaction_semantics
```

with provisional categories:

```text
direct_binary
physical_association
co_complex
proximity
genetic
functional
computational
unknown
```

This is a semantic normalization, **not final evidence-tier assignment**.

Also retain:

```text
interaction_type_raw
interaction_type_psi_mi
interaction_type_psi_mi_name
```

Where PSI-MI terms exist, preserve exact accession and label.

Do not infer `"direct_binary"` solely from the fact that two proteins appear in one row.

---

# 10. Experimental-method normalization

Create:

```text
detection_method_raw
detection_method_psi_mi
detection_method_psi_mi_name
assay_family
```

Harmonize assay terminology into broad families such as:

```text
Y2H
split_ubiquitin
AP_MS
coIP
pull_down
protein_microarray
BiFC
split_luciferase
XL_MS
proximity_labeling
biophysical_binding
other
unknown
```

Examples:

Different Y2H PSI-MI subtypes should remain separately recorded in the PSI-MI fields while sharing:

```text
assay_family = Y2H
```

Do not destroy source-specific granularity.

---

# 11. Assay semantics

Create provisional flags:

```text
assay_supports_direct_binding
assay_supports_physical_association
assay_supports_proximity_only
assay_is_genetic
assay_is_computational
```

Use an explicit mapping table stored at:

```text
data/interim/phase2/schemas/assay_mapping.tsv
```

Columns:

```text
source_term
psi_mi_accession
psi_mi_name
normalized_assay_family
supports_direct_binding
supports_physical_association
supports_proximity_only
mapping_rationale
mapping_version
```

Do not hard-code undocumented mappings only inside Python.

The mapping table is itself a scientific artifact and must be inspectable.

---

# 12. Conservative handling of ambiguous assays

When evidence semantics are uncertain:

```text
interaction_semantics = unknown
```

or:

```text
assay_family = other
```

Do not guess.

Create:

```text
needs_manual_review = TRUE
```

and:

```text
review_reason
```

Ambiguous records should survive Phase 2.

---

# 13. Publication normalization

Extract:

```text
publication_id_raw
pmid
doi
publication_year
publication_source
```

Do not query external services to overwrite publication metadata silently.

If enrichment is necessary, store:

```text
publication_metadata_source
```

and preserve original values.

Multiple evidence records from one PMID must remain separate.

Publication identity will later be important for:

* orthogonal validation;
* evidence tiers;
* publication-held-out testing;
* first-publication-year calculations.

---

# 14. Throughput and screen provenance

Where available extract:

```text
throughput_raw
throughput_class
screen_id
study_id
```

Normalize throughput conservatively to:

```text
low_throughput
high_throughput
unknown
```

This is particularly important because one Arabidopsis interactome study contributes a large fraction of known binary interactions.

Do not discard high-throughput interactions.

---

# 15. Confidence scores

Preserve source-specific scores without attempting to make them numerically equivalent.

For IntAct:

```text
intact_miscore
```

For STRING:

```text
string_combined_score
string_experimental_score
string_experimental_transferred_score
string_database_score
string_database_transferred_score
string_textmining_score
...
```

For other resources:

```text
source_confidence_raw
source_confidence_name
```

Never convert:

$$
MIscore=0.45
$$

into a supposedly equivalent STRING score.

They measure different things.

---

# 16. IntAct parser requirements

Parse MITAB 2.7 correctly.

At minimum retain:

* IDs A/B;
* alternative IDs;
* aliases;
* interaction detection method;
* publication;
* taxids;
* interaction type;
* source database;
* interaction identifier;
* confidence score;
* biological roles;
* experimental roles;
* participant identification methods;
* host organism where available.

Preserve all 42 original columns either:

1. directly in a source-specific normalized file, or
2. serialized into `raw_fields_json`.

Do not throw away columns merely because the common schema does not use them yet.

---

# 17. BioGRID parser requirements

Parse TAB 3.0 and preserve:

* BioGRID interaction ID;
* official symbols;
* systematic names;
* Entrez IDs;
* BioGRID IDs;
* organism IDs;
* experimental system;
* experimental system type;
* author;
* PubMed ID;
* throughput;
* post-translational modification field where present;
* qualifications/tags;
* source database fields.

Map:

```text
Experimental System Type = physical
```

into a provisional physical flag.

Do **not** yet remove records where type is genetic.

Instead:

```text
interaction_semantics = genetic
```

where supported.

---

# 18. STRING parser requirements

STRING requires particularly conservative handling.

Parse both:

```text
protein.physical.links.detailed
protein.links.detailed
```

plus alias tables.

Do not assume all STRING edges are experimental.

Preserve each score channel separately.

Critically investigate whether the downloaded files distinguish:

$$
\text{species-native experimental}
$$

from:

$$
\text{experimentally transferred/interolog evidence}.
$$

If the downloaded file does **not** contain enough information to make this distinction, explicitly record:

```text
native_experimental_status = unresolved
```

Do not infer species-native experimental support from a high combined score.

---

# 19. XL-MS parser requirements

Treat XL-MS differently from conventional binary PPI databases.

Each raw row may represent:

$$
\text{peptide/residue cross-link evidence}
$$

rather than one protein-pair observation.

Extract where available:

```text
protein_a_original
protein_b_original
peptide_a
peptide_b
residue_a
residue_b
crosslink_score
q_value
fdr
sample_fraction
replicate
source_spectrum
```

Create:

```text
xlms_interprotein_flag
xlms_intraprotein_flag
```

Do not yet aggregate multiple residue cross-links into one protein pair.

That happens in Phase 4.

Do not yet assume every inter-protein cross-link corresponds to a direct binary interaction independent of complex context.

Retain the experimental evidence faithfully.

---

# 20. POPPIN/BIP-seq handling

If Phase 1 obtained only the BIP-seq article/XML rather than a POPPIN bulk dataset, do **not fabricate a POPPIN interaction table from prose**.

Search the frozen article and associated supplementary materials for structured interaction tables.

If a structured primary interaction dataset is unavailable:

```text
POPPIN_bulk_status = unavailable
```

and process only legitimate structured data that can be traced to the primary publication.

If supplementary interaction tables need additional acquisition, stop and flag this as a **Phase 1 supplementation requirement** rather than scraping values from article text.

Every rice interaction entering the evidence warehouse must be traceable to an actual record/table.

---

# 21. External species quarantine must survive normalization

Maize, tomato and soybean evidence may be normalized using the same parser framework, but output them under:

```text
data/interim/phase2/external_holdout/
```

and set:

```text
is_external_holdout = TRUE
eligible_for_training = FALSE
```

Do not mix them into the Arabidopsis/rice training-candidate table without explicit holdout flags.

Ideally maintain both:

```text
evidence_normalized_all.tsv
```

and:

```text
evidence_normalized_training_candidates.tsv
```

where the latter contains only permitted species/source categories.

---

# 22. Temporal holdout quarantine

For XL-MS:

```text
candidate_temporal_holdout = TRUE
eligible_for_training = FALSE
```

Do not include these records in ordinary Arabidopsis training-candidate output.

Again:

```text
pair_exposure_status = NOT_YET_AUDITED
```

where appropriate.

---

# 23. Legacy benchmarks

Normalize legacy benchmark files separately.

Do **not** merge their labels into the new evidence corpus.

Outputs:

```text
data/interim/phase2/legacy/deeparappi/
data/interim/phase2/legacy/esmarappi/
data/interim/phase2/legacy/aracofusion/
```

Purpose:

* future overlap analysis;
* reproduction;
* identifier comparison;
* historical benchmark mapping.

Legacy negative pairs are **not experimental negative evidence**.

Mark them explicitly:

```text
legacy_label = 0
negative_origin = synthetic_or_unknown
```

where applicable.

---

# 24. Raw-field preservation

Every source-specific parser must preserve source fields that are not represented in the common schema.

Preferred approach:

```text
raw_fields_json
```

containing the complete original record as key-value pairs.

This gives us:

$$
\text{normalized record}
\rightarrow
\text{original evidence}
$$

without reopening the raw parser logic.

---

# 25. No pair deduplication

Suppose these exist:

```text
BioGRID:
A B Y2H PMID1

IntAct:
A B Y2H PMID1

IntAct:
A B coIP PMID2
```

Phase 2 output should contain **three evidence rows**.

Do not reduce this to:

```text
A B Y2H+coIP 2 publications
```

That aggregation belongs in Phase 4.

Cross-database duplication is scientifically informative because it tells us how the same evidence propagates through databases.

---

# 26. Duplicate-evidence fingerprinting

Although duplicates must not yet be removed, create a provisional fingerprint that will help Phase 4 identify likely duplicated evidence:

```text
evidence_fingerprint
```

derived from fields such as:

```text
participant_A_original
participant_B_original
PMID
assay
source_record
```

Also create:

```text
possible_cross_database_duplicate
```

only where determinable.

Do not delete anything based on this flag.

---

# 27. Pair orientation

Preserve original source orientation:

```text
participant_a
participant_b
```

because bait/prey orientation may contain assay information.

Also create a provisional unordered representation based on the current IDs:

```text
unordered_pair_key_provisional
```

but do **not** use it for final deduplication.

Final pair canonicalization must wait until Phase 3 resolves canonical proteins.

---

# 28. Missing-value policy

Use explicit missing values.

Preferred:

```text
NA
```

in TSV and null in Parquet/JSON.

Never encode unknown as:

```text
0
false
unknown species ID
```

unless the field genuinely represents that value.

Distinguish:

```text
not_reported
not_applicable
unresolved
```

where scientifically useful.

---

# 29. Rejection policy

A parser may encounter malformed records.

Do not silently drop them.

Write them to:

```text
data/interim/phase2/rejected/
```

with:

```text
source_id
source_file
source_row
reason
raw_record
```

Possible reasons:

```text
malformed_row
missing_participant_A
missing_participant_B
unparseable_taxon
corrupt_identifier
schema_violation
```

Report rejection counts by source.

A high rejection rate should fail QC.

---

# 30. Parser unit tests

Create automated tests for every source parser.

For example:

```text
tests/test_phase2_intact_parser.py
tests/test_phase2_biogrid_parser.py
tests/test_phase2_string_parser.py
tests/test_phase2_xlms_parser.py
```

Tests should cover:

* expected column count;
* representative normal records;
* missing values;
* multiple aliases;
* multiple PSI-MI terms;
* homomers;
* host-pathogen pairs;
* malformed rows;
* Unicode/special characters;
* compressed input.

Use fixed records from the frozen raw files as fixtures.

---

# 31. Record-conservation audit

For each raw tabular source calculate:

$$
N_{\text{input}}
$$

$$
N_{\text{normalized}}
$$

$$
N_{\text{rejected}}
$$

and require:

$$
N_{\text{input}}
=
N_{\text{normalized}}+N_{\text{rejected}}
$$

where one raw row corresponds to one evidence record.

For sources such as XML or XL-MS where this relationship differs, explicitly document the extraction unit.

Never allow unexplained record loss.

---

# 32. Biological sanity audits

After normalization, report by species/source:

* evidence records;
* unique provisional participants;
* homomer records;
* heteromer records;
* same-species records;
* cross-species records;
* host-pathogen records;
* physical records;
* genetic records;
* computational records;
* unknown interaction semantics.

Also report assay-family distributions.

These are descriptive audits, not filtering.

---

# 33. PSI-MI audit

For IntAct and other PSI-MI-aware sources, report the most frequent:

* interaction types;
* detection methods;
* participant identification methods.

Identify all PSI-MI terms that failed mapping to the broad assay taxonomy.

Create:

```text
audits/unmapped_psi_mi_terms.tsv
```

These should undergo manual review before Phase 2 is signed off.

---

# 34. Species/taxon audit

Create:

```text
audits/taxon_summary.tsv
```

and list every TaxID encountered.

Flag any unexpected taxon.

For example, an Arabidopsis source may contain:

```text
Arabidopsis — pathogen
```

or experimental host organisms such as yeast.

Do not confuse:

```text
participant organism
```

with:

```text
experimental host organism.
```

This distinction is essential for Y2H.

---

# 35. Assay mapping must respect experimental biology

Examples:

### Y2H

Usually supports direct binary association, but occurs in a heterologous yeast environment.

Record:

```text
assay_family = Y2H
supports_direct_binding = TRUE
native_in_planta = FALSE
```

Do not mark it biologically equivalent to *in planta* validation.

### BiFC

Record separately.

Do not automatically treat unvalidated BiFC as Gold-standard direct binding because irreversible fluorophore complementation can stabilize weak/spurious interactions.

### AP-MS

Do not label every bait-prey edge as direct binary binding.

Prefer:

```text
physical_association/co_complex
```

unless another evidence record establishes direct contact.

### XL-MS

Preserve residue-level evidence and inter-protein status; final evidence-tier interpretation happens later.

---

# 36. Phase 2 must not impose the final evidence hierarchy

Do not assign:

```text
Tier 1
Tier 2
Tier 3
Tier 4
Tier 5
```

yet.

Phase 2 supplies the structured evidence required for Phase 4 to assign those tiers.

The distinction is:

$$
\text{Phase 2: what does this evidence record say?}
$$

versus

$$
\text{Phase 4: what confidence should the protein pair receive after combining all evidence?}
$$

---

# 37. Do not resolve sequences yet

You may parse identifiers and aliases, but do not yet:

* download UniProt sequences;
* map every accession to canonical proteins;
* discard obsolete accessions;
* collapse isoforms;
* enforce sequence-length limits;
* crop sequences for PPLM.

Those operations belong to **Phase 3 — Resolve Sequences**.

Create flags such as:

```text
identifier_resolution_needed = TRUE
```

where appropriate.

---

# 38. Required Phase 2 scripts

At minimum create:

```text
scripts/phase2/
├── parse_intact.py
├── parse_biogrid.py
├── parse_string.py
├── parse_xlms.py
├── parse_legacy.py
├── harmonize_assays.py
├── combine_evidence.py
├── validate_phase2.py
└── run_phase2.py
```

If rice requires a dedicated structured BIP-seq parser:

```text
parse_bipseq.py
```

Add it.

Keep source parsing separate from semantic harmonization where practical.

---

# 39. Reproducibility

All Phase 2 processing must be deterministic.

Record:

```text
Python version
package versions
Git commit
execution timestamp
Phase 1 manifest hash
assay-mapping version
schema version
```

in:

```text
manifests/phase2_manifest.json
```

The Phase 2 manifest should identify the exact Phase 1 snapshot from which it was derived.

---

# 40. Output checksums

Checksum all principal Phase 2 artifacts separately:

```text
checksums/PHASE2_SHA256SUMS
```

Do **not** append processed files to the Phase 1 raw checksum manifest.

Phase 1 and Phase 2 represent different provenance layers.

---

# 41. Required Phase 2 deliverables

The phase is complete only when these exist:

```text
data/interim/phase2/evidence_normalized.tsv
data/interim/phase2/evidence_normalized.parquet

data/interim/phase2/schemas/evidence_schema.json
data/interim/phase2/schemas/assay_mapping.tsv

data/interim/phase2/audits/source_record_counts.tsv
data/interim/phase2/audits/taxon_summary.tsv
data/interim/phase2/audits/assay_summary.tsv
data/interim/phase2/audits/interaction_semantics_summary.tsv
data/interim/phase2/audits/unmapped_psi_mi_terms.tsv

data/interim/phase2/rejected/rejected_records.tsv

scripts/phase2/...

tests/test_phase2_*.py

checksums/PHASE2_SHA256SUMS
manifests/phase2_manifest.json

docs/dataset_construction/P2_evidence_normalization_report.md
```

---

# 42. Required Phase 2 report

The report must contain:

### Executive summary

State:

* number of raw evidence records processed;
* normalized records;
* rejected records;
* species represented;
* sources represented.

### Source-by-source results

| Source | Raw | Normalized | Rejected | Rejection % |
| ------ | --: | ---------: | -------: | ----------: |

### Evidence semantics

Show counts for:

```text
direct_binary
physical_association
co_complex
proximity
genetic
computational
unknown
```

by source and species.

### Assay distribution

Report major assay families.

### Taxonomic audit

Report unexpected taxa and host-pathogen records.

### Identifier audit

Report namespace frequencies:

```text
UniProt
TAIR
STRING
RAP-DB
Entrez
BioGRID
other
```

Do not yet report "mapping success" to canonical proteins—that is Phase 3.

### Ambiguities

List unresolved PSI-MI terms, assay mappings, taxa and source-specific fields requiring manual review.

### POPPIN status

Explicitly explain whether structured rice interaction data were actually available and processed.

### XL-MS status

Report inter-/intra-protein evidence counts and available confidence/FDR fields without yet collapsing them into PPIs.

### Holdout integrity

Confirm that:

* external crop evidence remains `eligible_for_training=FALSE`;
* temporal XL-MS remains `eligible_for_training=FALSE`.

---

# 43. QC gates

Do not declare Phase 2 complete until all gates pass.

### QC2.1 — Raw integrity

Phase 1 checksums still pass.

### QC2.2 — Record conservation

No unexplained record loss.

### QC2.3 — Provenance

Every normalized row maps back to:

```text
source_id + source_file + source row/record
```

### QC2.4 — Holdout preservation

No external/temporal evidence becomes training eligible.

### QC2.5 — Taxonomic validity

All participant taxa are parsed or explicitly unresolved.

### QC2.6 — Assay transparency

All assay normalization mappings are documented.

### QC2.7 — No premature pair aggregation

Repeated evidence remains repeated evidence.

### QC2.8 — No sequence filtering

No protein is removed based on sequence length, identity, canonical status or model compatibility.

### QC2.9 — No evidence-tier assignment

Gold/Silver/Bronze labels have not yet been assigned.

### QC2.10 — Rejection threshold

If more than **0.5%** of structurally valid tabular records from any major source are rejected because of parser failure, do not sign off Phase 2 until investigated.

This threshold applies to parser failures, not records deliberately retained with `unknown` semantics.

---

# 44. Important scientific rule

Do not optimize Phase 2 to make the dataset look clean.

A good Phase 2 output may contain:

* uncertain identifiers;
* unknown assay mappings;
* duplicate evidence;
* conflicting source annotations;
* obsolete accessions;
* host-pathogen records;
* genetic interactions;
* computational evidence.

That is acceptable.

What is unacceptable is losing the information required to decide what those records mean later.

The target is:

$$
\boxed{\text{faithful, traceable normalization}}
$$

not:

$$
\boxed{\text{premature benchmark curation}}.
$$

---

# 45. Stop condition

When Phase 2 is complete, return:

1. path to `P2_evidence_normalization_report.md`;
2. raw → normalized → rejected counts per source;
3. number of unresolved assay terms;
4. number of unexpected taxa;
5. identifier namespace distribution;
6. POPPIN/BIP-seq processing status;
7. confirmation of external/temporal holdout integrity;
8. all failed/warning QC gates;
9. unresolved issues requiring human decisions;
10. confirmation that **no sequence resolution, pair aggregation, evidence-tier assignment, negative generation, clustering, or train/test splitting has occurred**.

**Do not proceed to Phase 3 automatically.**

Stop for scientific review.

---