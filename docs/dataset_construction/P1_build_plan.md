# Phase 1: Freeze and Register Raw Plant PPI Data Sources

## Role

You are acting as a **bioinformatics data engineer and scientific data curator** constructing a reproducible benchmark for plant protein–protein interaction prediction.

You are executing **Phase 1 — Freeze Sources** of a larger plant PPI benchmark construction pipeline.

Your task is to acquire, freeze, verify, checksum, document, and register the exact raw source datasets that will later be processed into the benchmark.

**Do not normalize identifiers, filter interactions, deduplicate pairs, assign evidence tiers, generate negatives, retrieve canonical protein sequences, or construct train/test splits during this phase.**

Those operations belong to later phases.

The governing principle is:

> Preserve the original evidence exactly as obtained from the source before making any biological or computational transformations.

The output of Phase 1 must make it possible for another researcher to determine exactly **what data were downloaded, from where, when, under which release/version, and whether the local copy is byte-identical to the frozen source used for benchmark construction.**

---

# 1. Scientific context

The benchmark will ultimately support two supervised configurations:

1. **Arabidopsis-only benchmark**, providing continuity with DeepAraPPI and ESMAraPPI.
2. **Arabidopsis + rice benchmark**, intended to evaluate whether supervision from both a dicot and monocot improves plant-general PPI prediction.

Arabidopsis currently has the strongest evidence base, while rice has become sufficiently data-rich to warrant supervised inclusion. Maize, tomato, and soybean are intended primarily as **strictly held-out external species** rather than training species.

The benchmark must preserve distinctions between:

* direct binary interactions;
* physical associations/co-complex interactions;
* proximity evidence;
* computational/interolog predictions;
* multiple evidence records supporting the same protein pair.

Do **not** collapse these distinctions in Phase 1.

---

# 2. Primary objective

Create an **immutable raw-data snapshot** of all primary data sources required for subsequent benchmark construction.

For every source:

1. identify the exact release/version;
2. retrieve the original raw files;
3. retain the original files unchanged;
4. compute cryptographic checksums;
5. record acquisition metadata;
6. document licensing/access restrictions;
7. verify that the downloaded files are readable and plausibly correspond to the intended species/resource;
8. produce a machine-readable source registry;
9. produce a human-readable acquisition report.

No raw file may be modified after freezing.

---

# 3. Required source categories

Acquire sources in four groups.

## Group A — Arabidopsis core PPI sources

### A1. IntAct / IMEx

Target organism:

* *Arabidopsis thaliana*
* NCBI Taxonomy ID: **3702**

Preferred release:

* **IntAct Release 252, January 2026**, if the exact release remains retrievable.

Retrieve the richest available interaction format, preferably:

* PSI-MI TAB / MITAB;
* PSI-MI XML if available;
* associated release metadata/documentation.

Do **not** filter to MIscore ≥0.45 yet.

Do **not** filter to `"direct interaction"` or `"physical association"` yet.

Acquire the complete relevant raw records so filtering remains reproducible later.

Record:

* release;
* source URL;
* retrieval date/time;
* format;
* taxonomic scope;
* file size;
* checksum;
* license;
* official documentation URL.

---

### A2. BioGRID

Target:

* *Arabidopsis thaliana*
* TaxID 3702.

Preferred release:

* **BioGRID 5.0.261 (2026)**.

Retrieve an organism-specific file if the exact release provides one.

Otherwise retrieve the appropriate complete release file from which Arabidopsis records can later be extracted.

Prefer BioGRID's most information-rich tabular format.

Do **not** remove:

* genetic interactions;
* self-interactions;
* duplicates;
* non-physical records;

during Phase 1.

Preserve them because Phase 2/4 must document precisely why each record is retained or excluded.

Also retrieve the BioGRID release README/documentation describing the columns and experimental-system terminology.

---

### A3. STRING

Target:

* *Arabidopsis thaliana*, TaxID 3702.

Target version:

* **STRING v12.0**.

Acquire the raw files necessary to reconstruct:

* physical interaction links;
* experimental-channel scores;
* individual evidence channels;
* protein aliases/identifier mappings where provided.

Do **not** simply download a prefiltered list with:

$$
experimental\_score \ge 700
$$

because the benchmark needs to retain the original evidence values and perform filtering later.

Critically, retain enough STRING information to distinguish where possible:

* species-native experimental evidence;
* transferred/interolog evidence;
* database evidence;
* text-mining evidence;
* combined score.

STRING must **not** be treated as experimental ground truth merely because an edge has a nonzero `"experimental"` score.

---

### A4. 2026 Arabidopsis PhoX XL-MS dataset

Target study:

**Trinh et al., PNAS 2026 — Arabidopsis proteome-wide cross-linking mass spectrometry.**

Expected ProteomeXchange/PRIDE accessions include:

* `PXD066234`
* `PXD066291`

Verify these against the primary publication/repository rather than assuming they are correct.

Acquire, where available:

* final processed cross-link identification tables;
* inter-protein cross-link tables;
* protein-level interaction tables;
* peptide-level/residue-level cross-link results;
* search-result tables;
* metadata describing FDR;
* supplementary tables from the publication.

**Raw mass-spectrometry instrument files do not need to be downloaded unless required to reconstruct the published interaction list.**

The goal is to freeze the processed primary data needed to reproduce the published PPI evidence, while avoiding hundreds of GB of unnecessary raw spectra.

Document exactly which files were acquired and which were intentionally omitted.

This dataset is intended to become a potential **temporal/structural holdout**, so provenance is especially important.

---

# 4. Group B — Rice training sources

Target organism:

* *Oryza sativa*
* distinguish TaxID **4530** and *O. sativa* subsp. *japonica* TaxID **39947**.

Do not silently merge taxa during Phase 1.

## B1. POPPIN

Target:

**POPPIN — Port of Protein-Protein Interactomes**, 2026.

Because this resource is very recent, first verify:

* official publication/preprint;
* official database/repository;
* version/release date;
* downloadable datasets;
* field definitions.

Retrieve the **rawest available downloadable interaction/evidence tables**, not a manually cleaned subset.

It is particularly important to preserve fields distinguishing:

* experimentally measured PPIs;
* literature-derived interactions;
* text-mined/LLM-derived interaction clues;
* computational predictions;
* assay type;
* source publication;
* source database.

Do **not** assume all ~150k POPPIN entries are physical PPIs.

The benchmark will later include only appropriately supported experimental interactions.

If no bulk download is available, document this explicitly and investigate whether:

* an API exists;
* supplementary data accompany the preprint;
* an official GitHub/Zenodo repository exists.

Do not scrape the website in a manner prohibited by its terms of service.

---

## B2. IntAct rice

Retrieve the current/frozen IntAct records relevant to:

* *Oryza sativa*;
* *O. sativa* japonica where separately represented.

Prefer using the same IntAct release frozen for Arabidopsis.

Do not filter the interactions during acquisition.

---

## B3. BioGRID rice

Retrieve rice records from the same BioGRID release used for Arabidopsis.

Again, retain all raw evidence.

---

## B4. STRING rice

Acquire equivalent STRING v12.0 files for rice sufficient to reconstruct:

* physical links;
* experimental evidence;
* transferred evidence;
* aliases.

Do not apply score thresholds yet.

---

# 5. Group C — external zero-shot species

Acquire experimental PPI sources for:

### Maize

*Zea mays* — TaxID 4577.

### Tomato

*Solanum lycopersicum* — TaxID 4081.

### Soybean

*Glycine max* — TaxID 3847.

These species are intended to remain **completely outside model training and model selection**.

At minimum acquire:

* BioGRID;
* IntAct;
* relevant organism-specific curated resource identified in the evidence report;
* STRING raw files needed to distinguish species-native evidence from transferred interolog evidence.

Do not treat the very large STRING `"experimental"` channel counts for these species as direct experimental PPIs without verifying provenance.

Maintain these files under a clearly separated directory such as:

```text
raw/external_holdout/
```

to reduce the chance that they are accidentally included in training.

---

# 6. Group D — legacy benchmark datasets

Freeze the published datasets from:

### DeepAraPPI

Zheng et al., 2023.

Acquire:

* original Arabidopsis positive/negative datasets;
* C1;
* C2;
* C3;
* rice dataset;
* any relevant source-data files.

### ESMAraPPI

Zhou et al., 2023.

Acquire:

* C1;
* C2;
* C3;
* original supporting dataset files if available.

### ARACoFusion

Acquire any released datasets necessary to establish whether they are byte-/record-equivalent to ESMAraPPI and DeepAraPPI sources.

These are **not inputs to the new positive evidence corpus by default**.

They are frozen because we need:

1. historical benchmark comparability;
2. overlap analyses;
3. leakage analyses;
4. reproduction of published results.

Keep them under:

```text
raw/legacy_benchmarks/
```

rather than mixing them with contemporary source evidence.

---

# 7. Data-provenance requirements

This is mandatory.

Every downloaded file must have an associated provenance record.

Create:

```text
metadata/source_registry.tsv
```

and preferably:

```text
metadata/source_registry.json
```

with at least these fields:

```text
source_id
resource_name
resource_category
species
taxonomy_id
release_version
release_date
retrieval_date_utc
retrieval_url
landing_page_url
accession
original_filename
local_relative_path
file_format
compression
file_size_bytes
sha256
md5_if_source_provides_it
license_name
license_url
citation
doi
pmid
data_role
expected_evidence_type
is_training_candidate
is_external_holdout
is_temporal_holdout
download_method
download_command_or_script
notes
```

Use controlled values where practical.

Example:

```text
data_role =
    core_positive_source
    supplemental_source
    temporal_holdout
    external_species_holdout
    legacy_benchmark
    documentation
```

---

# 8. Cryptographic integrity

For every downloaded file compute:

```bash
sha256sum <file>
```

SHA-256 is mandatory.

If the original source publishes MD5/SHA checksums, retain and verify those too.

Create:

```text
checksums/SHA256SUMS
```

containing every immutable raw file.

After freezing, run:

```bash
sha256sum -c checksums/SHA256SUMS
```

and require:

```text
100% PASS
```

before Phase 1 can be considered complete.

Do not use timestamps alone as evidence of file identity.

---

# 9. Raw-data immutability

Use this directory model:

```text
data/
├── raw/
│   ├── arabidopsis/
│   │   ├── intact/
│   │   ├── biogrid/
│   │   ├── string/
│   │   └── xlms_2026/
│   │
│   ├── rice/
│   │   ├── poppin/
│   │   ├── intact/
│   │   ├── biogrid/
│   │   └── string/
│   │
│   ├── external_holdout/
│   │   ├── maize/
│   │   ├── tomato/
│   │   └── soybean/
│   │
│   └── legacy_benchmarks/
│       ├── deeparappi/
│       ├── esmarappi/
│       └── aracofusion/
│
├── metadata/
├── checksums/
├── manifests/
└── logs/
```

**Never write transformed data into `data/raw/`.**

Later phases should write to directories such as:

```text
data/interim/
data/processed/
```

Raw files must remain byte-identical to their downloaded versions.

---

# 10. Preserve documentation as data

For every database release, also save relevant:

* README;
* column specification;
* ontology description;
* release notes;
* license;
* evidence-code documentation.

For example, knowing that a BioGRID field contains `"Affinity Capture-MS"` is insufficient unless we also freeze the documentation defining the field.

Store these alongside the relevant source or under:

```text
metadata/source_documentation/
```

and checksum them as well.

---

# 11. Acquisition logging

All downloads must be scripted.

Create:

```text
scripts/phase1_download_sources.sh
```

or an equivalent reproducible Python downloader where necessary.

Do not make the final acquisition depend solely on manual browser downloads.

For every operation log:

```text
timestamp
source
requested URL
HTTP status
downloaded filename
bytes
checksum
success/failure
```

to:

```text
logs/phase1_acquisition.log
```

If authentication or manual acquisition is unavoidable, record that explicitly.

---

# 12. Source authenticity rules

Use only:

1. official database download servers;
2. official institutional repositories;
3. journal supplementary-data repositories;
4. ProteomeXchange/PRIDE;
5. author-designated GitHub/Zenodo repositories where appropriate.

Do not download benchmark data from:

* random mirrors;
* Kaggle reposts;
* unofficial GitHub forks;
* third-party file-sharing sites,

unless the primary source is unavailable and the exception is explicitly documented.

---

# 13. Do not silently substitute database versions

If the requested historical version is unavailable, **do not replace it silently with the latest version**.

For example, if IntAct release 252 cannot be downloaded:

record:

```text
requested_release: 252
status: unavailable
available_release: XXX
decision: pending
```

Then explain the problem in the Phase 1 report.

Similarly, if BioGRID 5.0.261 has been superseded, try to retrieve **5.0.261 from the release archive** because reproducibility requires an exact snapshot.

---

# 14. Taxonomic provenance

Do not collapse species/subspecies during acquisition.

For example, retain separately:

```text
Oryza sativa
TaxID 4530
```

and:

```text
Oryza sativa Japonica Group
TaxID 39947
```

if the source distinguishes them.

The later normalization phase will decide how these map to the chosen rice reference proteome.

Likewise, explicitly detect host-pathogen records where one participant belongs to another taxon.

Do not remove them yet.

---

# 15. No biological filtering in Phase 1

This constraint is important.

Do **not** remove records because they are:

* genetic interactions;
* self-interactions;
* proximity interactions;
* computational predictions;
* transferred interologs;
* low-confidence;
* duplicated;
* old identifiers;
* host-pathogen interactions;
* non-canonical isoforms.

The purpose of Phase 1 is preservation.

If a source allows downloading a sufficiently complete raw dataset, freeze that.

Filtering decisions must be explicit transformations in later phases.

---

# 16. Preliminary validation is allowed — modification is not

You should inspect the frozen files to verify that they appear valid.

For every tabular source report:

```text
number_of_rows
number_of_columns
column_names
first_3_rows
last_3_rows
parse_errors
```

For compressed archives:

* verify decompression succeeds;
* list archive contents;
* do not replace the original archive with only extracted data.

For FASTA/XML/JSON:

* perform syntax/readability checks.

These are **validation operations only**.

Do not rewrite or "clean" the original files.

---

# 17. Expected-count sanity checks

Compare raw downloaded data against approximate counts from the evidence-mapping stage, but treat those counts only as **sanity checks**, not targets.

For example, the evidence mapping estimated approximately:

* Arabidopsis: ~18,500 unique physical PPIs after future harmonization/filtering;
* rice: ~2,527 direct PPIs and ~14,200 broader associations;
* maize: ~850 usable physical PPIs;
* tomato: ~740;
* soybean: ~620.

These are **not expected raw row counts**.

Raw files may contain far more records because they include:

* repeated experimental observations;
* genetic interactions;
* multiple publications;
* multiple assay records;
* predictions;
* duplicated pair evidence.

Do not manipulate data to reproduce these estimates.

The final counts must emerge from the actual pipeline.

---

# 18. Licensing and redistribution audit

For each source determine:

* license;
* whether redistribution is allowed;
* whether modified redistribution is allowed;
* whether attribution is required;
* whether users should instead download the data themselves.

Create:

```text
metadata/licenses.tsv
```

with:

```text
resource
version
license
license_url
redistribution_allowed
attribution_required
notes
```

If the license is unclear, mark:

```text
UNKNOWN
```

rather than guessing.

This matters because the final benchmark may eventually be publicly released.

---

# 19. Citation registry

Create:

```text
metadata/citations.bib
```

containing the canonical citation for every major data resource.

Also create:

```text
metadata/citations.tsv
```

with:

```text
source_id
resource
citation
doi
pmid
publication_year
```

Dataset provenance should remain traceable all the way from a future benchmark pair back to the originating database/publication.

---

# 20. Machine-readable manifest

Create:

```text
manifests/phase1_manifest.json
```

with a structure conceptually equivalent to:

```json
{
  "phase": 1,
  "name": "freeze_sources",
  "created_utc": "...",
  "sources": [
    {
      "source_id": "...",
      "resource": "...",
      "release": "...",
      "species": "...",
      "taxid": "...",
      "files": [
        {
          "path": "...",
          "sha256": "...",
          "bytes": 123
        }
      ]
    }
  ]
}
```

The exact schema may be improved, but it must be deterministic and documented.

---

# 21. Special handling of external test data

Maize, tomato, soybean, and the intended 2026 temporal dataset must be visibly marked as **restricted benchmark holdouts**.

The source registry must contain:

```text
is_external_holdout = TRUE
```

or:

```text
is_temporal_holdout = TRUE
```

as appropriate.

This does not technically prevent later leakage, but creates a provenance-level safeguard.

Also generate:

```text
manifests/DO_NOT_TRAIN_ON_THESE_SOURCES.txt
```

listing the relevant source IDs.

---

# 22. Special handling of the 2026 XL-MS data

Do not yet assume that every reported XL-MS pair is a novel post-2024 interaction.

Phase 1 should establish only:

> This dataset was generated/published/deposited in 2026.

Later phases must audit whether individual pairs already existed in pre-2024 PDB/STRING or other sources.

Therefore do **not** label every pair:

```text
PPLM_exposure = E1
```

during Phase 1.

Instead record:

```text
candidate_temporal_holdout = TRUE
pair_exposure_status = NOT_YET_AUDITED
```

This prevents an unsupported leakage claim.

---

# 23. Special handling of POPPIN

Because POPPIN is a 2026 resource and may combine experimental evidence with text-mining/LLM-derived information, **do not use the total POPPIN entry count as the number of PPIs**.

The Phase 1 report must explicitly document:

* what constitutes one POPPIN record;
* which fields indicate experimental evidence;
* which fields identify source publication/database;
* whether raw evidence records are available;
* whether interaction confidence is supplied;
* whether downloadable records distinguish direct binding from physical association.

If these cannot be determined, mark them as unresolved questions for Phase 2.

---

# 24. Required deliverables

Phase 1 is complete only when the following exist.

### Deliverable 1 — Frozen raw datasets

All successfully acquired source files under:

```text
data/raw/
```

with no transformations.

### Deliverable 2 — Source registry

```text
metadata/source_registry.tsv
metadata/source_registry.json
```

### Deliverable 3 — Cryptographic manifest

```text
checksums/SHA256SUMS
```

and verification output showing all files pass.

### Deliverable 4 — Download pipeline

```text
scripts/phase1_download_sources.sh
```

and/or reproducible source-specific acquisition scripts.

### Deliverable 5 — Acquisition log

```text
logs/phase1_acquisition.log
```

### Deliverable 6 — Licensing table

```text
metadata/licenses.tsv
```

### Deliverable 7 — Citation registry

```text
metadata/citations.tsv
metadata/citations.bib
```

### Deliverable 8 — Machine-readable phase manifest

```text
manifests/phase1_manifest.json
```

### Deliverable 9 — Holdout warning manifest

```text
manifests/DO_NOT_TRAIN_ON_THESE_SOURCES.txt
```

### Deliverable 10 — Human-readable Phase 1 report

Create:

```text
docs/dataset_construction/phase1_source_freeze_report.md
```

---

# 25. Required contents of the Phase 1 report

The report must contain:

## Executive summary

State:

* what was successfully acquired;
* what could not be acquired;
* whether requested versions matched actual versions;
* whether any licensing/access issues occurred.

## Source inventory

Table:

| Source ID | Resource | Version | Species | Files | Size | SHA256 verified | Role |
| --------- | -------- | ------- | ------- | ----: | ---: | --------------- | ---- |

## Acquisition provenance

For every resource explain:

* official source;
* exact URL/accession;
* release;
* retrieval date;
* download method.

## Raw-data sanity checks

Report:

* row counts;
* schemas;
* parse status;
* obvious taxonomic issues;
* archive integrity.

Do not report these as final benchmark counts.

## Version deviations

Explicitly document cases such as:

```text
Requested: IntAct 252
Retrieved: IntAct 253
Reason: release 252 archive unavailable
```

Do not conceal substitutions.

## Licensing

Summarize redistribution implications.

## Unresolved issues

Examples:

* POPPIN lacks bulk download;
* uncertain license;
* ambiguous release metadata;
* inaccessible historical release;
* XL-MS supplementary table missing;
* source documentation inconsistent.

## Phase 1 completion status

For each source:

```text
COMPLETE
PARTIAL
BLOCKED
```

with an explanation.

---

# 26. QC gates

Do **not** declare Phase 1 complete unless:

### QC1 — Provenance completeness

Every raw file has:

$$
\text{source}+\text{version}+\text{retrieval date}+\text{URL/accession}.
$$

### QC2 — Integrity

Every raw file has SHA-256.

### QC3 — Registry consistency

Every file in:

```text
data/raw/
```

appears in the source registry.

Every source-registry raw file exists on disk.

### QC4 — No accidental transformation

No normalized/filtered/deduplicated interaction files are stored as though they were source files.

### QC5 — Holdout separation

External-species and temporal sources are explicitly marked.

### QC6 — Reproducibility

A second researcher can determine how every file was obtained.

### QC7 — Licensing

Every resource has either:

```text
known license
```

or explicitly:

```text
UNKNOWN
```

No license should be inferred.

---

# 27. Git/data-versioning rules

Commit to Git:

```text
scripts/
metadata/
manifests/
checksums/
docs/
```

Do not commit multi-GB raw datasets unless the repository is explicitly configured for Git LFS/DVC.

Add appropriate raw-data paths to `.gitignore`.

Where possible, use DVC or an equivalent data-versioning system to associate immutable files with the repository while preserving checksums.

The Git commit associated with completion of Phase 1 should be recorded in:

```text
manifests/phase1_manifest.json
```

if available.

---

# 28. Error-handling rules

Never silently skip a failed source.

If acquisition fails:

1. retry only where appropriate;
2. investigate official archive/repository alternatives;
3. record the failure;
4. preserve partial files only if clearly marked `.partial`;
5. do not include incomplete files in the frozen checksum manifest;
6. mark the source `PARTIAL` or `BLOCKED`.

Do not substitute data from a secondary source without explicitly documenting the change.

---

# 29. What you must NOT do

During Phase 1, do not:

* normalize UniProt IDs;
* map TAIR/RAP IDs;
* fetch canonical sequences for interaction proteins;
* collapse isoforms;
* deduplicate PPIs;
* remove self-interactions;
* remove host-pathogen interactions;
* apply MIscore thresholds;
* apply STRING thresholds;
* classify evidence tiers;
* remove genetic interactions;
* generate negatives;
* cluster sequences;
* construct C1/C2/C3;
* construct train/validation/test files;
* calculate model performance.

If you find something biologically questionable, **record it in the report** rather than fixing it.

---

# 30. Final response required from the agent

When Phase 1 is complete, return:

1. the path to `phase1_source_freeze_report.md`;
2. a concise status table for every requested source;
3. total number and storage size of frozen files;
4. checksum verification status;
5. blocked/partial sources;
6. unresolved questions requiring human decisions;
7. confirmation that **no biological filtering or normalization has yet been performed**;
8. whether the pipeline is ready to proceed to **Phase 2 — Normalize Evidence**.

Do not proceed to Phase 2 automatically.

Stop and await review.

---