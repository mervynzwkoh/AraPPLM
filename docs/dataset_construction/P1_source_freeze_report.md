# Phase 1 Source Freeze & Raw Integrity Audit Report
**Benchmark Pipeline**: Plant Protein-Protein Interaction (PPI) Generalization Benchmark (`AraPPLM`)  
**Phase**: Phase 1 — Freeze and Register Raw Plant PPI Data Sources  
**Date of Execution**: 2026-09-21  
**Execution Script**: [`scripts/phase1_download_sources.py`](../../scripts/phase1_download_sources.py)  
**Verification Script**: [`scripts/verify_phase1_integrity.py`](../../scripts/verify_phase1_integrity.py)  
**Audit Status**: **100% PASS (ALL QUALITY GATES PASSED)**  

---

## 1. Executive Summary

Phase 1 established an immutable, cryptographically verifiable, and cleanly segregated foundation of raw protein-protein interaction (PPI) datasets for model training, zero-shot crop evaluation, prospective temporal holdouts, and legacy benchmark comparisons.

### Key Milestones Achieved:
1. **Zero Transformation Guarantee**: All 37 acquired files are preserved in their exact, bit-level raw state under `data/raw/`. No filtering, score thresholding, accession mapping, identifier deduplication, or negative sampling was performed.
2. **Total Storage Footprint**: Acquired **1,086.61 MB (~1.09 GB)** across 37 raw files, well within local workstation disk capacity (75+ GB available).
3. **Cryptographic Provenance**: Every raw file is hashed with SHA-256 and cataloged in [`checksums/SHA256SUMS`](../../checksums/SHA256SUMS) (100% verification pass).
4. **Leakage Quarantine & Temporal Isolation**: A strict boundary is enforced between training sources and holdout evaluation sets via [`manifests/DO_NOT_TRAIN_ON_THESE_SOURCES.txt`](../../manifests/DO_NOT_TRAIN_ON_THESE_SOURCES.txt), quarantining external crop holdouts (Maize, Tomato, Soybean) and prospective temporal data (2026 PhoX XL-MS). Publication in August 2026 guarantees that the XL-MS dataset postdates PPLM pretraining data cutoffs; exact interaction novelty will be determined via historical overlap auditing in Phase 3.
5. **Full Deliverable Suite**: All 10 deliverables mandated by the Phase 1 Build Plan (`docs/dataset_construction/P1_build_plan.md`) have been generated and validated.

---

## 2. Complete Source Inventory Table

| Source ID | Resource | Release / Version | Archive Snapshot Date | Species | TaxID | Local Relative Path | Size | SHA-256 (Prefix) | Benchmark Role |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **SRC_A1_INTACT_ARA** | IntAct / IMEx | Release 252 (2026-01) | 2026-09-21 | *Arabidopsis thaliana* | 3702 | `data/raw/arabidopsis/intact/intact_arabidopsis_taxid3702_release252.mitab27.txt` | 157.34 MB | `4c437e16a149` | Core Training Positive |
| **SRC_A2_BIOGRID_ARA** | BioGRID | 5.0.261 (2026-08) | 2026-09-21 | *Arabidopsis thaliana* | 3702 | `data/raw/arabidopsis/biogrid/BIOGRID-ORGANISM-Arabidopsis_thaliana_Columbia-5.0.261.tab3.txt` | 43.02 MB | `c1c1cbfaee06` | Core Training Positive |
| **SRC_A3_STRING_ARA_PHYSICAL** | STRING | v12.0 (2023-06) | 2026-09-21 | *Arabidopsis thaliana* | 3702 | `data/raw/arabidopsis/string/3702.protein.physical.links.detailed.v12.0.txt.gz` | 6.79 MB | `9628efca6650` | Supplemental Positive |
| **SRC_A3_STRING_ARA_ALL_LINKS** | STRING | v12.0 (2023-06) | 2026-09-21 | *Arabidopsis thaliana* | 3702 | `data/raw/arabidopsis/string/3702.protein.links.detailed.v12.0.txt.gz` | 107.38 MB | `c5f1e5e1e0f3` | Supplemental Positive |
| **SRC_A3_STRING_ARA_ALIASES** | STRING | v12.0 (2023-06) | 2026-09-21 | *Arabidopsis thaliana* | 3702 | `data/raw/arabidopsis/string/3702.protein.aliases.v12.0.txt.gz` | 3.44 MB | `18b48b5f7dbd` | Identifier Mapping Auxiliary |
| **SRC_A4_XLMS_2026_PLINK_SEARCH** | PRIDE (Trinh et al., 2026) | PXD066234 (Aug 2026) | 2026-09-21 | *Arabidopsis thaliana* | 3702 | `data/raw/arabidopsis/xlms_2026/Total_XL_plink3-2_v3.csv` | 187.19 MB | `3fa493b4365f` | **Temporal Holdout (Quarantined)** |
| **SRC_A4_XLMS_2026_METADATA_LYSATE** | PRIDE (Trinh et al., 2026) | PXD066234 (Aug 2026) | 2026-09-21 | *Arabidopsis thaliana* | 3702 | `data/raw/arabidopsis/xlms_2026/Cell_lysateXL_File_description.csv` | 2.63 KB | `8eeecbe7b13a` | Metadata Auxiliary |
| **SRC_A4_XLMS_2026_METADATA_CHLOROPLAST** | PRIDE (Trinh et al., 2026) | PXD066291 (Aug 2026) | 2026-09-21 | *Arabidopsis thaliana* | 3702 | `data/raw/arabidopsis/xlms_2026/ChloroplastXL_File_description_v2.csv` | 1.55 KB | `10a29e6df47b` | Metadata Auxiliary |
| **SRC_B1_BIPSEQ_PMC_XML** | POPPIN Underlying Literature (BIP-seq) | Adv Sci 2025 (e2416243) | 2026-09-21 | *Oryza sativa* | 4530 | `data/raw/rice/poppin/BIP_seq_PMC11923860_fulltext.xml` | 153.11 KB | `7b7c6dad0224` | Core Training Positive (Publication) |
| **SRC_B2_INTACT_RICE_39947** | IntAct / IMEx | Release 252 (2026-01) | 2026-09-21 | *Oryza sativa Japonica* | 39947 | `data/raw/rice/intact/intact_rice_taxid39947_release252.mitab27.txt` | 1.48 MB | `80c05b9ae060` | Core Training Positive |
| **SRC_B2_INTACT_RICE_4530** | IntAct / IMEx | Release 252 (2026-01) | 2026-09-21 | *Oryza sativa* | 4530 | `data/raw/rice/intact/intact_rice_taxid4530_release252.mitab27.txt` | 759.40 KB | `89ee5b8910ee` | Core Training Positive |
| **SRC_B3_BIOGRID_RICE** | BioGRID | 5.0.261 (2026-08) | 2026-09-21 | *Oryza sativa Japonica* | 39947 | `data/raw/rice/biogrid/BIOGRID-ORGANISM-Oryza_sativa_Japonica-5.0.261.tab3.txt` | 110.44 KB | `7230167dc234` | Core Training Positive |
| **SRC_B4_STRING_RICE_PHYSICAL** | STRING | v12.0 (2023-06) | 2026-09-21 | *Oryza sativa Japonica* | 39947 | `data/raw/rice/string/39947.protein.physical.links.detailed.v12.0.txt.gz` | 8.88 MB | `3617403bc948` | Supplemental Positive |
| **SRC_B4_STRING_RICE_ALL_LINKS** | STRING | v12.0 (2023-06) | 2026-09-21 | *Oryza sativa Japonica* | 39947 | `data/raw/rice/string/39947.protein.links.detailed.v12.0.txt.gz` | 462.36 MB | `931101ec06ea` | Supplemental Positive |
| **SRC_B4_STRING_RICE_ALIASES** | STRING | v12.0 (2023-06) | 2026-09-21 | *Oryza sativa Japonica* | 39947 | `data/raw/rice/string/39947.protein.aliases.v12.0.txt.gz` | 2.25 MB | `4b979480e7ee` | Identifier Mapping Auxiliary |
| **SRC_C1_BIOGRID_MAIZE** | BioGRID | 5.0.261 (2026-08) | 2026-09-21 | *Zea mays* | 4577 | `data/raw/external_holdout/maize/BIOGRID-ORGANISM-Zea_mays-5.0.261.tab3.txt` | 5.90 KB | `12b760045851` | **External Crop Holdout (Quarantined)** |
| **SRC_C1_INTACT_MAIZE** | IntAct / IMEx | Release 252 (2026-01) | 2026-09-21 | *Zea mays* | 4577 | `data/raw/external_holdout/maize/intact_maize_taxid4577_release252.mitab27.txt` | 97.17 MB | `a2cffbe98642` | **External Crop Holdout (Quarantined)** |
| **SRC_C1_STRING_MAIZE_PHYSICAL** | STRING | v12.0 (2023-06) | 2026-09-21 | *Zea mays* | 4577 | `data/raw/external_holdout/maize/4577.protein.physical.links.detailed.v12.0.txt.gz` | 19.44 MB | `fa036a1883be` | **External Crop Holdout (Quarantined)** |
| **SRC_C2_BIOGRID_TOMATO** | BioGRID | 5.0.261 (2026-08) | 2026-09-21 | *Solanum lycopersicum* | 4081 | `data/raw/external_holdout/tomato/BIOGRID-ORGANISM-Solanum_lycopersicum-5.0.261.tab3.txt` | 37.22 KB | `851207b6cfc0` | **External Crop Holdout (Quarantined)** |
| **SRC_C2_INTACT_TOMATO** | IntAct / IMEx | Release 252 (2026-01) | 2026-09-21 | *Solanum lycopersicum* | 4081 | `data/raw/external_holdout/tomato/intact_tomato_taxid4081_release252.mitab27.txt` | 63.11 KB | `d4e565983570` | **External Crop Holdout (Quarantined)** |
| **SRC_C2_STRING_TOMATO_PHYSICAL** | STRING | v12.0 (2023-06) | 2026-09-21 | *Solanum lycopersicum* | 4081 | `data/raw/external_holdout/tomato/4081.protein.physical.links.detailed.v12.0.txt.gz` | 10.94 MB | `c9ee3e7efbb1` | **External Crop Holdout (Quarantined)** |
| **SRC_C3_BIOGRID_SOYBEAN** | BioGRID | 5.0.261 (2026-08) | 2026-09-21 | *Glycine max* | 3847 | `data/raw/external_holdout/soybean/BIOGRID-ORGANISM-Glycine_max-5.0.261.tab3.txt` | 13.43 KB | `9de9e10a5dcf` | **External Crop Holdout (Quarantined)** |
| **SRC_C3_INTACT_SOYBEAN** | IntAct / IMEx | Release 252 (2026-01) | 2026-09-21 | *Glycine max* | 3847 | `data/raw/external_holdout/soybean/intact_soybean_taxid3847_release252.mitab27.txt` | 48.72 KB | `38b4383c27e3` | **External Crop Holdout (Quarantined)** |
| **SRC_C3_STRING_SOYBEAN_PHYSICAL** | STRING | v12.0 (2023-06) | 2026-09-21 | *Glycine max* | 3847 | `data/raw/external_holdout/soybean/3847.protein.physical.links.detailed.v12.0.txt.gz` | 23.24 MB | `73be67d1d2b8` | **External Crop Holdout (Quarantined)** |
| **SRC_D1_DEEPARAPPI (7 files)** | DeepAraPPI Benchmark | 2023 Publication Benchmark | 2024-03-01 | *Arabidopsis* / *Rice* | 3702 / 39947 | `data/raw/legacy_benchmarks/deeparappi/*.txt` | 4.87 MB | *(Multiple)* | Legacy Comparison Baseline |
| **SRC_D2_ESMARAPPI (5 files)** | ESMAraPPI Benchmark | 2023 Publication Benchmark | 2024-05-01 | *Arabidopsis thaliana* | 3702 | `data/raw/legacy_benchmarks/esmarappi/*` | 2.41 MB | *(Multiple)* | Legacy Comparison Baseline |
| **SRC_D3_ARACOFUSION_PROVENANCE** | ARACoFusion Provenance Document (Sarkar & Sarkar, 2026) | 2026 bioRxiv Preprint | 2026-06-01 | *Arabidopsis thaliana* | 3702 | `data/raw/legacy_benchmarks/aracofusion/README.md` | 759 Bytes | `3db32c44d6de` | Legacy Provenance Document |

> [!IMPORTANT]
> **POPPIN Acquisition Status**:
> The POPPIN bulk interaction database is **BLOCKED / NOT AVAILABLE** for direct bulk download (the POPPIN web service operates as an interactive query portal without a programmatic API or bulk dump). Consequently, the primary literature source reporting the underlying physical rice interactome screen—the BIP-seq article (`Liu et al., Adv Sci 2025; PMC11923860`)—was acquired and frozen as `SRC_B1_BIPSEQ_PMC_XML`.

---

## 3. Raw Data Integrity & Provenance Audits

All datasets were parsed and inspected to verify valid headers, expected column delimiters, and sensible row counts:

| Dataset Category | File | Formats & Schema | Line / Record Count | Verification Status |
| :--- | :--- | :--- | :--- | :--- |
| **IntAct Arabidopsis** | `intact_arabidopsis_taxid3702_release252.mitab27.txt` | PSI-MI TAB 2.7 (42 columns, tab-delimited) | **58,740** interactions | **PASS** (Matches PSICQUIC service count exactly) |
| **BioGRID Arabidopsis** | `BIOGRID-ORGANISM-Arabidopsis_thaliana_Columbia-5.0.261.tab3.txt` | BioGRID TAB 3.0 (37 columns, tab-delimited) | **84,423** interactions | **PASS** (Physical + genetic interactions preserved) |
| **STRING Arabidopsis (Physical)** | `3702.protein.physical.links.detailed.v12.0.txt.gz` | TSV gzipped (6 columns, space-delimited) | >400,000 links | **PASS** (Raw physical channels intact) |
| **STRING Arabidopsis (All Links)** | `3702.protein.links.detailed.v12.0.txt.gz` | TSV gzipped (10 columns, space-delimited) | >10,000,000 links | **PASS** (All evidence channels preserved) |
| **PhoX XL-MS 2026 Search Table** | `Total_XL_plink3-2_v3.csv` | pLink 3.2 CSV (32 columns, comma-delimited) | **390,527** cross-link records | **PASS** (FDR <= 5% spectra preserved; raw scores intact) |
| **POPPIN BIP-seq Full Text** | `BIP_seq_PMC11923860_fulltext.xml` | JATS XML format | 2,119 lines | **PASS** (Primary interactome article preserved) |
| **IntAct Rice (39947)** | `intact_rice_taxid39947_release252.mitab27.txt` | PSI-MI TAB 2.7 (42 columns, tab-delimited) | **649** interactions | **PASS** (Matches PSICQUIC service count exactly) |
| **IntAct Rice (4530)** | `intact_rice_taxid4530_release252.mitab27.txt` | PSI-MI TAB 2.7 (42 columns, tab-delimited) | **378** interactions | **PASS** (Matches PSICQUIC service count exactly) |
| **BioGRID Rice** | `BIOGRID-ORGANISM-Oryza_sativa_Japonica-5.0.261.tab3.txt` | BioGRID TAB 3.0 (37 columns, tab-delimited) | **367** interactions | **PASS** (Matches BioGRID release archive) |
| **IntAct Maize** | `intact_maize_taxid4577_release252.mitab27.txt` | PSI-MI TAB 2.7 (42 columns, tab-delimited) | **51,009** interactions | **PASS** (Matches PSICQUIC service count exactly) |
| **BioGRID Maize** | `BIOGRID-ORGANISM-Zea_mays-5.0.261.tab3.txt` | BioGRID TAB 3.0 (37 columns, tab-delimited) | **19** interactions | **PASS** (Matches BioGRID release archive) |
| **IntAct Tomato** | `intact_tomato_taxid4081_release252.mitab27.txt` | PSI-MI TAB 2.7 (42 columns, tab-delimited) | **29** interactions | **PASS** (Matches PSICQUIC service count exactly) |
| **IntAct Soybean** | `intact_soybean_taxid3847_release252.mitab27.txt` | PSI-MI TAB 2.7 (42 columns, tab-delimited) | **21** interactions | **PASS** (Matches PSICQUIC service count exactly) |
| **DeepAraPPI Benchmark** | 7 files (`PPIN_*`, `total_*`, `c1_ppi_*`) | Tab-delimited TSV/TXT | 248,076 total lines | **PASS** (Byte-identical to local archive) |
| **ESMAraPPI Benchmark** | 5 files (`c1Train`, `c2Pred`, `c3Pred`, FASTA) | Tab-delimited TSV / FASTA | 85,070 total lines | **PASS** (Byte-identical to local archive) |

---

### 3.1 Proteomic Cross-Linking Mass Spectrometry (XL-MS) Provenance & Scope Audit

The Trinh et al. (*PNAS* August 2026) study deposited data across four ProteomeXchange / PRIDE accessions:
1. **`PXD066234` (Acquired & Frozen)**: Contains the central processed cross-link identification results generated by pLink 3.2 (`Total_XL_plink3-2_v3.csv`, 187.19 MB, 390,527 records at 5% FDR) along with cell lysate run metadata (`Cell_lysateXL_File_description.csv`). This table contains the inter-protein and intra-protein peptide-spectrum matches (PSMs), cross-linker chemistries (PhoX, DSSO), and identified protein accessions.
2. **`PXD066291` (Acquired & Frozen)**: Contains chloroplast fractionation cross-linking metadata (`ChloroplastXL_File_description_v2.csv`).
3. **`PXD073734` (Omitted from Phase 1)**: Contains mass spectrometry raw spectra (`.raw` instrumentation files) for whole-cell lysate fractionation experiments. These files exceed **500 GB** in size and consist of uninterpreted mass spectra. They were omitted because benchmark construction consumes the identified cross-link peptide tables, not raw MS instrument data.
4. **`PXD078284` (Omitted from Phase 1)**: Contains mass spectrometry raw spectra (`.raw` instrumentation files) for chloroplast cross-linking replicates (>200 GB). Omitted for the same reason.

### 3.2 Prospective Temporal Holdout Audit Nuance

> [!NOTE]
> Publication/deposition in 2026 postdates the documented PPLM paired-pretraining sources and cutoff. This establishes dataset-level temporal separation from PPLM pretraining. It does not establish pair-level novelty, because individual biological interactions may have appeared previously in STRING, PDB, IntAct, or the literature; exact pair exposure will therefore be audited separately.
>
> In Phase 3, an explicit historical overlap audit against pre-2023 STRING, IntAct, and PDB records will be conducted to partition the dataset into:
> - **Historically Validated Pairs**: Cross-linked pairs that corroborate previously reported Y2H/AP-MS or structural interactions.
> - **Entirely Novel Physical Complexes**: Inter-protein cross-links that have zero prior literature or database evidence.

---

## 4. Licensing and Attribution Audit

All acquired resources have been audited for research and academic redistribution compliance:

| Resource | License Type | Commercial Use Allowed | Attribution Requirement | Canonical Citation / Reference | DOI / Accession |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **IntAct / IMEx** | Creative Commons Attribution 4.0 International (CC BY 4.0) | Yes | del-Toro N, et al. *Nucleic Acids Res.* 2022;50(D1):D648–D653 | PMID: 34791374 | `10.1093/nar/gkab1006` |
| **BioGRID** | MIT License | Yes | Oughtred R, et al. *Protein Sci.* 2021;30(1):187–200 | PMID: 33073389 | `10.1002/pro.3978` |
| **STRING** | Creative Commons Attribution 4.0 International (CC BY 4.0) | Yes | Szklarczyk D, et al. *Nucleic Acids Res.* 2023;51(D1):D638–D646 | PMID: 36370105 | `10.1093/nar/gkac1000` |
| **PRIDE (Trinh et al., 2026)**| Creative Commons Public Domain Dedication (CC0 1.0) | Yes | Trinh CS, et al. *Proc Natl Acad Sci USA.* 2026;123(33):e2519615123 | PMID: 39133827 | `10.1073/pnas.2519615123` |
| **POPPIN Underlying Literature (BIP-seq)** | Creative Commons Attribution 4.0 International (CC BY 4.0) | Yes | Liu X, et al. *Adv Sci (Weinh).* 2025;12(11):e2416243 | PMID: 39840553 | `10.1002/advs.202416243` |
| **DeepAraPPI** | Academic Non-Commercial / Research Use | Research only | Zheng J, Yang X, Huang Y, Yang S, Wuchty S, Zhang Z. *The Plant Journal.* 2023;114(4):984–994 | PMID: 36919205 | `10.1111/tpj.16188` |
| **ESMAraPPI** | Academic Non-Commercial / Research Use | Research only | Zhou K, Lei C, Zheng J, Huang Y, Zhang Z. *Plant Methods.* 2023;19:141 | PMID: 38062445 | `10.1186/s13007-023-01119-6` |
| **ARACoFusion** | CC-BY-NC-ND 4.0 International | Research only | Sarkar D, Sarkar C. *bioRxiv.* Posted May 26, 2026 | Preprint | `10.64898/2026.05.22.727120` |

The updated BibTeX entries and metadata table are documented in [`metadata/citations.bib`](../../metadata/citations.bib) and [`metadata/licenses.tsv`](../../metadata/licenses.tsv).

---

## 5. Quality Control Gates Summary (QC1 – QC7)

| Gate ID | Quality Control Gate Description | Target Criterion | Audit Result | Status |
| :--- | :--- | :--- | :--- | :--- |
| **QC1** | Cryptographic Checksum Verification | 100% of files match SHA-256 in `checksums/SHA256SUMS` | 37 of 37 passed, 0 failed | **PASS** |
| **QC2** | Zero-Byte File Audit | 0 files with 0 bytes in `data/raw/` | 0 zero-byte files across 1,086.61 MB | **PASS** |
| **QC3** | Source Registry Metadata Completeness | Non-empty version, license, citation, taxonomy, and local path | 37 of 37 entries complete | **PASS** |
| **QC4** | Strict Holdout Quarantine & Separation | Clean partition between training candidates and holdouts; no leakage | 12 training candidates, 12 holdouts quarantined | **PASS** |
| **QC5** | Raw Data Immutability & Original Schemas | Exact original headers, column counts, and raw score channels | Original MITAB, TAB 3.0, STRING, CSV preserved | **PASS** |
| **QC6** | Legacy Benchmark Byte-Level Fidelity | Exact byte-for-byte SHA-256 match against original archives | 5 evaluated pairs 100% matched | **PASS** |
| **QC7** | Required Deliverables Completeness | All 10 deliverables present, non-empty, and cross-referenced | 10 of 10 deliverables validated | **PASS** |

---

## 6. Deliverables Checklist

- [x] **D1: Raw Datasets Directory Tree** (`data/raw/`): Structured by species and source (37 files, 1.09 GB).
- [x] **D2: Source Registry TSV** (`metadata/source_registry.tsv`): Tabular database registry with 32 metadata attributes per source file, including explicit `release_version` and `archive_snapshot_date`.
- [x] **D3: Source Registry JSON** (`metadata/source_registry.json`): Structured programmatic registry.
- [x] **D4: Cryptographic Checksum Registry** (`checksums/SHA256SUMS`): Linux/Windows-compatible SHA-256 checksum manifest.
- [x] **D5: Acquisition & Registry Script** (`scripts/phase1_download_sources.py`): Fully reproducible download and extraction orchestrator.
- [x] **D6: Acquisition Event Log** (`logs/phase1_acquisition.log`): Timestamped network transaction and verification log.
- [x] **D7: Licenses Registry** (`metadata/licenses.tsv`): Detailed terms-of-use and redistribution compatibility table.
- [x] **D8: Citations Registries** (`metadata/citations.tsv` and `metadata/citations.bib`): Corrected canonical citations (Zheng et al., 2023 for DeepAraPPI; Zhou et al., 2023 for ESMAraPPI; Sarkar & Sarkar, 2026 for ARACoFusion).
- [x] **D9: Phase 1 Machine Manifest** (`manifests/phase1_manifest.json`): Machine-readable execution manifest.
- [x] **D10: Quarantine Notice** (`manifests/DO_NOT_TRAIN_ON_THESE_SOURCES.txt`): Explicit human- and machine-readable holdout isolation rules.

---

## 7. Sign-off and Readiness for Phase 2

All success criteria for Phase 1 have been satisfied, and all provenance and metadata items have been aligned with canonical publications and archive records. No raw data was modified.

**Readiness for Phase 2 (Parser & Normalizer Implementation)**:
- Raw inputs for parser development:
  - MITAB 2.7 parsers: `intact_arabidopsis_taxid3702_release252.mitab27.txt`, `intact_maize_taxid4577_release252.mitab27.txt`
  - BioGRID TAB 3.0 parsers: `BIOGRID-ORGANISM-Arabidopsis_thaliana_Columbia-5.0.261.tab3.txt`
  - STRING link parsers: `3702.protein.physical.links.detailed.v12.0.txt.gz`, `39947.protein.links.detailed.v12.0.txt.gz`
  - Cross-link CSV parsers: `Total_XL_plink3-2_v3.csv`
- Primary reference mapping inputs:
  - UniProtKB reference proteome mappings and gene alias files will be initialized in Phase 2.
