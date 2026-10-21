# AraPPLM Data Directory Architecture

This directory houses all data assets for the AraPPLM project, organized to maintain a strict separation between raw downloads, verified positive interaction pairs, legacy evaluation benchmarks, and sequence preparation databases.

---

## Directory Structure Overview

```
data/
├── dataset_construction/          # [DATASET CONSTRUCTION] Progressive pipeline artifacts starting from verified positives
│   ├── README.md                  # Dataset construction pipeline data documentation
│   └── 01_positive_pairs/         # [STAGE 1: STARTING POINT] Verified positive interaction datasets
│       ├── arabidopsis/           # Arabidopsis thaliana (87,673 deduplicated pairs; 100% pure PPIs)
│       ├── maize/                 # Zea mays (50,984 deduplicated pairs; 50,983 pure PPIs)
│       ├── rice/                  # Oryza sativa (667 pure PPI pairs; 998 with nucleic acid complexes)
│       ├── soybean/               # Glycine max (48 deduplicated pairs)
│       └── tomato/                # Solanum lycopersicum (123 deduplicated pairs)
│
├── raw/                           # [IMMUTABLE] Raw, unmodified downloads from public databases
│   ├── arabidopsis/               # BioGRID v5.0.261, IntAct r252, STRING v12.0, PRIDE PhoX XL-MS
│   ├── rice/                      # BioGRID, IntAct (taxid 39947 + 4530), STRING, PoPPIN (PMC XML)
│   ├── external_holdout/          # Maize, Tomato, and Soybean raw downloads (BioGRID, IntAct, STRING)
│   └── legacy_benchmarks/         # Raw source files from DeepAraPPI, ESMAraPPI, and AraCoFusion
│
├── DeepAraPPI/                    # [LEGACY BENCHMARK] DeepAraPPI C1 (train), C2 (test), C3 (test), Rice test
├── ESMAraPPI/                     # [LEGACY BENCHMARK] ESMAraPPI C1 (train), C2 (test), C3 (test)
│
├── arabidopsis/                   # [DATA PREP] Arabidopsis sequence DB mapping & missing UniParc FASTA
└── rice/                          # [DATA PREP] Rice sequence DB mapping & missing UniParc FASTA
```

---

## 1. Verified Positive Pairs (`data/dataset_construction/01_positive_pairs/`)

This directory is the **authoritative starting point** for the fresh dataset construction pipeline. All records have been independently verified against the raw database dumps.

Each species subdirectory contains:
- `<species>_positive_evidence.tsv` & `.parquet`: Full-evidence table capturing provenance across all retained interactions.
- `<species>_positive_pairs_deduplicated.tsv` & `.parquet`: Master deduplicated canonical positive pairs (`participant_a <= participant_b`).
- `<species>_positive_audit.json`: Provenance metadata, exact filter counts, topology breakdown (heteromers vs. homomers), and literature metrics.

### Verified Dataset Summary

| Species | Taxonomy ID | Master Unique Pairs | Pure Protein-Protein Pairs | Heteromers | Homomers | Unique Proteins | Sources Included |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Arabidopsis** | `3702` | **87,673** | **87,673** | 81,900 | 5,773 | 18,437 | BioGRID, IntAct, STRING Full, PhoX XL-MS |
| **Maize** | `4577` | **50,984** | **50,983** | 50,409 | 575 | 14,888 | IntAct, BioGRID |
| **Rice** | `39947` / `4530` | **998** | **667** | 651 | 16 | 609 *(892 all)* | BioGRID, IntAct, STRING Full |
| **Tomato** | `4081` | **123** | **123** | 108 | 15 | 61 | BioGRID, IntAct |
| **Soybean** | `3847` | **48** | **48** | 47 | 1 | 51 | BioGRID, IntAct |

---

## 2. Raw Database Dumps (`data/raw/`)

Contains immutable raw archives downloaded directly from primary repositories:
- **BioGRID**: Release 5.0.261 TAB 3.0 organism files.
- **IntAct**: Release 252 PSI-MI TAB 2.7 organism dumps.
- **STRING**: Version 12.0 physical links (`full` and `detailed`).
- **PRIDE PhoX XL-MS**: High-density cross-linking mass spectrometry CSV (`Total_XL_plink3-2_v3.csv`, Trinh et al., 2026, Nature Communications).
- **PoPPIN**: PMC JATS XML full-text protocol article.
- **Legacy Benchmarks**: Raw benchmark sample files from prior publications.

---

## 3. Legacy Benchmark Datasets (`data/DeepAraPPI/` & `data/ESMAraPPI/`)

Preserved exclusively for model benchmarking and comparative evaluation:
- **DeepAraPPI**: C1 train, C2 test, C3 test, total samples, and independent rice holdout.
- **ESMAraPPI**: C1 train, C2 pred, C3 pred, and missing ID FASTA sequences.

---

## 4. Sequence Database Assets (`data/arabidopsis/` & `data/rice/`)

Contains sequence mapping files and supplementary FASTA files to ensure 100% sequence coverage for all interactors during PPLM embedding extraction.
