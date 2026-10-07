# External Species Interactome Curation & Physical PPI Breakdown

## Executive Summary

This report collates all interaction records available for external plant species (i.e., non-*Arabidopsis* and non-*Rice*) preserved in the AraPPLM benchmark warehouse under `data/raw/external_holdout/`. These external datasets serve as zero-shot generalization evaluation holdouts.

Three external crop species are available:
1. **Maize (*Zea mays*)** — NCBI Taxonomy ID `4577`
2. **Tomato (*Solanum lycopersicum*)** — NCBI Taxonomy ID `4081`
3. **Soybean (*Glycine max*)** — NCBI Taxonomy ID `3847`

All raw interaction records have been ingested, hashed, and audited with zero transformations. In addition, physical PPI filtering was applied according to the standardized AraPPLM interaction semantic rules (distinguishing direct binary biophysical and physical co-complex associations from non-physical genetic interactions and spatial proximity/colocalization).

---

## 1. Total Records Collation by Species and Data Resource

Across all three external species, a grand total of **8,568,960 raw records** are present across IntAct, BioGRID, and STRING Physical Links v12.0:

| Species (Common / Scientific) | NCBI TaxID | IntAct (Release 252) | BioGRID (5.0.261) | STRING Physical Links v12.0 | Total Records Available |
| :--- | :---: | --: | --: | --: | --: |
| **Maize** (*Zea mays*) | `4577` | 51,009 | 18 | 3,025,182 | **3,076,209** |
| **Tomato** (*Solanum lycopersicum*) | `4081` | 29 | 141 | 1,710,092 | **1,710,262** |
| **Soybean** (*Glycine max*) | `3847` | 21 | 46 | 3,782,422 | **3,782,489** |
| **TOTAL** | — | **51,059** | **205** | **8,517,696** | **8,568,960** |

> [!NOTE]
> In addition to the primary species files, **65 cross-species/heterologous interactions** involving these three taxa (54 in BioGRID, 11 in IntAct) were quarantined during Arabidopsis ingest to strictly prevent data leakage into training pools.

---

## 2. Physical PPI Filtering Methodology

Interaction records were classified and filtered according to PSI-MI controlled vocabularies and BioGRID experimental systems:

1. **Retained as Physical PPIs**:
   - `direct_binary`: Confirmed pairwise biophysical interaction (e.g., Yeast Two-Hybrid MI:0018, BiFC, split-ubiquitin, biophysical binding MI:0407 / MI:0107).
   - `physical_association`: Confirmed presence within shared complexes (e.g., AP-MS, co-immunoprecipitation MI:0006/MI:0007, biochemical pull-down MI:0096, kinase phosphorylation MI:0217).
2. **Excluded as Non-Physical**:
   - `genetic`: Phenotypic enhancement, synthetic lethality, or dosage rescue (BioGRID genetic systems).
   - `proximity`: Proximity labeling or colocalization assays without proof of physical binding (e.g., fluorescence microscopy colocalization MI:0403, PCA).
3. **STRING Physical Links Stratification**:
   - STRING v12.0 provides dedicated physical interaction link files (`protein.physical.links.detailed.v12.0.txt.gz`).
   - In STRING, links are bidirectional/symmetric (ratio 2.0 lines per unique pair).
   - We report both:
     - **All STRING Physical Links** (entire physical interactome, including orthology transfers and text mining).
     - **Experimental Physical Links** (`experimental > 0`, retaining interactions with direct experimental evidence).

---

## 3. Physical PPI Breakdown by Species

### 3.1 Curated Experimental Physical PPIs (IntAct & BioGRID)

Filtering out genetic and proximity assays yields **51,257 curated experimental physical records** (corresponding to **51,098 unique physical protein pairs** across IntAct and BioGRID):

| Species | Source | Total Raw | Excluded Records | Excluded Reasons | Physical Records | Unique Physical Pairs (Hetero / Homo) | Unique Proteins |
| :--- | :--- | --: | --: | :--- | --: | :---: | :---: |
| **Maize** | IntAct | 51,009 | 1 | Colocalization (MI:0403) | 51,008 | 50,972 (50,399 / 573) | 14,871 |
| | BioGRID | 18 | 1 | Proximity PCA | 17 | 12 (10 / 2) | 20 |
| | **Maize Subtotal** | **51,027** | **2** | — | **51,025** | **50,984** | **14,875** |
| **Tomato** | IntAct | 29 | 0 | None | 29 | 18 (16 / 2) | 21 |
| | BioGRID | 141 | 4 | Genetic (Dosage Rescue) | 137 | 107 (94 / 13) | 44 |
| | **Tomato Subtotal** | **170** | **4** | — | **166** | **125** | **58** |
| **Soybean** | IntAct | 21 | 1 | Colocalization (MI:0403) | 20 | 13 (13 / 0) | 12 |
| | BioGRID | 46 | 0 | None | 46 | 40 (39 / 1) | 45 |
| | **Soybean Subtotal** | **67** | **1** | — | **66** | **53** | **53** |
| **TOTAL** | **Curated Experimental** | **51,264** | **7** | **4 Genetic, 3 Proximity** | **51,257** | **51,162** | — |

*Key Findings on Experimental Datasets*:
- **Maize** possesses a massive curated physical dataset (**51,025 physical records**), largely derived from the high-throughput Y2H maize interactome (Han et al., PMID 36581701: 50,966 interactions).
- **Tomato** contains **166 experimental physical records** (137 BioGRID + 29 IntAct) covering 125 unique pairs across 58 proteins.
- **Soybean** contains **66 experimental physical records** (46 BioGRID + 20 IntAct) covering 53 unique pairs across 53 proteins.

---

### 3.2 STRING Physical Links Breakdown

Each STRING file contains pairwise undirected interactions duplicated in both directions ($A \to B$ and $B \to A$):

| Species | Total Physical Records | Unique Canonical Pairs | Experimental > 0 Records | Experimental > 0 Canonical Pairs | Computational / Transferred Only |
| :--- | --: | --: | --: | --: | --: |
| **Maize** (*Zea mays*) | 3,025,182 | 1,512,591 | 2,569,178 | 1,284,589 | 456,004 |
| **Tomato** (*Solanum lycopersicum*) | 1,710,092 | 855,046 | 1,501,486 | 750,743 | 208,606 |
| **Soybean** (*Glycine max*) | 3,782,422 | 1,891,211 | 3,267,484 | 1,633,742 | 514,938 |
| **TOTAL** | **8,517,696** | **4,258,848** | **7,338,148** | **3,669,074** | **1,179,548** |

---

### 3.3 Comprehensive Physical PPI Summary Matrix

Depending on the benchmark evaluation criteria, the physical PPI counts per species are:

| Metric Level | Maize (*Zea mays*) | Tomato (*Solanum lycopersicum*) | Soybean (*Glycine max*) | Total |
| :--- | --: | --: | --: | --: |
| **1. Curated Experimental Physical PPIs** (IntAct + BioGRID) | **51,025** | **166** | **66** | **51,257** |
| — Unique Experimental Physical Pairs | 50,984 | 125 | 53 | 51,162 |
| **2. STRING Direct Experimental Physical Links** (`exp > 0`) | **2,569,178** | **1,501,486** | **3,267,484** | **7,338,148** |
| — Unique STRING Experimental Pairs | 1,284,589 | 750,743 | 1,633,742 | 3,669,074 |
| **3. Combined Experimental Physical Records** (Levels 1 + 2) | **2,620,203** | **1,501,652** | **3,267,550** | **7,389,405** |
| **4. All Available Physical PPI Records** (Curated + All STRING Physical) | **3,076,207** | **1,710,258** | **3,782,488** | **8,568,953** |

*Note: Across all sources, exactly 7 records are non-physical (4 genetic, 3 proximity/colocalization).*

---

## 4. Audit Artefact Coordinates

- **Audit JSON**: `data/interim/phase2/audits/external_source_breakdown.json`
- **Analysis Script**: `scripts/analyze_external_breakdown.py`
- **Holdout Quarantine Manifest**: `manifests/DO_NOT_TRAIN_ON_THESE_SOURCES.txt`
