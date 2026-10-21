# Independent Positive Pair Verification: Rice, Tomato, and Soybean

**Verification Date**: 2026-10-16  
**Auditor**: Independent Data Audit Pipeline  
**Execution Mode**: Pure raw data parsing (zero reliance on historical documentation)  
**Target Taxonomies**:
- *Oryza sativa* (Rice, NCBI Taxonomy IDs: `39947` [Japonica Group] and `4530` [Species-Level])
- *Solanum lycopersicum* (Tomato, NCBI Taxonomy ID: `4081`)
- *Glycine max* (Soybean, NCBI Taxonomy ID: `3847`)

---

## 1. Executive Summary & Cross-Species Verification Dashboard

This document provides the definitive verification of positive protein-protein interaction (PPI) data for Rice, Tomato, and Soybean computed directly from the downloaded raw database dumps. All numbers are calculated by stream-parsing the raw files on disk.

### Master Cross-Species Verification Table

| Species (Common / Scientific) | Database / Source | Raw Downloaded File Path | File Size (Bytes) | 1. Raw Downloaded Records | 1b. Raw Canonical Pairs | 2a. Filter: Physical PPI Criteria | Records Passing Physical Filter | 2b. Filter: Protein-Only Criteria | Records Passing Protein Filter | 3. Unique Deduplicated Pairs (Database-Level) | Heteromeric Pairs | Homomeric Pairs | Unique Participating Proteins |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Rice** (*Oryza sativa*) | **BioGRID** (v5.0.261) | `data/raw/rice/biogrid/BIOGRID-ORGANISM-Oryza_sativa_Japonica-5.0.261.tab3.txt` | 110,438 | **366** | 343 | `sys_type == 'physical'` & exclude `Proximity Label-MS` (10), `PCA` (2), & exclude genetic (2) | 352 | All 352 represent protein products | 352 | **331** | 315 | 16 | 454 |
| | **IntAct** (r252: 39947 + 4530) | `data/raw/rice/intact/intact_rice_taxid39947_release252.mitab27.txt` + `...taxid4530...` | 2,240,301 *(sum)* | **1,027** *(649 + 378)* | 595 | Retain `MI:0915` (365) & `MI:0914` (661); exclude ubiquitination `MI:0220` (1) | 1,026 | Retain `MI:0326` (protein); exclude RNA (`MI:0320`: 507) & DNA (`MI:0319`: 154) | **365** *(1,026 with nucleic acids)* | **263** *(144 UniProt canonical; 594 with nucleic acids)* | 247 *(128)* | 16 *(16)* | 352 *(198)* |
| | **STRING** (v12.0 Physical Full) | `data/raw/rice/string/39947.protein.physical.links.full.v12.0.txt.gz` | 9,277,044 | **1,294,058** *(directed rows)* | 647,029 | Direct `experiments > 0`; strictly exclude transferred interologs (1,074,164) & non-exp (219,032) | 862 *(directed)* | 100% protein accessions | 862 *(directed)* | **431** *(undirected canonical pairs)* | 431 | 0 | 459 |
| | **PoPPIN** (BIP-seq) | `data/raw/rice/poppin/BIP_seq_PMC11923860_fulltext.xml` | 153,108 | 0 *(XML paper)* | 0 | Methodological paper; bulk interactome web portal only | 0 | 0 | 0 | **0** | 0 | 0 | 0 |
| | **Rice Combined Master** | — | — | **1,295,451** | — | **Direct In-Species Physical Filter** | **2,240** | **Pure Protein-Protein Strict Exclusion** | **1,579** *(2,241 in warehouse)* | **667** *(812 raw ID; 998 with nucleic acids)* | **651** *(796)* | **16** *(16)* | **609** *(892 with nucleic acids)* |
| **Tomato** (*Solanum lycopersicum*) | **BioGRID** (v5.0.261) | `data/raw/external_holdout/tomato/BIOGRID-ORGANISM-Solanum_lycopersicum-5.0.261.tab3.txt` | 37,222 | **141** | 109 | `sys_type == 'physical'` (137); exclude genetic `Dosage Rescue` (4) & 0 proximity | 137 | All 137 represent protein products | 137 | **107** *(105 canonical)* | 92 | 15 | 58 |
| | **IntAct** (Release 252) | `data/raw/external_holdout/tomato/intact_tomato_taxid4081_release252.mitab27.txt` | 63,108 | **29** | 18 | Physical association (`MI:0915`: 27), direct interaction (`MI:0407`: 2); 0 excluded | 29 | 100% `MI:0326` (protein) | 29 | **18** *(16 canonical)* | 18 | 0 | 17 |
| | **STRING** (v12.0 Physical Detailed) | `data/raw/external_holdout/tomato/4081.protein.physical.links.detailed.v12.0.txt.gz` | 10,942,272 | **1,710,092** *(directed rows)* | 855,046 | Physical links; `experimental > 0` (1,501,486 directed); exclude `experimental == 0` (208,606 directed) | 1,501,486 *(directed)* | 100% protein accessions | 1,501,486 *(directed)* | **750,743** *(undirected canonical pairs)* | 750,743 | 0 | 15,829 |
| | **Tomato Curated (IntAct + BioGRID)** | — | — | **170** | — | **Direct Curated Physical Filter** | **166** | **Pure Protein-Protein Strict Exclusion** | **166** | **123** *(125 raw ID)* | **108** | **15** | **61** |
| **Soybean** (*Glycine max*) | **BioGRID** (v5.0.261) | `data/raw/external_holdout/soybean/BIOGRID-ORGANISM-Glycine_max-5.0.261.tab3.txt` | 13,433 | **46** | 40 | `sys_type == 'physical'` (46); 0 genetic, 0 proximity | 46 | All 46 represent protein products | 46 | **40** | 39 | 1 | 42 |
| | **IntAct** (Release 252) | `data/raw/external_holdout/soybean/intact_soybean_taxid3847_release252.mitab27.txt` | 48,719 | **21** | 13 | Physical association (`MI:0915`: 13), direct interaction (`MI:0407`: 7); exclude colocalization (`MI:0403`: 1) | 20 | 100% `MI:0326` (protein) | 20 | **13** | 13 | 0 | 18 |
| | **STRING** (v12.0 Physical Detailed) | `data/raw/external_holdout/soybean/3847.protein.physical.links.detailed.v12.0.txt.gz` | 23,240,128 | **3,782,422** *(directed rows)* | 1,891,211 | Physical links; `experimental > 0` (3,267,484 directed); exclude `experimental == 0` (514,938 directed) | 3,267,484 *(directed)* | 100% protein accessions | 3,267,484 *(directed)* | **1,633,742** *(undirected canonical pairs)* | 1,633,742 | 0 | 25,540 |
| | **Soybean Curated (IntAct + BioGRID)** | — | — | **67** | — | **Direct Curated Physical Filter** | **66** | **Pure Protein-Protein Strict Exclusion** | **66** | **48** | **47** | **1** | **51** |

---

## 2. Species 1: Rice (*Oryza sativa*)

### 2.1 Raw Downloaded Data
- **BioGRID**: `data/raw/rice/biogrid/BIOGRID-ORGANISM-Oryza_sativa_Japonica-5.0.261.tab3.txt` (110,438 bytes, **366 records**, 343 raw pairs).
- **IntAct (TaxID 39947)**: `data/raw/rice/intact/intact_rice_taxid39947_release252.mitab27.txt` (1,480,902 bytes, **649 records**).
- **IntAct (TaxID 4530)**: `data/raw/rice/intact/intact_rice_taxid4530_release252.mitab27.txt` (759,399 bytes, **378 records**).
  - *Combined IntAct Rice*: **1,027 records**, 595 raw pairs.
- **STRING (Physical Links Full)**: `data/raw/rice/string/39947.protein.physical.links.full.v12.0.txt.gz` (9,277,044 bytes, **1,294,058 directed rows**, 647,029 canonical pairs).
- **PoPPIN (BIP-seq)**: `data/raw/rice/poppin/BIP_seq_PMC11923860_fulltext.xml` (153,108 bytes, **0 bulk pairs**). Audit confirmed that the full text PMC XML contains protocol descriptions; the bulk interactome is queryable solely via the web portal.

### 2.2 Physical PPI Filtering (Rice)
1. **BioGRID Filtering**:
   - `sys_type == 'physical'`: 364 records (99.45%)
   - `sys_type == 'genetic'` (`Synthetic Rescue`): **2 records** — **EXCLUDED**.
   - Proximity assays: `Proximity Label-MS` (**10 records**) and `PCA` (**2 records**) — **EXCLUDED**.
   - **Records Passing Physical Filter**: 366 - 2 - 12 = **352 records**.
2. **IntAct Filtering**:
   - `psi-mi:"MI:0915"(physical association)`: 365 records (328 in 39947, 37 in 4530)
   - `psi-mi:"MI:0914"(association)`: 661 records (320 in 39947, 341 in 4530)
   - `psi-mi:"MI:0220"(ubiquitination reaction)`: **1 record** (Line 338 in 39947) — **EXCLUDED** (enzymatic reaction).
   - **Records Passing Physical Filter**: 1,027 - 1 = **1,026 records**.
3. **STRING Filtering**:
   - Direct `experiments > 0`: **862 directed records** (**431 unique canonical pairs**).
   - `experiments == 0 and experiments_transferred > 0`: **1,074,164 directed records** — **EXCLUDED** (transferred interologs).
   - `experiments == 0 and experiments_transferred == 0`: **219,032 directed records** — **EXCLUDED** (textmining/database only).

### 2.3 Protein-Only Complex Filtering (Critical Verification Finding for Rice)
Analysis of Columns 20 and 21 in Rice IntAct reveals that over 64% of the downloaded records are **nucleic acid (RNA/DNA) binding interactions**, NOT protein-protein interactions:

| Participant Entity Type | TaxID 39947 (Interactor A) | TaxID 4530 (Interactor A) | Combined Interactor A | Interactor B (Both TaxIDs) | Classification | Action |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `psi-mi:"MI:0326"(protein)` | 329 | 37 | **366** | 1,027 | Pure Protein | **Retained** |
| `psi-mi:"MI:0320"(ribonucleic acid)` / RNA | 245 | 262 | **507** | 0 | RNA-Protein Complex | **Excluded** from pure PPI |
| `psi-mi:"MI:0319"(deoxyribonucleic acid)` / DNA | 75 | 79 | **154** | 0 | DNA-Protein Complex | **Excluded** from pure PPI |

- **Exact Impact**:
  - Out of 1,026 physical IntAct records, only **365 physical records** (328 in 39947, 37 in 4530) are pure protein-protein interactions.
  - Exactly **661 records** represent nucleic acid binding complexes (507 RNA binding, 154 DNA binding).
- BioGRID and STRING Rice records are 100% protein-protein.

### 2.4 Deduplication (Rice)
- **IntAct Deduplication**:
  - Pure protein-protein unique pairs: **263 pairs** (raw IDs) / **144 pairs** (UniProt canonical accessions).
  - All physical interactions (including RNA/DNA): **594 pairs** (raw IDs) / **468 pairs** (UniProt canonical).
- **BioGRID Deduplication**: **331 unique pairs** (raw IDs) / **318 pairs** (UniProt canonical).
- **STRING Deduplication**: **431 unique canonical pairs**.
- **Combined Master Rice Interactome**:
  - **All Physical Molecules (Protein + Nucleic Acids)**: **998 unique pairs** (468 IntAct-only, 220 BioGRID+STRING, 124 IntAct+STRING, 98 BioGRID-only, 85 STRING-only, 3 multi-overlap).
  - **Pure Protein-Protein Complexes Only**: Exactly **667 unique pairs** (**651 heteromers, 16 homomers**) spanning **609 unique proteins**.

---

## 3. Species 2: Tomato (*Solanum lycopersicum*, TaxID 4081)

### 3.1 Raw Downloaded Data
- **BioGRID**: `data/raw/external_holdout/tomato/BIOGRID-ORGANISM-Solanum_lycopersicum-5.0.261.tab3.txt` (37,222 bytes, **141 records**, 109 raw pairs).
- **IntAct**: `data/raw/external_holdout/tomato/intact_tomato_taxid4081_release252.mitab27.txt` (63,108 bytes, **29 records**, 18 raw pairs).
- **STRING**: `data/raw/external_holdout/tomato/4081.protein.physical.links.detailed.v12.0.txt.gz` (10,942,272 bytes, **1,710,092 directed rows**, 855,046 canonical pairs).

### 3.2 Physical PPI Filtering (Tomato)
1. **BioGRID Filtering**:
   - `sys_type == 'physical'`: **137 records** (97.16%)
     - Assays: `Two-hybrid` (118), `Reconstituted Complex` (13), `Affinity Capture-Western` (3), `FRET` (3).
     - Proximity assays: 0 PCA, 0 Proximity Label-MS.
   - `sys_type == 'genetic'` (`Dosage Rescue`): **4 records** (2.84%) — **EXCLUDED**.
   - **Records Passing Physical Filter**: 141 - 4 = **137 records**.
2. **IntAct Filtering**:
   - `psi-mi:"MI:0915"(physical association)`: **27 records**.
   - `psi-mi:"MI:0407"(direct interaction)`: **2 records**.
   - Exclusions: 0 colocalization, 0 proximity, 0 negative.
   - **Records Passing Physical Filter**: All **29 records**.
3. **STRING Filtering**:
   - `experimental > 0`: **1,501,486 directed records** (**750,743 canonical pairs**).
   - `experimental == 0`: **208,606 directed records** (**104,303 canonical pairs**) — non-experimental.

### 3.3 Protein-Only Complex Filtering (Tomato)
- All 137 retained BioGRID records and all 29 IntAct records are confirmed protein molecules (`MI:0326`). 0 DNA, 0 RNA, 0 small molecules.

### 3.4 Deduplication (Tomato Curated Benchmark)
- **BioGRID Unique Pairs**: 107 pairs (raw IDs) / 105 pairs (canonical gene symbols).
- **IntAct Unique Pairs**: 18 pairs (raw IDs) / 16 pairs (canonical gene symbols).
- **Combined Curated Benchmark (IntAct + BioGRID)**:
  - Total Evidence Records: 137 + 29 = **166 records**.
  - Total Unique Deduplicated Pairs: **123 unique pairs**.
  - Pair Topology: **108 heteromeric pairs** (87.80%), **15 homomeric pairs** (12.20%).
  - Source Overlap: BioGRID only = 105, IntAct only = 16, Both BioGRID + IntAct = 2.
  - Unique Participating Proteins: **61 proteins**.

---

## 4. Species 3: Soybean (*Glycine max*, TaxID 3847)

### 4.1 Raw Downloaded Data
- **BioGRID**: `data/raw/external_holdout/soybean/BIOGRID-ORGANISM-Glycine_max-5.0.261.tab3.txt` (13,433 bytes, **46 records**, 40 raw pairs).
- **IntAct**: `data/raw/external_holdout/soybean/intact_soybean_taxid3847_release252.mitab27.txt` (48,719 bytes, **21 records**, 13 raw pairs).
- **STRING**: `data/raw/external_holdout/soybean/3847.protein.physical.links.detailed.v12.0.txt.gz` (23,240,128 bytes, **3,782,422 directed rows**, 1,891,211 canonical pairs).

### 4.2 Physical PPI Filtering (Soybean)
1. **BioGRID Filtering**:
   - `sys_type == 'physical'`: All **46 records** (100.0%).
     - Assays: `Two-hybrid` (30), `Affinity Capture-Western` (8), `FRET` (5), `Co-fractionation` (1), `Biochemical Activity` (1), `Co-crystal Structure` (1).
     - Proximity assays: 0 PCA, 0 Proximity Label-MS.
   - `sys_type == 'genetic'`: 0 records.
   - **Records Passing Physical Filter**: All **46 records**.
2. **IntAct Filtering**:
   - `psi-mi:"MI:0915"(physical association)`: **13 records**.
   - `psi-mi:"MI:0407"(direct interaction)`: **7 records**.
   - `psi-mi:"MI:0403"(colocalization)`: **1 record** — **EXCLUDED**.
     - Line 19: Interactor A = `uniprotkb:O64471`, Interactor B = `uniprotkb:Q07185`.
     - Interaction Type: `colocalization`. Citation: `pubmed:12604780`. Excluded as spatial co-localization rather than physical complex.
   - **Records Passing Physical Filter**: 21 - 1 = **20 records**.
3. **STRING Filtering**:
   - `experimental > 0`: **3,267,484 directed records** (**1,633,742 canonical pairs**).
   - `experimental == 0`: **514,938 directed records** (**257,469 canonical pairs**) — non-experimental.

### 4.3 Protein-Only Complex Filtering (Soybean)
- All 46 retained BioGRID records and all 20 retained IntAct records are verified protein interactors (`MI:0326`). 0 DNA, 0 RNA, 0 small molecules.

### 4.4 Deduplication (Soybean Curated Benchmark)
- **BioGRID Unique Pairs**: 40 unique pairs.
- **IntAct Unique Pairs**: 13 unique pairs.
- **Combined Curated Benchmark (IntAct + BioGRID)**:
  - Total Evidence Records: 46 + 20 = **66 records**.
  - Total Unique Deduplicated Pairs: Exactly **48 unique pairs**.
  - Pair Topology: **47 heteromeric pairs** (97.92%), **1 homomeric pair** (2.08%).
  - Source Overlap: BioGRID only = 35, IntAct only = 8, Both BioGRID + IntAct = 5.
  - Unique Participating Proteins: **51 proteins**.

---

## 5. Comparative Cross-Species Summary

| Metric | Rice (*O. sativa*) | Tomato (*S. lycopersicum*) | Soybean (*G. max*) | Total External Holdout |
| :--- | :--- | :--- | :--- | :--- |
| **Raw Curated Records Downloaded** | 1,393 *(1,027 IntAct + 366 BioGRID)* | 170 *(29 IntAct + 141 BioGRID)* | 67 *(21 IntAct + 46 BioGRID)* | **1,630 records** |
| **Curated Physical Records Retained** | 1,378 *(1,026 IntAct + 352 BioGRID)* | 166 *(29 IntAct + 137 BioGRID)* | 66 *(20 IntAct + 46 BioGRID)* | **1,610 records** |
| **Curated Pure Protein-Protein Records** | **717** *(365 IntAct + 352 BioGRID)* | **166** *(29 IntAct + 137 BioGRID)* | **66** *(20 IntAct + 46 BioGRID)* | **949 records** |
| **Nucleic Acid Interactions Excluded** | **661** *(507 RNA + 154 DNA)* | **0** | **0** | **661 records** |
| **Proximity / Genetic Excluded** | **14** *(12 proximity + 2 genetic)* | **4** *(4 genetic)* | **1** *(1 colocalization)* | **19 records** |
| **STRING Direct Experimental Canonical Pairs** | **431 pairs** *(862 directed)* | **750,743 pairs** *(detailed dump)* | **1,633,742 pairs** *(detailed dump)* | — |
| **Master Curated Deduplicated Pairs** | **667 pure PPIs** *(998 with nucleic acids)* | **123 pure PPIs** | **48 pure PPIs** | **838 pure PPIs** |
| **Heteromers / Homomers** | 651 / 16 | 108 / 15 | 47 / 1 | 806 / 32 |
| **Unique Participating Proteins** | **609** *(892 with nucleic acids)* | **61** | **51** | **721 proteins** |

---

## 6. Verification Summary & Audit Conclusion

1. **Independent Verification**: All figures were computed by stream-parsing the physical source dumps without importing historical counts.
2. **Key Discovery on Rice Data**: Over 64% of Rice IntAct interactions (661 records) represent nucleic acid (RNA/DNA) binding events rather than protein-protein interactions. Filtering for pure protein-protein complexes is essential to avoid contaminating PPI benchmarks with protein-nucleic acid complexes.
3. **Exact Holdout Counts**:
   - Rice yields **667 pure protein-protein positive pairs** (or 998 when including nucleic acid binding).
   - Tomato yields **123 pure protein-protein positive pairs**.
   - Soybean yields **48 pure protein-protein positive pairs**.
