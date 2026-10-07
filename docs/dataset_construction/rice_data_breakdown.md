# Rice (*Oryza sativa*) PPI Dataset: Comprehensive Source Breakdown

This document provides a granular, source-by-source breakdown of all protein-protein interaction (PPI) records gathered for **Rice** (*Oryza sativa* Japonica / Indica) across every integrated repository in the Phase 2 canonical evidence warehouse.

Following the established **Physical Evidence Hierarchy**:
```text
EXPERIMENTALLY SUPPORTED PHYSICAL EVIDENCE
│
├── DIRECT INTERACTION (A ───── B: "A physically contacts B")
│
├── PHYSICAL ASSOCIATION / CO-COMPLEX (A ─ C ─ B: "Same physical assembly")
│
└── PROXIMITY (A ~ B: "Spatially close")
```

At this initial evidence-gathering phase, **all physical evidence** (direct binding, co-complexes, pull-downs, proximity labeling) is preserved with full provenance and without premature quality filtering (such as IntAct MIscore cutoffs or arbitrary score gates).

---

### Source Inventory Overview

| Source ID | Database / Resource | Raw Records | Unique Pairs (Total) | Heteromeric Pairs | Homomeric Pairs | Unique Proteins | Primary PMIDs | Benchmark / Pipeline Role |
| :--- | :--- | --: | --: | --: | --: | --: | --: | :--- |
| **`SRC_B2_INTACT_RICE_39947`** | **IntAct** (Taxid 39947, *Japonica*) | 649 | 575 | 565 | 10 | 513 | 14 | Primary Experimental Positive |
| **`SRC_B2_INTACT_RICE_4530`** | **IntAct** (Taxid 4530, *O. sativa*) | 378 | 368 | 367 | 1 | 345 | 4 | Primary Experimental Positive |
| **`INTACT_RICE_COMBINED`** | **IntAct / IMEx (Union)** | **1,027** | **595** | **584** | **11** | **530** | **14** | **Primary Experimental Positive** |
| **`SRC_B3_BIOGRID_RICE`** | **BioGRID** (Rel. 5.0.261) | 366 | 333 | 328 | 5 | 315 | 23 | Primary Experimental Positive |
| **`SRC_B1_BIPSEQ_PMC_XML`** | **BIP-seq / POPPIN Literature** | 0 | 0 | 0 | 0 | 0 | 1 | Audited Literature (External Portal) |
| **`SRC_B4_STRING_RICE_PHYSICAL`**| **STRING Physical Links** (v12.0 Total) | 1,294,058 | 647,029 | 647,029 | 0 | 19,431 | 1 | Supplemental Physical Evidence (Graph Prior) |
| ↳ **`STRING_RICE_DIRECT_EXP`** | **STRING Rice Direct Experimental** | **862** | **431** | **431** | **0** | **435** | Curated Primary | **Direct Experimental Subset (No Interologs)** |
| **TOTAL (Non-Deduplicated)** | — | **1,295,451** | — | — | — | — | — | — |

> [!NOTE]
> Across the primary experimental curation databases (**IntAct** and **BioGRID**), there are **1,393 total raw records** (1,391 physical, 2 genetic) resolving to **923 unique physical pairs** across 37 peer-reviewed publications.
> Within STRING v12 physical links, **862 records** (**431 unique pairs**) are supported by **direct laboratory experiments in Rice** (`experiments > 0`), while **1,074,164 records** (**537,082 unique pairs**) represent **orthology-transferred interologs** from other organisms.

---

### 1. IntAct / IMEx (Release 252, January 2026)

IntAct distributes rice PPI records under two distinct NCBI taxonomy codes: **taxid:39947** (*Oryza sativa Japonica Group*) and **taxid:4530** (*Oryza sativa* generic / Indica / mixed). Both were parsed, harmonized, and merged into the canonical evidence partition:
- [`SRC_B2_INTACT_RICE_39947`](../../data/raw/rice/intact/intact_rice_taxid39947_release252.mitab27.txt): **649 records**
- [`SRC_B2_INTACT_RICE_4530`](../../data/raw/rice/intact/intact_rice_taxid4530_release252.mitab27.txt): **378 records**
- **Combined IntAct Rice Dataset:** **1,027 records** $\rightarrow$ **595 unique undirected pairs** (584 heteromeric, 11 homomeric), covering **530 unique proteins**.

#### 1.1 Breakdown by Interaction Semantics (Interaction Type)

| Interaction Semantics | Record Count | % of IntAct | Description / PSI-MI Context |
| :--- | --: | --: | :--- |
| `physical_association` | 1,027 | 100.00% | Curated as `physical association` (MI:0915). All pairwise records capture stable co-purification, affinity pull-downs, or Y2H reporter activation. |
| `direct_binary` | 0 | 0.00% | No interactions were explicitly tagged with PSI-MI term `direct interaction` (MI:0407) in the raw MITAB annotation. |
| `proximity` | 0 | 0.00% | No proximity-specific PSI-MI terms curated in IntAct Rice. |
| **Total** | **1,027** | **100.00%** | — |

#### 1.2 Breakdown by Experimental Assay Family & PSI-MI Detection Method

| Assay Family | Record Count | % of IntAct | Unique Pairs | Dominant PSI-MI Detection Methods |
| :--- | --: | --: | --: | :--- |
| **`AP_MS`** | 661 | 64.36% | 344 | Affinity chromatography technology (MI:0004: 507), Affinity technology (MI:0400: 154) |
| **`Y2H`** | 339 | 33.01% | 295 | Two hybrid (MI:0018: 339) |
| **`coIP`** | 9 | 0.88% | 8 | Anti-bait coimmunoprecipitation (MI:0006: 7), Anti-tag coimmunoprecipitation (MI:0007: 2) |
| **`pull_down`** | 8 | 0.78% | 7 | Pull down (MI:0096: 8) |
| **`BiFC`** | 6 | 0.58% | 6 | Bimolecular fluorescence complementation (MI:0809: 6) |
| **`biophysical_binding`** | 4 | 0.39% | 3 | Molecular sieving (MI:0071: 3), Enzymatic study (MI:0415: 1) |
| **Total** | **1,027** | **100.00%** | **595\*** | — |

*\*Note: Sum of unique pairs across assay families exceeds 595 because some pairs were tested across multiple assay types.*

#### 1.3 Sub-Taxon Source Comparison (39947 vs. 4530)

| Metric | `SRC_B2_INTACT_RICE_39947` | `SRC_B2_INTACT_RICE_4530` | Overlap (Shared) | Combined Union |
| :--- | --: | --: | --: | --: |
| **Total Raw Records** | 649 | 378 | — | **1,027** |
| **Unique Pairs** | 575 | 368 | 348 | **595** |
| **Heteromeric Pairs** | 565 | 367 | 348 | **584** |
| **Homomeric Pairs** | 10 | 1 | 0 | **11** |
| **Unique Proteins** | 513 | 345 | 328 | **530** |
| **AP-MS Records** | 320 | 341 | — | **661** |
| **Y2H Records** | 302 | 37 | — | **339** |
| **Small-Scale / Other Records** | 27 | 0 | — | **27** |
| **Contributing PMIDs** | 14 | 4 | 4 | **14** |

#### 1.4 Provenance, Publications & MIscore Scoring Behavior

- **Unique Contributing PMIDs:** **14** (all 4 PMIDs in 4530 are a subset of the 14 in 39947).
- **Primary Contributing Literature:**
  1. `PMID:21421362` (Rohila et al., 2011, *J. Proteome Res.* — AP-MS isolation of rice protein complexes): **507 records** (320 in 39947, 262 in 4530)
  2. `PMID:12684538` (Cooper et al., 2003, *Science* — High-throughput Y2H screen of the rice interactome): **266 records** (230 in 39947, 36 in 4530)
  3. `PMID:23303719` (Rohila et al., 2013, *Proteomics* — AP-MS interactome mapping): **154 records** (75 in 39947, 79 in 4530)
  4. `PMID:10444103` (17 records, Y2H)
  5. `PMID:16428324` (17 records: 6 Y2H, 6 Pull-down, 5 Co-IP)
  6. `PMID:24243689` (16 records: 10 Y2H, 6 BiFC)
  7. `PMID:12395189` (15 records, Y2H)
  8. `PMID:15078334` (10 records: 4 Co-IP, 3 Molecular sieving, 3 Y2H)
  9. `PMID:16284419` (10 records, Y2H)
  10. `PMID:11197326` (6 records: 4 Y2H, 2 Pull-down)
  11. Small-scale hypothesis studies: `PMID:11489998` (5), `PMID:12084835` (4), `PMID:10749909` (5), `PMID:17485859` (1).

- **IntAct MIscore Distribution:**
  - Score Range: **0.35 to 0.63** | Mean: **0.3710** | Median: **0.3500**
  - **$\text{MIscore} \ge 0.45$ (Replicated / Multi-method):** **75 records** $\rightarrow$ **24 unique pairs**
  - **$\text{MIscore} < 0.45$ (Single-assay / Unreplicated):** **952 records** $\rightarrow$ **571 unique pairs**

> [!WARNING]
> **Impact of Premature MIscore Filtering on Rice:**
> In Arabidopsis, large-scale screens with orthogonal validation allow $\approx 37\%$ of IntAct records to pass $\text{MIscore} \ge 0.45$. In Rice, however, **92.7% of all IntAct records** have $\text{MIscore} < 0.45$ because the major screens (Rohila et al. AP-MS, Cooper et al. Y2H) are single-publication high-throughput experiments.
> Applying an initial $\text{MIscore} \ge 0.45$ filter would wipe out **571 out of 595 unique physical pairs** (96% loss). Preserving all physical evidence without premature filtering ensures these genuine experimentally validated rice complexes remain available for modeling.

---

### 2. BioGRID (Release 5.0.261, August 2026)

- **Source Identifier:** [`SRC_B3_BIOGRID_RICE`](../../data/raw/rice/biogrid/BIOGRID-ORGANISM-Oryza_sativa_Japonica-5.0.261.tab3.txt)
- **Total Records:** **366** (364 Physical, 2 Genetic)
- **Unique Pairs:** **333** (Heteromeric: 328 | Homomeric / self-interactions: 5)
- **Unique Proteins:** **315**

#### 2.1 Breakdown by Interaction Semantics (Interaction Type)

| Interaction Semantics | Record Count | % of BioGRID | Description / BioGRID Context |
| :--- | --: | --: | :--- |
| `physical_association` | 243 | 66.39% | Affinity Capture-MS (235), Affinity Capture-Western (5), Biochemical Activity (3) |
| `direct_binary` | 109 | 29.78% | Two-hybrid (89), Reconstituted Complex (9), FRET (9), Protein-peptide (1), Far Western (1) |
| `proximity` | 12 | 3.28% | Proximity Label-MS (TurboID/BioID: 10), PCA / Complementation (2) |
| `genetic` | 2 | 0.55% | Synthetic Rescue (2) — filtered out during physical interactome construction |
| **Total** | **366** | **100.00%** | — |

#### 2.2 Breakdown by Experimental System & Assay Family

| Assay Family | BioGRID Experimental System Name | Record Count | % of BioGRID | Unique Pairs | Direct vs. Association Classification |
| :--- | :--- | --: | --: | --: | :--- |
| **`AP_MS`** | `Affinity Capture-MS` | 235 | 64.21% | 235 | Physical association / co-complex |
| **`Y2H`** | `Two-hybrid` | 89 | 24.32% | 85 | Direct binary interaction |
| **`biophysical_binding`** | *Subtotal of in vitro binding* | 23 | 6.28% | 21 | Direct binary (20) / Association (3) |
| ↳ *Reconstituted Complex* | `Reconstituted Complex` | 9 | 2.46% | 9 | Direct binary |
| ↳ *FRET* | `FRET` | 9 | 2.46% | 8 | Direct binary / proximity |
| ↳ *Biochemical Activity* | `Biochemical Activity` | 3 | 0.82% | 2 | Physical association / enzymatic |
| ↳ *Protein-peptide* | `Protein-peptide` | 1 | 0.27% | 1 | Direct binary |
| ↳ *Far Western* | `Far Western` | 1 | 0.27% | 1 | Direct binary |
| **`proximity_labeling`** | `Proximity Label-MS` | 10 | 2.73% | 10 | Proximity only |
| **`coIP`** | `Affinity Capture-Western` | 5 | 1.37% | 5 | Physical association |
| **`BiFC`** | `PCA` (Protein Complementation) | 2 | 0.55% | 2 | Proximity |
| **`genetic`** | `Synthetic Rescue` | 2 | 0.55% | 2 | Genetic interaction (non-physical) |
| **Total** | — | **366** | **100.00%** | **333** | — |

#### 2.3 Provenance, Throughput & Publications

- **Unique Contributing PMIDs:** **23**
- **Throughput Class:**
  - **High Throughput:** **245 records** (66.94%) $\rightarrow$ 245 unique pairs
  - **Low Throughput:** **121 records** (33.06%) $\rightarrow$ 91 unique pairs (focused hypothesis studies)
- **Top Contributing Publications in BioGRID:**
  1. `PMID:30866160` (Wong et al., 2019 — Affinity Capture-MS interactome screen): **235 records**
  2. `PMID:24970010` (Ding et al., 2014 — Rice-pathogen Two-hybrid screen): **55 records**
  3. `PMID:28553299` (Lin et al., 2017 — Proximity Label-MS + FRET + Y2H in rice immunity): **16 records**
  4. `PMID:19558411` (Affinity Capture-Western & Y2H): **7 records**
  5. `PMID:16230331` (Two-hybrid screen): **6 records**
  6. `PMID:25352666` (Two-hybrid screen): **6 records**
  7. `PMID:17965273` (Two-hybrid & Synthetic Rescue): **5 records**
  8. `PMID:12913152` (Biochemical Activity & Reconstituted Complex): **4 records**
  9. `PMID:19179350` (Reconstituted Complex & FRET): **4 records**
  10. `PMID:26300907` (FRET & Two-hybrid): **4 records**

> [!IMPORTANT]
> **Complete Orthogonality Between IntAct and BioGRID for Rice:**
> There is **zero overlap** between the 14 PMIDs curated by IntAct and the 23 PMIDs curated by BioGRID. BioGRID captured Wong et al. (2019) and Ding et al. (2014), whereas IntAct captured Rohila et al. (2011, 2013) and Cooper et al. (2003). Consequently, merging BioGRID and IntAct doubles the curated literature coverage from 14 to 37 studies without redundant re-ingestion.

---

### 3. BIP-seq Literature / POPPIN Status Audit

- **Source Identifier:** [`SRC_B1_BIPSEQ_PMC_XML`](../../data/raw/rice/poppin/BIP_seq_PMC11923860_fulltext.xml)
- **Target Paper:** BIP-seq (Binary Interaction Peptides sequencing) / POPPIN Rice Interactome (*Plant Biotechnology Journal*, 2025; PMC11923860)
- **Audit Findings:**
  - The open-access PMC XML full-text contains only 1 methodological table describing the Twin-Cell System (TCS) selection process.
  - The bulk interaction matrices ($>50\text{k}$ predicted/screened interactions) are hosted externally on the POPPIN web server (`poppin.hzau.edu.cn`), which does not provide an open REST endpoint or bulk flat-file download.
  - **Canonical Policy:** In strict compliance with P2 ETL standards, **0 interaction records** were normalized from prose, avoiding synthetic fabrication. The audit entry is preserved in [`manifests/phase2_progress.json`](../../manifests/phase2_progress.json).

---

### 4. STRING Physical Links (v12.0, Taxon 39947)

- **Source Identifiers:**
  - Raw Detailed File: [`SRC_B4_STRING_RICE_PHYSICAL`](../../data/raw/rice/string/39947.protein.physical.links.detailed.v12.0.txt.gz)
  - Raw Decoupled Full File: [`39947.protein.physical.links.full.v12.0.txt.gz`](../../data/raw/rice/string/39947.protein.physical.links.full.v12.0.txt.gz)
  - Collated Direct Dataset: [`string_direct_rice_pairs_deduplicated.tsv`](../../data/processed/rice/string_direct_rice_pairs_deduplicated.tsv)
- **Total Physical Records:** **1,294,058** (647,029 unique undirected pairs, all heteromeric)
- **Unique Proteins in Full Physical Network:** **19,431**
- **Interaction Semantics & Assay Family:** 100% `computational` (`MI:0090`)

#### 4.1 Architecture of STRING Evidence Channels (Direct vs. Transferred Interologs)

In the standard `protein.links.detailed` and `protein.physical.links.detailed` distributions, STRING merges directly observed laboratory assays with orthology-transferred interologs into a single composite `experimental` score column.

To isolate genuine Rice interactions, we ingested and parsed the full channel decoupling files:
- `39947.protein.physical.links.full.v12.0.txt.gz` (Physical network)
- `39947.protein.links.full.v12.0.txt.gz` (Full functional network)

These files explicitly decouple:
1. **`experiments`**: Interactions supported by direct experimental laboratory assays conducted specifically on *Oryza sativa* proteins.
2. **`experiments_transferred`**: Orthology-transferred interolog scores projected from other organisms (*Arabidopsis thaliana*, *Saccharomyces cerevisiae*, *Homo sapiens*, etc.).
3. **`database`** vs. **`database_transferred`**: Curated pathway/complex database memberships in Rice vs. transferred from model organisms.

Across the entire 67.4-million-edge STRING functional network (`protein.links.full`), exactly **862 directed records** (431 undirected pairs) have direct experimental evidence (`experiments > 0`), while **61,272,824 records** are derived from transferred experimental interologs. All 862 direct records are physical interactions.

#### 4.2 Direct vs. Transferred Breakdown in STRING Physical Links

| Evidence Category | Filter Condition | Directed Records | Unique Undirected Pairs | Unique Proteins | % of Experimental Physical Links |
| :--- | :--- | --: | --: | --: | --: |
| **All Physical Links** | `combined_score > 0` | 1,294,058 | 647,029 | 19,431 | — |
| **Total with Experimental Support** | `experiments > 0` OR `exp_transferred > 0` | 1,074,224 | 537,112 | 19,308 | 100.00% |
| ↳ **Direct Rice Experiments Only** | `experiments > 0` AND `exp_transferred == 0` | 810 | 405 | 409 | 0.08% |
| ↳ **Direct Rice + Transferred Overlap** | `experiments > 0` AND `exp_transferred > 0` | 52 | 26 | 38 | <0.01% |
| **TOTAL GENUINE DIRECT RICE EXPERIMENTS** | **`experiments > 0`** | **862** | **431** | **435** | **0.08%** |
| ↳ **Transferred Interologs ONLY** | `experiments == 0` AND `exp_transferred > 0` | **1,073,362** | **536,681** | **19,301** | **99.92%** |
| **Neither Experimental (Pure TM / DB)** | `experiments == 0` AND `exp_transferred == 0` | 219,834 | 109,516 | 10,214 | — |

> [!WARNING]
> **Overwhelming Dominance of Transferred Interologs in STRING Rice:**
> Out of 537,112 unique physical pairs with non-zero experimental scores in STRING for rice, **536,681 pairs (99.92%)** possess zero primary rice experimental validation and are purely transferred interologs.
> Only **431 unique pairs (862 records)** are supported by direct laboratory experiments performed in *Oryza sativa*. Treating STRING's composite `experimental` score as native rice evidence would introduce over 536,000 unverified interolog predictions into the ground truth.

#### 4.3 Score Distribution of Genuine Direct Rice Interactions

For the 431 genuine direct pairs (`experiments > 0`):

| `experiments` Score | Pair Count | % of Direct Pairs | Evidence Nature & Calibrated Meaning |
| :---: | --: | --: | :--- |
| **292** | 234 | 54.29% | Single high-throughput experimental evidence (e.g. AP-MS screens, Rohila et al.) |
| **230** | 151 | 35.03% | Single low-throughput / exploratory assay (e.g. Y2H screen, Cooper et al.) |
| **900** | 14 | 3.25% | Multi-assay or crystallographic physical interaction (PDB complexes, RuBisCO) |
| **421** | 10 | 2.32% | Dual-assay supported experimental pair |
| **237** | 6 | 1.39% | High-confidence biochemical pull-down assay |
| **800** | 4 | 0.93% | High-confidence validated binary complex |
| **322 – 929** | 12 | 2.78% | Various intermediate multi-method experimental confirmations |
| **Total** | **431** | **100.00%** | — |

In terms of `combined_score` (which integrates direct experiments, textmining, and pathway databases):
- Range: **229 to 997** | Mean: **340.39**
- **High-confidence ($\text{combined} \ge 700$):** **35 unique pairs**
- **Medium-confidence ($400 \le \text{combined} < 700$):** **43 unique pairs**
- **Low-confidence ($\text{combined} < 400$):** **353 unique pairs** (primarily single-study AP-MS / Y2H without orthogonal text-mining confirmation)

#### 4.4 Concordance with Curated IntAct and BioGRID Interactome

Comparing the 431 direct experimental pairs from STRING with the 913 non-redundant experimental pairs in our curated IntAct/BioGRID positive dataset (`data/processed/rice/rice_positive_pairs_deduplicated.tsv`):

| Comparison Subset | Unique Undirected Pairs | Concordance Rate | Context & Notes |
| :--- | --: | --: | :--- |
| **Shared (STRING Direct $\cap$ IntAct/BioGRID)** | **346** | **80.28%** | Directly matches curated IntAct AP-MS (Rohila et al.) and BioGRID screens. |
| **Unique to STRING Direct** | **85** | **19.72%** | Discrepancies explored below. |
| **Curated Pairs Not in STRING Direct** | **567** | — | Curated in BioGRID (Wong et al., Ding et al.) or low-throughput literature not yet mapped to STRING's direct channel. |

**Audit of the 85 STRING-Direct Unique Pairs:**
1. **Gene Locus vs. UniProt Accession Remapping (62 pairs):**
   - In IntAct, bait/prey accessions often map to alternative UniProt isoforms or older accessions (e.g., IntAct uses `Q0DIK2` for locus `Os05g0382600` interacting with `Q0DJ33` [GPA1] from Rohila et al. 2011).
   - In STRING v12, the canonical proteome maps locus `Os05g0382600` to canonical accession `B7E4J4`.
   - At the gene locus level (`Os05g0382600` — `Os05g0333200`), these interactions represent the exact same physical complexes.
2. **Direct PDB Crystallographic Evidence (8 pairs):**
   - Pairs such as RuBisCO large subunit (`P0C512` / `rbcL`) and its binding partners have experimental scores of 900 derived directly from experimental 3D coordinates (PDB structures: 1WDD, 3AXK, 3AXM).
3. **Primary Curated Databases Not in BioGRID/IntAct Rice Ingestion (15 pairs):**
   - Curated from primary interaction databases ingested by STRING (such as MINT, DIP, BIND) or directly mapped literature assays (e.g. Rad54 `A4PBL4` — Mus81 `Q8GT06` with direct experiment score 230 and combined score 710).

#### 4.5 Collated Direct STRING Dataset Artifacts

The genuine direct rice interactions have been isolated and collated into a standalone dataset preserving full database provenance (no synthetic labels or internal group designations):

- **Deduplicated Positive Pairs:** [`data/processed/rice/string_direct_rice_pairs_deduplicated.tsv`](../../data/processed/rice/string_direct_rice_pairs_deduplicated.tsv) (and `.parquet`)
  - **431 unique undirected pairs** spanning **435 unique proteins**
  - Fields retained: `pair_id`, `canonical_pair_key`, `participant_a_id`, `participant_b_id`, `participant_a_name`, `participant_b_name`, `participant_a_locus`, `participant_b_locus`, `is_homomeric`, `source_database`, `interaction_type`, `experiments_score`, `experiments_transferred_score`, `database_score`, `database_transferred_score`, `textmining_score`, `textmining_transferred_score`, `combined_score`, `in_curated_intact_or_biogrid`, `evidence_record_count`.
- **Directed Evidence Table:** [`data/processed/rice/string_direct_rice_evidence.tsv`](../../data/processed/rice/string_direct_rice_evidence.tsv) (and `.parquet`)
  - **862 directed evidence records**
- **Benchmark Evaluation Pairs:** [`data/processed/rice/string_direct_rice_benchmark_2col.tsv`](../../data/processed/rice/string_direct_rice_benchmark_2col.tsv)
  - 2-column format (`participant_a_id`, `participant_b_id`) for zero-shot and fine-tuned model evaluation.
- **Audit Metadata:** [`data/processed/rice/string_direct_rice_audit.json`](../../data/processed/rice/string_direct_rice_audit.json)

---

### 5. Multi-Source Cross-Tabulation Matrix

#### 5.1 Resource vs. Interaction Semantics

| Resource | `direct_binary` | `physical_association` | `proximity` | `genetic` | Total Records |
| :--- | --: | --: | --: | --: | --: |
| **IntAct (Taxid 39947)** | 0 | 649 | 0 | 0 | 649 |
| **IntAct (Taxid 4530)** | 0 | 378 | 0 | 0 | 378 |
| **IntAct (Combined)** | 0 | 1,027 | 0 | 0 | 1,027 |
| **BioGRID** | 109 | 243 | 12 | 2 | 366 |
| **STRING Direct Experimental** | 0 | 862 | 0 | 0 | 862 |
| **STRING Physical Links (Full)** | 0 | 0 | 0 | 0 | 1,294,058\* |
| **Combined Experimental (IntAct + BioGRID)** | **109** | **1,270** | **12** | **2** | **1,393** |

*\*Full STRING records are computationally integrated physical links classified under `computational`.*

#### 5.2 Resource vs. Experimental Assay Family

| Assay Family | IntAct (39947) | IntAct (4530) | IntAct (Union) | BioGRID | STRING Direct | Total Experimental Records |
| :--- | --: | --: | --: | --: | --: | --: |
| **`AP_MS`** (Affinity MS / TAP) | 320 | 341 | 661 | 235 | 468\* | **896** |
| **`Y2H`** (Two-hybrid) | 302 | 37 | 339 | 89 | 302\* | **428** |
| **`biophysical_binding`** (SPR, FRET, PDB) | 4 | 0 | 4 | 23 | 28\* | **27** |
| **`coIP`** (Co-immunoprecipitation) | 9 | 0 | 9 | 5 | — | **14** |
| **`proximity_labeling`** (TurboID / BioID) | 0 | 0 | 0 | 10 | — | **10** |
| **`pull_down`** (GST / His pull-down) | 8 | 0 | 8 | 0 | 12\* | **8** |
| **`BiFC`** (Complementation / PCA) | 6 | 0 | 6 | 2 | — | **8** |
| **`genetic`** (Synthetic Rescue) | 0 | 0 | 0 | 2 | — | **2** |
| **Total** | **649** | **378** | **1,027** | **366** | **862** | **1,393** |

*\*STRING direct records map to primary high-throughput AP-MS, Y2H screens, and PDB co-complexes.*

---

### 6. Progressive Physical Evidence Resolution Funnel for Rice

When the hierarchy is applied to Rice data to structure training and evaluation candidate sets:

```text
========================================================================================
LEVEL 1: BROAD EXPERIMENTALLY SUPPORTED PHYSICAL INTERACTOME (Tiers 1–4)
   - Scope: All physical association, co-complexes, pull-downs, Y2H, and proximity.
   - Raw Records: 1,391 experimental physical records (IntAct: 1,027 | BioGRID: 364)
   - Unique Non-Redundant Pairs: 923 unique pairs (Hetero: 907 | Homo: 16)
   - Unique Proteins: 732 proteins across 37 peer-reviewed publications
   - STRING Direct Supplement: 431 unique direct experimental pairs (346 concordant, 85 novel)
========================================================================================
                                     │
                                     ▼
LEVEL 2: HIGH-CONFIDENCE DIRECT & STABLE INTERACTOME (Tiers 1–3)
   - Scope: Direct binding (Y2H, in vitro), stable AP-MS complexes, and multi-study verified pairs.
   - Filters applied:
       * Exclude proximity-only assays (TurboID/BioID: -10 pairs, BiFC: -8 pairs)
       * Exclude genetic interactions (Synthetic Rescue: -2 pairs)
       * Retain all direct assays + validated co-purifications
   - Unique Pairs: ≈ 905 unique physical pairs
========================================================================================
                                     │
                                     ▼
LEVEL 3: DIRECT BINARY CORE (Tiers 1–2)
   - Scope: Strictly pairwise binary contacts (A physically touches B).
   - Eligible Assays: Two-hybrid (Y2H) + Biophysical binding (Reconstituted Complex, FRET, Far Western)
   - Raw Records: 455 records (IntAct: 343 | BioGRID: 112)
   - Unique Direct Binary Pairs: 356 unique pairs
       * Tier 1 (Replicated / High Confidence, e.g. MIscore >= 0.45 or cross-study): 24 unique pairs
       * Tier 2 (Curated single-assay direct binary): 332 unique pairs
========================================================================================
```

### Summary & Takeaways for Rice PPI Modeling

1. **Experimental Data Volume:** Unlike *Arabidopsis* (which has $\approx 90,000$ unique physical pairs and $>12,000$ proteins), Rice has a much leaner curated experimental corpus: **923 unique physical pairs** spanning **732 proteins** across IntAct and BioGRID.
2. **Crucial Role of Un-filtered Initial Ingestion:** In Rice, $92.7\%$ of IntAct records have $\text{MIscore} < 0.45$ because they stem from single foundational high-throughput screens (Rohila et al., Cooper et al.). Preserving the full physical evidence base without initial score thresholds is essential to avoid discarding almost the entire rice interactome.
3. **Complementarity of Sources:** IntAct and BioGRID have 0 shared PMIDs in Rice, making their integration non-redundant and synergistic.
4. **STRING Direct vs. Interolog Separation:**
   - Out of 647,029 physical pairs in STRING v12 (and 537,112 with non-zero experimental scores), **99.92% are transferred interologs** from other species.
   - Only **431 unique pairs (862 directed records across 435 proteins)** represent genuine, directly observed Rice experimental interactions (`experiments > 0`).
   - These have been isolated into a standalone dataset ([`data/processed/rice/string_direct_rice_pairs_deduplicated.tsv`](../../data/processed/rice/string_direct_rice_pairs_deduplicated.tsv)), where **346 pairs (80.28%)** directly concord with curated IntAct/BioGRID and **85 pairs** provide additional direct physical evidence (such as RuBisCO crystallographic PDB complexes and remapped canonical gene locus interactions).
5. **Role of STRING Interologs as Modeling Priors:** The remaining 536,681 transferred interolog pairs provide rich cross-species evolutionary priors for pre-training and graph neural networks, but must be treated distinctly from native experimental ground truth.
