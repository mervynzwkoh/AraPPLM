### Source Inventory Overview

| Source ID | Database / Resource | Raw Records | Unique Pairs (Total) | Heteromeric Pairs | Homomeric Pairs | Unique Proteins | Primary PMIDs | Benchmark Role |
| :--- | :--- | --: | --: | --: | --: | --: | --: | :--- |
| **`SRC_A1_INTACT_ARA`** | **IntAct / IMEx** (Rel. 252) | 58,740 | 42,388 | 41,886 | 502 | 9,136 | 643 | Core Training Positive |
| **`SRC_A2_BIOGRID_ARA`** | **BioGRID** (Rel. 5.0.261) | 84,422 | 73,995 | 73,060 | 935 | 12,082 | 2,457 | Core Training Positive |
| **`SRC_A4_XLMS_2026`** | **PRIDE PhoX XL-MS** (2026) | 390,526 | 7,167 | 2,385 | 4,782 | 5,064 | 1 | **Temporal Holdout (Quarantined)** |
| **`SRC_A3_STRING_ARA`** | **STRING Physical Links** (v12.0) | 986,540 | 493,270 | 493,270 | 0 | 18,154 | 1 | Supplemental Positive |
| **TOTAL (Non-Deduplicated)** | — | **1,520,228** | — | — | — | — | — | — |

---

### 1. IntAct / IMEx (Release 252, January 2026)
- **Source Identifier:** [`SRC_A1_INTACT_ARA`](../../data/raw/arabidopsis/intact/intact_arabidopsis_taxid3702_release252.mitab27.txt)
- **Total Records:** **58,740**
- **Unique Pairs:** **42,388** (Heteromeric: 41,886 | Homomeric / self-interactions: 502)
- **Unique Proteins:** **9,136**

#### 1.1 Breakdown by Interaction Semantics (Interaction Type)
| Interaction Semantics | Record Count | % of IntAct | Description / PSI-MI Context |
| :--- | --: | --: | :--- |
| `physical_association` | 56,496 | 96.18% | Includes `physical association` (MI:0915) and general `association` (MI:0914) |
| `direct_binary` | 1,487 | 2.53% | Explicitly curated pairwise `direct interaction` (MI:0407) and enzymatic reactions |
| `proximity` | 757 | 1.29% | `proximity` (MI:2364) and spatial `colocalization` (MI:0403) |
| **Total** | **58,740** | **100.00%** | — |

#### 1.2 Breakdown by Experimental Assay Family
| Assay Family | Record Count | % of IntAct | Dominant PSI-MI Detection Methods |
| :--- | --: | --: | :--- |
| `Y2H` | 33,035 | 56.24% | Two-hybrid array (10,760), Cr-two hybrid (9,089), Validated two hybrid (5,542), Two hybrid (5,135), Prey pooling (2,508) |
| `split_ubiquitin` | 12,441 | 21.18% | Ubiquitin reconstruction (MI:0112; membrane yeast two-hybrid) |
| `biophysical_binding` | 4,683 | 7.97% | Solid-phase assay (3,720), Surface plasmon resonance (SPR), ITC, Fluorescence polarization |
| `AP_MS` | 2,430 | 4.14% | Tandem affinity purification (TAP: 1,441), Affinity chromatography (987) |
| `coIP` | 2,293 | 3.90% | Anti-tag coimmunoprecipitation (2,037), Anti-bait coimmunoprecipitation (253) |
| `pull_down` | 1,015 | 1.73% | GST / His pull-down (MI:0096) |
| `protein_microarray` | 799 | 1.36% | Protein array binding (MI:0089) |
| `proximity_labeling` | 760 | 1.29% | Proximity-dependent biotin identification (BioID/TurboID: 598) |
| `other` | 478 | 0.81% | Enzymatic assays (kinase, phosphatase, cleavage) |
| `BiFC` | 475 | 0.81% | Bimolecular fluorescence complementation (MI:0809) |
| `XL_MS` | 257 | 0.44% | Chemical cross-linking mass spectrometry |
| `split_luciferase` | 73 | 0.12% | Luciferase complementation assay |
| `unknown` | 1 | <0.01% | Unspecified legacy record |
| **Total** | **58,740** | **100.00%** | — |

#### 1.3 Provenance, Publications & Confidence Scoring
- **Unique Contributing PMIDs:** **643**
- **Top 5 Publications in IntAct:**
  1. `PMID:28650476` (Altmann et al., 2020 — Plant-pathogen Y2H interactome screen): **13,209 records**
  2. `PMID:24833385` (Jones et al., 2014 — MIND1 membrane split-ubiquitin screen): **12,185 records**
  3. `PMID:32612234` (Klopffleisch et al., 2011 — Heterotrimeric G-protein interactome): **6,942 records**
  4. `PMID:21798944` (Arabidopsis Interactome Mapping Consortium, 2011 — AI-1MAIN Y2H): **6,498 records**
  5. `PMID:30806640` (Vandereyken et al., 2018): **2,586 records**
- **IntAct MIscore Distribution:**
  - Mean MIscore: **0.4498** | Median: **0.3700** | Range: **0.27 to 0.97**
  - **$\text{MIscore} \ge 0.45$ (Medium/High Confidence):** **21,800 records** $\rightarrow$ **7,142 unique pairs**
  - **$\text{MIscore} < 0.45$ (Single-assay / Low Confidence):** **36,934 records** $\rightarrow$ **35,246 unique pairs**

---

### 2. BioGRID (Release 5.0.261, August 2026)
- **Source Identifier:** [`SRC_A2_BIOGRID_ARA`](../../data/raw/arabidopsis/biogrid/BIOGRID-ORGANISM-Arabidopsis_thaliana_Columbia-5.0.261.tab3.txt)
- **Total Records:** **84,422**
- **Unique Pairs:** **73,995** (Heteromeric: 73,060 | Homomeric / self-interactions: 935)
- **Unique Proteins:** **12,082**

#### 2.1 Breakdown by Interaction Semantics (Interaction Type)
| Interaction Semantics | Record Count | % of BioGRID | Description / BioGRID Context |
| :--- | --: | --: | :--- |
| `direct_binary` | 35,067 | 41.54% | Two-hybrid, biophysical binding, reconstituted complexes, cross-linking |
| `physical_association` | 32,997 | 39.09% | Co-fractionation (SEC-MS), AP-MS, co-purification, native Co-IP |
| `proximity` | 15,991 | 18.94% | Protein complementation assay (PCA / BiFC), Proximity Label-MS |
| `genetic` | 367 | 0.43% | Synthetic lethality, phenotypic enhancement, phenotypic suppression |
| **Total** | **84,422** | **100.00%** | — |

#### 2.2 Breakdown by Experimental System & Assay Family
| Assay Family | BioGRID Experimental System Name | Record Count | % of BioGRID | Unique Pairs |
| :--- | :--- | --: | --: | --: |
| **`AP_MS`** (Total) | *Subtotal of Co-fractionation + AP-MS* | **28,956** | 34.30% | **28,629** |
| ↳ *Co-fractionation* | `Co-fractionation` (SEC-MS co-elution) | 21,422 | 25.38% | 21,396 |
| ↳ *AP-MS* | `Affinity Capture-MS` | 7,451 | 8.83% | 7,238 |
| ↳ *Co-purification* | `Co-purification` | 83 | 0.10% | 79 |
| **`Y2H`** | `Two-hybrid` | **27,214** | 32.24% | **24,373** |
| **`BiFC`** | `PCA` (Protein Complementation Assay) | **15,042** | 17.82% | **14,855** |
| **`biophysical_binding`** | `Reconstituted Complex` (5,547), `Biochemical Activity` (1,922), `FRET` (673), `Far Western` (41), `Co-crystal Structure` (90), `Protein-peptide` (44) | **8,331** | 9.87% | **7,842** |
| **`coIP`** | `Affinity Capture-Western` (2,035), `Affinity Capture-Luminescence` (50) | **2,085** | 2.47% | **1,767** |
| **`XL_MS`** | `Cross-Linking-MS (XL-MS)` | **1,444** | 1.71% | **1,444** |
| **`proximity_labeling`** | `Proximity Label-MS` | **949** | 1.12% | **949** |
| **`genetic`** | `Phenotypic Enhancement/Suppression`, `Synthetic Lethality` | **367** | 0.43% | **354** |
| **`other`** | `Co-localization` (156), `RNA/Protein interactions` (34) | **34** | 0.04% | **34** |
| **Total** | — | **84,422** | **100.00%** | **73,995** |

#### 2.3 Provenance, Throughput & Publications
- **Unique Contributing PMIDs:** **2,457**
- **Throughput Class:**
  - **High Throughput:** **67,543 records** (80.01%) $\rightarrow$ 65,994 unique pairs
  - **Low Throughput:** **16,879 records** (19.99%) $\rightarrow$ 10,457 unique pairs (hypothesis-driven literature)
- **Top 5 Publications in BioGRID:**
  1. `PMID:32191846` (McWhite et al., 2020 *Cell* — SEC-MS Co-fractionation roadmap): **21,288 records**
  2. `PMID:24833385` (Jones et al., 2014 *Science* — MIND1 split-ubiquitin PCA screen): **12,097 records**
  3. `PMID:28650476` (Altmann et al., 2020 — Plant-pathogen Y2H screen): **8,327 records**
  4. `PMID:21798944` (AI-1MAIN, 2011 *Science* — Arabidopsis Interactome Consortium): **5,641 records**
  5. `PMID:29320478` (BioID proximity labeling): **2,586 records**

---

### 3. PhoX Cross-Linking Mass Spectrometry (Trinh et al., August 2026)
- **Source Identifier:** [`SRC_A4_XLMS_2026_PLINK_SEARCH`](../../data/raw/arabidopsis/xlms_2026/Total_XL_plink3-2_v3.csv)
- **Primary Publication:** Trinh et al., *Nature Communications* (August 2026); `PMID:39133827` (PRIDE PXD066234 / PXD066291)
- **Total Records:** **390,526 identified cross-linked peptide pairs**
- **Unique Pairs:** **7,167** (Inter-protein: 2,385 | Intra-protein: 4,782)
- **Unique Proteins Identified:** **5,064**

#### 3.1 Breakdown by Interaction Semantics & Cross-Link Topology
| Cross-Link Topology | Interaction Semantics | Record Count | % of XL-MS | Unique Locus Pairs | Biological Interpretation |
| :--- | :--- | --: | --: | --: | :--- |
| **Inter-Protein Cross-Links** | `direct_binary` | **39,765** | 10.18% | **2,385** | Pairwise cross-links bridging two distinct protein loci ($<35$Å contact) |
| **Intra-Protein Cross-Links** | `direct_binary` (structural) | **350,761** | 89.82% | **4,782** | Monolinks / homomeric cross-links within the same protein locus |
| **Total** | — | **390,526** | **100.00%** | **7,167** | — |

#### 3.2 Provenance, Assay Details & Confidence Metrics
- **Assay Family:** `XL_MS` (100.0%)
- **Cross-Linker:** PhoX (bis(sulfosuccinimidyl)phosphosuberic acid), reactive towards primary amine groups (lysine $\text{C}_\alpha–\text{C}_\alpha$ Euclidean distance constraint $<35$Å).
- **Search & FDR:** Identified via pLink 3.2 under a strict **$<1\%$ peptide-spectrum match (PSM) false discovery rate**.
- **Lysate Fractions:** Whole cell native tissue lysate (`PXD066234`) and enriched chloroplast fraction (`PXD066291`).
- **Benchmark Role:** **Held-Out Prospective Temporal Test Set** (`is_temporal_holdout = TRUE`, `eligible_for_training = FALSE`).

---

### 4. STRING Physical Links (v12.0, Stable Freeze)
- **Source Identifier:** [`SRC_A3_STRING_ARA_PHYSICAL`](../../data/raw/arabidopsis/string/3702.protein.physical.links.detailed.v12.0.txt.gz)
- **Total Records:** **986,540 directed links**
- **Unique Pairs (Undirected):** **493,270** (Heteromeric: 493,270 | Homomeric: 0)
- **Unique Proteins:** **18,154**

#### 4.1 Breakdown by Interaction Semantics & Assay Family
| Dimension | Category | Record Count | % of STRING |
| :--- | :--- | --: | --: |
| **Interaction Semantics** | `computational` | 986,540 | 100.00% |
| **Assay Family** | `computational` (PSI-MI: `MI:0090` computational combinatorial) | 986,540 | 100.00% |

#### 4.2 Score Channel Distribution & Quality Breakdown
STRING physical links report subscores ($0–1000$) across distinct evidence channels:
| Score Filter / Evidence Tier | Directed Records | Unique Undirected Pairs | % of Total Pairs | Benchmark Role |
| :--- | --: | --: | --: | :--- |
| **All Physical Links** | 986,540 | 493,270 | 100.00% | Unfiltered predicted physical universe |
| **Experimental Score $> 0$** | 811,900 | 405,950 | 82.30% | Possesses laboratory experimental evidence (transferred or direct) |
| **Combined Score $\ge 700$** | 83,908 | 41,954 | 8.51% | High-confidence aggregate physical links |
| **Experimental Score $\ge 700$** | **16,686** | **8,343** | **1.69%** | **High-Confidence Laboratory Experimental Core** |

- **Publication:** Szklarczyk et al., 2023 *Nucleic Acids Res* (PMID 36370105).
- **Benchmark Role:** Supplemental high-confidence filter (`experimental >= 700`).

---

### 5. Multi-Source Cross-Tabulation Matrix

Cross-tabulating **Interaction Semantics** against **Assay Family** across all 1,520,228 non-deduplicated Arabidopsis records:

| Interaction Semantics | Y2H | Split-Ubiquitin | AP-MS / Co-fract. | PCA / BiFC | Biophysical / Reconstituted | Co-IP | Proximity Labeling | Pull-down | Protein Array | XL-MS | Genetic | Computational / Other | TOTAL |
| :--- | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: |
| **`direct_binary`** | 27,215 | 0 | 11 | 17 | 6,700 | 87 | 1 | 271 | 799 | 391,974 | 0 | 5 | **427,080** |
| **`physical_association`** | 33,034 | 12,441 | 31,375 | 455 | 6,314 | 4,291 | 5 | 744 | 0 | 253 | 0 | 508 | **89,420** |
| **`proximity`** | 0 | 0 | 0 | 15,045 | 0 | 0 | 1,703 | 0 | 0 | 0 | 0 | 0 | **16,748** |
| **`genetic`** | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 367 | 0 | **367** |
| **`computational`** | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 986,540 | **986,540** |
| **TOTAL** | **60,249** | **12,441** | **31,386** | **15,517** | **13,014** | **4,378** | **1,709** | **1,015** | **799** | **392,227** | **367** | **987,053** | **1,520,228** |