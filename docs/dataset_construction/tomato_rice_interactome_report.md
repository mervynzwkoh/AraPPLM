# Tomato & Rice (with Soybean) Curated Physical Interactome & HIPPIE Scoring Report

**Dataset Version**: 1.0.0  
**Generated Date**: 2026-10-07  
**Species Covered**: 
- **Tomato** (*Solanum lycopersicum*, TaxID: `4081`) — External Zero-Shot Holdout
- **Rice** (*Oryza sativa Japonica*, TaxID: `39947` / `4530`) — Core Crop Training / Transfer Candidate
- **Soybean** (*Glycine max*, TaxID: `3847`) — External Zero-Shot Holdout  
**Pipeline Locations**: 
- `scripts/data_preparation/collate_rice_positives.py`
- `scripts/data_preparation/collate_external_positives.py`
- `scripts/data_preparation/compute_hippie_scores.py`  
**Output Directories**: `data/processed/tomato/`, `data/processed/rice/`, `data/processed/soybean/`

---

## Executive Summary

This report provides a unified, comparative analysis of the curated physical protein-protein interactomes and HIPPIE quality scores for **Tomato (*Solanum lycopersicum*)** and **Rice (*Oryza sativa*)**, alongside **Soybean (*Glycine max*)**.

Unlike Maize, which is dominated by a single ultra-high-throughput Y2H screen, Tomato, Rice, and Soybean represent curated literature interactomes derived from targeted individual studies, low-to-medium-throughput assays, and (for Rice) directly verified experimental STRING links.

```
┌────────────────────────────────────────────────────────────────────────────────────────────────┐
│                                 SPECIES COMPARISON SUMMARY                                     │
├─────────────────────┬──────────────┬──────────────┬──────────────┬─────────────────────────────┤
│ Metric              │ Tomato       │ Rice         │ Soybean      │ Primary Role in Benchmark   │
├─────────────────────┼──────────────┼──────────────┼──────────────┼─────────────────────────────┤
│ NCBI Taxonomy ID    │ 4081         │ 39947 / 4530 │ 3847         │ —                           │
│ Physical Evidence   │ 166 records  │ 2,241 records│ 66 records   │ Evidence-level observations │
│ Deduplicated Pairs  │ 123 pairs    │ 998 pairs    │ 48 pairs     │ Canonical unique PPIs       │
│ Unique Proteins     │ 61 proteins  │ 892 proteins │ 51 proteins  │ Participating interactors   │
│ Mean HIPPIE Score   │ 0.6559       │ 0.7544       │ 0.6571       │ Quality score [0.0, 1.0]    │
│ High Conf (S ≥ 0.72)│ 14 (11.4%)   │ 381 (38.2%)  │ 6 (12.5%)    │ Multi-validated PPIs        │
│ Benchmark Role      │ Zero-Shot    │ Training /   │ Zero-Shot    │ Quarantined holdout vs.     │
│                     │ Holdout      │ Transfer     │ Holdout      │ Training candidate pool     │
└─────────────────────┴──────────────┴──────────────┴──────────────┴─────────────────────────────┘
```

---

## 1. Primary Source Ingestion & Physical Filtering

For all three species, raw interaction records were parsed from IntAct Release 252 and BioGRID 5.0.261. For Rice, directly observed experimental records from STRING v12.0 were also incorporated:

| Species | Primary Sources | Ingested Raw | Excluded Records | Exclusion Rationale | Retained Physical Records |
| :--- | :--- | --: | --: | :--- | --: |
| **Tomato** | IntAct (taxid 4081) | 29 | 0 | None (all physical) | 29 |
| | BioGRID (*S. lycopersicum*) | 141 | 4 | Genetic: Dosage Rescue | 137 |
| | **Tomato Subtotal** | **170** | **4** | — | **166** |
| **Rice** | IntAct (taxid 39947 & 4530) | 1,027 | 0 | None (all physical association) | 1,027 |
| | BioGRID (*O. sativa Japonica*) | 366 | 14 | Genetic & Proximity assays | 352 |
| | STRING v12.0 (Physical Full) | 1,294,058 | 1,293,196 | Excluded transferred interologs | 862 |
| | **Rice Subtotal** | **1,295,451** | **1,293,210** | — | **2,241** |
| **Soybean** | IntAct (taxid 3847) | 21 | 1 | Proximity: Colocalization (`MI:0403`) | 20 |
| | BioGRID (*G. max*) | 46 | 0 | None (all physical) | 46 |
| | **Soybean Subtotal** | **67** | **1** | — | **66** |

### Crucial Filtering Distinctions:
1. **Genetic Exclusions in Tomato**: BioGRID contains 4 genetic interactions for Tomato under the "Dosage Rescue" experimental system. These were strictly removed to maintain purely biophysical interactions.
2. **STRING Interolog Decoupling in Rice**: In Rice, STRING contains 1.29 million physical links, but **1,074,164 were transferred interologs** from *Arabidopsis* and other species (`experiments_transferred > 0`, `experiments == 0`). Only the **862 directly observed experimental records** were retained.
3. **Colocalization Exclusion in Soybean**: A single fluorescence microscopy colocalization record (`MI:0403`) in Soybean IntAct was excluded because spatial proximity does not confirm physical binding.

---

## 2. Canonical Deduplication & Interaction Types

After filtering, records were deduplicated into sorted canonical pairs `(min(A, B), max(A, B))` and classified into homomers and heteromers:

| Metric | Tomato (*S. lycopersicum*) | Rice (*O. sativa*) | Soybean (*G. max*) |
| :--- | --: | --: | --: |
| **Total Physical Evidence Observations** | **166** | **2,241** | **66** |
| **Canonical Deduplicated Pairs** | **123** | **998** | **48** |
| — Heteromeric Pairs ($A \ne B$) | 108 (87.8%) | 982 (98.4%) | 47 (97.9%) |
| — Homomeric Pairs ($A = B$) | 15 (12.2%) | 16 (1.6%) | 1 (2.1%) |
| — Pure Protein-Protein Pairs | 123 (100.0%) | 667 (66.8%)* | 48 (100.0%) |
| **Participating Unique Proteins** | **61** | **892** | **51** |
| **Independent Literature Studies (PMIDs)** | **15** | **37** | **11** |

*\*Note on Rice: 331 pairs in Rice involve UniProt entries annotated with nucleic acid-binding or ribonucleoprotein machinery in legacy databases, whereas Tomato and Soybean are 100% pure protein-protein pairs.*

---

## 3. HIPPIE Quality Score Comparison

HIPPIE quality scores were computed using identical mathematical formulations across all species:

$$S = 0.6 \cdot s_s(n_s) + 0.3 \cdot s_t(n_t) + 0.1 \cdot s_o(n_o)$$

$$s_s(n) = \frac{2}{1 + e^{-2.3 \cdot n}} - 1, \quad s_t(n) = \frac{2}{1 + e^{-0.2 \cdot n}} - 1, \quad s_o(n) = \frac{2}{1 + e^{-1.6 \cdot n}} - 1$$

### 3.1 Comparative Score Statistics

| Metric | Tomato | Rice | Soybean |
| :--- | :---: | :---: | :---: |
| **Total Evaluated Pairs** | **123** | **998** | **48** |
| **High Confidence ($S \ge 0.72$)** | **14** (11.4%) | **381** (38.2%) | **6** (12.5%) |
| **Medium Confidence ($0.45 \le S < 0.72$)** | **109** (88.6%) | **617** (61.8%) | **42** (87.5%) |
| **Low Confidence ($S < 0.45$)** | **0** (0.0%) | **0** (0.0%) | **0** (0.0%) |
| **Mean HIPPIE Score** | **0.6559** | **0.7544** | **0.6571** |
| **Median HIPPIE Score** | **0.6293** | **0.7191** | **0.6293** |
| **Score Range [Min, Max]** | **[0.5641, 0.8868]** | **[0.6293, 0.9579]** | **[0.5206, 0.7799]** |

### 3.2 Why Rice Scores Substantially Higher:
1. **Multi-Study Validation ($n_s > 1$)**: Rice has 381 pairs validated by multiple independent publications or dual-source entries (IntAct + BioGRID + STRING direct), driving the study subscore $s_s(n_s)$ toward saturation ($> 0.98$).
2. **Technique Weight Diversity**: Rice includes direct biophysical and crystallographic confirmations (weight 10.0), whereas Tomato and Soybean are predominantly standard Yeast Two-Hybrid screens (technique weight 5.0, resulting in the baseline score of 0.6293).
3. **High-Confidence Tomato Pairs ($S \ge 0.72$)**: The 14 high-confidence Tomato pairs represent well-studied signaling complexes (e.g., Pto kinase interacting with AvrPto / AvrPtoB, and 14-3-3 protein complexes validated by multiple assays including Y2H, pull-down, and co-IP).

---

## 4. Benchmark Roles and Governance

* **Rice (*Oryza sativa*)**: 
  * Classified as **Eligible for Training / Transfer Learning** (`eligible_for_training = TRUE`).
  * Used for monocot training or supervised cross-species domain adaptation alongside *Arabidopsis*.
* **Tomato (*Solanum lycopersicum*) and Soybean (*Glycine max*)**:
  * Classified as **Quarantined External Zero-Shot Holdouts** (`eligible_for_training = FALSE`).
  * Strictly reserved for evaluating how well representations trained on *Arabidopsis* and Rice generalize to completely unseen dicot crop species.

---

## 5. Artifact Manifest & File Coordinates

All generated datasets follow the canonical AraPPLM schema:

### Tomato Files (`data/processed/tomato/`):
* `tomato_positive_evidence.tsv` / `.parquet` (166 records, 42 columns)
* `tomato_positive_pairs_deduplicated.tsv` / `.parquet` (123 pairs, 24 columns)
* `tomato_positive_hippie_provenance.tsv` / `.parquet` (123 pairs, 22 provenance columns)
* `tomato_positive_evidence_hippie_provenance.tsv` / `.parquet` (166 records, 12 provenance columns)
* `tomato_positive_audit.json` (machine-readable audit)

### Rice Files (`data/processed/rice/`):
* `rice_positive_evidence.tsv` / `.parquet` (2,241 records, 42 columns)
* `rice_positive_pairs_deduplicated.tsv` / `.parquet` (998 pairs, 24 columns)
* `rice_positive_hippie_provenance.tsv` / `.parquet` (998 pairs, 22 provenance columns)
* `rice_positive_evidence_hippie_provenance.tsv` / `.parquet` (2,241 records, 12 provenance columns)
* `rice_positive_audit.json` (machine-readable audit)

### Soybean Files (`data/processed/soybean/`):
* `soybean_positive_evidence.tsv` / `.parquet` (66 records, 42 columns)
* `soybean_positive_pairs_deduplicated.tsv` / `.parquet` (48 pairs, 24 columns)
* `soybean_positive_hippie_provenance.tsv` / `.parquet` (48 pairs, 22 provenance columns)
* `soybean_positive_evidence_hippie_provenance.tsv` / `.parquet` (66 records, 12 provenance columns)
* `soybean_positive_audit.json` (machine-readable audit)
