# Maize (*Zea mays*) Positive Interactome Curation & HIPPIE Scoring Report

**Dataset Version**: 1.0.0  
**Generated Date**: 2026-10-07  
**Species**: *Zea mays* (Maize)  
**NCBI Taxonomy ID**: `4577`  
**Benchmark Role**: Quarantined External Zero-Shot Holdout (`eligible_for_training = FALSE`)  
**Pipeline Location**: `scripts/data_preparation/collate_external_positives.py`  
**Output Directory**: `data/processed/maize/`

---

## Executive Summary

This report documents the curation, canonical pair deduplication, and HIPPIE quality scoring of the physical protein-protein interaction (PPI) dataset for **Maize (*Zea mays*)**. 

Unlike other crop species that possess low-throughput or fragmented experimental datasets, Maize features a large, genome-scale physical interactome (**51,025 physical evidence records** resolving into **50,984 unique physical protein pairs** across **14,888 distinct proteins**). This dataset is primarily driven by the genome-scale recombination-based library-versus-library yeast two-hybrid screen (RLL-Y2H) published by Han et al. (*Nature Genetics* 2023, PMID: [36581701](https://pubmed.ncbi.nlm.nih.gov/36581701/)).

In accordance with the AraPPLM benchmark standards:
1. All records were filtered strictly for **physical protein-protein interactions** (excluding non-physical spatial colocalizations and genetic associations).
2. Data provenance was preserved using standard database coordinates (`source_id`, `source_file`, `source_row_number`, `source_record_id`, `pmid`, `authors`).
3. Pair-level and evidence-level **HIPPIE quality scores** were computed using saturating sigmoidal functions across studies ($w_s=0.6$), experimental techniques ($w_t=0.3$), and orthology ($w_o=0.1$).

---

## 1. Primary Source Ingestion & Physical Filtering

Maize physical interaction evidence was harmonized across two primary curated databases:

| Source Database | Release Version | Raw Records Ingested | Excluded Records | Exclusion Reason | Retained Physical Records | Unique Physical Pairs |
| :--- | :---: | --: | --: | :--- | --: | --: |
| **IntAct** | Release 252 | 51,009 | 1 | Colocalization (`MI:0403`) | 51,008 | 50,972 |
| **BioGRID** | Release 5.0.261 | 18 | 1 | Proximity PCA | 17 | 12 |
| **TOTAL** | — | **51,027** | **2** | — | **51,025** | **50,984** |

*(Note: In Maize, 12 protein pairs overlap between IntAct and BioGRID, yielding 50,984 net unique canonical pairs).*

### Filtering Rules Applied:
* **IntAct Filtering**:
  * Retained interaction semantics: `physical_association` (`MI:0915`), `direct_binary` (`MI:0407`), and enzymatic phosphorylation (`MI:0217`).
  * Excluded: Fluorescence microscopy spatial colocalization (`MI:0403`).
* **BioGRID Filtering**:
  * Retained experimental systems: Two-hybrid, Affinity Capture-Western, Reconstituted Complex, Biochemical Activity, Far Western, Co-purification.
  * Excluded: Protein complementation assay (`PCA` / proximity without direct binding).
* **STRING Separation**:
  * STRING Physical Links v12.0 for Maize (`4577.protein.physical.links.detailed.v12.0.txt.gz`) contains 3,025,182 records. As verified by audit, $> 98\%$ of these represent computationally transferred interologs rather than native maize experiments. STRING links are quarantined and separated from the gold-standard experimental master database.

---

## 2. Deduplication & Network Properties

The 51,025 master evidence observations collapse into **50,984 unique canonical protein pairs**:

```
51,025 Evidence Records  ──►  Deduplication (min(A,B) | max(A,B))  ──►  50,984 Canonical Pairs
                                                                           ├── 50,409 Heteromers (98.9%)
                                                                           └──    575 Homomers (1.1%)
```

### 2.1 Evidence Redundancy Profile
Only 36 records represent multiple observations of the same pair:
* **Single Observation (1x)**: 50,948 pairs
* **Dual Observations (2x)**: 33 pairs (reciprocal bait/prey screening or Y2H + pull-down)
* **Triple Observations (3x)**: 1 pair
* **Quadruple Observations (4x)**: 2 pairs

### 2.2 Network Topology & Protein Coverage
* **Total Distinct Proteins**: 14,888 proteins
* **Pure Protein-Protein Pairs**: 50,983 pairs (1 pair involves an RNA-associated chaperone)
* **Median Degree**: 3 interaction partners per protein (Mean degree: 6.9)
* **Singletons (Degree = 1)**: 4,611 proteins
* **Hub Proteins (Degree $\ge 20$)**: 1,196 proteins ($\sim 8\%$)

---

## 3. HIPPIE Quality Score Formulation & Results

HIPPIE (Human Integrated Protein-Protein Interaction rEference) scores were computed according to Schaefer et al. (2012):

$$S = w_s \cdot s_s(n_s) + w_t \cdot s_t(n_t) + w_o \cdot s_o(n_o)$$

where each subscore follows the saturating sigmoidal response:
$$s_i(n) = \frac{2}{1 + e^{-a_i \cdot n}} - 1$$

* **Hyperparameters**:
  * Study weight: $w_s = 0.6$, steepness: $a_s = 2.3$
  * Technique weight: $w_t = 0.3$, steepness: $a_t = 0.2$
  * Orthology weight: $w_o = 0.1$, steepness: $a_o = 1.6$
  * High-confidence cutoff: $S \ge 0.72$

### 3.1 Experimental Technique Weights in Maize
* **Two-hybrid (`MI:0018`)**: Weight = `5.0` ($s_t(5.0) = 0.4621$)
* **Biochemical kinase assay (`MI:0424`)**: Weight = `7.5` ($s_t(7.5) = 0.6351$)
* **In vitro pull-down (`MI:0096`)**: Weight = `2.5` ($s_t(2.5) = 0.2449$)
* **Reconstituted complex / X-ray crystallography (`MI:0114`)**: Weight = `10.0` ($s_t(10.0) = 0.7616$)

### 3.2 Maize HIPPIE Score Distribution
Because 50,966 interactions originate from a single high-throughput Y2H screen ($n_s = 1$, $n_t = 5.0$, $n_o = 0$):

$$s_s(1) = \frac{2}{1 + e^{-2.3}} - 1 = 0.8176$$
$$s_t(5) = \frac{2}{1 + e^{-1.0}} - 1 = 0.4621$$
$$S_{\text{baseline}} = (0.6 \times 0.8176) + (0.3 \times 0.4621) = 0.4906 + 0.1386 = \mathbf{0.6293}$$

| Metric | Value | Interpretation |
| :--- | :---: | :--- |
| **Total Evaluated Pairs** | **50,984** | 100% of deduplicated canonical pairs |
| **High Confidence ($S \ge 0.72$)** | **6** (0.01%) | Multivalidated pairs (Y2H + pull-down / kinase / multi-paper) |
| **Medium Confidence ($0.45 \le S < 0.72$)** | **50,978** (99.99%) | Single-study high-throughput Y2H pairs ($S \approx 0.6293$) |
| **Low Confidence ($S < 0.45$)** | **0** (0.00%) | None |
| **Mean Score** | **0.6294** | Uniform, robust baseline |
| **Score Range** | **[0.5641, 0.7866]** | Min: unconfirmed low-weight assay; Max: dual-study multi-assay |

---

## 4. Benchmark Governance & Quarantine Policy

> [!IMPORTANT]
> **Strict Evaluation Holdout Policy**:
> Although the 50,984 maize interactions are physical, protein-only, and experimental, they are **strictly quarantined from model training** under [`manifests/DO_NOT_TRAIN_ON_THESE_SOURCES.txt`](../../manifests/DO_NOT_TRAIN_ON_THESE_SOURCES.txt).
> 
> *Rationale*:
> 1. Preserves an unseen, genome-wide monocot plant evaluation holdout to benchmark zero-shot cross-species generalization.
> 2. Avoids single-screen overfitting to the specific biophysical library composition and false-positive characteristics of the RLL-Y2H platform.

---

## 5. Artifact Manifest & Coordinates

All files are located in `data/processed/maize/`:

| File Name | Format | Record Count | Description |
| :--- | :---: | --: | :--- |
| `maize_positive_evidence.tsv` | TSV | 51,025 | Master physical evidence records with database provenance |
| `maize_positive_evidence.parquet` | Parquet | 51,025 | Compressed columnar master evidence records |
| `maize_positive_pairs_deduplicated.tsv` | TSV | 50,984 | Canonical unique pairs with aggregated provenance & HIPPIE score |
| `maize_positive_pairs_deduplicated.parquet` | Parquet | 50,984 | Compressed columnar deduplicated pairs table |
| `maize_positive_hippie_provenance.tsv` | TSV | 50,984 | Pair-level HIPPIE audit ($n_s, s_s, n_t, s_t, n_o, s_o$, weights) |
| `maize_positive_evidence_hippie_provenance.tsv` | TSV | 51,025 | Evidence-level HIPPIE record scores |
| `maize_positive_audit.json` | JSON | — | Machine-readable audit summary |
