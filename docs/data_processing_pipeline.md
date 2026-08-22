# Technical Documentation: Sequence Database Preparation & Data Processing Pipeline

## 1. Overview & Objective

In deep learning models for Protein-Protein Interaction (PPI) prediction like **PPLM (Paired Sequence Language Model)**, pairs of protein accessions (e.g., `Q9LUI9` and `F4IDF1`) must be converted into their corresponding full-length amino acid sequences before tokenization and feature extraction.

The primary objective of this data processing pipeline is to construct fast, indexed, high-coverage sequence lookup databases (`data/arabidopsis/uniprot_final.pkl` and `data/rice/uniprot_rice_final.pkl`) that guarantee **100% sequence coverage** across all benchmark datasets in the DeepAraPPI suite.

```
+------------------------------------+
|  UniProt Proteome FASTA            |
|  (Swiss-Prot + TrEMBL, ~50k-54k)   |
+-----------------+------------------+
                  |
                  v  [Step 1: Regex Accession Indexing]
+-----------------+------------------+      +-------------------------------+
|  Initial Index (~90-97% Coverage)   |      |  UniParc Missing Sequences    |
|  (Unresolved / Subspecies IDs)     |      |  (Retrieved via UniParc/API)  |
+-----------------+------------------+      +---------------+---------------+
                  |                                          |
                  +-------------------+----------------------+
                                      |
                                      v  [Step 2: UniParc Rescue & Merge]
                      +---------------+---------------+
                      |  Final Sequence DB (100%)     |
                      |  `uniprot_final.pkl`          |
                      |  (~100k-110k Mappings)        |
                      +---------------+---------------+
                                      |
                                      v  [Step 3: Verification & Inference]
                      +---------------+---------------+
                      |  `batch_predict.py` (PPLM)    |
                      |  DeepAraPPI C1 / C2 / C3 /Rice|
                      +-------------------------------+
```

---

## 2. Key Challenges & Technical Solutions

### Challenge 1: Header & Accession Format Mismatches
* **Problem:** 
  PPI datasets (e.g., DeepAraPPI) reference proteins using bare UniProt accessions:
  ```tsv
  Protein1    Protein2    label
  Q9LUI9      F4IDF1      0
  A0A178VKB9  Q9M8S6      1
  ```
  However, standard UniProt FASTA headers contain compound headers with source prefixes and entry names:
  ```fasta
  >sp|Q9LUI9|PIAL1_ARATH Protein PIAL1 OS=Arabidopsis thaliana ...
  >tr|A0A178VKB9|A0A178VKB9_ARATH Uncharacterized protein OS=Arabidopsis thaliana ...
  ```
* **Solution:**
  The parser utilizes regular expression matching:
  $$\text{Regex: } \verb|^(?:sp|tr)\|([A-Z0-9]+)\||$$
  Both the raw header (`sp|Q9LUI9|PIAL1_ARATH`) and the extracted bare accession (`Q9LUI9`) are registered in the dictionary, allowing instantaneous $O(1)$ lookups regardless of query format.

---

### Challenge 2: The Arabidopsis "Missing 2.4%" Legacy/TrEMBL Discrepancy
* **Problem:**
  When querying the initial UniProt *Arabidopsis* proteome against the 130,478 pair DeepAraPPI dataset, **2.4% of protein pairs (1,348 pairs across 348 unique protein IDs)** failed to resolve.
  * *Root Cause:* DeepAraPPI (compiled from historical IntAct, BioGRID, DIP, TAIR, and MINT records) contains older TrEMBL accessions that have been retired, merged, or moved into UniParc (UniProt Archive).
* **Solution (UniParc Cross-Referencing & Sequence Rescue):**
  1. All unresolved protein IDs were isolated into `data/arabidopsis/missing_ids.txt`.
  2. The exact matching sequence records were retrieved from UniParc and saved to `data/arabidopsis/uniparc_348.fasta`.
  3. An ID mapping table `data/arabidopsis/uniparc_id_mapping.json` links original DeepAraPPI accessions directly to their UniParc identifier (`UPI...`).
  4. Merging maps all 348 rescued sequences into the main dictionary.
* **Result:** **0 missing sequences (100.0% coverage across all 130,478 pairs).**

---

### Challenge 3: The Rice (*Oryza sativa*) Subspecies Coverage Gap
* **Problem:**
  When querying the downloaded *Oryza sativa subsp. japonica* proteome (Taxonomy: `39947`, 49,985 sequences) against the 6,721 pair Rice benchmark, **10.15% of pairs (682 pairs across 163 unique protein IDs)** were missing.
  * *Root Cause:* DeepAraPPI includes interaction pairs gathered across multiple rice cultivars and subspecies, including *Oryza sativa subsp. indica* (Taxonomy: `39946`, e.g., `Q6YW54`, `Q8W0F9`) and general *Oryza sativa* (Taxonomy: `4530`).
* **Solution (Targeted UniParc Rescue for Rice):**
  1. [`scripts/data_preparation/verify_coverage.py`](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/scripts/data_preparation/verify_coverage.py) exported the exact 163 missing accessions to `data/rice/missing_rice_ids.txt`.
  2. Matching sequences were retrieved from UniParc into `data/rice/missing_rice.fasta`.
  3. [`scripts/data_preparation/merge_rice_uniparc.py`](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/scripts/data_preparation/merge_rice_uniparc.py) merged these sequences into `data/rice/uniprot_rice_final.pkl`.
* **Result:** **6,721 out of 6,721 pairs found (100.00% sequence coverage).**

---

## 3. Dataset Architecture & Benchmark Suites

All interaction ground truth datasets are located in [`data/DeepAraPPI/`](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/data/DeepAraPPI):

### 3.1 Arabidopsis thaliana Partitions (Park & Marcotte Scheme)

| Dataset File | Positive Pairs | Negative Pairs | Total Pairs | Class Ratio | Difficulty Level | What it Tests |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `c1_ppi_sample_DeepAraPPI.txt` | 1,000 | 10,000 | 11,000 | 1 : 10 | **Low (Task 1)** | Both proteins seen in training domain |
| `c2_ppi_sample_DeepAraPPI.txt` | 1,000 | 10,000 | 11,000 | 1 : 10 | **Medium (Task 2)** | **One protein unseen** (partially novel pairs) |
| `c3_ppi_sample_DeepAraPPI.txt` | 1,000 | 10,000 | 11,000 | 1 : 10 | **High (Task 3)** | **Both proteins unseen** (true zero-shot generalization) |
| `total_positive_negative_samples_DeepAraPPI.txt` | 11,858 | 118,620 | 130,478 | 1 : 10 | **Full Interactome** | HIPPIE score $\ge 0.72$ with subcellular exclusion |
| `total_positive_samples_DeepAraPPI.txt` | 11,858 | 0 | 11,858 | Positive only | Gold-standard physical interactions |

### 3.2 Cross-Species Evaluation (Monocot Transfer)

| Dataset File | Positive Pairs | Negative Pairs | Total Pairs | Class Ratio | Purpose |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `all_rice_positive_negative_DeepAraPPI.txt` | 611 | 6,110 | 6,721 | 1 : 10 | **Task 4:** *Oryza sativa* cross-species transfer benchmark |
| `all_rice_PPI_positive_DeepAraPPI.txt` | 611 | 0 | 611 | Positive only | Curated non-redundant physical rice interactions |

---

## 4. Pipeline Execution Reference

### 4.1 Building the Arabidopsis Sequence Database
```bash
python scripts/data_preparation/build_sequence_db.py \
    --fasta data/arabidopsis/uniprot_arabidopsis.fasta \
    --uniparc_fasta data/arabidopsis/uniparc_348.fasta \
    --mapping data/arabidopsis/uniparc_id_mapping.json \
    --output data/arabidopsis/uniprot_final.pkl \
    --verify data/DeepAraPPI/total_positive_negative_samples_DeepAraPPI.txt
```

### 4.2 Building & Merging the Rice Sequence Database
```bash
# 1. Build primary Rice database
python scripts/data_preparation/build_sequence_db.py \
    --fasta data/rice/uniprot_rice.fasta \
    --output data/rice/uniprot_rice_final.pkl

# 2. Merge missing UniParc sequences
python scripts/data_preparation/merge_rice_uniparc.py \
    --db data/rice/uniprot_rice_final.pkl \
    --missing_ids data/rice/missing_rice_ids.txt \
    --fasta data/rice/missing_rice.fasta \
    --output data/rice/uniprot_rice_final.pkl \
    --verify data/DeepAraPPI/all_rice_positive_negative_DeepAraPPI.txt
```

### 4.3 Verifying Coverage for Any Dataset
```bash
python scripts/data_preparation/verify_coverage.py \
    --dataset data/DeepAraPPI/all_rice_positive_negative_DeepAraPPI.txt \
    --seq_db data/rice/uniprot_rice_final.pkl \
    --output_missing data/rice/missing_rice_ids.txt
```

---

## 5. Summary of Sequence Databases

| Database File | Organism | Accessions / Mappings | Coverage on DeepAraPPI | Status |
| :--- | :--- | :--- | :--- | :--- |
| `data/arabidopsis/uniprot_final.pkl` | *Arabidopsis thaliana* | 109,640 entries | **100.00%** (130,478 / 130,478 pairs) | ✅ Complete |
| `data/rice/uniprot_rice_final.pkl` | *Oryza sativa* | 100,297 entries | **100.00%** (6,721 / 6,721 pairs) | ✅ Complete |
