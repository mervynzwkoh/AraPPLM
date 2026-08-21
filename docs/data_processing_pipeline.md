# Technical Documentation: Sequence Database Preparation & Data Processing Pipeline

## 1. Overview & Objective

In deep learning models for Protein-Protein Interaction (PPI) prediction like **PPLM (Paired Sequence Language Model)**, pairs of protein accessions (e.g., `Q9LUI9` and `F4IDF1`) must be converted into their corresponding full-length amino acid sequences before tokenization and feature extraction.

The primary objective of this data processing pipeline is to construct a fast, indexed, high-coverage sequence lookup database (`uniprot_final.pkl`) that guarantees **100% sequence coverage** across all benchmark datasets (such as DeepAraPPI *Arabidopsis thaliana* and *Oryza sativa* interaction sets).

```
+------------------------------------+
|  UniProt Proteome FASTA            |
|  (Swiss-Prot + TrEMBL, ~54k seqs)  |
+-----------------+------------------+
                  |
                  v  [Step 1: Regex Accession Indexing]
+-----------------+------------------+      +-------------------------------+
|  Initial Index (97.6% Coverage)     |      |  UniParc Missing Sequences    |
|  (Missing 348 legacy/deleted IDs)   |      |  (348 retrieved sequences)    |
+-----------------+------------------+      +---------------+---------------+
                  |                                          |
                  +-------------------+----------------------+
                                      |
                                      v  [Step 2: UniParc Rescue & Merge]
                      +---------------+---------------+
                      |  Final Sequence DB (100%)     |
                      |  `uniprot_final.pkl`          |
                      |  (~110k Key-Sequence Mappings)|
                      +---------------+---------------+
                                      |
                                      v  [Step 3: Verification & Inference]
                      +---------------+---------------+
                      |  `batch_predict.py` (PPLM)    |
                      |  DeepAraPPI C1 / C2 / C3      |
                      +-------------------------------+
```

---

## 2. Key Challenges & Solutions

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

### Challenge 2: The "Missing 2.4%" Legacy / TrEMBL Discrepancy
* **Problem:**
  When querying the initial UniProt *Arabidopsis* proteome against the 130,478 pair DeepAraPPI dataset, **2.4% of protein pairs (1,348 pairs across 348 unique protein IDs)** failed to resolve.
  * *Root Cause:* DeepAraPPI (compiled from historical IntAct, BioGRID, DIP, TAIR, and MINT records) contains older TrEMBL accessions that have been retired, merged, or moved into UniParc (UniProt Archive).
* **Solution (UniParc Cross-Referencing & Sequence Rescue):**
  1. All unresolved protein IDs were isolated into `missing_ids.txt`.
  2. The exact matching sequence records were retrieved from UniParc and saved to `data/uniparc_348.fasta`.
  3. An ID mapping table `data/uniparc_id_mapping.json` links original DeepAraPPI accessions directly to their UniParc identifier (`UPI...`).
  4. The merge phase maps all 348 rescued sequences into the main dictionary.
* **Result:** **0 missing sequences (100.0% coverage across all 130,478 pairs).**

---

## 3. Dataset Architecture & Benchmarks

The benchmark datasets are organized in [`data/DeepAraPPI/`](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/data/DeepAraPPI):

### 3.1 Arabidopsis thaliana Partitions (Park & Marcotte Evaluation Scheme)

| Dataset File | Positive Pairs | Negative Pairs | Total Pairs | Class Ratio | Purpose / Evaluation Level |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `c1_ppi_sample_DeepAraPPI.txt` | 1,000 | 10,000 | 11,000 | 1 : 10 | **C1 Benchmark:** Both interacting proteins may be present in the training distribution. Tests standard pattern recognition. |
| `c2_ppi_sample_DeepAraPPI.txt` | 1,000 | 10,000 | 11,000 | 1 : 10 | **C2 Benchmark:** Exactly one protein in each test pair is unseen in training. Tests single-partner generalization. |
| `c3_ppi_sample_DeepAraPPI.txt` | 1,000 | 10,000 | 11,000 | 1 : 10 | **C3 Benchmark:** Both proteins in the test pair are unseen in training. Hardest benchmark; tests true zero-shot representation. |
| `total_positive_negative_samples_DeepAraPPI.txt` | 11,858 | 118,620 | 130,478 | 1 : 10 | **Full Interactome:** Complete Arabidopsis dataset filtered by HIPPIE confidence score $\ge 0.72$ with subcellular exclusion negatives. |
| `total_positive_samples_DeepAraPPI.txt` | 11,858 | 0 | 11,858 | Positive only | Gold-standard physical interactions from IntAct, BioGRID, TAIR, DIP, and MINT. |

### 3.2 Cross-Species Evaluation (Generalization)

| Dataset File | Positive Pairs | Negative Pairs | Description |
| :--- | :--- | :--- | :--- |
| `all_rice_positive_negative_DeepAraPPI.txt` | 611 | 6,110 | *Oryza sativa* (Rice) cross-species generalization benchmark. Evaluates cross-species transfer from Arabidopsis/Animal models to monocot plants. |
| `all_rice_PPI_positive_DeepAraPPI.txt` | 611 | 0 | Curated positive physical interaction pairs for Rice. |

---

## 4. Pipeline Execution & Usage

The data preparation tools reside in [`scripts/data_preparation/`](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/scripts/data_preparation).

### 4.1 Step 1: Building the Complete Database

Use `build_sequence_db.py` to parse the raw UniProt FASTA and merge UniParc sequences in a single step:

```bash
python scripts/data_preparation/build_sequence_db.py \
    --fasta data/arabidopsis/uniprot_arabidopsis.fasta \
    --uniparc_fasta data/arabidopsis/uniparc_348.fasta \
    --mapping data/arabidopsis/uniparc_id_mapping.json \
    --output data/arabidopsis/uniprot_final.pkl \
    --verify data/DeepAraPPI/total_positive_negative_samples_DeepAraPPI.txt
```

#### CLI Parameters:
* `--fasta`: Path to the primary organism FASTA download from UniProt.
* `--uniparc_fasta` (optional): Path to rescued UniParc sequence FASTA.
* `--mapping` (optional): Path to original-to-UniParc ID JSON mapping.
* `--output`: Output path for the serialized Python dictionary (`.pkl`).
* `--verify` (optional): Path to a PPI dataset TSV to verify 100% coverage immediately after building.

---

### 4.2 Step 2: Standalone Coverage Verification

To verify any dataset against any sequence database before launching long GPU runs on the HPC:

```bash
python scripts/data_preparation/verify_coverage.py \
    --dataset data/DeepAraPPI/c1_ppi_sample_DeepAraPPI.txt \
    --seq_db data/arabidopsis/uniprot_final.pkl
```

#### Expected Output:
```text
Loading sequence database: data/arabidopsis/uniprot_final.pkl
Loaded 109,640 sequence/accession entries.

Checking dataset: data/DeepAraPPI/c1_ppi_sample_DeepAraPPI.txt
============================================================
Total pairs in dataset:     11,000
Pairs with BOTH sequences:  11,000 (100.00%)
Pairs with missing seqs:    0
Unique missing protein IDs: 0
============================================================

✅ Verification SUCCESS: 100% of protein pairs are present in the sequence database!
```

---

## 5. Downstream Integration with PPLM Inference

Once `uniprot_final.pkl` is built, [`scripts/batch_predict.py`](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/scripts/batch_predict.py) loads the sequence dictionary into memory:

```python
with open(args.seq_db, "rb") as f:
    seq_db = pickle.load(f)

# Instant O(1) sequence retrieval per pair
seqA = seq_db.get(protA, None)
seqB = seq_db.get(protB, None)
```

The sequences are subsequently tokenized, concatenated with inter-chain attention masks, and passed through PPLM's 33-layer Transformer backbone to extract mean/max pooling features for PPI classification.

---

## 6. Extending to Other Plant Species (Rice / Maize)

The pipeline is modular and reusable for future datasets (e.g. Maize PPIM, Rice PRIN):

1. Download the species proteome FASTA from UniProt into `data/<species>/uniprot_<species>.fasta`.
2. Run `build_sequence_db.py`:
   ```bash
   python scripts/data_preparation/build_sequence_db.py \
       --fasta data/rice/uniprot_rice.fasta \
       --output data/rice/uniprot_rice_final.pkl \
       --verify data/DeepAraPPI/all_rice_positive_negative_DeepAraPPI.txt
   ```
3. If any legacy IDs are missing, `verify_coverage.py` will print the exact missing ID list to be retrieved from UniParc.
