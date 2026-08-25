# Technical Documentation: Sequence Database Preparation & Data Processing Pipeline

## 1. Overview & Objective

In deep learning models for Protein-Protein Interaction (PPI) prediction like **PPLM (Paired Sequence Language Model)**, pairs of protein accessions (e.g., `Q9LUI9` and `F4IDF1`) must be converted into their corresponding full-length amino acid sequences before tokenization and feature extraction.

The primary objective of this data processing pipeline is to construct fast, indexed, high-coverage sequence lookup databases (`data/arabidopsis/uniprot_final.pkl` and `data/rice/uniprot_rice_final.pkl`) that guarantee **100% sequence coverage** across all plant benchmark suites (**DeepAraPPI** and **ESMAraPPI**), and to execute PPLM inference efficiently within HPC hardware constraints.

```
+-----------------------------------------------------------------------------------+
|                           PROTEOME SEQUENCE INGESTION                             |
|  UniProt Arabidopsis Reference (54k seqs) + UniProt Rice Proteome (50k seqs)      |
+----------------------------------------+------------------------------------------+
                                         |
                                         v  [Step 1: Regex Accession Indexing]
+----------------------------------------+------------------------------------------+
|  Initial Lookup Index (~90-97% Base Coverage across PPI Benchmark Sets)          |
+--------------------+-------------------+--------------------+---------------------+
                     |                                        |
                     v  [Step 2: Legacy TrEMBL Rescue]        v  [Step 3: sec_acc Resolution]
+--------------------+-------------------+  +-----------------+---------------------+
|  UniParc Sequence Retrieval            |  |  UniProtKB Secondary Accession Query  |
|  - Arabidopsis: 348 legacy IDs rescued |  |  - ESMAraPPI: 7 secondary IDs resolved|
|  - Rice: 163 subspecies IDs rescued    |  |  - FASTA generated & merged           |
+--------------------+-------------------+  +-----------------+---------------------+
                     |                                        |
                     +-------------------+--------------------+
                                         |
                                         v  [Step 4: Unified Indexed PKL Databases]
+----------------------------------------+------------------------------------------+
|  FINAL SEQUENCE DATABASES (100.00% GUARANTEED COVERAGE)                           |
|  - data/arabidopsis/uniprot_final.pkl (109,994 entries, 100% on DeepAra & ESM)   |
|  - data/rice/uniprot_rice_final.pkl   (100,297 entries, 100% on Rice Transfer)   |
+----------------------------------------+------------------------------------------+
                                         |
                                         v  [Step 5: HPC-Optimized Inference & Eval]
+----------------------------------------+------------------------------------------+
|  PPLM BATCH INFERENCE & BENCHMARK EVALUATION ENGINE                               |
|  - Combined Cropping (LA + LB <= 1020) | - PyTorch Expandable Segments Allocator  |
|  - Instant VRAM Deallocation           | - Resilient OOM Catch-and-Retry Handler  |
|  ================================================================================ |
|  DeepAraPPI Suite (Task 1, 2, 3, Rice) | ESMAraPPI Suite (Task C2, C3)           |
+-----------------------------------------------------------------------------------+
```

---

## 2. Key Sequence Database Challenges & Technical Solutions

### Challenge 1: Header & Accession Format Mismatches
* **The Problem:** 
  PPI benchmark datasets reference proteins using bare UniProt accessions (`Q9LUI9`, `F4IDF1`), whereas standard UniProt FASTA files contain compound headers (`>sp|Q9LUI9|PIAL1_ARATH`, `>tr|A0A178VKB9|...`).
* **The Solution:**
  The parser utilizes regex matching ($\verb|^(?:sp|tr)\|([A-Z0-9]+)\||$) to register both raw compound headers and stripped bare accessions simultaneously into a Python dictionary, enabling instantaneous $O(1)$ lookups regardless of query format.

---

### Challenge 2: The Arabidopsis "Missing 2.4%" Legacy/TrEMBL Discrepancy
* **The Problem:**
  When querying the initial UniProt *Arabidopsis* proteome against the 130,478 pair DeepAraPPI dataset, **2.4% of pairs (1,348 pairs across 348 unique protein IDs)** failed to resolve due to retired or merged TrEMBL accessions from older database releases (BioGRID, DIP, TAIR, MINT).
* **The Solution (UniParc Cross-Referencing & Rescue):**
  1. Unresolved accessions were isolated into `data/arabidopsis/missing_ids.txt`.
  2. Rescued matching sequence records from the UniProt Archive into `data/arabidopsis/uniparc_348.fasta`.
  3. Mapped accessions via `data/arabidopsis/uniparc_id_mapping.json` into the sequence dictionary.
* **Result:** **100.00% coverage (130,478 / 130,478 pairs).**

---

### Challenge 3: The Rice (*Oryza sativa*) Subspecies Coverage Gap
* **The Problem:**
  The initial *Oryza sativa subsp. japonica* proteome download (Taxonomy: `39947`) missed **10.15% of pairs (682 pairs across 163 unique protein IDs)** in the Rice benchmark because DeepAraPPI contains interactions across multiple rice cultivars and subspecies (including *indica* `39946` and general *Oryza sativa* `4530`).
* **The Solution (UniParc Cross-Species Rescue):**
  1. [`scripts/data_preparation/verify_coverage.py`](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/scripts/data_preparation/verify_coverage.py) extracted the 163 missing accessions to `data/rice/missing_rice_ids.txt`.
  2. Retrieved matching UniParc sequences into `data/rice/missing_rice.fasta`.
  3. Merged into `data/rice/uniprot_rice_final.pkl` using [`scripts/data_preparation/merge_rice_uniparc.py`](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/scripts/data_preparation/merge_rice_uniparc.py).
* **Result:** **100.00% coverage (6,721 / 6,721 pairs).**

---

### Challenge 4: ESMAraPPI Secondary Accession Resolution (`sec_acc`)
* **The Problem:**
  When querying `data/arabidopsis/uniprot_final.pkl` against the ESMAraPPI benchmark datasets ([`c2Pred.txt`](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/data/ESMAraPPI/c2Pred.txt) and [`c3Pred.txt`](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/data/ESMAraPPI/c3Pred.txt)), **7 unique protein IDs** (`P25069`, `P25854`, `P29512`, `P59263`, `Q9FDW0`, `Q9LUF3`, `Q9ZNT9`) failed to resolve directly because they were historical secondary accessions merged into newer primary entries in UniProtKB.
* **The Solution (Secondary Accession Mapping via REST API):**
  Queried the UniProtKB search API using the secondary accession operator (`https://rest.uniprot.org/uniprotkb/search?query=sec_acc:{acc}`) to map each ID to its active primary sequence:
  * `P25069` $\rightarrow$ `P0DH97` (Calmodulin-2, 149 aa)
  * `P25854` $\rightarrow$ `P0DH95` (Calmodulin-1, 149 aa)
  * `P29512` $\rightarrow$ `Q56YW9` (Tubulin beta-2 chain, 450 aa)
  * `P59263` $\rightarrow$ `B9DHA6` (Ubiquitin-ribosomal protein eL40z fusion, 128 aa)
  * `Q9FDW0` $\rightarrow$ `P0DH90` (Protein FRIGIDA, 609 aa)
  * `Q9LUF3` $\rightarrow$ `P0DI12` (SUMO-activating enzyme subunit 1B-1, 320 aa)
  * `Q9ZNT9` $\rightarrow$ `P0DKJ8` (Polycomb group protein FIS2, 755 aa)
  Sequences were saved to `data/ESMAraPPI/missing_esmarappi.fasta` and merged into `data/arabidopsis/uniprot_final.pkl`.
* **Result:** **100.00% coverage across all 37,444 pairs in C2 and all 8,866 pairs in C3.**

---

## 3. Benchmark Dataset Architectures & Partition Schemes

### 3.1 DeepAraPPI Benchmark Suite (*Zheng et al., The Plant Journal 2023*)
* **Positive Data Source:** BioGRID, DIP, IntAct, MINT, TAIR with HIPPIE quality score $\ge 0.72$.
* **Negative Data Generation:** Random pairing from reference proteome pool with non-overlapping subcellular localization.
* **Class Ratio:** 1:10 positive-to-negative ratio.
* **Partitions (Park & Marcotte 2012 scheme):**

| Dataset File | Positive Pairs | Negative Pairs | Total Pairs | Class Ratio | Difficulty Level & Evaluation Objective |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `c1_ppi_sample_DeepAraPPI.txt` | 2,844 | 28,440 | **31,284** | 1 : 10 | **Task 1 (Low Difficulty):** Standard random 80/20 test split (seen domain). |
| `c2_ppi_sample_DeepAraPPI.txt` | 6,005 | 60,050 | **66,055** | 1 : 10 | **Task 2 (Medium Difficulty):** One protein in each test pair is unseen. |
| `c3_ppi_sample_DeepAraPPI.txt` | 3,009 | 30,090 | **33,099** | 1 : 10 | **Task 3 (High Difficulty):** Both proteins in test pair are unseen (zero-shot). |
| `total_positive_negative_samples_DeepAraPPI.txt` | 11,858 | 118,620 | **130,478** | 1 : 10 | **Full Interactome:** Complete filtered Arabidopsis interactome. |
| `all_rice_positive_negative_DeepAraPPI.txt` | 611 | 6,110 | **6,721** | 1 : 10 | **Task 4 (Cross-Species Transfer):** Non-redundant curated *Oryza sativa* interactome. |

---

### 3.2 ESMAraPPI Benchmark Suite (*Zhou et al., Plant Methods 2023*)
* **Positive Data Source:** IntAct physical interactions with MIscore $\ge 0.45$.
* **Negative Data Generation & Redundancy Filtering:**
  * Candidate negative proteins filtered at a strict **40% sequence identity threshold** relative to positive proteins, plus 40% internal redundancy reduction, retaining 8,382 non-homologous background proteins.
  * Random pairing generated at controlled 1:10 ratio (77,290 negative pairs total).
* **Partitions (Park & Marcotte 2012 scheme):**

| Dataset File | Positive Pairs | Negative Pairs | Total Pairs | Class Ratio | Difficulty Level & Evaluation Objective |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `c1` (Training Set) | 3,519 | 35,190 | **38,709** | 1 : 10 | **Training Baseline:** Used for supervised training in ESMAraPPI study. |
| `c2Pred.txt` (Test Set) | 3,404 | 34,040 | **37,444** | 1 : 10 | **Task C2 (Medium Difficulty):** Exactly one protein in each test pair is unseen. |
| `c3Pred.txt` (Test Set) | 806 | 8,060 | **8,866** | 1 : 10 | **Task C3 (High Difficulty):** Both proteins in test pair are completely unseen. |

---

## 4. Summary of Sequence Databases & Multi-Dataset Coverage

| Sequence Database File | Organism | Accession / Mapping Entries | Verified Benchmark Datasets | Coverage Rate | Status |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `data/arabidopsis/uniprot_final.pkl` | *Arabidopsis thaliana* | **109,994** entries | • DeepAraPPI C1 (31,284 pairs)<br>• DeepAraPPI C2 (66,055 pairs)<br>• DeepAraPPI C3 (33,099 pairs)<br>• DeepAraPPI Total (130,478 pairs)<br>• ESMAraPPI C2 (37,444 pairs)<br>• ESMAraPPI C3 (8,866 pairs) | **100.00%** (0 missing pairs) | ✅ Complete & Verified |
| `data/rice/uniprot_rice_final.pkl` | *Oryza sativa* | **100,297** entries | • DeepAraPPI Rice (6,721 pairs) | **100.00%** (0 missing pairs) | ✅ Complete & Verified |

---

## 5. HPC Optimization & Model Execution Design Decisions

To evaluate PPLM (a 33-layer, 650M Transformer) across large plant interactomes (~106k held-out DeepAraPPI pairs + ~46k ESMAraPPI pairs) within the physical constraints of the NSCC Aspire 2A A100 GPU (40GB VRAM) and PBS walltime limits, the following engineering design decisions are implemented across [`scripts/batch_predict.py`](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/scripts/batch_predict.py) and PBS batch scripts:

### Decision 1: Combined Sequence Pair Length Cropping ($L_A + L_B \le 1020$)
* **The Problem:** 
  Plant proteins can exceed 2,000–3,000 residues. In PPLM, feature extraction requires retaining the full $[33, 20, L, L]$ attention tensor. For uncropped plant pairs where $L \approx 4,000$, a single pair requires:
  $$33 \text{ layers} \times 20 \text{ heads} \times (4000 \times 4000) \times 4 \text{ bytes} \approx \mathbf{42.2 \text{ GB VRAM}}$$
  This instantly caused `torch.OutOfMemoryError` on 40GB A100 GPUs.
* **The Justified Solution (Matching PPLM Pretraining):**
  The original PPLM paper explicitly specifies:
  > *"To accommodate GPU memory limitations, sequence pairs with full lengths exceeding 1024 residues are cropped."*
  We implemented **proportional budget cropping** that enforces $L_A + L_B \le 1020$ (leaving 4 tokens for `<cls>`, `<sep>`, and boundary tokens):
  $$\text{Budget}_A = \max\left(50, \left\lfloor 1020 \times \frac{\text{len}(\text{seqA})}{\text{len}(\text{seqA}) + \text{len}(\text{seqB})} \right\rfloor\right), \quad \text{Budget}_B = 1020 - \text{Budget}_A$$
  $$\text{seqA} = \text{seqA}[:\text{Budget}_A], \quad \text{seqB} = \text{seqB}[:\text{Budget}_B]$$
* **Impact:** Guarantees total input tokens $\le 1024$. Peak attention tensor VRAM dropped from $>42\text{ GB}$ to $\mathbf{\le 2.6\text{ GB}}$, ensuring $100\%$ zero-OOM execution across all evaluated pairs.

---

### Decision 2: Instant VRAM Deallocation (`del` + Periodic Garbage Collection)
* **The Problem:** 
  PyTorch caches intermediate tensor blocks in VRAM. Leaving large 33-layer attention tensors ($[33, 20, L, L]$) in memory across thousands of consecutive pairs causes severe memory fragmentation.
* **The Solution:**
  1. Immediately after extracting the 660-dimensional mean-pooled and max-pooled attention vectors, all raw bulky tensors are explicitly deleted:
     ```python
     del out, attns, attn_AA, attn_AB, attn_BA, attn_BB, inter_attn, embed_A, embed_B, tokens, inter_chain_mask
     ```
  2. `torch.cuda.empty_cache()` is called periodically every 200 pairs to flush freed blocks back to the GPU allocator.

---

### Decision 3: Autograd Suppression & Tensor Detachment (`torch.no_grad()` + `.detach()`)
* **The Problem:** 
  During inference, PyTorch builds dynamic computation graphs by default, doubling VRAM consumption and throwing `RuntimeError: Can't call numpy() on Tensor that requires grad`.
* **The Solution:**
  1. Wrapped all feature extraction and 10-fold classifier forward passes in `with torch.no_grad():`.
  2. Safely detached prediction tensors before converting to NumPy floats:
     ```python
     final_score = torch.mean(predictions, dim=0).squeeze().detach().cpu().numpy()
     ```

---

### Decision 4: CUDA Memory Allocator Optimization
* **The Problem:** 
  Evaluating tens of thousands of variable-length sequence pairs causes thousands of dynamic allocations of different sizes, fragmenting VRAM over hours of execution.
* **The Solution:**
  Added the following environment variable to all PBS scripts:
  ```bash
  export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
  ```
  This instructs the PyTorch Caching Allocator to create expandable virtual memory segments instead of failing when large contiguous physical blocks are temporarily unavailable.

---

### Decision 5: Resilient OOM Catch-and-Retry Fallback
* **The Problem:** 
  If an edge-case sequence pair caused a transient memory spike, a standard script would crash, aborting a multi-hour PBS job and wasting GPU allocation hours.
* **The Solution:**
  Wrapped the prediction call in a two-stage fallback:
  ```python
  try:
      features = get_features(model, batch_converter, seqA, seqB, device, max_pair_len=args.max_pair_len)
      score = predict_with_weights(ppi_model, features, ppi_weights)
  except torch.cuda.OutOfMemoryError:
      torch.cuda.empty_cache()
      try:
          # Retry with tighter sequence budget (800 residues)
          features = get_features(model, batch_converter, seqA, seqB, device, max_pair_len=800)
          score = predict_with_weights(ppi_model, features, ppi_weights)
      except Exception as e:
          print(f"WARNING: Skipping pair {protA}:{protB} ({e})")
  ```

---

### Decision 6: Unbuffered Real-Time HPC Streaming (`python -u`)
* **The Problem:** 
  PBS Pro and Python default to block-buffering stdout when output is redirected to a log file. Users could not see live progress via `tail -f *.log`.
* **The Solution:**
  Added the unbuffered flag `-u` to all Python executions in the PBS scripts (`python -u scripts/batch_predict.py ...`), enabling real-time streaming of pairs-per-second progress.

---

### Decision 7: Modular Per-Dataset & Per-Task Output Decoupling
* **The Problem:** 
  Standard PBS normal queues enforce walltime limits (e.g. 3 to 4 hours). A monolithic script that saves results only at the end risks losing all computed predictions if walltime expires during downstream tasks.
* **The Solution:**
  1. Outputs are strictly isolated into dataset subdirectories (`results/DeepAraPPI/` and `results/ESMAraPPI/`).
  2. Each task writes its predictions independently to CSV files as soon as that task completes, allowing modular resumption.

---

## 6. Pipeline Execution Reference

### 6.1 Building & Updating Sequence Databases
```bash
# 1. Build initial Arabidopsis database with UniParc legacy rescue
python scripts/data_preparation/build_sequence_db.py \
    --fasta data/arabidopsis/uniprot_arabidopsis.fasta \
    --uniparc_fasta data/arabidopsis/uniparc_348.fasta \
    --mapping data/arabidopsis/uniparc_id_mapping.json \
    --output data/arabidopsis/uniprot_final.pkl

# 2. Merge ESMAraPPI secondary accession sequences into Arabidopsis database
python scripts/data_preparation/merge_esmarappi_sequences.py \
    --db data/arabidopsis/uniprot_final.pkl \
    --fasta data/ESMAraPPI/missing_esmarappi.fasta \
    --output data/arabidopsis/uniprot_final.pkl

# 3. Build Rice database with UniParc rescue
python scripts/data_preparation/build_sequence_db.py \
    --fasta data/rice/uniprot_rice.fasta \
    --output data/rice/uniprot_rice_final.pkl

python scripts/data_preparation/merge_rice_uniparc.py \
    --db data/rice/uniprot_rice_final.pkl \
    --missing_ids data/rice/missing_rice_ids.txt \
    --fasta data/rice/missing_rice.fasta \
    --output data/rice/uniprot_rice_final.pkl
```

---

### 6.2 Verifying Sequence Coverage for Any Dataset
```bash
# Verify DeepAraPPI datasets
python scripts/data_preparation/verify_coverage.py \
    --dataset data/DeepAraPPI/c1_ppi_sample_DeepAraPPI.txt \
    --seq_db data/arabidopsis/uniprot_final.pkl

# Verify ESMAraPPI datasets
python scripts/data_preparation/verify_coverage.py \
    --dataset data/ESMAraPPI/c2Pred.txt \
    --seq_db data/arabidopsis/uniprot_final.pkl

python scripts/data_preparation/verify_coverage.py \
    --dataset data/ESMAraPPI/c3Pred.txt \
    --seq_db data/arabidopsis/uniprot_final.pkl
```

---

### 6.3 Executing Benchmarks on NSCC HPC

```bash
# 1. Submit DeepAraPPI Benchmark Suite (Tasks 1, 2, 3, and Rice)
qsub scripts/run_all_benchmarks_nscc.pbs

# 2. Submit ESMAraPPI Benchmark Suite (Tasks C2 and C3)
qsub scripts/run_esmarappi_benchmarks_nscc.pbs

# 3. Monitor live execution
tail -f pplm_esmarappi_benchmarks.log
```

---

### 6.4 Standalone Evaluation & Baseline Comparisons
```bash
# Evaluate DeepAraPPI results
python scripts/evaluate_pplm.py \
    --c1 results/DeepAraPPI/deepara_c1_scores.csv \
    --c2 results/DeepAraPPI/deepara_c2_scores.csv \
    --c3 results/DeepAraPPI/deepara_c3_scores.csv \
    --rice results/DeepAraPPI/deepara_rice_scores.csv \
    --output_dir results/DeepAraPPI

# Evaluate ESMAraPPI results against paper baselines (ESMAraPPI, TAGPPI, PIPR, D-SCRIPT)
python scripts/evaluate_pplm.py \
    --esm_c2 results/ESMAraPPI/esmarappi_c2_scores.csv \
    --esm_c3 results/ESMAraPPI/esmarappi_c3_scores.csv \
    --output_dir results/ESMAraPPI
```
