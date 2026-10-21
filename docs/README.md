# AraPPLM Documentation Architecture & Index

This directory serves as the centralized technical knowledge base for **AraPPLM (Arabidopsis Paired Protein Language Model)**. 

Following the codebase sanitization, legacy and unverified dataset construction documentation has been purged. The documentation is organized into clear functional pillars: authoritative raw data verification, sequence database preparation, PPI head training and benchmarking, and core project methodology.

---

## Documentation Structure

```
docs/
├── README.md                              # Master documentation index (this file)
│
├── dataset_construction/                  # [DATASET CONSTRUCTION] Methodology & pipeline starting from verified positives
│   ├── README.md                          # Dataset construction pipeline overview & roadmap
│   └── 01_positive_pairs/                 # [STAGE 1: STARTING POINT] Verified positive PPI numbers from raw data
│       ├── 01_arabidopsis_positive_pairs.md # Arabidopsis thaliana verification (BioGRID, IntAct, STRING, PhoX XL-MS)
│       ├── 02_maize_positive_pairs.md     # Zea mays verification (IntAct, BioGRID, STRING)
│       └── 03_rice_tomato_soybean_positive_pairs.md # Rice, Tomato, and Soybean verification and complex analysis
│
├── data_processing_pipeline.md            # Ingestion, UniParc rescue, and sequence DB creation (100% benchmark coverage)
├── plant-pplm-methodology.md              # End-to-end scientific methodology for plant PPI prediction via PPLM
├── project_overview.md                    # High-level architecture, research objectives, and pipeline stages
├── technical_summary.md                   # Computational optimizations, VRAM management, and benchmark performance
│
├── ppi_head_retraining/                   # Retraining the classification head on plant benchmarks
│   ├── ppi_head_retraining_methodology.md # Architecture, cross-entropy training, hyperparameter choices
│   ├── ppi_head_retraining_walkthrough.md # Step-by-step reproduction and training procedure
│   ├── benchmark_analysis_pretrained_pplm_deeparappi.md # DeepAraPPI benchmark head retraining evaluation
│   └── benchmark_analysis_pretrained_pplm_esmarappi.md  # ESMAraPPI benchmark head retraining evaluation
│
├── zero_shot_benchmarking/                # Zero-shot evaluation of pretrained base PPLM on plant datasets
│   ├── benchmark_analysis_deeparappi_vs_pplm.md # DeepAraPPI zero-shot baseline comparison
│   └── benchmark_analysis_esmarappi_vs_pplm.md  # ESMAraPPI zero-shot baseline comparison
│
└── lit_review/                            # Literature review and survey of existing plant PPI models
```

---

## 1. Verified Positive Pair Documentation (`docs/dataset_construction/01_positive_pairs/`)

These documents represent the **authoritative starting point (Stage 1)** for the fresh dataset construction pipeline. All metrics were audited directly from raw database archives without relying on past documentation.

| Document | Target Species | Verified Raw Unique Pairs | Pure PPI Pairs | Key Data Sources Audited |
| :--- | :--- | :--- | :--- | :--- |
| [`01_arabidopsis_positive_pairs.md`](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/AraPPLM/docs/dataset_construction/01_positive_pairs/01_arabidopsis_positive_pairs.md) | *Arabidopsis thaliana* (`3702`) | **87,673** | **87,673** | BioGRID v5.0.261, IntAct r252, STRING v12.0, PRIDE PhoX XL-MS |
| [`02_maize_positive_pairs.md`](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/AraPPLM/docs/dataset_construction/01_positive_pairs/02_maize_positive_pairs.md) | *Zea mays* (`4577`) | **50,984** | **50,983** | IntAct r252, BioGRID v5.0.261 |
| [`03_rice_tomato_soybean_positive_pairs.md`](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/AraPPLM/docs/dataset_construction/01_positive_pairs/03_rice_tomato_soybean_positive_pairs.md) | *Oryza sativa* (`39947` / `4530`)<br>*Solanum lycopersicum* (`4081`)<br>*Glycine max* (`3847`) | **998** (Rice)<br>**123** (Tomato)<br>**48** (Soybean) | **667** (Rice)<br>**123** (Tomato)<br>**48** (Soybean) | BioGRID v5.0.261, IntAct r252, STRING v12.0 |

### Key Filter Definitions
- **Experimental Physical Interactions**: Enforces MI-compliant physical interaction tags (e.g., `MI:0915 physical association`, `MI:0407 direct interaction`). Proximity labeling (`MI:2364`) and genetic interactions are strictly segregated or excluded.
- **Pure Protein Complexes**: Filters out non-protein participants (DNA `MI:0319`, RNA `MI:0320`, small molecules `MI:0328`, gene IDs).
- **Stringent STRING Direct Evidence**: Retains strictly experimental channel (`experiments > 0`), excluding transferred interologs (`experiments == 0 and experiments_transferred > 0`) and computational co-expression links.

---

## 2. Sequence Preparation & Data Pipeline (`docs/data_processing_pipeline.md`)

Covers the technical pipeline used to achieve **100.00% sequence coverage** across legacy benchmark datasets:
- **Identifier Harmonization**: Regex parsing across complex UniProt / NCBI headers.
- **Legacy ID Rescue**: Fetching retired and secondary accessions from UniParc and UniProtKB archives (`data/arabidopsis/uniprot_final.pkl` and `data/rice/uniprot_rice_final.pkl`).
- **Memory & VRAM Management**: Cropping sequences (`L_A + L_B <= 1020`), dynamically managing PyTorch memory segments, and resilient inference execution on NSCC HPC clusters.

---

## 3. PPI Head Retraining & Benchmarking (`docs/ppi_head_retraining/` & `docs/zero_shot_benchmarking/`)

Comprehensive technical walkthroughs and experimental analyses:
- **Head Architecture**: 2-layer MLP classifier on top of pooled paired language model embeddings.
- **Cross-Validation**: Performance evaluations on DeepAraPPI Tasks (C1, C2, C3, and Rice transfer) and ESMAraPPI Tasks (C2, C3).
- **Zero-Shot Baselines**: Establishing the performance boundaries of base PPLM before plant-specific fine-tuning.

---

## 4. Methodology & Project Architecture

- [`plant-pplm-methodology.md`](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/AraPPLM/docs/plant-pplm-methodology.md): Full scientific motivation, negative pair generation strategies, sequence homology clustering (CD-HIT / MMseqs2), and evaluation metrics.
- [`technical_summary.md`](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/AraPPLM/docs/technical_summary.md): Summary of architectural choices, compute infrastructure (NSCC PBS jobs), and key findings.
- [`project_overview.md`](file:///c:/Users/User/OneDrive/Desktop/NUS/AIS/AIS5281/AraPPLM/docs/project_overview.md): High-level system overview and scope.
