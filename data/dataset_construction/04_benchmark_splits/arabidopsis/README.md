# Arabidopsis PPI Benchmark Dataset Splits

This directory houses the finalized, graph-theoretically partitioned datasets for *Arabidopsis thaliana* protein-protein interaction (PPI) prediction at sequence identity threshold $\theta = 0.50$ and reciprocal alignment coverage $c_{\min} = 0.80$.

The datasets are organized into three dedicated subdirectories based on sample class and benchmarking usage:

```
data/dataset_construction/04_benchmark_splits/arabidopsis/
├── positive_data/        # Pure positive interaction pairs & partition summaries
├── negative_data/        # Non-redundant negative interaction pairs, audit logs, & summaries
└── benchmark_data/       # Full benchmark datasets containing both positive and negative samples
```

---

## 1. Directory Structure & File Catalog

### 1.1 Positive Data (`positive_data/`)
Contains non-redundant, high-confidence physical positive interactions ($S_{\text{HIPPIE}} \ge 0.72$) partitioned across C1, C2, C3, and C4.

| File Name | Role | Pair Count | File Size | Description |
| :--- | :--- | :---: | :---: | :--- |
| `ara_positive_c1_train_id0_50_cov0_80.tsv` | C1 Train | 19,462 | 576.6 KB | Both interactors in training set ($u, v \in V(C_1) \subseteq V_{\text{train}}$) |
| `ara_positive_c2_test_id0_50_cov0_80.tsv` | C2 Test | 1,145 | 34.1 KB | Cut boundary test ($u \in V_{\text{train}}, v \in V_{\text{test}}$) |
| `ara_positive_c3_test_id0_50_cov0_80.tsv` | C3 Test | 1,145 | 34.1 KB | Both unseen test ($u, v \in V_{\text{test}}$, zero homology leakage) |
| `ara_positive_c4_test_id0_50_cov0_80.tsv` | C4 Test | 1,145 | 33.9 KB | Both seen test ($u, v \in V(C_1)$, novel test edge) |
| `summary_positive_c1234_partitions.md` | Audit Doc | — | 3.7 KB | Markdown summary of C1/C2/C3/C4 positive partitioning |
| `summary_positive_c1234_partitions.tsv` | Audit Table | — | 419 B | Tabular metrics of positive partitioning |
| `summary_positive_c123_partitions.md` | Legacy Doc | — | 3.7 KB | Parameter sweep report (8 permutations) |
| `summary_positive_c123_partitions.tsv` | Legacy Table | — | 419 B | Tabular parameter sweep metrics |

**TSV Schema (`positive_data/`)**:
```tsv
participant_a_id	participant_b_id	cluster_id	hippie_score	evidence_count
P43273	P93002	1	0.9000	29
```

---

### 1.2 Negative Data (`negative_data/`)
Contains non-redundant negative interaction pairs sampled from disjoint subcellular compartments, filtered against positive interologs ($S \ge 0.50$), and partitioned at a 1:10 ratio with zero leakage.

| File Name | Role | Pair Count | File Size | Description |
| :--- | :--- | :---: | :---: | :--- |
| `ara_negative_c1_train_id0_50_cov0_80.tsv` | C1 Train | 194,620 | 17.5 MB | Both interactors in $V(C_1)$ (1:10 ratio to positive C1) |
| `ara_negative_c2_test_id0_50_cov0_80.tsv` | C2 Test | 11,450 | 1.0 MB | Cut boundary test ($u \in V_{\text{train}}, v \in V_{\text{test}}$) |
| `ara_negative_c3_test_id0_50_cov0_80.tsv` | C3 Test | 10,965 | 985.2 KB | Both unseen test ($u, v \in V_{\text{test}}$, capped at max non-redundant capacity) |
| `ara_negative_c4_test_id0_50_cov0_80.tsv` | C4 Test | 11,450 | 1.0 MB | Both seen test ($u, v \in V(C_1)$, edge-disjoint from C1) |
| `ara_negative_c1234_partitions_audit.json` | Audit JSON | — | 2.5 KB | Programmatic audit verification for all 9 graph integrity criteria |
| `summary_negative_c1234_partitions.md` | Audit Doc | — | 3.2 KB | Negative partition summary report with execution metrics |

**TSV Schema (`negative_data/`)**:
```tsv
participant_a_id	participant_b_id	locus_a	locus_b	compartments_a	compartments_b
F4HZB2	Q9XIE2	AT1G03060	AT1G59870	nucleus|cytosol_cytoplasm	chloroplast_plastid|mitochondrion
```

---

### 1.3 Benchmark Data (`benchmark_data/`)
Full datasets that contain **both positive (`label=1`) and negative (`label=0`) samples** unified into standard machine learning benchmark splits.

| File Name | Role | Positives | Negatives | Total Pairs | Pos:Neg Ratio | File Size |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| `ara_benchmark_c1_train_id0_50_cov0_80.tsv` | C1 Train | 19,462 | 194,620 | 214,082 | 1 : 10.00 | 9.4 MB |
| `ara_benchmark_c2_test_id0_50_cov0_80.tsv` | C2 Test | 1,145 | 11,450 | 12,595 | 1 : 10.00 | 554.1 KB |
| `ara_benchmark_c3_test_id0_50_cov0_80.tsv` | C3 Test | 1,145 | 10,965 | 12,110 | 1 : 9.58 | 532.3 KB |
| `ara_benchmark_c4_test_id0_50_cov0_80.tsv` | C4 Test | 1,145 | 11,450 | 12,595 | 1 : 10.00 | 553.9 KB |
| **Total** | **Full Benchmark** | **22,897** | **228,485** | **251,382** | **1 : 9.98** | **~11.0 MB** |

**TSV Schema (`benchmark_data/`)**:
```tsv
participant_a_id	participant_b_id	label	partition	locus_a	locus_b
P43273	P93002	1	C1_pos	AT5G06950	AT1G64280
F4HZB2	Q9XIE2	0	C1_neg	AT1G03060	AT1G59870
```

---

## 2. Usage Examples

### Loading Benchmark Splits in Python
```python
import pandas as pd

# Load C1 Training benchmark (Pos + Neg)
c1_df = pd.read_csv("data/dataset_construction/04_benchmark_splits/arabidopsis/benchmark_data/ara_benchmark_c1_train_id0_50_cov0_80.tsv", sep="\t")
print(f"C1 Benchmark: {len(c1_df):,} pairs ({sum(c1_df['label'] == 1):,} pos, {sum(c1_df['label'] == 0):,} neg)")

# Load C3 Zero-Shot Generalization test benchmark
c3_df = pd.read_csv("data/dataset_construction/04_benchmark_splits/arabidopsis/benchmark_data/ara_benchmark_c3_test_id0_50_cov0_80.tsv", sep="\t")
print(f"C3 Test: {len(c3_df):,} pairs")
```

### Loading Pure Positives or Negatives
```python
# Pure positive C1 training pairs
pos_c1 = pd.read_csv("data/dataset_construction/04_benchmark_splits/arabidopsis/positive_data/ara_positive_c1_train_id0_50_cov0_80.tsv", sep="\t")

# Pure negative C1 training pairs
neg_c1 = pd.read_csv("data/dataset_construction/04_benchmark_splits/arabidopsis/negative_data/ara_negative_c1_train_id0_50_cov0_80.tsv", sep="\t")
```
