# Arabidopsis Positive PPI C1/C2/C3 Partition Benchmark Report

## 1. Overview

This document summarizes the **cluster-level protein partitioning** of non-redundant positive protein-protein interactions across all **8 parameter sweeps** ($\theta \in \{0.40, 0.70\} \times c_{\min} \in \{0.50, 0.70, 0.80, 0.90\}$).

Partitions strictly satisfy the classical Park & Marcotte (2012) evaluation scheme:
- **C1 (Train, ~80%)**: Both interactors in $V_{\text{train}}$ ($u, v \in V_{\text{train}}$).
- **C2 (One Unseen Test, ~10%)**: Exactly one interactor unseen ($u \in V_{\text{train}}, v \in V_{\text{test}}$ or vice versa).
- **C3 (Both Unseen Test, ~10%)**: Both interactors unseen in training ($u, v \in V_{\text{test}}$).
- **Strict Zero-Homology / Zero-Protein Leakage**: $V_{\text{train}} \cap V_{\text{test}} = \emptyset \implies V_{\text{C1}} \cap V_{\text{C3}} = \emptyset$.

---

## 2. Partition Summary Table Across All 8 Permutations

| Identity ($\theta$) | Coverage ($c_{\min}$) | Total Pairs | Total Proteins | C1 Pairs | C1 % | C1 Proteins | C2 Pairs | C2 % | C2 Proteins | C3 Pairs | C3 % | C3 Proteins | C1/C3 Overlap | Verification |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **0.40** | **0.50** | 19,939 | 8,177 | 15,951 | **80.00%** | 6,762 | 1,994 | **10.00%** | 1,802 | 1,994 | **10.00%** | 1,298 | **0** | **PASS** |
| **0.40** | **0.70** | 20,798 | 8,256 | 16,639 | **80.00%** | 7,215 | 2,079 | **10.00%** | 1,751 | 2,080 | **10.00%** | 753 | **0** | **PASS** |
| **0.40** | **0.80** | 21,480 | 8,310 | 17,184 | **80.00%** | 6,858 | 2,148 | **10.00%** | 1,730 | 2,148 | **10.00%** | 1,385 | **0** | **PASS** |
| **0.40** | **0.90** | 22,546 | 8,406 | 18,037 | **80.00%** | 6,701 | 2,255 | **10.00%** | 2,030 | 2,254 | **10.00%** | 1,502 | **0** | **PASS** |
| **0.70** | **0.50** | 26,170 | 8,618 | 20,937 | **80.00%** | 7,233 | 2,616 | **10.00%** | 1,987 | 2,617 | **10.00%** | 1,236 | **0** | **PASS** |
| **0.70** | **0.70** | 26,307 | 8,623 | 21,045 | **80.00%** | 6,879 | 2,631 | **10.00%** | 2,058 | 2,631 | **10.00%** | 1,689 | **0** | **PASS** |
| **0.70** | **0.80** | 26,413 | 8,632 | 21,131 | **80.00%** | 7,294 | 2,641 | **10.00%** | 1,958 | 2,641 | **10.00%** | 1,120 | **0** | **PASS** |
| **0.70** | **0.90** | 26,807 | 8,653 | 21,445 | **80.00%** | 6,908 | 2,681 | **10.00%** | 2,116 | 2,681 | **10.00%** | 1,713 | **0** | **PASS** |

---

## 3. Detailed Protein Allocation & Subgraph Topology

| Identity ($\theta$) | Coverage ($c_{\min}$) | Total Nodes ($|V|$) | Train Nodes ($|V_{\text{train}}|$) | Train % | Test Nodes ($|V_{\text{test}}|$) | Test % | C2 Train Prots | C2 Test Prots | C3 Internal Density | Runtime |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **0.40** | **0.50** | 8,177 | 6,857 | **83.9%** | 1,320 | **16.1%** | 1,071 | 731 | 1.54 pairs/prot | 0.35s |
| **0.40** | **0.70** | 8,256 | 7,503 | **90.9%** | 753 | **9.1%** | 1,452 | 299 | 2.76 pairs/prot | 0.64s |
| **0.40** | **0.80** | 8,310 | 6,918 | **83.2%** | 1,392 | **16.8%** | 1,015 | 715 | 1.55 pairs/prot | 0.58s |
| **0.40** | **0.90** | 8,406 | 6,896 | **82.0%** | 1,510 | **18.0%** | 1,274 | 756 | 1.50 pairs/prot | 1.04s |
| **0.70** | **0.50** | 8,618 | 7,368 | **85.5%** | 1,250 | **14.5%** | 1,426 | 561 | 2.12 pairs/prot | 1.21s |
| **0.70** | **0.70** | 8,623 | 6,910 | **80.1%** | 1,713 | **19.9%** | 1,093 | 965 | 1.56 pairs/prot | 1.10s |
| **0.70** | **0.80** | 8,632 | 7,508 | **87.0%** | 1,124 | **13.0%** | 1,499 | 459 | 2.36 pairs/prot | 3.32s |
| **0.70** | **0.90** | 8,653 | 6,938 | **80.2%** | 1,715 | **19.8%** | 1,166 | 950 | 1.57 pairs/prot | 1.09s |

---

## 4. Verification Methodology

Each generated dataset underwent the following automated validation checks:

1. **Edge Exhaustiveness**: $N_{\text{C1}} + N_{\text{C2}} + N_{\text{C3}} = N_{\text{total}}$ (Verified: 100% of non-redundant pairs accounted for).
2. **Edge Disjointness**: $E_{\text{C1}} \cap E_{\text{C2}} = \emptyset$, $E_{\text{C1}} \cap E_{\text{C3}} = \emptyset$, $E_{\text{C2}} \cap E_{\text{C3}} = \emptyset$.
3. **Disjoint Protein Sets**: $V_{\text{train}} \cap V_{\text{test}} = \emptyset$ across all proteomic nodes.
4. **C1 Validity**: Every interaction in C1 connects two proteins strictly belonging to $V_{\text{train}}$.
5. **C2 Validity**: Every interaction in C2 connects exactly one protein in $V_{\text{train}}$ and one protein in $V_{\text{test}}$.
6. **C3 Validity**: Every interaction in C3 connects two proteins strictly belonging to $V_{\text{test}}$.
7. **Zero Protein Leakage**: The intersection of unique proteins in C1 (Train) and C3 (Test) is **strictly empty** ($|V_{\text{C1}} \cap V_{\text{C3}}| = 0$).

---

## 5. Output File Inventory

All partition TSV files have been exported to `data/dataset_construction/04_benchmark_splits/arabidopsis/`:

### Sweep $\theta = 0.40, c_{\min} = 0.50$
- C1 Train TSV: `ara_positive_c1_train_id0_40_cov0_50.tsv` (15,951 pairs, 6,762 proteins)
- C2 Test TSV:  `ara_positive_c2_test_id0_40_cov0_50.tsv` (1,994 pairs, 1,802 proteins)
- C3 Test TSV:  `ara_positive_c3_test_id0_40_cov0_50.tsv` (1,994 pairs, 1,298 proteins)

### Sweep $\theta = 0.40, c_{\min} = 0.70$
- C1 Train TSV: `ara_positive_c1_train_id0_40_cov0_70.tsv` (16,639 pairs, 7,215 proteins)
- C2 Test TSV:  `ara_positive_c2_test_id0_40_cov0_70.tsv` (2,079 pairs, 1,751 proteins)
- C3 Test TSV:  `ara_positive_c3_test_id0_40_cov0_70.tsv` (2,080 pairs, 753 proteins)

### Sweep $\theta = 0.40, c_{\min} = 0.80$
- C1 Train TSV: `ara_positive_c1_train_id0_40_cov0_80.tsv` (17,184 pairs, 6,858 proteins)
- C2 Test TSV:  `ara_positive_c2_test_id0_40_cov0_80.tsv` (2,148 pairs, 1,730 proteins)
- C3 Test TSV:  `ara_positive_c3_test_id0_40_cov0_80.tsv` (2,148 pairs, 1,385 proteins)

### Sweep $\theta = 0.40, c_{\min} = 0.90$
- C1 Train TSV: `ara_positive_c1_train_id0_40_cov0_90.tsv` (18,037 pairs, 6,701 proteins)
- C2 Test TSV:  `ara_positive_c2_test_id0_40_cov0_90.tsv` (2,255 pairs, 2,030 proteins)
- C3 Test TSV:  `ara_positive_c3_test_id0_40_cov0_90.tsv` (2,254 pairs, 1,502 proteins)

### Sweep $\theta = 0.70, c_{\min} = 0.50$
- C1 Train TSV: `ara_positive_c1_train_id0_70_cov0_50.tsv` (20,937 pairs, 7,233 proteins)
- C2 Test TSV:  `ara_positive_c2_test_id0_70_cov0_50.tsv` (2,616 pairs, 1,987 proteins)
- C3 Test TSV:  `ara_positive_c3_test_id0_70_cov0_50.tsv` (2,617 pairs, 1,236 proteins)

### Sweep $\theta = 0.70, c_{\min} = 0.70$
- C1 Train TSV: `ara_positive_c1_train_id0_70_cov0_70.tsv` (21,045 pairs, 6,879 proteins)
- C2 Test TSV:  `ara_positive_c2_test_id0_70_cov0_70.tsv` (2,631 pairs, 2,058 proteins)
- C3 Test TSV:  `ara_positive_c3_test_id0_70_cov0_70.tsv` (2,631 pairs, 1,689 proteins)

### Sweep $\theta = 0.70, c_{\min} = 0.80$
- C1 Train TSV: `ara_positive_c1_train_id0_70_cov0_80.tsv` (21,131 pairs, 7,294 proteins)
- C2 Test TSV:  `ara_positive_c2_test_id0_70_cov0_80.tsv` (2,641 pairs, 1,958 proteins)
- C3 Test TSV:  `ara_positive_c3_test_id0_70_cov0_80.tsv` (2,641 pairs, 1,120 proteins)

### Sweep $\theta = 0.70, c_{\min} = 0.90$
- C1 Train TSV: `ara_positive_c1_train_id0_70_cov0_90.tsv` (21,445 pairs, 6,908 proteins)
- C2 Test TSV:  `ara_positive_c2_test_id0_70_cov0_90.tsv` (2,681 pairs, 2,116 proteins)
- C3 Test TSV:  `ara_positive_c3_test_id0_70_cov0_90.tsv` (2,681 pairs, 1,713 proteins)

