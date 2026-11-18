# Arabidopsis Positive PPI Benchmark C1/C2/C3/C4 Partition Summary

This report documents the rigorous graph-theoretic benchmark partitioning of the high-confidence, non-redundant *Arabidopsis thaliana* positive protein-protein interactions (at sequence identity $\theta = 0.50$ and alignment coverage $c_{\min} = 0.80$) into **C1 (Train, 85%)**, **C2 (One Unseen Test, 5%)**, **C3 (Both Unseen Test, 5%)**, and **C4 (Both Seen Test, 5%)**.

---

## 1. Partition Statistics Table

| Partition | Role | Number of Pairs | % of Total Pairs | Unique Proteins | Vertex Membership Property |
|:---:|:---:|:---:|:---:|:---:|:---|
| **C1** | **Training Set** | **19,462** | **85.00%** | **7,844** | $u, v \in V_{\text{train}}$ (Both interactors seen in train) |
| **C2** | **One Unseen Test** | **1,145** | **5.00%** | **1,069** | $u \in V_{\text{train}}, v \in V_{\text{test}}$ (One seen, one unseen) |
| **C3** | **Both Unseen Test** | **1,145** | **5.00%** | **456** | $u, v \in V_{\text{test}}$ (Both interactors completely unseen; $V(C_3) \cap V(C_1) = \emptyset$) |
| **C4** | **Both Seen Test** | **1,145** | **5.00%** | **1,482** | $u, v \in V(C_1)$ (Both interactors seen in C1; novel test interaction pair) |
| **Total** | **Full NR Dataset** | **22,897** | **100.00%** | **8,426** | Disjoint graph union $V = V_{\text{train}} \cup V_{\text{test}}$ |

---

## 2. Vertex & Graph Partition Properties

- **Total Unique Proteins**: **8,426**
- **Disjoint Train Vertex Set ($|V_{\text{train}}|$)**: **7,970** proteins (94.6%)
- **Disjoint Test Vertex Set ($|V_{\text{test}}|$)**: **456** proteins (5.4%)
- **C2 Cross-Boundary Node Distribution**: 801 proteins from $V_{\text{train}}$ and 268 proteins from $V_{\text{test}}$
- **C4 Seen Node Coverage**: **100.0%** of proteins in C4 (1,482 / 1,482) strictly appear in the C1 training set.
- **C3 Unseen Node Isolation**: **0** proteins in C3 overlap with C1 or C4 ($V(C_3) \cap V(C_1) = \emptyset$).

---

## 3. Strict Verification & Integrity Audit

| Verification Check | Target Requirement | Evaluated Result | Status |
|:---|:---|:---:|:---:|
| **1. Edge Exhaustiveness** | $|C_1| + |C_2| + |C_3| + |C_4| == N_{\text{total}}$ | 19,462 + 1,145 + 1,145 + 1,145 = 22,897 | **PASS** |
| **2. Pairwise Edge Disjointness** | Zero shared interaction pairs across all 4 sets | $\max_{i \neq j} |C_i \cap C_j| = 0$ | **PASS** |
| **3. Vertex Disjointness** | $V_{\text{train}} \cap V_{\text{test}} = \emptyset$ | $|V_{\text{train}} \cap V_{\text{test}}| = 0$ | **PASS** |
| **4. C1 Vertex Containment** | $V(C_1) \subseteq V_{\text{train}}$ | All 7,844 nodes $\in V_{\text{train}}$ | **PASS** |
| **5. C4 Both Seen Requirement** | $\forall (u, v) \in C_4: u, v \in V(C_1)$ | All 1,482 nodes $\in V(C_1)$ (100%) | **PASS** |
| **6. C3 Vertex Containment** | $V(C_3) \subseteq V_{\text{test}}$ | All 456 nodes $\in V_{\text{test}}$ | **PASS** |
| **7. C1 / C3 Zero Protein Overlap** | $V(C_1) \cap V(C_3) = \emptyset$ | $|V(C_1) \cap V(C_3)| = 0$ | **PASS** |
| **8. C4 / C3 Zero Protein Overlap** | $V(C_4) \cap V(C_3) = \emptyset$ | $|V(C_4) \cap V(C_3)| = 0$ | **PASS** |
| **9. C2 Cut Boundary Property** | Exactly one protein $\in V_{\text{train}}$, one $\in V_{\text{test}}$ | 100% of C2 pairs cross cut boundary | **PASS** |

---

## 4. Output Partition Files

All partition files have been written to `data/dataset_construction/04_benchmark_splits/arabidopsis/`:
- `ara_positive_c1_train_id0_50_cov0_80.tsv` (19,462 pairs)
- `ara_positive_c2_test_id0_50_cov0_80.tsv` (1,145 pairs)
- `ara_positive_c3_test_id0_50_cov0_80.tsv` (1,145 pairs)
- `ara_positive_c4_test_id0_50_cov0_80.tsv` (1,145 pairs)
