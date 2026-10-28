# Stage 2.6: Non-Redundant Negative Candidate Pool Audit Report

**Execution Timestamp**: `2026-10-27T03:16:02Z`

## 1. Summary Overview

Pair-level redundancy reduction was executed across the clean negative candidate pool using greedy CD-HIT representative selection.

| Metric | Value |
| :--- | :--- |
| **Sequence Identity Threshold ($\theta$)** | `0.50` |
| **Reciprocal Coverage Threshold ($c_{\min}$)** | `0.80` |
| **E-value Cutoff** | `0.001` |
| **Total Clean Candidate Negatives Evaluated** | `101,851,542` |
| **Non-Redundant Representatives Retained** | `50,308,730` (`49.39%`) |
| **Redundant Paralogous Negatives Pruned** | `51,542,812` (`50.61%`) |
| **Unique Proteins Represented** | `15,909` |
| **Homology Component-Pair Signatures** | `45,785,498` |
| **Singleton Component Pairs** | `27,800,435` |
| **Multi-Pair Paralog Groups** | `17,985,063` |
| **Output Filtered Parquet Chunks** | `51` chunks |
| **Total Execution Time** | `627.64` seconds |

## 2. Methodology & Representative Priority Hierarchy

Pairs within each homology component pair are collapsed using exact CD-HIT greedy clustering, prioritized by:
1. **TAIR Subcellular Evidence Confidence**: Direct experimental codes (`EXP`, `IDA`, `IPI`, `HDA` = 3.0) > Curated codes (`TAS`, `NAS` = 2.0) > Computational predictions (`IEA`, `ISS` = 1.0).
2. **Combined Sequence Length**: $(L_A + L_B)$ descending (longer canonical proteins provide richer domain context).
3. **Lexicographical Tie-Breaking**: Deterministic sorting by `participant_a_id`, `participant_b_id`.

## 3. Mathematical Guarantee

All surviving negative pairs strictly satisfy:
$$\forall P_i, P_j \in \mathcal{P}_{\text{neg, NR}} \ (i \ne j): \quad S(P_i, P_j) < 0.50$$
where $S(P_i, P_j) = \max [ \min(s(A, C), s(B, D)), \min(s(A, D), s(B, C)) ]$.
