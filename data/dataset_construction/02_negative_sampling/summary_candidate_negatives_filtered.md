# Stage 2.5: Candidate Negative Homology Filtering Audit Report

**Execution Timestamp**: `2026-10-26T15:18:34Z`

## 1. Summary Overview

Candidate negative interactions were screened against all positive interactions in the latest benchmark splits to eliminate positive interologs and pseudo-negatives.

| Metric | Value |
| :--- | :--- |
| **Sequence Identity Threshold ($\theta$)** | `0.50` |
| **Reciprocal Coverage Threshold ($c_{\min}$)** | `0.80` |
| **E-value Cutoff** | `0.001` |
| **Positive Interactions Screened (C1/C2/C3/C4)** | `22,897` pairs |
| **Precomputed Forbidden Interolog Signatures** | `244,950` pairs |
| **Total Candidate Negatives Evaluated** | `101,884,330` |
| **Positive Interologs Discarded** | `32,788` (`0.03%`) |
| **Clean Non-Redundant Negatives Retained** | `101,851,542` (`99.97%`) |
| **Output Filtered Parquet Chunks** | `102` chunks |
| **Streaming Execution Time** | `252.94` seconds |

## 2. Positive Reference Breakdown

| Split | Purpose | Verified Positive Pairs |
| :--- | :--- | :--- |
| **C1** | Benchmark Positive Set | 19,462 |
| **C2** | Benchmark Positive Set | 1,145 |
| **C3** | Benchmark Positive Set | 1,145 |
| **C4** | Benchmark Positive Set | 1,145 |
| **Total** | Full Non-Redundant Positive Interactome | **22,897** |

## 3. Mathematical Equivalence

The pair similarity metric defined in `docs/dataset_construction/03_homology_clustering/01_pair_redundancy_reduction.md` is:
$$S(P_{\text{neg}}, P_{\text{pos}}) = \max \Big[ \min(s(A, C), s(B, D)), \min(s(A, D), s(B, C)) \Big]$$

A candidate negative pair $(A, B)$ satisfies $S \ge 0.50$ if and only if $(A, B) \in \mathcal{H}(C) \times \mathcal{H}(D)$ or $(A, B) \in \mathcal{H}(D) \times \mathcal{H}(C)$ for some verified positive pair $(C, D)$. All such pairs have been rigorously purged.
