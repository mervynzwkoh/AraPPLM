# Summary of Arabidopsis Pair Redundancy Clustering Parameter Sweep

The parameter sweep evaluates **5** configuration(s) across sequence identity ($\theta \in \{0.40, 0.50, 0.60, 0.70, 0.80\}$) and alignment coverage ($c_{\min} \in \{0.80\}$) on the **30,376** high-confidence heteromeric positive pairs (comprising **8,970** unique proteins):

| Identity Threshold ($\theta$) | Coverage ($c_{\min}$) | Non-Redundant Pairs | Unique Proteins | Redundancy Removed (%) | Total Clusters | Avg Cluster Size | Max Cluster Size | Min Cluster Size | Multi-Pair Clusters | Singletons | Runtime |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **0.40** | **0.80** | 21,480 | 8,310 | **29.29%** | 21,480 | 1.414 | 113 | 1 | 3,969 | 17,511 | 0.8s |
| **0.50** | **0.80** | 22,897 | 8,426 | **24.62%** | 22,897 | 1.327 | 65 | 1 | 3,866 | 19,031 | 0.5s |
| **0.60** | **0.80** | 24,771 | 8,540 | **18.45%** | 24,771 | 1.226 | 52 | 1 | 3,349 | 21,422 | 0.4s |
| **0.70** | **0.80** | 26,413 | 8,632 | **13.05%** | 26,413 | 1.150 | 52 | 1 | 2,797 | 23,616 | 0.4s |
| **0.80** | **0.80** | 28,139 | 8,737 | **7.36%** | 28,139 | 1.079 | 10 | 1 | 1,797 | 26,342 | 0.3s |


*Note: Redundancy removed is defined as $(N_{\text{total}} - N_{\text{rep}}) / N_{\text{total}} \times 100\%$, where $N_{\text{total}} = 30,376$.*

## Parameters Used in the MMseqs2 & Clustering Pipeline

### A. MMseqs2 All-vs-All Search Parameters (`mmseqs search`)
- **`-s 7.5` (Sensitivity)**: High sensitivity search threshold (MMseqs2 default is 5.7). Captures remote homologs down to the twilight zone with sensitivity matching BLASTP/Smith-Waterman.
- **`--min-seq-id 0.30`**: Minimum sequence identity fraction during alignment search (30%). Permissive cutoff allowing all downstream $\theta$ thresholds to be evaluated from one database.
- **`-c 0.30`**: Minimum alignment coverage fraction during search phase (30%).
- **`--cov-mode 0`**: Reciprocal / bidirectional coverage mode: $\min(\text{cov}_{\text{query}}, \text{cov}_{\text{target}}) \ge c$.
- **`--threads 16`**: 16 CPU cores utilized for multi-threaded alignment.

### B. Pair Similarity & Clustering Parameters
- **$\theta$ (Identity Threshold)**: Minimum sequence identity fraction required for single-protein homology ($s(X, Y) \ge \theta$). In this sensitivity sweep, $\theta \in [0.40, 0.50, 0.60, 0.70, 0.80]$ systematically evaluates the transition from structural fold conservation (0.40) up to near-duplicate paralogs (0.80).
- **$c_{\min}$ (Coverage Threshold)**: Minimum reciprocal alignment coverage required for homology ($\min(\text{qcov}, \text{tcov}) \ge c_{\min}$). Fixed at $0.80$ (80% reciprocal alignment coverage) to prevent short local motif matches (e.g. 20-30 residue ATP loops or zinc fingers) from falsely flagging full-length proteins as homologous.
- **`max_evalue 1e-3`**: Maximum E-value cutoff to ensure sequence alignment statistical significance.
- **Pair Similarity Metric**: $S(P_i, P_j) = \max[\min(s(A,C), s(B,D)), \min(s(A,D), s(B,C))]$. The `min` operator enforces that both interactors must be homologous; `max` accounts for interaction symmetry.
- **Representative Selection Priority**: Greedy priority hierarchy: HIPPIE confidence score (`score` desc) $\rightarrow$ total evidence count (`evidence_count` desc) $\rightarrow$ combined protein length ($L_A + L_B$ desc) $\rightarrow$ lexicographical IDs.
