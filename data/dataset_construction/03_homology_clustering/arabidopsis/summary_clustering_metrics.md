# Summary of Arabidopsis Pair Redundancy Clustering Parameter Sweep

## 1. Parameter Sweep Results Table

The grid search evaluates 12 combinations across sequence identity ($\theta \in \{0.40, 0.70, 0.90\}$) and alignment coverage ($c_{\min} \in \{0.50, 0.70, 0.80, 0.90\}$) on the **30,376** high-confidence heteromeric positive pairs (comprising **8,970** unique proteins):

| Identity Threshold ($\theta$) | Coverage ($c_{\min}$) | Non-Redundant Pairs | Unique Proteins | Redundancy Removed (%) | Total Clusters | Avg Cluster Size | Max Cluster Size | Min Cluster Size | Multi-Pair Clusters | Singletons | Runtime |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **0.40** | **0.50** | 19,939 | 8,177 | **34.36%** | 19,939 | 1.523 | 115 | 1 | 4,161 | 15,778 | 0.9s |
| **0.40** | **0.70** | 20,798 | 8,256 | **31.53%** | 20,798 | 1.461 | 115 | 1 | 4,091 | 16,707 | 0.7s |
| **0.40** | **0.80** | 21,480 | 8,310 | **29.29%** | 21,480 | 1.414 | 113 | 1 | 3,969 | 17,511 | 0.6s |
| **0.40** | **0.90** | 22,546 | 8,406 | **25.78%** | 22,546 | 1.347 | 62 | 1 | 3,805 | 18,741 | 0.5s |
| **0.70** | **0.50** | 26,170 | 8,618 | **13.85%** | 26,170 | 1.161 | 52 | 1 | 2,904 | 23,266 | 0.3s |
| **0.70** | **0.70** | 26,307 | 8,623 | **13.40%** | 26,307 | 1.155 | 52 | 1 | 2,846 | 23,461 | 0.3s |
| **0.70** | **0.80** | 26,413 | 8,632 | **13.05%** | 26,413 | 1.150 | 52 | 1 | 2,797 | 23,616 | 0.3s |
| **0.70** | **0.90** | 26,807 | 8,653 | **11.75%** | 26,807 | 1.133 | 43 | 1 | 2,618 | 24,189 | 0.3s |
| **0.90** | **0.50** | 29,290 | 8,852 | **3.58%** | 29,290 | 1.037 | 9 | 1 | 878 | 28,412 | 0.2s |
| **0.90** | **0.70** | 29,294 | 8,852 | **3.56%** | 29,294 | 1.037 | 9 | 1 | 874 | 28,420 | 0.2s |
| **0.90** | **0.80** | 29,327 | 8,853 | **3.45%** | 29,327 | 1.036 | 9 | 1 | 842 | 28,485 | 0.2s |
| **0.90** | **0.90** | 29,331 | 8,855 | **3.44%** | 29,331 | 1.036 | 9 | 1 | 839 | 28,492 | 0.2s |

*Note: Redundancy removed is defined as $(N_{\text{total}} - N_{\text{rep}}) / N_{\text{total}} \times 100\%$, where $N_{\text{total}} = 30,376$.*

---

## 2. Parameters Used in the MMseqs2 & Clustering Pipeline

The pipeline consists of two stages: (1) an initial high-sensitivity all-vs-all search on the 8,970 unique proteins in **MMseqs2**, followed by (2) an accelerated pair-level similarity calculation and greedy representative clustering in **Python**.

### A. MMseqs2 All-vs-All Search Parameters (`mmseqs search`)

| Parameter | Value | Definition & Purpose |
|---|---|---|
| `-s` (Sensitivity) | `7.5` | **Search sensitivity**. Controls the k-mer matching threshold and diagonal extension. The MMseqs2 default is `5.7`. Setting `-s 7.5` significantly increases sensitivity to detect remote sequence homologs down to the twilight zone, matching the sensitivity of standard Smith-Waterman / BLASTP while running orders of magnitude faster. |
| `--min-seq-id` | `0.30` | **Minimum sequence identity fraction**. Only protein pairs sharing $\ge 30\%$ sequence identity in their alignment are retained in the base alignment database. This permissive cutoff ensures that all downstream thresholds ($\theta = 0.40, 0.70, 0.90$) are captured in a single run without re-running the alignment search. |
| `-c` (Coverage) | `0.30` | **Minimum alignment coverage fraction**. At least 30% of the sequence must be covered by the local alignment during the initial search phase. |
| `--cov-mode` | `0` | **Coverage calculation mode**. Mode `0` enforces **bidirectional / reciprocal coverage**: $\min(\text{cov}_{\text{query}}, \text{cov}_{\text{target}}) \ge c$. This ensures that neither sequence is a tiny local fragment relative to the other. (Modes: `0` = bidirectional, `1` = target only, `2` = query only). |
| `--threads` | `16` | **Parallel processing**. Utilizes 16 CPU cores for multi-threaded k-mer matching and vectorized Smith-Waterman alignment. |

### B. Single-Protein Homology Filtering Parameters (Python Stage)

Before evaluating pair similarity, alignments from MMseqs2 are evaluated against the specific grid parameters:

| Parameter | Value(s) | Definition & Purpose |
|---|---|---|
| $\theta$ (Identity Threshold) | `0.40`, `0.70`, `0.90` | **Minimum sequence identity fraction ($s(X, Y) \ge \theta$)**: <br>• **`0.40` (Structural / Fold-Level Cutoff)**: Based on the classical Rost / Chothia curve; proteins sharing $\ge 40\%$ identity over their interacting domain share the same 3D tertiary fold and interface geometry. Pruning at 0.40 eliminates homology-based data leakage for ML benchmark evaluation.<br>• **`0.70` (Paralog / Duplicate Gene Cutoff)**: Targets recent gene duplications (e.g. from *Arabidopsis* $\alpha/\beta$ whole-genome duplications). Removes copy-paste duplicate interactions while preserving divergent family members.<br>• **`0.90` (Near-Identity Deduplication)**: Prunes only near-identical alleles, isoforms, or virtually identical duplicate copies ($<4\%$ of dataset). |
| $c_{\min}$ (Coverage Threshold) | `0.50`, `0.70`, `0.80`, `0.90` | **Minimum reciprocal alignment coverage ($\min(\text{qcov}, \text{tcov}) \ge c_{\min}$)**: Requires that the alignment covers at least $c_{\min}$ fraction of *both* proteins. This prevents false-positive homology calls caused by small, promiscuous 20–30 residue motifs (e.g., ATP-binding loops, single zinc-finger motifs) between proteins that have completely different structural scaffolds. |
| `max_evalue` | `1e-3` | **Statistical significance cutoff**. Alignments with an Expectation value $> 10^{-3}$ are rejected to eliminate false-positive matches arising from biased amino acid compositions or random sequence noise. |

### C. Pair-Level Similarity Function $S(P_i, P_j)$

For two interacting pairs $P_i = (A, B)$ and $P_j = (C, D)$, pair similarity is defined as:

$$S(P_i, P_j) = \max \Big[ \underbrace{\min(s(A, C), s(B, D))}_{\text{Orientation 1}}, \underbrace{\min(s(A, D), s(B, C))}_{\text{Orientation 2 (Symmetric)}} \Big]$$

- **The `min` operator** enforces the biological rule that **both** interactors must be homologous. If protein $A \approx C$ (90% identity) but protein $B$ is unrelated to $D$ (0% identity), $\min(0.90, 0.0) = 0.0$. Thus, $S(P_i, P_j) = 0.0$ and the pairs are **not** redundant.
- **The `max` operator** accounts for the undirected physical symmetry of heteromeric interactions ($(A, B) \equiv (B, A)$).

### D. Cluster Representative Selection Hierarchy

Within each redundant cluster, a single representative pair is selected using a deterministic greedy priority hierarchy:

1. **HIPPIE Confidence Score (`score`, descending)**: Experimentally validated, high-confidence interactions ($S \ge 0.72$) with higher confidence scores are retained over lower-scoring counterparts.
2. **Total Evidence Count (`evidence_count`, descending)**: Pairs supported by more independent publications and experimental detection methods are prioritized.
3. **Combined Protein Length ($L_A + L_B$, descending)**: Prefers full-length canonical protein isoforms over shorter fragments or truncated isoforms.
4. **Lexicographical Pair IDs (ascending)**: Deterministic tie-breaker ensuring 100% reproducible cluster representative selection across runs.

---

## 3. Key Observations from the 12-Point Grid

1. **Identity is the primary driver of redundancy removal**:
   - $\theta = 0.90$ collapses **~3.5%** of pairs (near-identical duplicates).
   - $\theta = 0.70$ collapses **~13.5%** of pairs (recent WGD paralog complexes).
   - $\theta = 0.40$ collapses **~30%–34%** of pairs (entire structural fold families).
2. **High Proteome Breadth Preserved**:
   - Even under the strictest setting ($\theta = 0.40, c = 0.50$), **8,177 out of 8,970 unique proteins (91.2%)** are retained. Because pair-level clustering only prunes pairs where *both* interactors are homologous, proteins that interact with distinct partners are preserved.
3. **Coverage Impact at $\theta = 0.40$**:
   - Raising coverage from $0.50 \to 0.70 \to 0.80$ gradually recovers ~700–800 pairs per step while keeping the maximum cluster size stable (~113–115).
   - At $c = 0.90$, the maximum cluster size drops sharply from 113 to 62, because requiring 90% full-length alignment fails to cluster plant paralogs that have diverged terminal intrinsically disordered regions (IDRs).
   - Hence, **$c_{\min} = 0.50$ or $0.70$** remains the most structurally sound domain-level threshold.
