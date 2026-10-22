# Summary of Arabidopsis Pair Redundancy Clustering Parameter Sweep

| Identity Threshold ($\theta$) | Coverage ($c_{\min}$) | Non-Redundant Pairs | Unique Proteins | Redundancy Removed (%) | Total Clusters | Avg Cluster Size | Max Cluster Size | Min Cluster Size | Multi-Pair Clusters | Singletons | Runtime |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **0.40** | **0.40** | 19,333 | 8,137 | **36.35%** | 19,333 | 1.571 | 155 | 1 | 4,232 | 15,101 | 1.0s |
| **0.40** | **0.50** | 19,939 | 8,177 | **34.36%** | 19,939 | 1.523 | 115 | 1 | 4,161 | 15,778 | 0.9s |
| **0.40** | **0.70** | 20,798 | 8,256 | **31.53%** | 20,798 | 1.461 | 115 | 1 | 4,091 | 16,707 | 0.7s |
| **0.70** | **0.40** | 26,088 | 8,611 | **14.12%** | 26,088 | 1.164 | 52 | 1 | 2,935 | 23,153 | 0.3s |
| **0.70** | **0.50** | 26,170 | 8,618 | **13.85%** | 26,170 | 1.161 | 52 | 1 | 2,904 | 23,266 | 0.3s |
| **0.70** | **0.70** | 26,307 | 8,623 | **13.40%** | 26,307 | 1.155 | 52 | 1 | 2,846 | 23,461 | 0.3s |


*Note: Redundancy removed is defined as (Total Input Pairs - Non-Redundant Pairs) / Total Input Pairs.*

---

## Compact Summary (Core Clustering Metrics)

| Identity Threshold ($\theta$) | Coverage ($c_{\min}$) | Non-Redundant Pairs | Unique Proteins | Redundancy Removed (%) | Total Clusters | Max Cluster Size | Multi-Pair Clusters | Singletons |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **0.40** | **0.40** | 19,333 | 8,137 | **36.35%** | 19,333 | 155 | 4,232 | 15,101 |
| **0.40** | **0.50** | 19,939 | 8,177 | **34.36%** | 19,939 | 115 | 4,161 | 15,778 |
| **0.40** | **0.70** | 20,798 | 8,256 | **31.53%** | 20,798 | 115 | 4,091 | 16,707 |
| **0.70** | **0.40** | 26,088 | 8,611 | **14.12%** | 26,088 | 52 | 2,935 | 23,153 |
| **0.70** | **0.50** | 26,170 | 8,618 | **13.85%** | 26,170 | 52 | 2,904 | 23,266 |
| **0.70** | **0.70** | 26,307 | 8,623 | **13.40%** | 26,307 | 52 | 2,846 | 23,461 |

