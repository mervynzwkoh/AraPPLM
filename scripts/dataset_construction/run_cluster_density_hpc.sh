#!/usr/bin/env bash
# ==============================================================================
# run_cluster_density_hpc.sh
# ==============================================================================
# Calculates the intra-cluster clique density metric f_k(theta), pair similarity,
# and representative-to-member cohesion for the 20 largest clusters across:
#   - Fixed sequence identity threshold: theta = 0.50
#   - Alignment coverage sweep: cmin in [0.50, 0.60, 0.70, 0.80, 0.90]
#
# If mapping files do not exist yet for theta=0.50, they will be auto-generated
# on-the-fly using the precomputed MMseqs2 alignments.
#
# Usage:
#   # Run with defaults (theta=0.50, cmins in 0.50 0.60 0.70 0.80 0.90, top 20 clusters):
#   bash scripts/dataset_construction/run_cluster_density_hpc.sh
#
#   # Override via environment variables if desired:
#   THETAS="0.50" CMINS="0.50 0.60 0.70 0.80 0.90" bash scripts/dataset_construction/run_cluster_density_hpc.sh
# ==============================================================================

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_ROOT}"

echo "========================================================================"
echo "    AraPPLM: Evaluating Intra-Cluster Pair Similarity & Clique Density  "
echo "========================================================================"

# Sourcing conda/python if available (disable nounset to avoid system /etc/bashrc unbound variable errors)
set +u 2>/dev/null || true
if [ -f "$HOME/.bashrc" ]; then
    source "$HOME/.bashrc" 2>/dev/null || true
fi
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh" 2>/dev/null || true
elif [ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniforge3/etc/profile.d/conda.sh" 2>/dev/null || true
elif [ -f "/app1/common/miniconda3/etc/profile.d/conda.sh" ]; then
    source "/app1/common/miniconda3/etc/profile.d/conda.sh" 2>/dev/null || true
fi

conda activate python11 2>/dev/null || true
for py_dir in \
    "$HOME/.conda/envs/python11/bin" \
    "$HOME/miniconda3/envs/python11/bin" \
    "$HOME/miniforge3/envs/python11/bin" \
    "/app1/common/miniconda3/envs/python11/bin"; do
    if [ -f "${py_dir}/python3" ]; then
        export PATH="${py_dir}:$PATH"
        break
    fi
done

PYTHON_CMD="python3"
if ! command -v python3 &> /dev/null; then
    if command -v python &> /dev/null; then
        PYTHON_CMD="python"
    else
        echo "ERROR: Neither 'python3' nor 'python' found in PATH." >&2
        exit 1
    fi
fi

# Configurations (Fixed theta=0.50, coverage sweep across 0.50 to 0.90)
THETAS="${THETAS:-0.50}"
CMINS="${CMINS:-0.50 0.60 0.70 0.80 0.90}"
TOP_N="${TOP_N:-20}"

ALN_TSV="data/dataset_construction/03_homology_clustering/arabidopsis/ara_prot_all_vs_all.tsv"
MAPPINGS_DIR="data/dataset_construction/03_homology_clustering/arabidopsis/mappings"
OUTPUT_MD="${OUTPUT_MD:-data/dataset_construction/03_homology_clustering/arabidopsis/top20_cluster_density_report.md}"
OUTPUT_TSV="${OUTPUT_TSV:-data/dataset_construction/03_homology_clustering/arabidopsis/top20_cluster_density_metrics.tsv}"

echo "Repository Root:     ${REPO_ROOT}"
echo "Date / Time:         $(date)"
echo "Identity Thresholds: ${THETAS}"
echo "Coverage Thresholds: ${CMINS}"
echo "Top N Clusters:      ${TOP_N}"
echo "Output Markdown:     ${OUTPUT_MD}"
echo "Output TSV:          ${OUTPUT_TSV}"
echo "========================================================================"

if [ ! -f "${ALN_TSV}" ]; then
    echo "ERROR: Alignment file missing: ${ALN_TSV}" >&2
    echo "Please ensure MMseqs2 all-vs-all search has completed first." >&2
    exit 1
fi

${PYTHON_CMD} scripts/dataset_construction/evaluate_cluster_density.py \
    --alignments_tsv "${ALN_TSV}" \
    --mappings_dir "${MAPPINGS_DIR}" \
    --output_md "${OUTPUT_MD}" \
    --output_tsv "${OUTPUT_TSV}" \
    --thetas ${THETAS} \
    --cmins ${CMINS} \
    --top_n "${TOP_N}"

echo ""
echo "========================================================================"
echo "Cluster density evaluation completed successfully!"
echo "Outputs generated in data directory:"
echo "  - Markdown Report: ${OUTPUT_MD}"
echo "  - Tabular TSV:     ${OUTPUT_TSV}"
echo "========================================================================"
