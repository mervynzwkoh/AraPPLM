#!/usr/bin/env bash
# ==============================================================================
# run_cluster_density_hpc.sh
# ==============================================================================
# Calculates the intra-cluster clique density metric f_k(theta) and minimum pair
# similarity for the 20 largest clusters at theta=0.40 and theta=0.70.
#
# Can be run directly on an HPC login/compute node (runs in ~2 seconds).
#
# Usage:
#   bash scripts/dataset_construction/run_cluster_density_hpc.sh
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_ROOT}"

echo "========================================================================"
echo "    AraPPLM: Evaluating Intra-Cluster Pair Similarity & Clique Density  "
echo "========================================================================"

# Sourcing conda/python if available
if [ -f "$HOME/.bashrc" ]; then
    source "$HOME/.bashrc"
fi
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniforge3/etc/profile.d/conda.sh"
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
    PYTHON_CMD="python"
fi

ALN_TSV="data/dataset_construction/03_homology_clustering/arabidopsis/ara_prot_all_vs_all.tsv"
MAPPINGS_DIR="data/dataset_construction/03_homology_clustering/arabidopsis/mappings"
SUMMARY_MD="data/dataset_construction/03_homology_clustering/arabidopsis/summary_clustering_metrics.md"
OUTPUT_TSV="data/dataset_construction/03_homology_clustering/arabidopsis/top20_cluster_density_metrics.tsv"

if [ ! -f "${ALN_TSV}" ]; then
    echo "ERROR: Alignment file missing: ${ALN_TSV}" >&2
    echo "Please ensure MMseqs2 all-vs-all search has completed first." >&2
    exit 1
fi

${PYTHON_CMD} scripts/dataset_construction/evaluate_cluster_density.py \
    --alignments_tsv "${ALN_TSV}" \
    --mappings_dir "${MAPPINGS_DIR}" \
    --summary_md "${SUMMARY_MD}" \
    --output_tsv "${OUTPUT_TSV}" \
    --cmin 0.50 \
    --top_n 20

echo ""
echo "========================================================================"
echo "Cluster density evaluation complete!"
echo "Results updated in:"
echo "  - ${SUMMARY_MD}"
echo "  - ${OUTPUT_TSV}"
echo "========================================================================"
