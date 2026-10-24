#!/usr/bin/env bash
# ==============================================================================
# run_pair_clustering_pipeline.sh
# ==============================================================================
# End-to-end runner for MMseqs2-accelerated Arabidopsis PPI pair clustering.
#
# Steps:
#   1. Runs MMseqs2 all-vs-all search on unique proteins (one-time permissive run,
#      skipped automatically if ara_prot_all_vs_all.tsv already exists).
#   2. Exports tabular alignment metrics (identity, coverage, evalue).
#   3. Runs Python multi-parameter sensitivity sweep across sequence identity
#      thresholds (default: theta in [0.40, 0.50, 0.60, 0.70, 0.80] at fixed cmin 0.80).
#   4. Generates representative datasets, cluster mappings, and summary metrics.
#
# Usage:
#   # Run with default identity sweep (0.40, 0.50, 0.60, 0.70, 0.80) at cmin=0.80:
#   bash scripts/dataset_construction/run_pair_clustering_pipeline.sh
#
#   # Override thresholds via environment variables if desired:
#   THETAS="0.40 0.50 0.60 0.70 0.80" CMINS="0.80" bash scripts/dataset_construction/run_pair_clustering_pipeline.sh
# ==============================================================================

set -eo pipefail

# Sourcing conda/environment if available (disable nounset during profile sourcing)
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

# Activate environments if present
conda activate mmseqs2_env 2>/dev/null || true
for env_dir in \
    "$HOME/.conda/envs/mmseqs2_env/bin" \
    "$HOME/miniconda3/envs/mmseqs2_env/bin" \
    "$HOME/miniforge3/envs/mmseqs2_env/bin" \
    "/app1/common/miniconda3/envs/mmseqs2_env/bin"; do
    if [ -f "${env_dir}/mmseqs" ]; then
        export PATH="${env_dir}:$PATH"
        break
    fi
done

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

set -u

# Dynamic root resolution
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_ROOT}"

# Parameter configurations (sweep theta across [0.40, 0.50, 0.60, 0.70, 0.80] with fixed cmin=0.80)
THETAS="${THETAS:-0.40 0.50 0.60 0.70 0.80}"
CMINS="${CMINS:-0.80}"
THREADS="${THREADS:-8}"

echo "========================================================================"
echo "    AraPPLM: MMseqs2 Pair Redundancy Clustering Pipeline                "
echo "========================================================================"
echo "Repository Root:     ${REPO_ROOT}"
echo "Date / Time:         $(date)"
echo "Identity Thresholds: ${THETAS}"
echo "Coverage Thresholds: ${CMINS}"
echo "Threads:             ${THREADS}"
echo "========================================================================"

# Python command check
PYTHON_CMD="python3"
if ! command -v python3 &> /dev/null; then
    if command -v python &> /dev/null; then
        PYTHON_CMD="python"
    else
        echo "ERROR: Neither 'python3' nor 'python' found in PATH." >&2
        exit 1
    fi
fi

# File and directory configurations
INPUT_FASTA="data/dataset_construction/01_positive_pairs/arabidopsis/ara_positive_pairs_high_only_hetero.fasta"
PAIRS_TSV="data/dataset_construction/01_positive_pairs/arabidopsis/ara_positive_pairs_high_only_hetero.tsv"
METADATA_TSV="data/dataset_construction/01_positive_pairs/arabidopsis/ara_positive_pairs_deduplicated.tsv"

OUT_DIR="data/dataset_construction/03_homology_clustering/arabidopsis"
DB_DIR="${OUT_DIR}/mmseqs_db"
TMP_DIR="${DB_DIR}/tmp"
ALN_TSV="${OUT_DIR}/ara_prot_all_vs_all.tsv"

mkdir -p "${DB_DIR}" "${TMP_DIR}" "${OUT_DIR}/mappings"

# Verify inputs exist
if [ ! -f "${INPUT_FASTA}" ]; then
    echo "ERROR: FASTA file missing: ${INPUT_FASTA}" >&2
    exit 1
fi
if [ ! -f "${PAIRS_TSV}" ]; then
    echo "ERROR: Pairs TSV missing: ${PAIRS_TSV}" >&2
    exit 1
fi

# MMseqs2 check: only required if ALN_TSV does not exist
if [ ! -s "${ALN_TSV}" ]; then
    if ! command -v mmseqs &> /dev/null; then
        echo "ERROR: 'mmseqs' command not found in PATH and ${ALN_TSV} does not exist." >&2
        echo "If on HPC, please load the module before running, e.g.:" >&2
        echo "    module load mmseqs2" >&2
        echo "Or activate your conda environment containing mmseqs2." >&2
        exit 1
    fi
fi

# 1. MMseqs2 All-vs-All Search (Skip if ALN_TSV already exists and is non-empty)
if [ -s "${ALN_TSV}" ]; then
    echo "Notice: Existing alignments file found at ${ALN_TSV}."
    echo "Skipping MMseqs2 database build and alignment search."
    echo "If you wish to force re-computation, delete or rename ${ALN_TSV}."
else
    echo "--- Step 1/3: Creating MMseqs2 Sequence Database ---"
    mmseqs createdb "${INPUT_FASTA}" "${DB_DIR}/prot_db"

    echo "--- Step 2/3: Running Permissive All-vs-All Search ---"
    # Permissive thresholds (min-seq-id 0.3, coverage 0.3) so all grid points can be swept in Python
    mmseqs search \
        "${DB_DIR}/prot_db" \
        "${DB_DIR}/prot_db" \
        "${DB_DIR}/res_db" \
        "${TMP_DIR}" \
        -s 7.5 \
        --min-seq-id 0.3 \
        -c 0.3 \
        --cov-mode 0 \
        --threads "${THREADS}"

    echo "--- Step 3/3: Exporting Alignments to TSV ---"
    mmseqs convertalis \
        "${DB_DIR}/prot_db" \
        "${DB_DIR}/prot_db" \
        "${DB_DIR}/res_db" \
        "${ALN_TSV}" \
        --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qcov,tcov" \
        --threads "${THREADS}"
        
    echo "Generated alignment file: ${ALN_TSV} ($(wc -l < "${ALN_TSV}") lines)"
fi

# 2. Run Python Parameter Sweep Engine
echo ""
echo "========================================================================"
echo "    Executing Python Pair Redundancy Grid Sweep                         "
echo "========================================================================"

${PYTHON_CMD} scripts/dataset_construction/cluster_positive_pairs_mmseqs.py \
    --pairs_tsv "${PAIRS_TSV}" \
    --alignments_tsv "${ALN_TSV}" \
    --metadata_tsv "${METADATA_TSV}" \
    --fasta "${INPUT_FASTA}" \
    --output_dir "${OUT_DIR}" \
    --thetas ${THETAS} \
    --cmins ${CMINS} \
    --sweep

echo ""
echo "========================================================================"
echo "Pipeline completed successfully!"
echo "Check summary results in:"
echo "  - ${OUT_DIR}/summary_clustering_metrics.tsv"
echo "  - ${OUT_DIR}/summary_clustering_metrics.md"
echo "========================================================================"
