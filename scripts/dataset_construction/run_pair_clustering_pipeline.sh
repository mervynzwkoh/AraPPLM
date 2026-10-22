#!/usr/bin/env bash
# ==============================================================================
# run_pair_clustering_pipeline.sh
# ==============================================================================
# End-to-end runner for MMseqs2-accelerated Arabidopsis PPI pair clustering.
#
# Steps:
#   1. Runs MMseqs2 all-vs-all search on unique proteins (one-time permissive run)
#   2. Exports tabular alignment metrics (identity, coverage, evalue)
#   3. Runs Python multi-parameter sweep across 12 combinations:
#      - theta in [0.40, 0.70, 0.90]
#      - cmin in [0.50, 0.70, 0.80, 0.90]
#   4. Generates representative datasets, cluster mappings, and summary metrics.
#
# Usage:
#   bash scripts/dataset_construction/run_pair_clustering_pipeline.sh
#   THREADS=16 bash scripts/dataset_construction/run_pair_clustering_pipeline.sh
# ==============================================================================

set -euo pipefail

# Dynamic root resolution
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_ROOT}"

echo "========================================================================"
echo "    AraPPLM: MMseqs2 Pair Redundancy Clustering Pipeline                "
echo "========================================================================"
echo "Repository Root: ${REPO_ROOT}"
echo "Date / Time:     $(date)"
echo "Threads:         ${THREADS:-8}"
echo "========================================================================"

# 1. Dependency checks
if ! command -v mmseqs &> /dev/null; then
    echo "ERROR: 'mmseqs' command not found in PATH." >&2
    echo "If on HPC, please load the module before running, e.g.:" >&2
    echo "    module load mmseqs2" >&2
    echo "Or activate your conda environment containing mmseqs2." >&2
    exit 1
fi

PYTHON_CMD="python3"
if ! command -v python3 &> /dev/null; then
    if command -v python &> /dev/null; then
        PYTHON_CMD="python"
    else
        echo "ERROR: Neither 'python3' nor 'python' found in PATH." >&2
        exit 1
    fi
fi

# 2. File and directory configurations
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

# 3. MMseqs2 All-vs-All Search (Skip if ALN_TSV already exists and is non-empty)
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
        --threads "${THREADS:-8}"

    echo "--- Step 3/3: Exporting Alignments to TSV ---"
    mmseqs convertalis \
        "${DB_DIR}/prot_db" \
        "${DB_DIR}/prot_db" \
        "${DB_DIR}/res_db" \
        "${ALN_TSV}" \
        --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qcov,tcov" \
        --threads "${THREADS:-8}"
        
    echo "Generated alignment file: ${ALN_TSV} ($(wc -l < "${ALN_TSV}") lines)"
fi

# 4. Run Python Parameter Sweep Engine
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
    --thetas 0.40 0.70 0.90 \
    --cmins 0.50 0.70 0.80 0.90 \
    --sweep

echo ""
echo "========================================================================"
echo "Pipeline completed successfully!"
echo "Check summary results in:"
echo "  - ${OUT_DIR}/summary_clustering_metrics.tsv"
echo "  - ${OUT_DIR}/summary_clustering_metrics.md"
echo "========================================================================"
