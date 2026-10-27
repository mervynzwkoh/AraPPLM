#!/usr/bin/env bash
# ==============================================================================
# run_cluster_candidate_negatives_hpc.sh
# ==============================================================================
# End-to-end execution script for Arabidopsis Candidate Negative Pair-Level Redundancy Reduction.
#
# Steps:
#   1. Resolves environment (Conda python11 / PPLM_NSCC_A100 / base).
#   2. Pre-flight check: validates MMseqs2 alignments and filtered candidate chunks.
#   3. Executes Python CD-HIT Streaming Clustering (04_cluster_candidate_negatives_mmseqs_hpc.py):
#      - Evaluates pair similarity S(P_i, P_j) >= 0.50 (cmin=0.80, evalue <= 1e-3).
#      - Decouples computation across homology component pairs.
#      - Retains non-redundant representatives prioritized by Evidence Quality -> Length -> Lexicographical.
#      - Outputs final non-redundant negative candidate Parquet chunks.
#
# Usage:
#   bash scripts/dataset_construction/run_cluster_candidate_negatives_hpc.sh
# ==============================================================================

set -eo pipefail

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
else
    module load miniforge3 2>/dev/null || module load miniconda3 2>/dev/null || true
fi

# Activate Conda environment
TARGET_ENV="${CONDA_ENV_NAME:-${CONDA_DEFAULT_ENV:-python11}}"
for env_name in "$TARGET_ENV" python11 PPLM_NSCC_A100 base; do
    if [ -n "$env_name" ] && conda activate "$env_name" 2>/dev/null; then
        echo "Activated Conda environment: $env_name"
        break
    fi
done

# Resolve python binary
PYTHON_BIN="$(which python3 2>/dev/null || which python 2>/dev/null || echo "python3")"
if ! "${PYTHON_BIN}" -c "import numpy; import pyarrow" 2>/dev/null; then
    for candidate in \
        ${CONDA_PREFIX:+"${CONDA_PREFIX}/bin/python3"} \
        "$HOME/.conda/envs/"*"/bin/python3" \
        "$HOME/miniconda3/envs/"*"/bin/python3" \
        "$HOME/miniforge3/envs/"*"/bin/python3" \
        "/app1/common/miniconda3/envs/"*"/bin/python3" \
        "$HOME/miniconda3/bin/python3"; do
        if [ -x "$candidate" ] && "$candidate" -c "import numpy; import pyarrow" 2>/dev/null; then
            PYTHON_BIN="$candidate"
            break
        fi
    done
fi

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_ROOT}"

THETA="${THETA:-0.50}"
CMIN="${CMIN:-0.80}"

echo "========================================================================"
echo "    AraPPLM: Negative Pair Redundancy Reduction Pipeline (HPC)          "
echo "========================================================================"
echo "Date / Time:     $(date)"
echo "Host:            $(hostname)"
echo "Repo Root:       ${REPO_ROOT}"
echo "Python:          ${PYTHON_BIN}"
echo "Theta (Identity): ${THETA}"
echo "Cmin (Coverage):  ${CMIN}"
echo "========================================================================"

STREAM_PATH="data/arabidopsis/uniprot_arabidopsis_stream.tsv.gz"
TAIR_GAF_PATH="data/raw/arabidopsis/tair/tair.gaf.gz"
OUT_BASE="data/dataset_construction/02_negative_sampling"
ALN_TSV="${OUT_BASE}/ara_canonical_plus_positives_alignments.tsv"
CANDIDATES_DIR="${OUT_BASE}/filtered_candidates"
FINAL_DIR="${OUT_BASE}/final_candidates_nr"
AUDIT_JSON="${OUT_BASE}/ara_candidate_negatives_redundancy_reduction_audit.json"
SUMMARY_MD="${OUT_BASE}/summary_candidate_negatives_nr.md"

# Pre-flight check
if [ ! -s "${ALN_TSV}" ]; then
    echo "ERROR: Alignments TSV missing: ${ALN_TSV}" >&2
    echo "Please ensure Stage 2.5 (interolog filtering) completed and generated this file." >&2
    exit 1
fi

if [ ! -d "${CANDIDATES_DIR}" ] || [ -z "$(ls -A "${CANDIDATES_DIR}" 2>/dev/null)" ]; then
    echo "ERROR: Clean filtered candidate chunks not found in: ${CANDIDATES_DIR}" >&2
    exit 1
fi

mkdir -p "${FINAL_DIR}"

echo ""
echo "========================================================================"
echo "Starting Pair Redundancy Reduction across Clean Candidates..."
echo "Started at: $(date)"
echo "========================================================================"

"${PYTHON_BIN}" scripts/dataset_construction/04_cluster_candidate_negatives_mmseqs_hpc.py \
    --stream_path "${STREAM_PATH}" \
    --tair_gaf_path "${TAIR_GAF_PATH}" \
    --alignments_tsv "${ALN_TSV}" \
    --candidates_dir "${CANDIDATES_DIR}" \
    --out_dir "${FINAL_DIR}" \
    --audit_path "${AUDIT_JSON}" \
    --summary_md "${SUMMARY_MD}" \
    --theta "${THETA}" \
    --cmin "${CMIN}" \
    --max_evalue 1e-3 \
    --chunk_size 1000000

echo ""
echo "========================================================================"
echo "Pipeline finished successfully at $(date)!"
echo "Audit summary:"
cat "${AUDIT_JSON}"
echo ""
echo "Final Non-Redundant Negative Chunks in: ${FINAL_DIR}/"
ls -lh "${FINAL_DIR}"
echo "========================================================================"
