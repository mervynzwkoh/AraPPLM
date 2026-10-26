#!/usr/bin/env bash
# ==============================================================================
# run_filter_negative_interologs_hpc.sh
# ==============================================================================
# End-to-end execution script for Arabidopsis Candidate Negative Homology Filtering.
#
# Pipeline Steps:
#   1. Generates unified FASTA (canonical proteome + positive interactors, ~30.5k proteins).
#   2. Executes MMseqs2 all-vs-all search on the unified proteome (if alignments file missing).
#   3. Executes Python Streaming Homology Filter (03_filter_negative_positive_interologs_hpc.py):
#      - Evaluates pair similarity S(P_neg, P_pos) >= 0.50 (cmin=0.80, evalue <= 1e-3).
#      - Discards positive interologs and streams clean candidate negative chunks.
#
# Usage:
#   bash scripts/dataset_construction/run_filter_negative_interologs_hpc.sh
#
# Environment variables:
#   THREADS=16 bash scripts/dataset_construction/run_filter_negative_interologs_hpc.sh
# ==============================================================================

set -eo pipefail

# Disable unbound variable checking during shell profile sourcing
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

# Activate Conda environment with numpy and pyarrow
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

# Ensure MMseqs2 is in PATH
if ! command -v mmseqs &> /dev/null; then
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
fi

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_ROOT}"

THREADS="${THREADS:-16}"
THETA="${THETA:-0.50}"
CMIN="${CMIN:-0.80}"

echo "========================================================================"
echo "    AraPPLM: Negative Interolog Homology Filtering Pipeline             "
echo "========================================================================"
echo "Date / Time:     $(date)"
echo "Host:            $(hostname)"
echo "Repo Root:       ${REPO_ROOT}"
echo "Python:          ${PYTHON_BIN}"
echo "MMseqs2:         $(command -v mmseqs || echo 'NOT FOUND')"
echo "Threads:         ${THREADS}"
echo "Theta (Identity): ${THETA}"
echo "Cmin (Coverage):  ${CMIN}"
echo "========================================================================"

STREAM_PATH="data/arabidopsis/uniprot_arabidopsis_stream.tsv.gz"
HETERO_FASTA="data/dataset_construction/01_positive_pairs/arabidopsis/ara_positive_pairs_high_only_hetero.fasta"
UNIFIED_FASTA="data/arabidopsis/arabidopsis_canonical_plus_positives.fasta"
SPLITS_DIR="data/dataset_construction/04_benchmark_splits/arabidopsis"
OUT_BASE="data/dataset_construction/02_negative_sampling"
MMSEQS_DIR="${OUT_BASE}/mmseqs_db"
MMSEQS_TMP="${MMSEQS_DIR}/tmp"
ALN_TSV="${OUT_BASE}/ara_canonical_plus_positives_alignments.tsv"
CANDIDATES_DIR="${OUT_BASE}/candidates"
FILTERED_DIR="${OUT_BASE}/filtered_candidates"
AUDIT_JSON="${OUT_BASE}/ara_candidate_negatives_interolog_removal_audit.json"
SUMMARY_MD="${OUT_BASE}/summary_candidate_negatives_filtered.md"

mkdir -p "${MMSEQS_DIR}" "${MMSEQS_TMP}" "${FILTERED_DIR}"

# ------------------------------------------------------------------------------
# STEP 1: Generate Unified FASTA (Canonical Proteome + Positive Interactors)
# ------------------------------------------------------------------------------
if [ ! -s "${UNIFIED_FASTA}" ]; then
    echo ""
    echo "========================================================================"
    echo "Step 1/3: Building Unified Proteome FASTA (~30.5k sequences)"
    echo "Started at: $(date)"
    echo "========================================================================"
    "${PYTHON_BIN}" scripts/dataset_construction/03_filter_negative_positive_interologs_hpc.py \
        --stream_path "${STREAM_PATH}" \
        --hetero_fasta "${HETERO_FASTA}" \
        --unified_fasta "${UNIFIED_FASTA}" \
        --build_fasta_only
    echo "Step 1 completed at $(date)."
else
    echo "Notice: Unified FASTA already exists at ${UNIFIED_FASTA}. Skipping generation."
fi

# ------------------------------------------------------------------------------
# STEP 2: MMseqs2 All-vs-All Search
# ------------------------------------------------------------------------------
if [ -s "${ALN_TSV}" ]; then
    echo "Notice: Existing MMseqs2 alignments file found at: ${ALN_TSV}"
    echo "Skipping MMseqs2 alignment search."
else
    echo ""
    echo "========================================================================"
    echo "Step 2/3: Executing MMseqs2 All-vs-All Search on Unified Proteome"
    echo "Started at: $(date)"
    echo "========================================================================"

    if ! command -v mmseqs &> /dev/null; then
        echo "ERROR: 'mmseqs' executable not found in PATH." >&2
        echo "Please load mmseqs2 module on HPC, e.g.: module load mmseqs2" >&2
        exit 1
    fi

    echo "Building MMseqs2 Sequence Database..."
    mmseqs createdb "${UNIFIED_FASTA}" "${MMSEQS_DIR}/prot_db"

    echo "Running Permissive MMseqs2 All-vs-All Search (-s 7.5, min-seq-id 0.40, cov 0.70)..."
    mmseqs search \
        "${MMSEQS_DIR}/prot_db" \
        "${MMSEQS_DIR}/prot_db" \
        "${MMSEQS_DIR}/res_db" \
        "${MMSEQS_TMP}" \
        -s 7.5 \
        --min-seq-id 0.40 \
        -c 0.70 \
        --cov-mode 0 \
        --threads "${THREADS}"

    echo "Exporting Tabular Alignments to: ${ALN_TSV}..."
    mmseqs convertalis \
        "${MMSEQS_DIR}/prot_db" \
        "${MMSEQS_DIR}/prot_db" \
        "${MMSEQS_DIR}/res_db" \
        "${ALN_TSV}" \
        --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qcov,tcov" \
        --threads "${THREADS}"

    echo "Step 2 completed at $(date). Exported $(wc -l < "${ALN_TSV}") alignments."
fi

# ------------------------------------------------------------------------------
# STEP 3: Streaming Candidate Negative Interolog Purge
# ------------------------------------------------------------------------------
echo ""
echo "========================================================================"
echo "Step 3/3: Running Streaming Negative Interolog Homology Filter"
echo "Started at: $(date)"
echo "========================================================================"

"${PYTHON_BIN}" scripts/dataset_construction/03_filter_negative_positive_interologs_hpc.py \
    --stream_path "${STREAM_PATH}" \
    --hetero_fasta "${HETERO_FASTA}" \
    --unified_fasta "${UNIFIED_FASTA}" \
    --splits_dir "${SPLITS_DIR}" \
    --alignments_tsv "${ALN_TSV}" \
    --candidates_dir "${CANDIDATES_DIR}" \
    --out_dir "${FILTERED_DIR}" \
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
echo "Clean filtered chunks in: ${FILTERED_DIR}/"
ls -lh "${FILTERED_DIR}"
echo "========================================================================"
