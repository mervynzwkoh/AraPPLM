#!/usr/bin/env bash
# ==============================================================================
# run_partition_negative_hpc.sh
# ==============================================================================
# End-to-end execution script for Arabidopsis Negative PPI Benchmark Partitioning.
#
# Partitions the final non-redundant negative candidate pool into:
#   - C1 (Train Set): 194,620 pairs (10x ratio of 19,462 positive train pairs)
#   - C2 (One Unseen Test): 11,450 pairs (10x ratio of 1,145 positive test pairs)
#   - C3 (Both Unseen Test): 11,450 pairs (10x ratio of 1,145 positive test pairs)
#   - C4 (Both Seen Test): 11,450 pairs (10x ratio of 1,145 positive test pairs)
#
# Total Negative Benchmark Pairs: 228,970 pairs
# Combined Benchmark Pairs (Pos + Neg): 251,867 pairs
#
# Usage:
#   bash scripts/dataset_construction/run_partition_negative_hpc.sh
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

THETA="${THETA:-0_50}"
CMIN="${CMIN:-0_80}"
RATIO="${RATIO:-10.0}"
SEED="${SEED:-42}"

echo "========================================================================"
echo "   AraPPLM: Negative Benchmark Partitioning Pipeline (HPC)              "
echo "========================================================================"
echo "Date / Time:       $(date)"
echo "Host:              $(hostname)"
echo "Repo Root:         ${REPO_ROOT}"
echo "Python:            ${PYTHON_BIN}"
echo "Theta (Identity):  ${THETA}"
echo "Cmin (Coverage):   ${CMIN}"
echo "Pos:Neg Ratio:     1:${RATIO}"
echo "Seed:              ${SEED}"
echo "========================================================================"

SPLITS_DIR="data/dataset_construction/04_benchmark_splits/arabidopsis"
CANDIDATES_DIR="data/dataset_construction/02_negative_sampling/final_candidates_nr"
STREAM_PATH="data/arabidopsis/uniprot_arabidopsis_stream.tsv.gz"
TAIR_GAF_PATH="data/raw/arabidopsis/tair/tair.gaf.gz"
AUDIT_JSON="${SPLITS_DIR}/ara_negative_c1234_partitions_audit.json"
SUMMARY_MD="${SPLITS_DIR}/summary_negative_c1234_partitions.md"

# Pre-flight check: positive splits exist
POS_C1="${SPLITS_DIR}/ara_positive_c1_train_id${THETA}_cov${CMIN}.tsv"
POS_C2="${SPLITS_DIR}/ara_positive_c2_test_id${THETA}_cov${CMIN}.tsv"
POS_C3="${SPLITS_DIR}/ara_positive_c3_test_id${THETA}_cov${CMIN}.tsv"
POS_C4="${SPLITS_DIR}/ara_positive_c4_test_id${THETA}_cov${CMIN}.tsv"

for f in "${POS_C1}" "${POS_C2}" "${POS_C3}" "${POS_C4}"; do
    if [ ! -s "${f}" ]; then
        echo "ERROR: Positive split missing: ${f}" >&2
        echo "Please verify that the positive C1/C2/C3/C4 partitions are generated." >&2
        exit 1
    fi
done

# Pre-flight check: candidate negative chunks exist
if [ ! -d "${CANDIDATES_DIR}" ] || [ -z "$(ls -A "${CANDIDATES_DIR}" 2>/dev/null)" ]; then
    echo "ERROR: Non-redundant candidate negative chunks not found in: ${CANDIDATES_DIR}" >&2
    echo "Please ensure Stage 2.6 (redundancy reduction) has completed successfully." >&2
    exit 1
fi

TOTAL_CHUNKS=$(ls -1 "${CANDIDATES_DIR}"/ara_candidate_negatives_nr_chunk_*.parquet 2>/dev/null | wc -l)
echo "Verified: Found ${TOTAL_CHUNKS} candidate negative Parquet chunks in ${CANDIDATES_DIR}."

echo ""
echo "[Step 1/1] Executing C1/C2/C3/C4 Negative Partitioning (1:10 Ratio)..."
"${PYTHON_BIN}" scripts/dataset_construction/05_partition_negative_c1234_hpc.py \
    --splits-dir "${SPLITS_DIR}" \
    --candidates-dir "${CANDIDATES_DIR}" \
    --stream-path "${STREAM_PATH}" \
    --tair-gaf "${TAIR_GAF_PATH}" \
    --theta "${THETA}" \
    --cmin "${CMIN}" \
    --ratio "${RATIO}" \
    --seed "${SEED}" \
    --out-dir "${SPLITS_DIR}" \
    --audit-path "${AUDIT_JSON}" \
    --summary-md "${SUMMARY_MD}"

echo ""
echo "========================================================================"
echo " Pipeline Complete: Negative Benchmark Partitions Generated Successfully!"
echo "========================================================================"
echo "Output Directory: ${SPLITS_DIR}"
echo ""
echo "Negative Partitions:"
ls -lh "${SPLITS_DIR}"/ara_negative_c*_id${THETA}_cov${CMIN}.tsv
echo ""
echo "Combined Benchmarks (Pos + Neg):"
ls -lh "${SPLITS_DIR}"/ara_benchmark_c*_id${THETA}_cov${CMIN}.tsv
echo ""
echo "Audit Report: ${SUMMARY_MD}"
echo "========================================================================"
