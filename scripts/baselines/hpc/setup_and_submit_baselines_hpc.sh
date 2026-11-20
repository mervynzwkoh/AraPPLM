#!/bin/bash
# ==============================================================================
# NSCC Login-Node Setup & Job Submission Script for Legacy Baselines on AraPPLM
# (ESMAraPPI, ARACoFusion, DeepAraPPI)
#
# Run this script on the NSCC login node (which has internet access):
#   bash scripts/baselines/hpc/setup_and_submit_baselines_hpc.sh
#
# What it does:
#   1. Clones zjy1125/DeepAraPPI and keiwo/ESMAraPPI into external_repos/
#   2. Ensures `fair-esm` and `scikit-learn` are installed in PPLM_NSCC_A100
#   3. Pre-caches `esm1b_t33_650M_UR50S.pt` on the login node (since NSCC compute
#      nodes do not have outbound internet access)
#   4. Merges `ara_positive_pairs_high_only_hetero.fasta` into `data/arabidopsis/uniprot_final.pkl`
#   5. Submits `scripts/baselines/hpc/run_baselines_araplm_nscc.pbs` via `qsub`
# ==============================================================================
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
cd "$PROJECT_ROOT"

echo "=== [1/5] Preparing External Repositories in external_repos/ ==="
mkdir -p external_repos
if [ ! -d "external_repos/DeepAraPPI/.git" ]; then
    echo "Cloning https://github.com/zjy1125/DeepAraPPI.git ..."
    git clone https://github.com/zjy1125/DeepAraPPI.git external_repos/DeepAraPPI
else
    echo "external_repos/DeepAraPPI already present."
fi

if [ ! -d "external_repos/ESMAraPPI/.git" ]; then
    echo "Cloning https://github.com/keiwo/ESMAraPPI.git ..."
    git clone https://github.com/keiwo/ESMAraPPI.git external_repos/ESMAraPPI
else
    echo "external_repos/ESMAraPPI already present."
fi

echo -e "\n=== [2/5] Activating Conda Environment (PPLM_NSCC_A100) ==="
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniforge3/etc/profile.d/conda.sh"
fi
conda activate PPLM_NSCC_A100

python -c "import esm" 2>/dev/null || {
    echo "Installing fair-esm..."
    pip install fair-esm
}

echo -e "\n=== [3/5] Pre-caching ESM-1b (esm1b_t33_650M_UR50S) Weights on Login Node ==="
python - << 'EOF'
import esm
print("Verifying ESM-1b checkpoint in torch hub cache...")
model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
del model
print("ESM-1b weights cached and verified!")
EOF

echo -e "\n=== [4/5] Ensuring 100% Sequence Coverage in data/arabidopsis/uniprot_final.pkl ==="
python - << 'EOF'
import os, pickle

pkl_path = "data/arabidopsis/uniprot_final.pkl"
fasta_path = "data/dataset_construction/01_positive_pairs/arabidopsis/ara_positive_pairs_high_only_hetero.fasta"

if os.path.exists(pkl_path) and os.path.exists(fasta_path):
    with open(pkl_path, "rb") as f:
        seq_db = pickle.load(f)
    before = len(seq_db)
    cur_id, cur_seq = None, []
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if cur_id:
                    seq_db[cur_id] = "".join(cur_seq)
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line)
        if cur_id:
            seq_db[cur_id] = "".join(cur_seq)
    if len(seq_db) > before:
        with open(pkl_path, "wb") as f:
            pickle.dump(seq_db, f, protocol=pickle.HIGHEST_PROTOCOL)
        print(f"Updated {pkl_path}: {before:,} -> {len(seq_db):,} IDs")
    else:
        print(f"{pkl_path} already up to date ({len(seq_db):,} IDs)")
EOF

echo -e "\n=== [5/5] Submitting PBS Job: run_baselines_araplm_nscc.pbs ==="
qsub scripts/baselines/hpc/run_baselines_araplm_nscc.pbs
