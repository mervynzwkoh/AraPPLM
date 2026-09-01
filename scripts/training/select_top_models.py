#!/usr/bin/env python3
"""
Select top-K PPI model checkpoints and package them into a single ensemble
weights file compatible with the PPLM-PPI inference pipeline.

After training completes (10 folds × 3 pooling modes), this script:
  1. Parses all recording files to find best epoch per fold
  2. Ranks all fold×epoch checkpoints by validation AUPRC
  3. Selects top-K per pooling mode (default: top-5)
  4. Packages state_dicts into ppi_plant_models.pkl
     Format: {'mean': [state_dict, ...], 'max': [state_dict, ...]}

This format is directly compatible with the existing batch_predict.py
ensemble inference pipeline.

Usage:
    python scripts/training/select_top_models.py \
        --model_dir models/DeepAraPPI/ \
        --output models/DeepAraPPI/ppi_plant_models.pkl \
        --top_k 5
"""

import os
import sys
import argparse
import re
import torch
import glob


def parse_recording_files(model_dir, mode):
    """
    Parse all recording files for a given pooling mode and extract
    per-fold, per-epoch AUPRC values.

    Returns: List of (auprc, checkpoint_path) tuples
    """
    pattern = os.path.join(model_dir, f"plant_ppi.{mode}.cv_*.recording")
    recording_files = sorted(glob.glob(pattern))

    if not recording_files:
        print(f"  No recording files found for mode={mode}")
        return []

    candidates = []

    for rec_path in recording_files:
        # Extract fold directory from recording path
        fold_dir = rec_path.replace(".recording", "")

        if not os.path.isdir(fold_dir):
            print(f"  WARNING: Fold directory not found: {fold_dir}")
            continue

        with open(rec_path, "r") as f:
            lines = f.readlines()

        # Parse accuracy lines to extract AUPRC and epoch
        for line in lines:
            line = line.strip()
            if not line.startswith("========== Accuracy"):
                continue

            # Extract epoch number
            epoch_match = re.search(r"Epoch:\s*(\d+)/", line)
            if not epoch_match:
                continue
            epoch = int(epoch_match.group(1))

            # Extract AUPRC
            auprc_match = re.search(r"AUPRC:\s*([\d.]+)", line)
            if not auprc_match:
                continue
            auprc = float(auprc_match.group(1))

            # Construct checkpoint path
            checkpoint_path = os.path.join(fold_dir, f"model_{epoch}.pkl")
            if os.path.isfile(checkpoint_path):
                candidates.append((auprc, checkpoint_path))

    return candidates


def main():
    parser = argparse.ArgumentParser(
        description="Select top-K PPI model checkpoints and package into ensemble weights"
    )
    parser.add_argument(
        "--model_dir",
        required=True,
        help="Directory containing trained model fold directories and recording files",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output path for ensemble weights file (ppi_plant_models.pkl)",
    )
    parser.add_argument(
        "--top_k",
        type=int,
        default=5,
        help="Number of top checkpoints to select per pooling mode (default: 5)",
    )
    parser.add_argument(
        "--modes",
        nargs="+",
        default=["mean", "max"],
        help="Pooling modes to include in ensemble (default: mean max)",
    )
    args = parser.parse_args()

    print(f"Selecting top-{args.top_k} models from: {args.model_dir}")
    print(f"Pooling modes: {args.modes}")

    ensemble_weights = {}

    for mode in args.modes:
        print(f"\n--- Mode: {mode} ---")

        candidates = parse_recording_files(args.model_dir, mode)
        if not candidates:
            print(f"  ERROR: No candidates found for mode={mode}")
            continue

        # Sort by AUPRC descending
        candidates.sort(key=lambda x: x[0], reverse=True)

        # Select top-K
        top_k = candidates[: args.top_k]

        print(f"  Found {len(candidates)} total checkpoints")
        print(f"  Selected top-{args.top_k}:")
        state_dicts = []
        for rank, (auprc, ckpt_path) in enumerate(top_k):
            print(f"    #{rank + 1}: AUPRC={auprc:.4f} | {os.path.basename(ckpt_path)}")
            checkpoint = torch.load(ckpt_path, map_location="cpu", weights_only=False)
            state_dicts.append(checkpoint["model_state_dict"])

        ensemble_weights[mode] = state_dicts

        # Also write the top-K list to a text file for reference
        list_path = os.path.join(
            args.model_dir, f"plant_ppi.{mode}.top{args.top_k}_list"
        )
        with open(list_path, "w") as f:
            for auprc, ckpt_path in top_k:
                f.write(f"{ckpt_path}\t{auprc:.6f}\n")
        print(f"  Model list saved: {list_path}")

    # Save ensemble weights
    if not ensemble_weights:
        print("\nERROR: No models selected. Check that training has completed.")
        sys.exit(1)

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    torch.save(ensemble_weights, args.output)

    total_models = sum(len(v) for v in ensemble_weights.values())
    print(f"\nEnsemble weights saved: {args.output}")
    print(f"Total models: {total_models} ({', '.join(f'{k}={len(v)}' for k, v in ensemble_weights.items())})")


if __name__ == "__main__":
    main()
