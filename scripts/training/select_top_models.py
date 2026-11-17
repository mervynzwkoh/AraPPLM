#!/usr/bin/env python3
"""
Select top-K PPI model checkpoints and package them into ensemble lists
and a consolidated state_dict archive.

Aligned with the original PPLM author's model selection step
(ppi_training_code/run_testing.sh):
  1. Parses per-fold .recording files for all pooling modes (mean, max, min)
  2. Extracts epoch and validation AUC_PR (column 14 of `Accuarcy Epoch:` lines)
  3. Ranks all fold x epoch checkpoints globally per mode by AUC_PR descending
  4. Writes author-compatible path-list files:
       {prefix}.{mode}.round_{seed_round}.top{top_k}_list
     (one checkpoint path per line, directly consumed by --mean_path/--max_path/--min_path)
  5. Packages state_dicts into ppi_plant_models.pkl:
       {'mean': [state_dict, ...], 'max': [state_dict, ...], 'min': [state_dict, ...]}

Usage:
    python scripts/training/select_top_models.py \
        --model_dir models/DeepAraPPI/ \
        --output models/DeepAraPPI/ppi_plant_models.pkl \
        --top_k 5 --modes mean max min --seed_round 1
"""

import argparse
import glob
import os
import re
import sys

import torch


def find_recording_files(model_dir, mode, prefix, seed_round):
    """Discover recording files for a pooling mode (author format first, then legacy)."""
    patterns = [
        os.path.join(model_dir, f"{prefix}.{mode}.round_{seed_round}.cv_*.recording"),
        os.path.join(model_dir, f"{prefix}.{mode}.round_*.cv_*.recording"),
        os.path.join(model_dir, f"plant_ppi.{mode}.cv_*.recording"),
    ]
    for pattern in patterns:
        matches = sorted(glob.glob(pattern))
        if matches:
            return matches
    return []


def parse_recording_files(model_dir, mode, prefix="PPLM-PPI.inter_intra_attn_embed", seed_round=1):
    """
    Parse all recording files for a given pooling mode and extract
    per-fold, per-epoch AUC_PR values.

    Returns:
        List of (auprc, checkpoint_path) tuples deduplicated by checkpoint_path.
    """
    recording_files = find_recording_files(model_dir, mode, prefix, seed_round)
    if not recording_files:
        print(f"  No recording files found for mode={mode} in {model_dir}")
        return []

    ckpt_scores = {}

    for rec_path in recording_files:
        fold_dir = rec_path.replace(".recording", "")
        if not os.path.isdir(fold_dir):
            print(f"  WARNING: Fold directory not found: {fold_dir}")
            continue

        with open(rec_path, "r") as f:
            for raw_line in f:
                line = raw_line.strip()
                if not (
                    line.startswith("========== Accuarcy")
                    or line.startswith("========== Accuracy")
                ):
                    continue

                epoch_match = re.search(r"Epoch:\s*(\d+)/", line)
                if not epoch_match:
                    continue
                epoch = int(epoch_match.group(1))

                # Legacy tagged format ("AUPRC: 0.8765") vs author positional format (col 14)
                auprc_match = re.search(r"AUPRC:\s*([\d.eE+-]+)", line)
                if auprc_match:
                    auprc = float(auprc_match.group(1))
                else:
                    parts = line.split()
                    if len(parts) < 14:
                        continue
                    try:
                        auprc = float(parts[13])
                    except ValueError:
                        continue

                checkpoint_path = os.path.join(fold_dir, f"model_{epoch}.pkl")
                if os.path.isfile(checkpoint_path):
                    ckpt_scores[checkpoint_path] = auprc

    return [(auprc, ckpt_path) for ckpt_path, auprc in ckpt_scores.items()]


def main():
    parser = argparse.ArgumentParser(
        description="Select top-K PPI model checkpoints and package into ensemble lists/weights"
    )
    parser.add_argument(
        "--model_dir",
        required=True,
        help="Directory containing trained model fold directories and recording files",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output path for packaged ensemble weights file (ppi_plant_models.pkl)",
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
        default=["mean", "max", "min"],
        help="Pooling modes to include in ensemble (default: mean max min)",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default="PPLM-PPI.inter_intra_attn_embed",
        help="Model folder/recording prefix (default: PPLM-PPI.inter_intra_attn_embed)",
    )
    parser.add_argument(
        "--seed_round",
        type=int,
        default=1,
        help="Seed round number matching training runs (default: 1)",
    )
    args = parser.parse_args()

    print(f"Selecting top-{args.top_k} models from: {args.model_dir}")
    print(f"Pooling modes: {args.modes} | seed_round={args.seed_round}")

    ensemble_weights = {}

    for mode in args.modes:
        print(f"\n--- Mode: {mode} ---")
        candidates = parse_recording_files(
            args.model_dir, mode, prefix=args.prefix, seed_round=args.seed_round
        )
        if not candidates:
            print(f"  ERROR: No candidates found for mode={mode}")
            continue

        candidates.sort(key=lambda x: x[0], reverse=True)
        top_k_items = candidates[: args.top_k]

        print(f"  Found {len(candidates)} total checkpoints")
        print(f"  Selected top-{len(top_k_items)}:")

        state_dicts = []
        for rank, (auprc, ckpt_path) in enumerate(top_k_items):
            fold_label = os.path.basename(os.path.dirname(ckpt_path))
            print(
                f"    #{rank + 1}: AUC_PR={auprc:.4f} | "
                f"{fold_label}/{os.path.basename(ckpt_path)}"
            )
            checkpoint = torch.load(ckpt_path, map_location="cpu", weights_only=False)
            state_dicts.append(checkpoint["model_state_dict"])

        ensemble_weights[mode] = state_dicts

        # 1. Author-format path list (one path per line for --mean_path/--max_path/--min_path)
        author_list_path = os.path.join(
            args.model_dir,
            f"{args.prefix}.{mode}.round_{args.seed_round}.top{args.top_k}_list",
        )
        with open(author_list_path, "w") as f:
            for _, ckpt_path in top_k_items:
                f.write(f"{ckpt_path}\n")
        print(f"  Author top-K path list saved: {author_list_path}")

        # 2. Scored reference list
        ref_list_path = os.path.join(
            args.model_dir, f"plant_ppi.{mode}.top{args.top_k}_list"
        )
        with open(ref_list_path, "w") as f:
            for auprc, ckpt_path in top_k_items:
                f.write(f"{ckpt_path}\t{auprc:.6f}\n")
        print(f"  Scored reference list saved:  {ref_list_path}")

    if not ensemble_weights:
        print("\nERROR: No models selected. Check that training has completed.")
        sys.exit(1)

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    torch.save(ensemble_weights, args.output)

    total_models = sum(len(v) for v in ensemble_weights.values())
    summary = ", ".join(f"{k}={len(v)}" for k, v in ensemble_weights.items())
    print(f"\nEnsemble weights saved: {args.output}")
    print(f"Total models: {total_models} ({summary})")


if __name__ == "__main__":
    main()
