#!/usr/bin/env python3
"""
Evaluate plant-trained PPLM-PPI models on held-out test sets (C2, C3, Rice).

This script:
  1. Loads the plant-trained ensemble weights (ppi_plant_models.pkl)
  2. Runs the full PPLM backbone on each test pair (same as benchmarking)
  3. Passes features through all ensemble classifiers with symmetric averaging
  4. Outputs prediction CSVs compatible with the existing evaluate_pplm.py

The inference procedure is identical to the zero-shot benchmarking pipeline
(scripts/benchmarking/batch_predict.py), except it loads plant-trained
weights instead of the original human-trained ppi_models.pkl.

Usage:
    python scripts/training/test_ppi_head.py \
        --input data/DeepAraPPI/c2_ppi_sample_DeepAraPPI.txt \
        --output results/DeepAraPPI/deepara_c2_plant_scores.csv \
        --model_weights models/DeepAraPPI/ppi_plant_models.pkl \
        --seq_db data/arabidopsis/uniprot_final.pkl \
        --gpu_id 0
"""

import os
import sys
import argparse
import csv
import pickle
import torch
import numpy as np
import time

# Path setup
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))
PPLM_DIR = os.path.join(PROJECT_ROOT, "PPLM")
sys.path.insert(0, PPLM_DIR)
sys.path.insert(0, SCRIPT_DIR)

from pplm import PPLM, Alphabet
from ppi_model import PPI_inter_intra_attn_embed_single_pooling


def load_models(device, model_weights_path):
    """Load PPLM backbone and plant-trained PPI classifier weights."""
    print("Loading PPLM backbone...")

    # Load PPLM backbone
    alphabet = Alphabet.from_architecture()
    batch_converter = alphabet.get_batch_converter()
    model_path = os.path.join(PPLM_DIR, "weights", "pplm_t33_650M.pt")

    if not os.path.exists(model_path):
        print(f"ERROR: PPLM weights not found at {model_path}")
        sys.exit(1)

    model_data = torch.load(model_path, map_location="cpu")
    model_param = model_data["param"]
    model_state = model_data["model"]

    pplm_model = PPLM(
        num_layers=model_param["encoder_layers"],
        embed_dim=model_param["encoder_embed_dim"],
        attention_heads=model_param["encoder_attention_heads"],
        token_dropout=False,
        alphabet=alphabet,
    )
    pplm_model.to(device)
    pplm_model.load_state_dict(model_state, strict=False)
    pplm_model.eval()

    # Load plant-trained PPI classifier weights
    print(f"Loading plant-trained PPI weights: {model_weights_path}")
    ppi_weights = torch.load(model_weights_path, map_location=device)

    ppi_model = PPI_inter_intra_attn_embed_single_pooling()
    ppi_model.to(device)
    ppi_model.eval()

    # Report ensemble composition
    total_models = sum(len(v) for v in ppi_weights.items())
    for mode, weights_list in ppi_weights.items():
        print(f"  {mode} models: {len(weights_list)}")

    print(f"Models loaded. Device: {device}")
    return pplm_model, batch_converter, ppi_model, ppi_weights


def get_features(model, batch_converter, seqA, seqB, device, max_pair_len=1020):
    """
    Extract PPLM features for a protein pair.
    Identical to scripts/benchmarking/batch_predict.py::get_features().
    """
    total_len = len(seqA) + len(seqB)
    if total_len > max_pair_len:
        ratio_A = len(seqA) / total_len
        budget_A = max(50, int(max_pair_len * ratio_A))
        budget_B = max_pair_len - budget_A
        seqA = seqA[:budget_A]
        seqB = seqB[:budget_B]

    _, _, seqA_tokens = batch_converter([("seqA", seqA)])
    _, _, seqB_tokens = batch_converter([("seqB", seqB)])
    tokens = torch.cat([seqA_tokens, seqB_tokens], dim=-1).to(device)

    lenA = len(seqA)
    lenB = len(seqB)
    total_tok_len = lenA + 2 + lenB + 2
    inter_chain_mask = torch.ones((total_tok_len, total_tok_len), device=device)
    inter_chain_mask[: lenA + 2, : lenA + 2] = 0
    inter_chain_mask[lenA + 2 :, lenA + 2 :] = 0

    with torch.no_grad():
        out = model(
            tokens, inter_chain_mask,
            repr_layers=[33], need_head_weights=True, return_contacts=False,
        )

        embed_A = out["representations"][33][0, 1 : lenA + 1, :]
        embed_B = out["representations"][33][0, -(lenB + 1) : -1, :]

        attns = out["attentions"].squeeze()
        n_heads = 33 * 20
        attn_AA = attns[:, :, 1 : lenA + 1, 1 : lenA + 1].reshape(n_heads, lenA, lenA)
        attn_AB = attns[:, :, 1 : lenA + 1, -(lenB + 1) : -1].reshape(n_heads, lenA, lenB)
        attn_BA = attns[:, :, -(lenB + 1) : -1, 1 : lenA + 1].reshape(n_heads, lenB, lenA)
        attn_BB = attns[:, :, -(lenB + 1) : -1, -(lenB + 1) : -1].reshape(n_heads, lenB, lenB)
        inter_attn = (attn_AB + attn_BA.transpose(1, 2)) / 2

        features = {
            "mean_inter_attn": inter_attn.mean(dim=[1, 2]).unsqueeze(0),
            "mean_attn_AA": attn_AA.mean(dim=[1, 2]).unsqueeze(0),
            "mean_attn_BB": attn_BB.mean(dim=[1, 2]).unsqueeze(0),
            "mean_embed_A": embed_A.mean(dim=0).unsqueeze(0),
            "mean_embed_B": embed_B.mean(dim=0).unsqueeze(0),
            "max_inter_attn": torch.amax(inter_attn, dim=(1, 2)).unsqueeze(0),
            "max_attn_AA": torch.amax(attn_AA, dim=(1, 2)).unsqueeze(0),
            "max_attn_BB": torch.amax(attn_BB, dim=(1, 2)).unsqueeze(0),
            "max_embed_A": torch.amax(embed_A, dim=0).unsqueeze(0),
            "max_embed_B": torch.amax(embed_B, dim=0).unsqueeze(0),
        }

        del out, attns, attn_AA, attn_AB, attn_BA, attn_BB
        del inter_attn, embed_A, embed_B, tokens, inter_chain_mask

    return features


def predict_with_plant_weights(ppi_model, features, ppi_weights):
    """
    Predict using the plant-trained ensemble with symmetric pair averaging.
    Supports any combination of pooling modes present in ppi_weights.
    """
    with torch.no_grad():
        predictions_list = []

        # Iterate over all pooling modes in the ensemble
        for mode in ["mean", "max", "min"]:
            if mode not in ppi_weights:
                continue

            prefix = mode
            inter_key = f"{prefix}_inter_attn"
            aa_key = f"{prefix}_attn_AA"
            bb_key = f"{prefix}_attn_BB"
            ea_key = f"{prefix}_embed_A"
            eb_key = f"{prefix}_embed_B"

            # Skip if features for this mode are not available
            if inter_key not in features:
                continue

            for model_weight in ppi_weights[mode]:
                ppi_model.load_state_dict(model_weight)

                # Forward pass
                pred = ppi_model(
                    features[inter_key], features[aa_key], features[bb_key],
                    features[ea_key], features[eb_key],
                )
                # Swapped order (symmetric averaging)
                pred_swap = ppi_model(
                    features[inter_key], features[bb_key], features[aa_key],
                    features[eb_key], features[ea_key],
                )
                pred_sym = (pred + pred_swap) / 2
                predictions_list.append(pred_sym)

        # Average across all ensemble models
        predictions = torch.stack(predictions_list)
        final_score = torch.mean(predictions, dim=0).squeeze().detach().cpu().numpy()
        return float(final_score)


def main():
    parser = argparse.ArgumentParser(
        description="Test plant-trained PPLM-PPI on held-out datasets"
    )
    parser.add_argument(
        "--input", required=True,
        help="Input TSV file (Protein1\\tProtein2\\tlabel)",
    )
    parser.add_argument(
        "--output", required=True, help="Output CSV file"
    )
    parser.add_argument(
        "--model_weights", required=True,
        help="Path to plant-trained ensemble weights (ppi_plant_models.pkl)",
    )
    parser.add_argument(
        "--seq_db", required=True, help="Path to pickled sequence dictionary"
    )
    parser.add_argument(
        "--gpu_id", type=int, default=0, help="GPU device ID"
    )
    parser.add_argument(
        "--max_pair_len", type=int, default=1020,
        help="Maximum combined sequence pair length",
    )
    args = parser.parse_args()

    # Setup
    device = torch.device(
        f"cuda:{args.gpu_id}" if torch.cuda.is_available() else "cpu"
    )
    print(f"Device: {device}")

    # Load models
    pplm_model, batch_converter, ppi_model, ppi_weights = load_models(
        device, args.model_weights
    )

    # Load sequence database
    print(f"Loading sequence database: {args.seq_db}")
    with open(args.seq_db, "rb") as f:
        seq_db = pickle.load(f)
    print(f"Loaded {len(seq_db)} sequences")

    # Read input
    print(f"Reading input: {args.input}")
    rows = []
    with open(args.input, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("Protein"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            rows.append((parts[0], parts[1], int(parts[2])))
    print(f"Total pairs: {len(rows)}")

    # Predict
    print("Starting prediction...")
    start_time = time.time()
    results = []
    skipped = 0

    for i, (protA, protB, label) in enumerate(rows):
        seqA = seq_db.get(protA)
        seqB = seq_db.get(protB)

        if seqA is None or seqB is None:
            print(f"WARNING: Missing sequence for {protA} or {protB}, skipping")
            skipped += 1
            continue

        try:
            features = get_features(
                pplm_model, batch_converter, seqA, seqB, device, args.max_pair_len
            )
            score = predict_with_plant_weights(ppi_model, features, ppi_weights)
            results.append((f"{protA}:{protB}", score, label))
        except torch.cuda.OutOfMemoryError:
            torch.cuda.empty_cache()
            try:
                features = get_features(
                    pplm_model, batch_converter, seqA, seqB, device, max_pair_len=800
                )
                score = predict_with_plant_weights(ppi_model, features, ppi_weights)
                results.append((f"{protA}:{protB}", score, label))
            except Exception as e:
                print(f"WARNING: Error {protA}:{protB} ({e}), skipping")
                skipped += 1
                torch.cuda.empty_cache()

        # Periodic cleanup and progress
        if (i + 1) % 200 == 0 and torch.cuda.is_available():
            torch.cuda.empty_cache()

        if (i + 1) % 1000 == 0 or (i + 1) >= len(rows):
            elapsed = time.time() - start_time
            rate = (i + 1) / max(0.1, elapsed)
            print(f"  Progress: {i + 1}/{len(rows)} ({rate:.1f} pairs/sec)")

    # Save results
    print(f"Saving results to: {args.output}")
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["pair_id", "pred_score", "true_label"])
        writer.writerows(results)

    elapsed = time.time() - start_time
    print(
        f"Done! Predicted {len(results)} pairs in {elapsed:.1f}s "
        f"({len(results) / elapsed:.1f} pairs/sec)"
    )
    print(f"Skipped: {skipped}")


if __name__ == "__main__":
    main()
