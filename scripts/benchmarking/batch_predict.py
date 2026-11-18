#!/usr/bin/env python3
"""
PPLM Batch PPI Prediction Script.
Predicts protein-protein interactions for a dataset of protein pairs.

Supports:
  - Zero-shot prediction (default: PPLM/weights/ppi_models.pkl)
  - Plant-retrained ensemble weights (--model_weights models/.../ppi_plant_models.pkl)
  - Any pooling combination via --modes (e.g. mean max min)
  - Optional pre-extracted feature directory (--feat_dir) for fast batched inference

Usage:
    python scripts/benchmarking/batch_predict.py \
        --input data/legacy_benchmarks/DeepAraPPI/c2_ppi_sample_DeepAraPPI.txt \
        --output results/DeepAraPPI/zero-shot/deepara_c2_scores.csv \
        --batch_size 4 --gpu_id 0 \
        --seq_db data/arabidopsis/uniprot_final.pkl
"""

import argparse
import csv
import os
import pickle
import sys
import time

import numpy as np
import torch
from torch.utils.data import DataLoader

# Path setup: PROJECT_ROOT is two levels up from scripts/benchmarking/
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))
PPLM_DIR = os.path.join(PROJECT_ROOT, "PPLM")
TRAINING_DIR = os.path.join(PROJECT_ROOT, "scripts", "training")
DATA_DIR = os.path.join(PROJECT_ROOT, "data")

sys.path.insert(0, PPLM_DIR)
sys.path.insert(0, TRAINING_DIR)

print(f"Script directory: {SCRIPT_DIR}")
print(f"Project root:     {PROJECT_ROOT}")
print(f"PPLM directory:   {PPLM_DIR}")
print(f"Data directory:   {DATA_DIR}")


def load_ppi_head(device, model_weights_path=None):
    """Load PPI classifier model and ensemble weights dictionary."""
    from pplm_ppi import PPLM_PPI

    ppi_path = model_weights_path or os.path.join(PPLM_DIR, "weights", "ppi_models.pkl")
    if not os.path.exists(ppi_path):
        print(f"ERROR: PPI weights not found at {ppi_path}")
        sys.exit(1)

    print(f"Loading PPI classifier weights: {ppi_path}")
    ppi_weights = torch.load(ppi_path, map_location=device, weights_only=False)

    ppi_model = PPLM_PPI()
    ppi_model.to(device)
    ppi_model.eval()
    return ppi_model, ppi_weights


def load_models(device, model_weights_path=None):
    """Load frozen PPLM backbone and PPI classifier weights."""
    from pplm import Alphabet, PPLM

    print("Loading PPLM backbone...")
    alphabet = Alphabet.from_architecture()
    batch_converter = alphabet.get_batch_converter()
    model_path = os.path.join(PPLM_DIR, "weights", "pplm_t33_650M.pt")
    if not os.path.exists(model_path):
        print(f"ERROR: PPLM backbone weights not found at {model_path}")
        sys.exit(1)

    model_data = torch.load(model_path, map_location="cpu", weights_only=False)
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

    ppi_model, ppi_weights = load_ppi_head(device, model_weights_path)
    print(f"Models loaded. Device: {device}")
    return pplm_model, batch_converter, ppi_model, ppi_weights


def get_features(model, batch_converter, seqA, seqB, device, max_pair_len=1020):
    """
    Extract mean, max, and min pooled PPLM features for a protein pair.
    Enforces (len(seqA) + len(seqB) <= max_pair_len) using proportional cropping.
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
    inter_chain_mask = torch.ones(
        (lenA + 2 + lenB + 2, lenA + 2 + lenB + 2),
        device=device,
    )
    inter_chain_mask[: lenA + 2, : lenA + 2] = 0
    inter_chain_mask[lenA + 2 :, lenA + 2 :] = 0

    with torch.no_grad():
        out = model(
            tokens,
            inter_chain_mask,
            repr_layers=[33],
            need_head_weights=True,
            return_contacts=False,
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
            # Mean pooling
            "mean_inter_attn": inter_attn.mean(dim=[1, 2]).unsqueeze(0),
            "mean_attn_AA": attn_AA.mean(dim=[1, 2]).unsqueeze(0),
            "mean_attn_BB": attn_BB.mean(dim=[1, 2]).unsqueeze(0),
            "mean_embed_A": embed_A.mean(dim=0).unsqueeze(0),
            "mean_embed_B": embed_B.mean(dim=0).unsqueeze(0),
            # Max pooling
            "max_inter_attn": torch.amax(inter_attn, dim=(1, 2)).unsqueeze(0),
            "max_attn_AA": torch.amax(attn_AA, dim=(1, 2)).unsqueeze(0),
            "max_attn_BB": torch.amax(attn_BB, dim=(1, 2)).unsqueeze(0),
            "max_embed_A": torch.amax(embed_A, dim=0).unsqueeze(0),
            "max_embed_B": torch.amax(embed_B, dim=0).unsqueeze(0),
            # Min pooling
            "min_inter_attn": torch.amin(inter_attn, dim=(1, 2)).unsqueeze(0),
            "min_attn_AA": torch.amin(attn_AA, dim=(1, 2)).unsqueeze(0),
            "min_attn_BB": torch.amin(attn_BB, dim=(1, 2)).unsqueeze(0),
            "min_embed_A": torch.amin(embed_A, dim=0).unsqueeze(0),
            "min_embed_B": torch.amin(embed_B, dim=0).unsqueeze(0),
        }

        del (
            out,
            attns,
            attn_AA,
            attn_AB,
            attn_BA,
            attn_BB,
            inter_attn,
            embed_A,
            embed_B,
            tokens,
            inter_chain_mask,
        )

    return features


def predict_with_weights(ppi_model, features, ppi_weights, modes=("mean", "max"), device=None):
    """Predict using the ensemble across the requested pooling modes with symmetric averaging."""
    with torch.no_grad():
        predictions_list = []

        for mode in modes:
            if mode not in ppi_weights:
                continue

            inter_attn = features[f"{mode}_inter_attn"]
            attn_aa = features[f"{mode}_attn_AA"]
            attn_bb = features[f"{mode}_attn_BB"]
            embed_a = features[f"{mode}_embed_A"]
            embed_b = features[f"{mode}_embed_B"]

            if device is not None:
                inter_attn = inter_attn.to(device)
                attn_aa = attn_aa.to(device)
                attn_bb = attn_bb.to(device)
                embed_a = embed_a.to(device)
                embed_b = embed_b.to(device)

            for model_weight in ppi_weights[mode]:
                ppi_model.load_state_dict(model_weight)
                pred = ppi_model(inter_attn, attn_aa, attn_bb, embed_a, embed_b)
                pred_swap = ppi_model(inter_attn, attn_bb, attn_aa, embed_b, embed_a)
                predictions_list.append((pred + pred_swap) / 2)

        if not predictions_list:
            raise ValueError(f"No ensemble weights found for requested modes: {modes}")

        stacked = torch.stack(predictions_list)
        return torch.mean(stacked, dim=0).view(-1)


def run_batched_from_feat_dir(args, device):
    """Fast batched prediction using pre-extracted .pkl test features."""
    from dataset import Dataset_PPI_Interact

    ppi_model, ppi_weights = load_ppi_head(device, args.model_weights)
    dataset = Dataset_PPI_Interact(
        args.input,
        esm_dir=args.feat_dir,
        order_mode="both",
        testing=True,
    )
    loader_kwargs = {
        "batch_size": args.batch_size,
        "shuffle": False,
        "num_workers": args.num_workers,
        "pin_memory": torch.cuda.is_available(),
    }
    if args.num_workers > 0:
        loader_kwargs["prefetch_factor"] = 4
    dataloader = DataLoader(dataset, **loader_kwargs)

    results = []
    start_time = time.time()
    for data in dataloader:
        scores = (
            predict_with_weights(ppi_model, data, ppi_weights, modes=args.modes, device=device)
            .detach()
            .cpu()
            .numpy()
            .tolist()
        )
        labels = data["label"].numpy().astype(int).tolist()
        for sample_name, score, label in zip(data["sample"], scores, labels):
            results.append((sample_name.replace("@", ":"), float(score), int(label)))

    elapsed = time.time() - start_time
    return results, 0, elapsed


def main():
    parser = argparse.ArgumentParser(description="PPLM Batch PPI Prediction")
    parser.add_argument(
        "--input", required=True, help="Input TSV file (Protein1\\tProtein2\\tlabel)"
    )
    parser.add_argument("--output", required=True, help="Output CSV file")
    parser.add_argument(
        "--feat_dir",
        default=None,
        help="Optional directory of pre-extracted .pkl features for fast batched inference",
    )
    parser.add_argument(
        "--seq_db",
        default=None,
        help="Path to pickled sequence dictionary (required if --feat_dir is not provided)",
    )
    parser.add_argument(
        "--model_weights",
        default=None,
        help="Path to PPI head weights (defaults to PPLM/weights/ppi_models.pkl)",
    )
    parser.add_argument(
        "--modes",
        nargs="+",
        default=["mean", "max"],
        choices=["mean", "max", "min"],
        help="Pooling modes to include from model_weights (default: mean max)",
    )
    parser.add_argument(
        "--batch_size", type=int, default=4, help="Batch size for inference (use 512 with --feat_dir)"
    )
    parser.add_argument("--num_workers", type=int, default=4, help="DataLoader workers for --feat_dir")
    parser.add_argument("--gpu_id", type=int, default=0, help="GPU device ID")
    parser.add_argument(
        "--max_pair_len",
        type=int,
        default=1020,
        help="Maximum combined sequence pair length (LA + LB <= 1020) matching PPLM pretraining",
    )
    args = parser.parse_args()

    if not args.feat_dir and not args.seq_db:
        parser.error("Must specify either --feat_dir or --seq_db.")

    device = torch.device(f"cuda:{args.gpu_id}" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    if args.feat_dir:
        results, skipped, elapsed = run_batched_from_feat_dir(args, device)
    else:
        pplm_model, batch_converter, ppi_model, ppi_weights = load_models(
            device, args.model_weights
        )

        print(f"Loading sequence database: {args.seq_db}")
        with open(args.seq_db, "rb") as f:
            seq_db = pickle.load(f)
        print(f"Loaded {len(seq_db)} sequences")

        print(f"Reading input: {args.input}")
        rows = []
        with open(args.input, "r") as f:
            for line in f:
                line = line.strip()
                if (
                    not line
                    or line.startswith("#")
                    or line.startswith("Protein")
                    or line.startswith("participant_a_id")
                ):
                    continue
                parts = line.split("\t")
                if len(parts) < 3 or parts[2] == "label":
                    continue
                rows.append((parts[0], parts[1], int(parts[2])))
        print(f"Total rows: {len(rows)}")

        print("Starting prediction...")
        start_time = time.time()
        results = []
        skipped = 0

        for i in range(0, len(rows), args.batch_size):
            batch = rows[i : i + args.batch_size]
            batch_scores = []

            for protA, protB, label in batch:
                seqA = seq_db.get(protA, None)
                seqB = seq_db.get(protB, None)

                if seqA is None or seqB is None:
                    print(f"WARNING: Missing sequence for {protA} or {protB}, skipping")
                    skipped += 1
                    continue

                try:
                    features = get_features(
                        pplm_model,
                        batch_converter,
                        seqA,
                        seqB,
                        device,
                        max_pair_len=args.max_pair_len,
                    )
                    score = float(
                        predict_with_weights(
                            ppi_model, features, ppi_weights, modes=args.modes
                        )
                        .squeeze()
                        .cpu()
                        .item()
                    )
                    batch_scores.append((f"{protA}:{protB}", score, label))
                except torch.cuda.OutOfMemoryError:
                    torch.cuda.empty_cache()
                    try:
                        features = get_features(
                            pplm_model, batch_converter, seqA, seqB, device, max_pair_len=800
                        )
                        score = float(
                            predict_with_weights(
                                ppi_model, features, ppi_weights, modes=args.modes
                            )
                            .squeeze()
                            .cpu()
                            .item()
                        )
                        batch_scores.append((f"{protA}:{protB}", score, label))
                    except Exception as e:
                        print(f"WARNING: Error predicting pair {protA}:{protB} ({e}), skipping")
                        skipped += 1
                        torch.cuda.empty_cache()

            results.extend(batch_scores)

            if i % 200 == 0 and torch.cuda.is_available():
                torch.cuda.empty_cache()

            if (i + len(batch)) % 1000 == 0 or (i + len(batch)) >= len(rows):
                elapsed_now = time.time() - start_time
                rate = (i + len(batch)) / max(0.1, elapsed_now)
                print(f"  Progress: {i + len(batch)}/{len(rows)} ({rate:.1f} pairs/sec)")

        elapsed = time.time() - start_time

    print(f"Saving results to: {args.output}")
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["pair_id", "pred_score", "true_label"])
        writer.writerows(results)

    print(
        f"Done! Predicted {len(results)} pairs in {elapsed:.1f} seconds "
        f"({len(results) / max(0.01, elapsed):.1f} pairs/sec)"
    )
    print(f"Skipped {skipped} pairs due to missing sequences")


if __name__ == "__main__":
    main()
