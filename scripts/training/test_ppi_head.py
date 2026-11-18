#!/usr/bin/env python3
"""
Evaluate plant-trained PPLM-PPI models on held-out test sets (C2, C3, Rice).

Aligned with the original PPLM author's multi-model multi-pooling testing script
(ppi_training_code/testing.PPI.single_gpu.multi-model.multi_pooling.py):
  - Supports loading ensembles via --mean_path / --max_path / --min_path
    (top-5 checkpoint path lists) OR packaged --model_weights (ppi_plant_models.pkl)
  - Primary fast batched mode (--feat_dir): evaluates pre-extracted .pkl test
    features with DataLoader(batch_size=512, num_workers=8) in seconds
  - Fallback on-the-fly mode (--seq_db): extracts mean/max/min features directly
    from the frozen PPLM backbone if --feat_dir is omitted
  - Supports --order_mode (both, normal, reverse, random; default: both)
  - Outputs both pair_id,pred_score,true_label CSV (--output) and optional
    author-format .preds.csv / .preds.statistic files (--out_path_prefix)

Usage:
    # Fast batched evaluation (matching run_testing.sh):
    python scripts/training/test_ppi_head.py \
        --input data/legacy_benchmarks/DeepAraPPI/c2_ppi_sample_DeepAraPPI.txt \
        --feat_dir features/DeepAraPPI_C2/ \
        --mean_path models/DeepAraPPI/PPLM-PPI.inter_intra_attn_embed.mean.round_1.top5_list \
        --max_path  models/DeepAraPPI/PPLM-PPI.inter_intra_attn_embed.max.round_1.top5_list \
        --min_path  models/DeepAraPPI/PPLM-PPI.inter_intra_attn_embed.min.round_1.top5_list \
        --output results/DeepAraPPI/pretrained/deepara_c2_combine_mean_max_min_scores.csv \
        --batch_size 512 --order_mode both --gpu_id 0

    # Or using packaged weights:
    python scripts/training/test_ppi_head.py \
        --input data/legacy_benchmarks/DeepAraPPI/c2_ppi_sample_DeepAraPPI.txt \
        --feat_dir features/DeepAraPPI_C2/ \
        --model_weights models/DeepAraPPI/ppi_plant_models.pkl \
        --modes mean max min \
        --output results/DeepAraPPI/pretrained/deepara_c2_plant_scores.csv
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

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))
PPLM_DIR = os.path.join(PROJECT_ROOT, "PPLM")
sys.path.insert(0, PPLM_DIR)
sys.path.insert(0, SCRIPT_DIR)

from dataset import Dataset_PPI_Interact, POOLING_MODES
from ppi_model import PPI_inter_intra_attn_embed_single_pooling, evaluate


def load_ensemble_from_paths(device, mean_path="", max_path="", min_path=""):
    """Load top-K checkpoint models per pooling mode from path-list files."""
    path_map = {"mean": mean_path, "max": max_path, "min": min_path}
    model_lists = {"mean": [], "max": [], "min": []}

    for mode, list_file in path_map.items():
        if not list_file:
            continue
        with open(list_file, "r") as f:
            ckpt_paths = [line.strip().split("\t")[0] for line in f if line.strip()]

        print(f"Loading {len(ckpt_paths)} {mode} models from: {list_file}")
        for ckpt_path in ckpt_paths:
            model = PPI_inter_intra_attn_embed_single_pooling()
            checkpoint = torch.load(ckpt_path, map_location=device, weights_only=False)
            state_dict = checkpoint.get("model_state_dict", checkpoint)
            model.load_state_dict(state_dict)
            model.to(device)
            model.eval()
            model_lists[mode].append(model)

    return model_lists


def load_ensemble_from_weights(device, model_weights_path, modes):
    """Load ensemble models from a packaged ppi_plant_models.pkl dictionary."""
    print(f"Loading packaged PPI weights: {model_weights_path}")
    ppi_weights = torch.load(model_weights_path, map_location=device, weights_only=False)
    model_lists = {"mean": [], "max": [], "min": []}

    for mode in modes:
        if mode not in ppi_weights:
            continue
        for state_dict in ppi_weights[mode]:
            model = PPI_inter_intra_attn_embed_single_pooling()
            model.load_state_dict(state_dict)
            model.to(device)
            model.eval()
            model_lists[mode].append(model)
        print(f"  {mode} models loaded: {len(model_lists[mode])}")

    return model_lists


def load_pplm_backbone(device):
    """Load frozen PPLM backbone for fallback on-the-fly feature extraction."""
    from pplm import Alphabet, PPLM

    print("Loading PPLM backbone for on-the-fly feature extraction...")
    alphabet = Alphabet.from_architecture()
    batch_converter = alphabet.get_batch_converter()
    model_path = os.path.join(PPLM_DIR, "weights", "pplm_t33_650M.pt")

    if not os.path.exists(model_path):
        print(f"ERROR: PPLM weights not found at {model_path}")
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
    return pplm_model, batch_converter


def get_features(model, batch_converter, seqA, seqB, device, max_pair_len=1020):
    """Extract mean/max/min PPLM features for a single pair on-the-fly."""
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

        del out, attns, attn_AA, attn_AB, attn_BA, attn_BB
        del inter_attn, embed_A, embed_B, tokens, inter_chain_mask

    return features


def predict_batch(model_lists, batch_data, device, order_mode="both"):
    """
    Run all loaded models across pooling modes on a batch (or single pair).
    Applies symmetric pair averaging when order_mode == 'both'.
    """
    with torch.no_grad():
        predictions_list = []

        for mode in POOLING_MODES:
            models = model_lists.get(mode, [])
            if not models:
                continue

            inter_attn = batch_data[f"{mode}_inter_attn"].to(device)
            attn_aa = batch_data[f"{mode}_attn_AA"].to(device)
            attn_bb = batch_data[f"{mode}_attn_BB"].to(device)
            embed_a = batch_data[f"{mode}_embed_A"].to(device)
            embed_b = batch_data[f"{mode}_embed_B"].to(device)

            for model in models:
                pred = model(inter_attn, attn_aa, attn_bb, embed_a, embed_b)
                if order_mode == "both":
                    pred_swap = model(inter_attn, attn_bb, attn_aa, embed_b, embed_a)
                    pred = (pred + pred_swap) / 2
                predictions_list.append(pred)

        stacked = torch.stack(predictions_list)
        return torch.mean(stacked, dim=0).view(-1)


def run_batched_feature_testing(args, model_lists, device):
    """Fast batched evaluation over pre-extracted .pkl test features (author method)."""
    test_samples = Dataset_PPI_Interact(
        args.input,
        esm_dir=args.feat_dir,
        esm_flag=False,
        order_mode=args.order_mode,
        testing=True,
        debug=args.debug,
        n_debug_target=args.n_debug_target,
    )
    loader_kwargs = {
        "batch_size": args.batch_size,
        "shuffle": False,
        "num_workers": args.num_workers,
        "pin_memory": torch.cuda.is_available(),
    }
    if args.num_workers > 0:
        loader_kwargs["prefetch_factor"] = 4

    test_dataloader = DataLoader(test_samples, **loader_kwargs)

    results = []
    batch_losses = []
    all_preds = []
    all_labels = []

    with torch.no_grad():
        for data in test_dataloader:
            labels = data["label"].to(device).view(-1).float()
            samples = data["sample"]

            preds = predict_batch(model_lists, data, device, order_mode=args.order_mode)
            loss = torch.nn.functional.binary_cross_entropy(preds, labels)
            batch_losses.append(loss.item())

            preds_np = preds.detach().cpu().numpy().tolist()
            labels_np = labels.detach().cpu().numpy().astype(int).tolist()

            all_preds.extend(preds_np)
            all_labels.extend(labels_np)

            for sample_name, score, lbl in zip(samples, preds_np, labels_np):
                results.append((sample_name, float(score), int(lbl)))

    avg_loss = float(np.mean(batch_losses)) if batch_losses else float("nan")
    return results, all_preds, all_labels, avg_loss


def run_onthefly_testing(args, model_lists, device):
    """Fallback pair-by-pair evaluation using the frozen PPLM backbone."""
    pplm_model, batch_converter = load_pplm_backbone(device)

    print(f"Loading sequence database: {args.seq_db}")
    with open(args.seq_db, "rb") as f:
        seq_db = pickle.load(f)
    print(f"Loaded {len(seq_db)} sequences")

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
            if args.debug and len(rows) >= args.n_debug_target:
                break

    print(f"Total pairs to evaluate on-the-fly: {len(rows)}")
    start_time = time.time()
    results = []
    all_preds = []
    all_labels = []
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
            score = float(
                predict_batch(model_lists, features, device, order_mode=args.order_mode)
                .squeeze()
                .cpu()
                .item()
            )
            results.append((f"{protA}@{protB}", score, label))
            all_preds.append(score)
            all_labels.append(label)
        except torch.cuda.OutOfMemoryError:
            torch.cuda.empty_cache()
            try:
                features = get_features(
                    pplm_model, batch_converter, seqA, seqB, device, max_pair_len=800
                )
                score = float(
                    predict_batch(
                        model_lists, features, device, order_mode=args.order_mode
                    )
                    .squeeze()
                    .cpu()
                    .item()
                )
                results.append((f"{protA}@{protB}", score, label))
                all_preds.append(score)
                all_labels.append(label)
            except Exception as e:
                print(f"WARNING: Error {protA}:{protB} ({e}), skipping")
                skipped += 1
                torch.cuda.empty_cache()

        if (i + 1) % 200 == 0 and torch.cuda.is_available():
            torch.cuda.empty_cache()
        if (i + 1) % 1000 == 0 or (i + 1) >= len(rows):
            elapsed = time.time() - start_time
            rate = (i + 1) / max(0.1, elapsed)
            print(f"  Progress: {i + 1}/{len(rows)} ({rate:.1f} pairs/sec)")

    if skipped > 0:
        print(f"Skipped {skipped} pairs due to missing sequences/errors.")

    avg_loss = float("nan")
    if all_preds:
        preds_t = torch.tensor(all_preds, dtype=torch.float32)
        labels_t = torch.tensor(all_labels, dtype=torch.float32)
        avg_loss = float(torch.nn.functional.binary_cross_entropy(preds_t, labels_t).item())

    return results, all_preds, all_labels, avg_loss


def main():
    parser = argparse.ArgumentParser(
        description="Test plant-trained PPLM-PPI models (aligned with author's multi-model multi-pooling testing)"
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input test file (TSV Protein1\\tProtein2\\tlabel or sample ID list)",
    )
    parser.add_argument(
        "--feat_dir",
        default=None,
        help="Directory of pre-extracted .pkl features for fast batched inference (recommended)",
    )
    parser.add_argument(
        "--seq_db",
        default=None,
        help="Path to pickled sequence dictionary (fallback on-the-fly mode if --feat_dir omitted)",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output CSV path (pair_id,pred_score,true_label)",
    )
    parser.add_argument(
        "--out_path_prefix",
        default=None,
        help="Optional prefix to write author-format .preds.csv and .preds.statistic",
    )
    parser.add_argument(
        "--test_name",
        default=None,
        help="Optional test set label for author-format outputs",
    )
    parser.add_argument(
        "--mean_path",
        "-mean_path",
        type=str,
        default="",
        help="Path to top-K checkpoint list for mean pooling",
    )
    parser.add_argument(
        "--max_path",
        "-max_path",
        type=str,
        default="",
        help="Path to top-K checkpoint list for max pooling",
    )
    parser.add_argument(
        "--min_path",
        "-min_path",
        type=str,
        default="",
        help="Path to top-K checkpoint list for min pooling",
    )
    parser.add_argument(
        "--model_weights",
        default=None,
        help="Path to packaged ensemble weights (ppi_plant_models.pkl)",
    )
    parser.add_argument(
        "--modes",
        nargs="+",
        default=["mean", "max", "min"],
        choices=["mean", "max", "min"],
        help="Pooling modes to use when loading --model_weights (default: mean max min)",
    )
    parser.add_argument(
        "--order_mode",
        "-order_mode",
        type=str,
        default="both",
        choices=["both", "normal", "reverse", "random"],
        help="Order averaging mode (default: both)",
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        default=512,
        help="Batch size for pre-extracted feature testing (default: 512)",
    )
    parser.add_argument(
        "--num_workers",
        type=int,
        default=8,
        help="DataLoader worker count for batched testing (default: 8)",
    )
    parser.add_argument(
        "--gpu_id",
        "-gpu",
        type=int,
        default=0,
        help="GPU device ID (default: 0)",
    )
    parser.add_argument(
        "--max_pair_len",
        type=int,
        default=1020,
        help="Maximum combined sequence pair length for on-the-fly extraction",
    )
    parser.add_argument(
        "--debug",
        "-debug",
        action="store_true",
        default=False,
        help="Debug mode: limit test samples",
    )
    parser.add_argument(
        "--n_debug_target",
        type=int,
        default=1000,
        help="Max samples in debug mode (default: 1000)",
    )
    args = parser.parse_args()

    if not args.output and not args.out_path_prefix:
        parser.error("Must specify at least one of --output or --out_path_prefix.")
    if not args.feat_dir and not args.seq_db:
        parser.error("Must specify either --feat_dir (fast batched) or --seq_db (on-the-fly).")

    has_path_lists = any([args.mean_path, args.max_path, args.min_path])
    if not has_path_lists and not args.model_weights:
        parser.error("Must specify --mean_path/--max_path/--min_path or --model_weights.")

    device = torch.device(f"cuda:{args.gpu_id}" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device} | order_mode={args.order_mode}")

    if has_path_lists:
        model_lists = load_ensemble_from_paths(
            device, args.mean_path, args.max_path, args.min_path
        )
    else:
        model_lists = load_ensemble_from_weights(device, args.model_weights, args.modes)

    total_models = sum(len(v) for v in model_lists.values())
    if total_models == 0:
        print("ERROR: No ensemble models were loaded.")
        sys.exit(1)
    print(
        f"Loaded {total_models} ensemble models ("
        + ", ".join(f"{k}={len(v)}" for k, v in model_lists.items())
        + ")"
    )

    start_time = time.time()
    if args.feat_dir:
        results, all_preds, all_labels, avg_loss = run_batched_feature_testing(
            args, model_lists, device
        )
    else:
        results, all_preds, all_labels, avg_loss = run_onthefly_testing(
            args, model_lists, device
        )
    elapsed = time.time() - start_time

    if not results:
        print("ERROR: Zero pairs predicted.")
        sys.exit(1)

    # Save evaluate_pplm.py-compatible CSV
    if args.output:
        os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
        with open(args.output, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["pair_id", "pred_score", "true_label"])
            for sample_name, score, label in results:
                writer.writerow([sample_name.replace("@", ":"), score, label])
        print(f"Saved prediction CSV to: {args.output}")

    # Compute author evaluation metrics
    (
        precision,
        recall,
        accuracy,
        F1_score,
        Top10,
        Top20,
        Top50,
        AUC_ROC,
        AUC_PR,
    ) = evaluate(all_labels, all_preds)

    test_set_label = args.test_name or os.path.splitext(os.path.basename(args.input))[0]
    stat_line = (
        f"========== {test_set_label} | Loss: {avg_loss} | Accuarcy: "
        f"{precision} {recall} {accuracy} {F1_score} "
        f"{Top10} {Top20} {Top50} {AUC_ROC} {AUC_PR}"
    )
    print(stat_line)
    print(
        f"Completed {len(results)} pairs in {elapsed:.1f}s "
        f"({len(results) / max(0.01, elapsed):.1f} pairs/sec) | "
        f"AUC_PR={AUC_PR:.4f}, AUC_ROC={AUC_ROC:.4f}"
    )

    # Optional author-format .preds.csv and .preds.statistic
    if args.out_path_prefix:
        os.makedirs(os.path.dirname(args.out_path_prefix) or ".", exist_ok=True)
        preds_csv_path = f"{args.out_path_prefix}.{test_set_label}.preds.csv"
        stat_path = f"{args.out_path_prefix}.preds.statistic"

        with open(preds_csv_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["name", "label", "prediction"])
            for sample_name, score, label in results:
                writer.writerow([sample_name, label, score])

        with open(stat_path, "w") as f:
            f.write(f"output_statistic: {stat_path}\n")
            f.write(f"output_csv: {preds_csv_path}\n")
            f.write(stat_line + "\n")

        print(f"Saved author-format CSV:  {preds_csv_path}")
        print(f"Saved author statistics:  {stat_path}")


if __name__ == "__main__":
    main()
