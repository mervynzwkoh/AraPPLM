#!/usr/bin/env python3
"""
Train the PPLM-PPI MLP head on plant C1 training data using K-fold
stratified cross-validation.

Aligned with the original PPLM author's training script
(ppi_training_code/training.PPI.single_gpu.cross_validation.py):
  - Pre-generated per-fold .trainset / .testset split files (Stratified 10-fold)
  - Supports per-fold execution (--cross_round 0..9) or sequential (--all_folds)
  - Supports all 3 pooling modes (mean, max, min)
  - Author seed table with --seed_round 1 (seed = 73591024) and cudnn.enabled = False
  - AdamW (lr=5e-5, weight_decay=0.01), BCE loss, 15 epochs, batch_size=512
  - Checkpoint resumption from model.pkl and exact .recording column format
    (AUC_PR at column 14 for `sort -rn -k14` compatibility)

Usage:
    # Single fold matching original author's run_training.sh:
    python scripts/training/train_ppi_head.py \
        models/DeepAraPPI/PPLM-PPI.inter_intra_attn_embed.mean \
        --feat_dir features/DeepAraPPI_C1/ \
        --split_dir models/DeepAraPPI/splits \
        --pair_list data/legacy_benchmarks/DeepAraPPI/c1_ppi_sample_DeepAraPPI.txt \
        --dataset_prefix deeparappi_c1 \
        --mode PPI_inter_intra_attn_embed_mean \
        --cross_round 0 --seed_round 1

    # Or run all 10 folds sequentially:
    python scripts/training/train_ppi_head.py \
        --feat_dir features/DeepAraPPI_C1/ \
        --output_dir models/DeepAraPPI/ \
        --mode mean --all_folds --seed_round 1
"""

import argparse
import os
import random
import sys
from os.path import isdir, isfile, join
from time import time

import numpy as np
import torch
from torch.utils.data import DataLoader

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

from dataset import Dataset_PPI_Interact, generate_cv_split_files
from ppi_model import PPI_inter_intra_attn_embed_single_pooling, evaluate


SEEDS = [
    26240761,
    73591024,
    48136297,
    95317406,
    62780953,
    58302741,
    94716038,
    12649570,
    80421359,
    39287564,
]

MODE_ALIASES = {
    "mean": "mean",
    "max": "max",
    "min": "min",
    "PPI_inter_intra_attn_embed_mean": "mean",
    "PPI_inter_intra_attn_embed_max": "max",
    "PPI_inter_intra_attn_embed_min": "min",
}


def seed_everything(seed=73591024):
    """Set random seeds and deterministic cuDNN flags matching author's script."""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.enabled = False


def make_dataloader(dataset, batch_size, shuffle, num_workers):
    """Build DataLoader using author's prefetch/pin_memory settings safely."""
    kwargs = {
        "batch_size": batch_size,
        "shuffle": shuffle,
        "num_workers": num_workers,
        "pin_memory": torch.cuda.is_available(),
    }
    if num_workers > 0:
        kwargs["prefetch_factor"] = 4
    return DataLoader(dataset, **kwargs)


def train_single_fold(
    fold_idx,
    trainset_path,
    validset_path,
    folder_prefix,
    short_mode,
    canonical_mode,
    seed,
    device,
    args,
):
    """Train and validate a single cross-validation fold with checkpoint resumption."""
    loss_name = f"{folder_prefix}.round_{args.seed_round}.cv_{fold_idx}"
    if args.linear_all:
        loss_name = loss_name + ".linear_all"
    out_recording_file = loss_name + ".recording"

    print(f"\n{'=' * 60}")
    print(f"FOLD {fold_idx}/{args.n_folds - 1} | mode={canonical_mode} | seed={seed}")
    print(f"Checkpoint dir: {loss_name}")
    print(f"{'=' * 60}")

    model = PPI_inter_intra_attn_embed_single_pooling()
    model.to(device)

    total_num = sum(p.numel() for p in model.parameters())
    trainable_num = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"Number of Network Parameters | Total: {total_num} Trainable: {trainable_num}")

    if isfile(join(loss_name, "model.pkl")):
        print("Loading a checkpoint...")
        checkpoint = torch.load(
            join(loss_name, "model.pkl"), map_location=device, weights_only=False
        )
        model.load_state_dict(checkpoint["model_state_dict"])
        start_epoch = checkpoint["epoch"] + 1
        train_loss = checkpoint["train_loss"]
        valid_loss = checkpoint["valid_loss"]
        best_models = checkpoint.get("best_models", [])
        if not args.silent:
            print("Restarting at epoch", start_epoch)
    else:
        print("Training a new model...")
        start_epoch = 1
        train_loss = {"total": []}
        valid_loss = {"total": []}
        best_models = []
        if not isdir(loss_name):
            if not args.silent:
                print("Creating a new dir at", loss_name)
            os.makedirs(loss_name, exist_ok=True)

    if start_epoch > args.epochs:
        print(f"Fold {fold_idx} already completed ({args.epochs} epochs). Skipping.")
        return -1.0, args.epochs

    os.makedirs(os.path.dirname(out_recording_file) or ".", exist_ok=True)
    out_recording = open(out_recording_file, "a+")
    out_recording.write(loss_name + "\n")
    out_recording.write(f"PPLM: PPLM_PPI | mode:{canonical_mode}\n")
    out_recording.write(f"random.seed: {seed}\n")
    out_recording.write(f"trainset_path: {trainset_path}\n")
    out_recording.write(f"validset_path: {validset_path}\n")
    out_recording.write(
        f"Number of Network Parameters | Total:{total_num} Trainable: {trainable_num}\n"
    )
    out_recording.write(f"learning rate: {args.learning_rate}\n")
    out_recording.write(f"batch_size: {args.batch_size}\n")

    if not args.silent:
        print("Loading dataset...")
    print("train_pre_load_path:", trainset_path)
    print("valid_pre_load_path:", validset_path)

    train_samples = Dataset_PPI_Interact(
        trainset_path,
        esm_dir=args.feat_dir,
        esm_flag=False,
        order_mode="random",
        debug=args.debug,
        n_debug_target=args.n_debug_train,
        validation=False,
    )
    train_dataloader = make_dataloader(
        train_samples, batch_size=args.batch_size, shuffle=True, num_workers=args.num_workers
    )

    valid_samples = Dataset_PPI_Interact(
        validset_path,
        esm_dir=args.feat_dir,
        esm_flag=False,
        order_mode="random",
        debug=args.debug,
        n_debug_target=args.n_debug_valid,
        validation=True,
    )
    valid_dataloader = make_dataloader(
        valid_samples, batch_size=args.batch_size, shuffle=False, num_workers=args.num_workers
    )

    optimizer = torch.optim.AdamW(
        model.parameters(), lr=args.learning_rate, weight_decay=0.01
    )

    embed_a_key = f"{short_mode}_embed_A"
    embed_b_key = f"{short_mode}_embed_B"
    attn_aa_key = f"{short_mode}_attn_AA"
    attn_bb_key = f"{short_mode}_attn_BB"
    inter_key = f"{short_mode}_inter_attn"

    best_auprc = -1.0
    best_epoch = -1

    for epoch in range(start_epoch, args.epochs + 1):
        begin_time_0 = time()

        # ==================== Training ====================
        model.train(True)
        temp_loss = {"total": []}

        for data in train_dataloader:
            optimizer.zero_grad()

            labels = data["label"].to(device).view(-1).float()
            embed_a = data[embed_a_key].to(device)
            embed_b = data[embed_b_key].to(device)
            attn_aa = data[attn_aa_key].to(device)
            attn_bb = data[attn_bb_key].to(device)
            inter_attn = data[inter_key].to(device)

            predictions = model(inter_attn, attn_aa, attn_bb, embed_a, embed_b).view(-1)
            loss = torch.nn.functional.binary_cross_entropy(predictions, labels)

            loss.backward()
            optimizer.step()
            temp_loss["total"].append(loss.item())

        train_loss["total"].append(float(np.mean(temp_loss["total"])))

        # ==================== Validation ====================
        temp_loss = {"total": []}
        model.eval()
        with torch.no_grad():
            total_preds_ave = torch.tensor([], device=device)
            total_labels = torch.tensor([], device=device)

            for data in valid_dataloader:
                labels = data["label"].to(device).view(-1).float()
                embed_a = data[embed_a_key].to(device)
                embed_b = data[embed_b_key].to(device)
                attn_aa = data[attn_aa_key].to(device)
                attn_bb = data[attn_bb_key].to(device)
                inter_attn = data[inter_key].to(device)

                predictions = model(inter_attn, attn_aa, attn_bb, embed_a, embed_b).view(-1)
                loss = torch.nn.functional.binary_cross_entropy(predictions, labels)
                temp_loss["total"].append(loss.item())

                total_preds_ave = torch.cat((total_preds_ave, predictions), 0)
                total_labels = torch.cat((total_labels, labels), 0)

            valid_loss["total"].append(float(np.mean(temp_loss["total"])))

            total_preds_np = total_preds_ave.cpu().numpy().flatten()
            total_labels_np = total_labels.cpu().numpy().flatten()

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
            ) = evaluate(total_labels_np, total_preds_np)

        if AUC_PR > best_auprc:
            best_auprc = AUC_PR
            best_epoch = epoch

        ckpt_payload = {
            "epoch": epoch,
            "model_state_dict": model.state_dict(),
            "train_loss": train_loss,
            "valid_loss": valid_loss,
            "best_models": best_models,
            "AUC_PR": AUC_PR,
            "AUC_ROC": AUC_ROC,
        }
        torch.save(ckpt_payload, join(loss_name, f"model_{epoch}.pkl"))
        torch.save(ckpt_payload, join(loss_name, "model.pkl"))

        elapsed = time() - begin_time_0
        lr_val = optimizer.param_groups[0]["lr"]
        loss_record = (
            f"========== Loss Epoch: {epoch}/{args.epochs} lr: {lr_val} time: {elapsed} | "
            f"Train_loss: {train_loss['total'][-1]} Valid_loss: {valid_loss['total'][-1]}"
        )
        accu_record = (
            f"========== Accuarcy Epoch: {epoch}/{args.epochs} | "
            f"{precision} {recall} {accuracy} {F1_score} "
            f"{Top10} {Top20} {Top50} {AUC_ROC} {AUC_PR}"
        )

        out_recording.write(loss_record + "\n")
        out_recording.write(accu_record + "\n")
        out_recording.flush()

        print(loss_record)
        print(accu_record)

    out_recording.close()
    return best_auprc, best_epoch


def main():
    parser = argparse.ArgumentParser(
        description="Train PPLM-PPI MLP head on plant C1 data (aligned with original PPLM protocol)"
    )
    parser.add_argument(
        "folder",
        nargs="?",
        default=None,
        help="Folder prefix for checkpoints/recordings (e.g. models/DeepAraPPI/PPLM-PPI.inter_intra_attn_embed.mean)",
    )
    parser.add_argument(
        "--output_dir",
        default=None,
        help="Base output directory if positional folder prefix is omitted",
    )
    parser.add_argument(
        "--feat_dir",
        required=True,
        help="Directory containing pre-extracted .pkl feature files",
    )
    parser.add_argument(
        "--split_dir",
        default=None,
        help="Directory for pre-generated .trainset / .testset fold files (default: <output_dir>/splits)",
    )
    parser.add_argument(
        "--pair_list",
        default=None,
        help="Optional C1 TSV file for fast stratified split file generation",
    )
    parser.add_argument(
        "--dataset_prefix",
        default="plant",
        help="Filename prefix for generated .trainset / .testset files (default: plant)",
    )
    parser.add_argument(
        "--cross_round",
        "-cross_round",
        type=int,
        default=0,
        help="Cross-validation fold index (0..n_folds-1), or -1 to run all folds",
    )
    parser.add_argument(
        "--all_folds",
        action="store_true",
        default=False,
        help="Run all CV folds (0..n_folds-1) sequentially in one process",
    )
    parser.add_argument(
        "--n_folds",
        type=int,
        default=10,
        help="Number of stratified CV folds (default: 10)",
    )
    parser.add_argument(
        "--seed_round",
        "-seed_round",
        type=int,
        default=1,
        help="Index into author's SEEDS table (default: 1 -> 73591024, matching run_training.sh)",
    )
    parser.add_argument(
        "--seed_add",
        "-seed",
        type=int,
        default=0,
        help="Integer offset added to SEEDS[seed_round] (default: 0)",
    )
    parser.add_argument(
        "--linear_all",
        "-linear_all",
        action="store_true",
        default=False,
        help="Kept for CLI parity with original training script",
    )
    parser.add_argument(
        "--mode",
        "-mode",
        type=str,
        default="PPI_inter_intra_attn_embed_mean",
        choices=list(MODE_ALIASES.keys()),
        help="Pooling mode (e.g. PPI_inter_intra_attn_embed_mean or mean/max/min)",
    )
    parser.add_argument(
        "--batch_size",
        "-batch_size",
        type=int,
        default=512,
        help="Batch size (default: 512)",
    )
    parser.add_argument(
        "--learning_rate",
        "-lr",
        type=float,
        default=0.00005,
        help="Learning rate (default: 5e-5)",
    )
    parser.add_argument(
        "--gpu_id",
        "-gpu",
        type=int,
        default=0,
        help="GPU device ID (default: 0)",
    )
    parser.add_argument(
        "--epochs",
        "-e",
        type=int,
        default=15,
        help="Training epochs per fold (default: 15)",
    )
    parser.add_argument(
        "--num_workers",
        type=int,
        default=8,
        help="DataLoader worker count (default: 8)",
    )
    parser.add_argument(
        "--debug",
        "-debug",
        action="store_true",
        default=False,
        help="Debug mode: limit sample counts (default: False)",
    )
    parser.add_argument(
        "--n_debug_train",
        type=int,
        default=5000,
        help="Max training samples in debug mode (default: 5000)",
    )
    parser.add_argument(
        "--n_debug_valid",
        type=int,
        default=1000,
        help="Max validation samples in debug mode (default: 1000)",
    )
    parser.add_argument(
        "--silent",
        "-s",
        action="store_true",
        default=False,
        help="Run in silent mode (default: False)",
    )
    args = parser.parse_args()

    short_mode = MODE_ALIASES[args.mode]
    canonical_mode = f"PPI_inter_intra_attn_embed_{short_mode}"

    if args.folder:
        folder_prefix = args.folder
        base_dir = os.path.dirname(folder_prefix) or "."
    elif args.output_dir:
        base_dir = args.output_dir
        folder_prefix = join(base_dir, f"PPLM-PPI.inter_intra_attn_embed.{short_mode}")
    else:
        parser.error("Must provide either positional `folder` prefix or `--output_dir`.")

    if not (0 <= args.seed_round < len(SEEDS)):
        parser.error(f"--seed_round must be between 0 and {len(SEEDS) - 1}")

    seed = SEEDS[args.seed_round] + args.seed_add
    seed_everything(seed)

    print(f"PPLM: PPLM_PPI | mode: {canonical_mode} | seed_round={args.seed_round} | seed={seed}")

    device = torch.device(f"cuda:{args.gpu_id}" if torch.cuda.is_available() else "cpu")
    print("   device:", device)

    split_dir = args.split_dir or join(base_dir, "splits")
    split_paths = generate_cv_split_files(
        feat_dir=args.feat_dir,
        split_dir=split_dir,
        dataset_prefix=args.dataset_prefix,
        pair_list_path=args.pair_list,
        n_folds=args.n_folds,
        seed=seed,
    )

    if args.all_folds or args.cross_round < 0:
        folds_to_run = list(range(args.n_folds))
    else:
        if not (0 <= args.cross_round < args.n_folds):
            parser.error(f"--cross_round must be between 0 and {args.n_folds - 1} (or -1)")
        folds_to_run = [args.cross_round]

    fold_results = []
    for fold_idx in folds_to_run:
        trainset_path, validset_path = split_paths[fold_idx]
        best_auprc, best_epoch = train_single_fold(
            fold_idx=fold_idx,
            trainset_path=trainset_path,
            validset_path=validset_path,
            folder_prefix=folder_prefix,
            short_mode=short_mode,
            canonical_mode=canonical_mode,
            seed=seed,
            device=device,
            args=args,
        )
        fold_results.append((fold_idx, best_auprc, best_epoch))

    if len(fold_results) > 1:
        print(f"\n{'=' * 60}")
        print(f"TRAINING COMPLETE | mode={canonical_mode}")
        print(f"{'=' * 60}")
        for fold_idx, auprc, best_ep in fold_results:
            print(f"  Fold {fold_idx}: AUPRC={auprc:.4f} @ epoch {best_ep}")


if __name__ == "__main__":
    main()
