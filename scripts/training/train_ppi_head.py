#!/usr/bin/env python3
"""
Train the PPLM-PPI MLP head on plant C1 training data using K-fold
stratified cross-validation.

This script replicates the original PPLM-PPI training protocol
(10-fold CV, AdamW, BCE loss, 15 epochs, batch 512) but operates on
plant Arabidopsis C1 data rather than human D-SCRIPT data.

Prerequisites:
  - Pre-extracted PPLM features in --feat_dir (from extract_features.py)
  - GPU with ≥4GB VRAM (MLP head is only ~3.9M params)

Usage:
    # Full training (10 folds × 3 pooling modes)
    for mode in mean max min; do
        python scripts/training/train_ppi_head.py \
            --feat_dir features/DeepAraPPI_C1/ \
            --output_dir models/DeepAraPPI/ \
            --mode $mode --n_folds 10 --epochs 15 --batch_size 512
    done

    # Debug mode (2 folds, 2 epochs, limited samples)
    python scripts/training/train_ppi_head.py \
        --feat_dir features/DeepAraPPI_C1/ \
        --output_dir models/DeepAraPPI/ \
        --mode mean --n_folds 2 --epochs 2 --debug
"""

import sys
import os
import argparse
import random
import numpy as np
import torch
from torch.utils.data import DataLoader
from time import time
from sklearn.model_selection import StratifiedKFold

# Add script directory to path
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

from ppi_model import PPI_inter_intra_attn_embed_single_pooling, evaluate
from dataset import PlantPPIDataset


def seed_everything(seed):
    """Set all random seeds for reproducibility."""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def get_all_samples(feat_dir):
    """
    Discover all available .pkl feature files in the feature directory
    and return (pair_id, label) tuples.
    """
    import pickle

    samples = []
    pkl_files = [f for f in os.listdir(feat_dir) if f.endswith(".pkl")]
    print(f"Found {len(pkl_files)} .pkl files in {feat_dir}")

    for fname in pkl_files:
        pair_id = fname[:-4]  # strip .pkl
        pkl_path = os.path.join(feat_dir, fname)
        try:
            with open(pkl_path, "rb") as f:
                data = pickle.load(f)
            label = int(data["label"])
            samples.append((pair_id, label))
        except Exception as e:
            print(f"WARNING: Could not read {fname}: {e}")

    print(f"Loaded {len(samples)} valid samples")
    return samples


def train_one_fold(
    model, train_loader, valid_loader, device, args, fold_name, recording_file
):
    """
    Train one fold of cross-validation.

    Returns:
        best_auprc: Best validation AUPRC across epochs
        best_epoch: Epoch number of best AUPRC
    """
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=args.learning_rate, weight_decay=0.01
    )

    best_auprc = -1
    best_epoch = -1

    for epoch in range(1, args.epochs + 1):
        epoch_start = time()

        # ==================== Training ====================
        model.train()
        train_losses = []

        for batch_data in train_loader:
            optimizer.zero_grad()

            labels = batch_data["label"].to(device).float()
            inter_attn = batch_data["inter_attn"].to(device)
            attn_AA = batch_data["attn_AA"].to(device)
            attn_BB = batch_data["attn_BB"].to(device)
            embed_A = batch_data["embed_A"].to(device)
            embed_B = batch_data["embed_B"].to(device)

            predictions = model(inter_attn, attn_AA, attn_BB, embed_A, embed_B)
            loss = torch.nn.functional.binary_cross_entropy(
                predictions.squeeze(), labels
            )

            loss.backward()
            optimizer.step()
            train_losses.append(loss.item())

        avg_train_loss = np.mean(train_losses)

        # ==================== Validation ====================
        model.eval()
        valid_losses = []
        all_preds = []
        all_labels = []

        with torch.no_grad():
            for batch_data in valid_loader:
                labels = batch_data["label"].to(device).float()
                inter_attn = batch_data["inter_attn"].to(device)
                attn_AA = batch_data["attn_AA"].to(device)
                attn_BB = batch_data["attn_BB"].to(device)
                embed_A = batch_data["embed_A"].to(device)
                embed_B = batch_data["embed_B"].to(device)

                predictions = model(inter_attn, attn_AA, attn_BB, embed_A, embed_B)
                loss = torch.nn.functional.binary_cross_entropy(
                    predictions.squeeze(), labels
                )
                valid_losses.append(loss.item())

                all_preds.extend(predictions.squeeze().cpu().numpy().tolist())
                all_labels.extend(labels.cpu().numpy().tolist())

        avg_valid_loss = np.mean(valid_losses)

        # Compute metrics
        precision, recall, accuracy, F1, specificity, MCC, \
            Top10, Top20, Top50, AUC_ROC, AUC_PR = evaluate(all_labels, all_preds)

        elapsed = time() - epoch_start

        # Save checkpoint
        checkpoint_path = os.path.join(fold_name, f"model_{epoch}.pkl")
        torch.save(
            {
                "epoch": epoch,
                "model_state_dict": model.state_dict(),
                "train_loss": avg_train_loss,
                "valid_loss": avg_valid_loss,
                "AUC_PR": AUC_PR,
                "AUC_ROC": AUC_ROC,
            },
            checkpoint_path,
        )

        # Also save as "model.pkl" (latest)
        torch.save(
            {
                "epoch": epoch,
                "model_state_dict": model.state_dict(),
                "train_loss": avg_train_loss,
                "valid_loss": avg_valid_loss,
                "AUC_PR": AUC_PR,
                "AUC_ROC": AUC_ROC,
            },
            os.path.join(fold_name, "model.pkl"),
        )

        # Track best
        if AUC_PR > best_auprc:
            best_auprc = AUC_PR
            best_epoch = epoch

        # Logging
        loss_record = (
            f"========== Loss Epoch: {epoch}/{args.epochs} "
            f"lr: {optimizer.param_groups[0]['lr']:.6f} "
            f"time: {elapsed:.1f}s | "
            f"Train_loss: {avg_train_loss:.6f} "
            f"Valid_loss: {avg_valid_loss:.6f}"
        )
        metrics_record = (
            f"========== Accuracy Epoch: {epoch}/{args.epochs} | "
            f"Prec: {precision:.4f} Rec: {recall:.4f} "
            f"Acc: {accuracy:.4f} F1: {F1:.4f} "
            f"Spec: {specificity:.4f} MCC: {MCC:.4f} "
            f"AUROC: {AUC_ROC:.4f} AUPRC: {AUC_PR:.4f}"
        )

        print(loss_record)
        print(metrics_record)

        recording_file.write(loss_record + "\n")
        recording_file.write(metrics_record + "\n")
        recording_file.flush()

    return best_auprc, best_epoch


def main():
    parser = argparse.ArgumentParser(
        description="Train PPLM-PPI MLP head on plant C1 data (K-fold CV)"
    )
    parser.add_argument(
        "--feat_dir",
        required=True,
        help="Directory containing pre-extracted .pkl feature files",
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Output directory for checkpoints and recordings",
    )
    parser.add_argument(
        "--mode",
        type=str,
        default="mean",
        choices=["mean", "max", "min"],
        help="Pooling mode (default: mean)",
    )
    parser.add_argument(
        "--n_folds", type=int, default=10, help="Number of CV folds (default: 10)"
    )
    parser.add_argument(
        "--epochs", type=int, default=15, help="Training epochs per fold (default: 15)"
    )
    parser.add_argument(
        "--batch_size", type=int, default=512, help="Batch size (default: 512)"
    )
    parser.add_argument(
        "--learning_rate",
        "-lr",
        type=float,
        default=5e-5,
        help="Learning rate (default: 5e-5)",
    )
    parser.add_argument(
        "--gpu_id", type=int, default=0, help="GPU device ID (default: 0)"
    )
    parser.add_argument(
        "--seed", type=int, default=26240761, help="Random seed (default: 26240761)"
    )
    parser.add_argument(
        "--num_workers", type=int, default=4, help="DataLoader workers (default: 4)"
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        default=False,
        help="Debug mode: limit samples (default: False)",
    )
    args = parser.parse_args()

    # Seed
    seed_everything(args.seed)
    print(f"Plant-PPLM PPI Training | mode={args.mode} | seed={args.seed}")

    # Device
    device = torch.device(
        f"cuda:{args.gpu_id}" if torch.cuda.is_available() else "cpu"
    )
    print(f"Device: {device}")

    # Load all available samples from feature directory
    all_samples = get_all_samples(args.feat_dir)
    if args.debug:
        all_samples = all_samples[:5000]
        print(f"DEBUG: Limited to {len(all_samples)} samples")

    pair_ids = [s[0] for s in all_samples]
    labels = np.array([s[1] for s in all_samples])

    print(
        f"Dataset: {len(all_samples)} pairs | "
        f"Positives: {np.sum(labels)} ({100 * np.mean(labels):.1f}%) | "
        f"Negatives: {len(labels) - np.sum(labels)} ({100 * (1 - np.mean(labels)):.1f}%)"
    )

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Stratified K-fold
    skf = StratifiedKFold(n_splits=args.n_folds, shuffle=True, random_state=args.seed)

    fold_results = []

    for fold_idx, (train_indices, valid_indices) in enumerate(skf.split(pair_ids, labels)):
        print(f"\n{'='*60}")
        print(f"FOLD {fold_idx}/{args.n_folds - 1} | mode={args.mode}")
        print(f"{'='*60}")

        # Split data
        train_pairs = [(pair_ids[i], int(labels[i])) for i in train_indices]
        valid_pairs = [(pair_ids[i], int(labels[i])) for i in valid_indices]
        train_labels_fold = labels[train_indices]

        print(
            f"  Train: {len(train_pairs)} "
            f"(pos={np.sum(train_labels_fold)}, "
            f"neg={len(train_labels_fold) - np.sum(train_labels_fold)})"
        )
        print(f"  Valid: {len(valid_pairs)}")

        # Create datasets
        train_dataset = PlantPPIDataset(
            train_pairs, args.feat_dir, mode=args.mode, shuffle=True, augment=True
        )
        valid_dataset = PlantPPIDataset(
            valid_pairs, args.feat_dir, mode=args.mode, shuffle=False, augment=False
        )

        train_loader = DataLoader(
            train_dataset,
            batch_size=args.batch_size,
            shuffle=True,
            num_workers=args.num_workers,
            prefetch_factor=4,
            pin_memory=True,
        )
        valid_loader = DataLoader(
            valid_dataset,
            batch_size=args.batch_size,
            shuffle=False,
            num_workers=args.num_workers,
            prefetch_factor=4,
            pin_memory=True,
        )

        # Initialize model (fresh for each fold)
        model = PPI_inter_intra_attn_embed_single_pooling()
        model.to(device)

        total_params = sum(p.numel() for p in model.parameters())
        trainable_params = sum(
            p.numel() for p in model.parameters() if p.requires_grad
        )
        print(f"  Parameters: Total={total_params:,} Trainable={trainable_params:,}")

        # Create fold output directory
        fold_name = os.path.join(
            args.output_dir, f"plant_ppi.{args.mode}.cv_{fold_idx}"
        )
        os.makedirs(fold_name, exist_ok=True)

        # Recording file
        recording_path = fold_name + ".recording"
        recording = open(recording_path, "a+")
        recording.write(f"fold: {fold_name}\n")
        recording.write(f"mode: {args.mode}\n")
        recording.write(f"seed: {args.seed}\n")
        recording.write(f"train_samples: {len(train_pairs)}\n")
        recording.write(f"valid_samples: {len(valid_pairs)}\n")
        recording.write(f"learning_rate: {args.learning_rate}\n")
        recording.write(f"batch_size: {args.batch_size}\n")
        recording.write(f"parameters: {trainable_params}\n")

        # Train
        best_auprc, best_epoch = train_one_fold(
            model, train_loader, valid_loader, device, args, fold_name, recording
        )

        recording.write(f"best_epoch: {best_epoch}\n")
        recording.write(f"best_auprc: {best_auprc:.6f}\n")
        recording.close()

        fold_results.append((fold_idx, best_auprc, best_epoch))
        print(f"\n  Fold {fold_idx} best: AUPRC={best_auprc:.4f} @ epoch {best_epoch}")

    # Summary
    print(f"\n{'='*60}")
    print(f"TRAINING COMPLETE | mode={args.mode}")
    print(f"{'='*60}")
    for fold_idx, auprc, epoch in fold_results:
        print(f"  Fold {fold_idx}: AUPRC={auprc:.4f} @ epoch {epoch}")
    avg_auprc = np.mean([r[1] for r in fold_results])
    print(f"  Average AUPRC: {avg_auprc:.4f}")


if __name__ == "__main__":
    main()
