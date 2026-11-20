#!/usr/bin/env python3
"""
Train and Evaluate ESMAraPPI (Zhou et al., Plant Methods 2023) on AraPPLM (C1, C2, C3, C4).

Adheres strictly to the official `MLP.py` implementation from https://github.com/keiwo/ESMAraPPI:
  - Random seed: `torch.manual_seed(0)`
  - Input: ESM-1b (`esm1b_t33_650M_UR50S`) layer-33 mean-pooled 1280-d embeddings
  - Pair combination: Hadamard product `torch.mul(input1, input2)`
  - MLP architecture:
      Linear(1280, 1024) -> F.relu
      Linear(1024, 512)  -> F.relu
      Linear(512, 128)   -> F.relu
      Linear(128, 16)    -> F.relu
      Linear(16, 1)      -> nn.Sigmoid()
  - Training regime:
      DataLoader(batch_size=128, shuffle=True, num_workers=0)
      Adam(lr=0.001)
      nn.BCELoss()
      EPOCH = 40
  - Evaluates on C1, C2, C3, and C4 and invokes `scripts/benchmarking/evaluate_pplm.py`.
"""

import argparse
import csv
import os
import pickle
import subprocess
import sys
import time

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import average_precision_score, roc_auc_score
from torch.utils.data import DataLoader, Dataset


DEFAULT_SPLIT_DIR = (
    "data/dataset_construction/04_benchmark_splits/arabidopsis/benchmark_data"
)
DEFAULT_EMB_PKL = "features/baselines/esm1b_t33_650M_UR50S_ara.pkl"
DEFAULT_OUTPUT_DIR = "results/AraPPLM/baselines/ESMAraPPI"
DEFAULT_MODEL_PATH = "models/AraPPLM/baselines/esmarappi_mlp.pt"


class ESMAraPPI_MLP(nn.Module):
    """Exact MLP architecture from keiwo/ESMAraPPI/MLP.py."""

    def __init__(self):
        super().__init__()
        self.l1 = nn.Linear(1280, 1024)
        self.l2 = nn.Linear(1024, 512)
        self.l3 = nn.Linear(512, 128)
        self.l4 = nn.Linear(128, 16)
        self.l5 = nn.Linear(16, 1)
        self.sig = nn.Sigmoid()

    def forward(self, input1, input2):
        x = torch.mul(input1.reshape(-1, 1280), input2.reshape(-1, 1280))
        x = F.relu(self.l1(x))
        x = F.relu(self.l2(x))
        x = F.relu(self.l3(x))
        x = F.relu(self.l4(x))
        x = self.sig(self.l5(x))
        return x


class PairEmbeddingDataset(Dataset):
    """Memory-efficient Dataset indexing pre-loaded 1280-d protein tensors by pair."""

    def __init__(self, pairs, embedding_dict):
        self.pairs = pairs
        self.embedding_dict = embedding_dict

    def __len__(self):
        return len(self.pairs)

    def __getitem__(self, idx):
        p1, p2, label = self.pairs[idx]
        e1 = self.embedding_dict[p1]
        e2 = self.embedding_dict[p2]
        return e1, e2, float(label)


def load_pairs_from_tsv(tsv_path, embedding_dict):
    """Load (prot1, prot2, label) pairs from TSV, skipping header lines."""
    pairs = []
    skipped = 0
    with open(tsv_path, "r", encoding="utf-8") as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            if i == 0 and (
                parts[0] in ("Protein1", "participant_a_id", "protein1")
                or parts[2] == "label"
            ):
                continue
            p1, p2, lbl = parts[0], parts[1], int(parts[2])
            if p1 in embedding_dict and p2 in embedding_dict:
                pairs.append((p1, p2, lbl))
            else:
                skipped += 1
    if skipped > 0:
        print(f"  WARNING: Skipped {skipped} pairs missing from embeddingDict in {tsv_path}")
    return pairs


@torch.no_grad()
def evaluate_split(model, pairs, embedding_dict, device, output_csv, batch_size=2048):
    """Run batched inference on a split and write standardized prediction CSV."""
    model.eval()
    dataset = PairEmbeddingDataset(pairs, embedding_dict)
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=False, num_workers=0)

    preds = []
    labels = []
    for e1, e2, y in loader:
        e1 = e1.to(device)
        e2 = e2.to(device)
        out = model(e1, e2).view(-1).detach().cpu().numpy().tolist()
        preds.extend(out)
        labels.extend(y.numpy().astype(int).tolist())

    os.makedirs(os.path.dirname(os.path.abspath(output_csv)), exist_ok=True)
    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Protein1", "Protein2", "pred_score", "true_label"])
        for (p1, p2, _), score, lbl in zip(pairs, preds, labels):
            writer.writerow([p1, p2, f"{score:.6f}", int(lbl)])

    auprc = average_precision_score(labels, preds)
    auroc = roc_auc_score(labels, preds)
    return auprc, auroc


def main():
    parser = argparse.ArgumentParser(
        description="Train and evaluate ESMAraPPI (MLP.py) on AraPPLM (C1, C2, C3, C4)"
    )
    parser.add_argument(
        "--c1",
        default=os.path.join(
            DEFAULT_SPLIT_DIR, "ara_benchmark_c1_train_id0_50_cov0_80.tsv"
        ),
        help="Path to C1 train TSV",
    )
    parser.add_argument(
        "--c2",
        default=os.path.join(
            DEFAULT_SPLIT_DIR, "ara_benchmark_c2_test_id0_50_cov0_80.tsv"
        ),
        help="Path to C2 test TSV",
    )
    parser.add_argument(
        "--c3",
        default=os.path.join(
            DEFAULT_SPLIT_DIR, "ara_benchmark_c3_test_id0_50_cov0_80.tsv"
        ),
        help="Path to C3 test TSV",
    )
    parser.add_argument(
        "--c4",
        default=os.path.join(
            DEFAULT_SPLIT_DIR, "ara_benchmark_c4_test_id0_50_cov0_80.tsv"
        ),
        help="Path to C4 test TSV",
    )
    parser.add_argument(
        "--embeddings",
        default=DEFAULT_EMB_PKL,
        help="Path to ESM-1b 1280-d embedding dictionary (.pkl)",
    )
    parser.add_argument(
        "--output_dir",
        default=DEFAULT_OUTPUT_DIR,
        help="Output directory for prediction CSVs and benchmark_summary.csv",
    )
    parser.add_argument(
        "--model_out",
        default=DEFAULT_MODEL_PATH,
        help="Output path to save trained ESMAraPPI MLP checkpoint",
    )
    parser.add_argument("--epochs", type=int, default=40, help="Epochs (MLP.py: 40)")
    parser.add_argument(
        "--batch_size", type=int, default=128, help="Batch size (MLP.py: 128)"
    )
    parser.add_argument(
        "--lr", type=float, default=0.001, help="Learning rate (MLP.py: 0.001)"
    )
    parser.add_argument("--seed", type=int, default=0, help="Random seed (MLP.py: 0)")
    parser.add_argument("--gpu_id", type=int, default=0, help="GPU device ID")
    args = parser.parse_args()

    # Exact seed from MLP.py: torch.manual_seed(0)
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    device = torch.device(
        f"cuda:{args.gpu_id}" if torch.cuda.is_available() else "cpu"
    )
    print(f"=== ESMAraPPI Training & Evaluation (Device: {device}) ===")

    print(f"Loading ESM-1b embeddings: {args.embeddings}")
    with open(args.embeddings, "rb") as f:
        raw_dict = pickle.load(f)
    embedding_dict = {
        k: (
            v.detach().cpu().float()
            if torch.is_tensor(v)
            else torch.as_tensor(v, dtype=torch.float32)
        )
        for k, v in raw_dict.items()
    }
    print(f"Loaded {len(embedding_dict):,} protein embeddings")

    c1_pairs = load_pairs_from_tsv(args.c1, embedding_dict)
    c2_pairs = load_pairs_from_tsv(args.c2, embedding_dict)
    c3_pairs = load_pairs_from_tsv(args.c3, embedding_dict)
    c4_pairs = load_pairs_from_tsv(args.c4, embedding_dict)
    print(
        f"Pairs loaded | C1: {len(c1_pairs):,} | C2: {len(c2_pairs):,} | "
        f"C3: {len(c3_pairs):,} | C4: {len(c4_pairs):,}"
    )

    train_dataset = PairEmbeddingDataset(c1_pairs, embedding_dict)
    train_loader = DataLoader(
        dataset=train_dataset,
        batch_size=args.batch_size,
        shuffle=True,
        num_workers=0,
    )

    model = ESMAraPPI_MLP().to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)
    loss_func = nn.BCELoss().to(device)

    print(
        f"Training Starts (epochs={args.epochs}, batch_size={args.batch_size}, lr={args.lr})"
    )
    start_total = time.time()
    for epoch in range(args.epochs):
        model.train()
        start_ep = time.time()
        running_loss = 0.0
        for e1, e2, batch_y in train_loader:
            e1 = e1.to(device)
            e2 = e2.to(device)
            batch_y = batch_y.to(device).float()

            output = model(e1, e2)
            loss = loss_func(torch.flatten(output).float(), batch_y)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            running_loss += loss.item() * len(batch_y)

        ep_loss = running_loss / max(1, len(c1_pairs))
        if (epoch + 1) % 5 == 0 or epoch == 0 or (epoch + 1) == args.epochs:
            print(
                f"  Epoch {epoch + 1:02d}/{args.epochs:02d} | "
                f"BCE Loss: {ep_loss:.5f} | "
                f"Time: {time.time() - start_ep:.1f}s (Total: {time.time() - start_total:.0f}s)"
            )

    if args.model_out:
        os.makedirs(os.path.dirname(os.path.abspath(args.model_out)), exist_ok=True)
        torch.save(model.state_dict(), args.model_out)
        print(f"Saved trained ESMAraPPI model to: {args.model_out}")

    os.makedirs(args.output_dir, exist_ok=True)
    c1_csv = os.path.join(args.output_dir, "esmarappi_c1_scores.csv")
    c2_csv = os.path.join(args.output_dir, "esmarappi_c2_scores.csv")
    c3_csv = os.path.join(args.output_dir, "esmarappi_c3_scores.csv")
    c4_csv = os.path.join(args.output_dir, "esmarappi_c4_scores.csv")

    for tag, split_pairs, out_csv in [
        ("C1 (Train)", c1_pairs, c1_csv),
        ("C2 (Test)", c2_pairs, c2_csv),
        ("C3 (Test)", c3_pairs, c3_csv),
        ("C4 (Test)", c4_pairs, c4_csv),
    ]:
        pr, roc = evaluate_split(model, split_pairs, embedding_dict, device, out_csv)
        print(f"  {tag:<12} | AUPRC: {pr:.4f} | AUROC: {roc:.4f} -> {out_csv}")

    eval_script = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "benchmarking",
        "evaluate_pplm.py",
    )
    cmd = [
        sys.executable,
        eval_script,
        "--c1",
        c1_csv,
        "--c2",
        c2_csv,
        "--c3",
        c3_csv,
        "--c4",
        c4_csv,
        "--suite",
        "AraPPLM",
        "--output_dir",
        args.output_dir,
    ]
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
