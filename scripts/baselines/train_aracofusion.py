#!/usr/bin/env python3
"""
Train and Evaluate ARACoFusion (Sarkar & Sarkar, bioRxiv 2026) on AraPPLM (C1, C2, C3, C4).

Implements the full architecture and Optuna-tuned training regime from:
  "ARACoFusion: Uncertainty-aware calibrated deep learning for protein-protein
   interaction network prediction in Arabidopsis thaliana" (DOI: 10.64898/2026.05.22.727120)

Architecture (Section 2.3, Equations 1-14 & Figure 2):
  1. Input: ESM-1b (esm1b_t33_650M_UR50S) layer-33 mean-pooled vectors p1, p2 in R^{1280}.
  2. Reciprocal Cross-Attention Encoder (Section 2.3.2, Eqs. 3-6):
     - Outer projections: Q1 = W_q1 p1, (K1, V1) = chunk(W_kv1 p1);
                          Q2 = W_q2 p2, (K2, V2) = chunk(W_kv2 p2).
     - Multi-head cross-attention (H = 8 heads) over reshaped sub-vector tokens [B, 8, 160]:
       c1 = MultiHeadAttention(Q1, K2, V2) in R^{1280}
       c2 = MultiHeadAttention(Q2, K1, V1) in R^{1280}
  3. Feature Adaptation & Latent Interaction Projector (Section 2.3.3, Eqs. 7-10 & Fig. 2):
     - P_prod = p1 * p2 in R^{1280}
     - P_diff = |p1 - p2| in R^{1280}
     - a_i = W2 * Dropout(GELU(W1 p_i + b1)) + b2 in R^{512}, i in {1, 2}
  4. Multi-Stage Fusion Classifier Head (Section 2.3.4, Eqs. 11-14):
     - Z = [p1; p2; c1; c2; P_prod; P_diff; a1; a2] in R^{8704}
     - h1 = Dropout(GELU(Linear(BatchNorm1d(Z))))
     - h2 = Dropout(GELU(Linear(h1)))
     - u  = Linear(h2), y_hat = sigmoid(u / T)

Training Regime (Section 2.4 & Table 2 Optuna Optimal Parameters):
  - Attention heads H = 8
  - AdamW optimizer, lr = 0.000681383
  - Focal Loss gamma = 3.892573, Label smoothing = 0.0514500
  - Dropout rate = 0.118479
  - MC-Dropout uncertainty variance weight lambda = 0.579057 across M = 3 forward passes (Eq. 15)
  - Post-hoc Temperature Scaling T = 0.58 (Eq. 16 & Section 3.2.3)
  - Evaluates on C1, C2, C3, and C4 and invokes `scripts/benchmarking/evaluate_pplm.py`.
"""

import argparse
import csv
import math
import os
import pickle
import random
import subprocess
import sys
import time

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedShuffleSplit
from torch.utils.data import DataLoader, Dataset


DEFAULT_SPLIT_DIR = (
    "data/dataset_construction/04_benchmark_splits/arabidopsis/benchmark_data"
)
DEFAULT_EMB_PKL = "features/baselines/esm1b_t33_650M_UR50S_ara.pkl"
DEFAULT_OUTPUT_DIR = "results/AraPPLM/baselines/ARACoFusion"
DEFAULT_MODEL_PATH = "models/AraPPLM/baselines/aracofusion_best.pt"


class SubVectorMultiHeadCrossAttention(nn.Module):
    """
    Multi-Head Cross-Attention over 1280-d protein vectors partitioned into H sub-vector
    tokens of dimension d_head = d / H (Equations 3-6), ensuring non-trivial H x H
    attention weights and full gradient flow through Q, K, and V.
    """

    def __init__(self, embed_dim=1280, num_heads=8, dropout=0.118479):
        super().__init__()
        if embed_dim % num_heads != 0:
            raise ValueError(
                f"embed_dim ({embed_dim}) must be divisible by num_heads ({num_heads})"
            )
        self.embed_dim = embed_dim
        self.num_heads = num_heads
        self.head_dim = embed_dim // num_heads
        self.scale = math.sqrt(self.head_dim)

        self.w_q = nn.Linear(self.head_dim, self.head_dim)
        self.w_k = nn.Linear(self.head_dim, self.head_dim)
        self.w_v = nn.Linear(self.head_dim, self.head_dim)
        self.w_o = nn.Linear(embed_dim, embed_dim)
        self.attn_drop = nn.Dropout(dropout)

    def forward(self, q_vec, k_vec, v_vec):
        bsz = q_vec.size(0)
        # Reshape [B, 1280] -> [B, H, d_head]
        q = self.w_q(q_vec.view(bsz, self.num_heads, self.head_dim))
        k = self.w_k(k_vec.view(bsz, self.num_heads, self.head_dim))
        v = self.w_v(v_vec.view(bsz, self.num_heads, self.head_dim))

        # Scaled dot-product cross-attention across H sub-vector tokens: [B, H, H]
        scores = torch.bmm(q, k.transpose(1, 2)) / self.scale
        attn = self.attn_drop(F.softmax(scores, dim=-1))
        context = torch.bmm(attn, v).reshape(bsz, self.embed_dim)
        return self.w_o(context)


class ARACoFusion(nn.Module):
    """Full ARACoFusion architecture (Section 2.3, Equations 2-14 & Figure 2)."""

    def __init__(
        self,
        embed_dim=1280,
        num_heads=8,
        latent_hidden_dim=1024,
        latent_out_dim=512,
        clf_h1_dim=1024,
        clf_h2_dim=256,
        dropout=0.118479,
    ):
        super().__init__()
        # Section 2.3.2: Outer Query and Key-Value transformations for p1 and p2
        self.w_q1 = nn.Linear(embed_dim, embed_dim)
        self.w_kv1 = nn.Linear(embed_dim, 2 * embed_dim)
        self.w_q2 = nn.Linear(embed_dim, embed_dim)
        self.w_kv2 = nn.Linear(embed_dim, 2 * embed_dim)

        self.cross_attn_1to2 = SubVectorMultiHeadCrossAttention(
            embed_dim=embed_dim, num_heads=num_heads, dropout=dropout
        )
        self.cross_attn_2to1 = SubVectorMultiHeadCrossAttention(
            embed_dim=embed_dim, num_heads=num_heads, dropout=dropout
        )

        # Section 2.3.3 & Fig. 2: Feature adaptation & Latent interaction projector
        self.latent_proj = nn.Sequential(
            nn.Linear(embed_dim, latent_hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(latent_hidden_dim, latent_out_dim),
        )

        # Section 2.3.4: Z = [p1; p2; c1; c2; P_prod; P_diff; a1; a2] -> 6*1280 + 2*512 = 8704
        fusion_dim = 6 * embed_dim + 2 * latent_out_dim
        self.bn = nn.BatchNorm1d(fusion_dim)
        self.fc1 = nn.Linear(fusion_dim, clf_h1_dim)
        self.drop1 = nn.Dropout(dropout)
        self.fc2 = nn.Linear(clf_h1_dim, clf_h2_dim)
        self.drop2 = nn.Dropout(dropout)
        self.out = nn.Linear(clf_h2_dim, 1)

    def forward_logits(self, p1, p2):
        q1 = self.w_q1(p1)
        k1, v1 = torch.chunk(self.w_kv1(p1), 2, dim=-1)

        q2 = self.w_q2(p2)
        k2, v2 = torch.chunk(self.w_kv2(p2), 2, dim=-1)

        c1 = self.cross_attn_1to2(q1, k2, v2)
        c2 = self.cross_attn_2to1(q2, k1, v1)

        p_prod = p1 * p2
        p_diff = torch.abs(p1 - p2)

        a1 = self.latent_proj(p1)
        a2 = self.latent_proj(p2)

        z = torch.cat([p1, p2, c1, c2, p_prod, p_diff, a1, a2], dim=-1)
        h1 = self.drop1(F.gelu(self.fc1(self.bn(z))))
        h2 = self.drop2(F.gelu(self.fc2(h1)))
        logits = self.out(h2).view(-1)
        return logits

    def forward(self, p1, p2, temperature=1.0):
        logits = self.forward_logits(p1, p2)
        return torch.sigmoid(logits / temperature)


class FocalLossWithLabelSmoothing(nn.Module):
    """Focal Loss with Label Smoothing (Section 2.4.1, Table 2)."""

    def __init__(self, gamma=3.892573, smoothing=0.0514500, eps=1e-7):
        super().__init__()
        self.gamma = gamma
        self.smoothing = smoothing
        self.eps = eps

    def forward(self, probs, targets):
        y_smooth = targets * (1.0 - self.smoothing) + 0.5 * self.smoothing
        p = torch.clamp(probs, self.eps, 1.0 - self.eps)
        loss_pos = -y_smooth * torch.pow(1.0 - p, self.gamma) * torch.log(p)
        loss_neg = -(1.0 - y_smooth) * torch.pow(p, self.gamma) * torch.log(1.0 - p)
        return torch.mean(loss_pos + loss_neg)


class ARACoFusionDataset(Dataset):
    """Dataset supporting 50% random A<->B pair swapping during training."""

    def __init__(self, pairs, embedding_dict, random_swap=False):
        self.pairs = pairs
        self.embedding_dict = embedding_dict
        self.random_swap = random_swap

    def __len__(self):
        return len(self.pairs)

    def __getitem__(self, idx):
        p1, p2, label = self.pairs[idx]
        if self.random_swap and random.random() >= 0.5:
            p1, p2 = p2, p1
        return self.embedding_dict[p1], self.embedding_dict[p2], float(label)


def load_pairs_from_tsv(tsv_path, embedding_dict):
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
def predict_pairs(
    model, pairs, embedding_dict, device, temperature=0.58, batch_size=1024, symmetric=True
):
    model.eval()
    dataset = ARACoFusionDataset(pairs, embedding_dict, random_swap=False)
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=False, num_workers=0)

    preds = []
    labels = []
    for e1, e2, y in loader:
        e1 = e1.to(device)
        e2 = e2.to(device)
        p_ab = model(e1, e2, temperature=temperature)
        if symmetric:
            p_ba = model(e2, e1, temperature=temperature)
            p = 0.5 * (p_ab + p_ba)
        else:
            p = p_ab
        preds.extend(p.detach().cpu().numpy().tolist())
        labels.extend(y.numpy().astype(int).tolist())
    return np.asarray(preds, dtype=np.float64), np.asarray(labels, dtype=np.int64)


def save_and_score_split(
    model, pairs, embedding_dict, device, output_csv, temperature=0.58
):
    preds, labels = predict_pairs(
        model, pairs, embedding_dict, device, temperature=temperature, symmetric=True
    )
    os.makedirs(os.path.dirname(os.path.abspath(output_csv)), exist_ok=True)
    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Protein1", "Protein2", "pred_score", "true_label"])
        for (p1, p2, _), score, lbl in zip(pairs, preds, labels):
            writer.writerow([p1, p2, f"{float(score):.6f}", int(lbl)])

    auprc = average_precision_score(labels, preds)
    auroc = roc_auc_score(labels, preds)
    return auprc, auroc


def main():
    parser = argparse.ArgumentParser(
        description="Train and evaluate ARACoFusion on AraPPLM (C1, C2, C3, C4)"
    )
    parser.add_argument(
        "--c1",
        default=os.path.join(
            DEFAULT_SPLIT_DIR, "ara_benchmark_c1_train_id0_50_cov0_80.tsv"
        ),
    )
    parser.add_argument(
        "--c2",
        default=os.path.join(
            DEFAULT_SPLIT_DIR, "ara_benchmark_c2_test_id0_50_cov0_80.tsv"
        ),
    )
    parser.add_argument(
        "--c3",
        default=os.path.join(
            DEFAULT_SPLIT_DIR, "ara_benchmark_c3_test_id0_50_cov0_80.tsv"
        ),
    )
    parser.add_argument(
        "--c4",
        default=os.path.join(
            DEFAULT_SPLIT_DIR, "ara_benchmark_c4_test_id0_50_cov0_80.tsv"
        ),
    )
    parser.add_argument("--embeddings", default=DEFAULT_EMB_PKL)
    parser.add_argument("--output_dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--model_out", default=DEFAULT_MODEL_PATH)
    # Optuna-optimal parameters from Table 2 of ARACoFusion.pdf
    parser.add_argument("--num_heads", type=int, default=8)
    parser.add_argument("--lr", type=float, default=0.000681383)
    parser.add_argument("--gamma", type=float, default=3.892573)
    parser.add_argument("--dropout", type=float, default=0.118479)
    parser.add_argument("--smoothing", type=float, default=0.0514500)
    parser.add_argument("--uncertainty_weight", type=float, default=0.579057)
    parser.add_argument("--mc_passes", type=int, default=3)
    parser.add_argument("--temperature", type=float, default=0.58)
    parser.add_argument("--epochs", type=int, default=40)
    parser.add_argument("--batch_size", type=int, default=128)
    parser.add_argument("--patience", type=int, default=8)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--gpu_id", type=int, default=0)
    args = parser.parse_args()

    random.seed(args.seed)
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    device = torch.device(
        f"cuda:{args.gpu_id}" if torch.cuda.is_available() else "cpu"
    )
    print(f"=== ARACoFusion Training & Evaluation (Device: {device}) ===")

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
    print(f"Loaded {len(embedding_dict):,} ESM-1b protein embeddings")

    c1_pairs = load_pairs_from_tsv(args.c1, embedding_dict)
    c2_pairs = load_pairs_from_tsv(args.c2, embedding_dict)
    c3_pairs = load_pairs_from_tsv(args.c3, embedding_dict)
    c4_pairs = load_pairs_from_tsv(args.c4, embedding_dict)

    # Stratified 90/10 split of C1 for early stopping / checkpoint selection
    c1_labels = np.array([lbl for _, _, lbl in c1_pairs])
    sss = StratifiedShuffleSplit(n_splits=1, test_size=0.10, random_state=args.seed)
    train_idx, val_idx = next(sss.split(np.zeros(len(c1_labels)), c1_labels))
    train_sub = [c1_pairs[i] for i in train_idx]
    val_sub = [c1_pairs[i] for i in val_idx]
    print(
        f"C1 Split: {len(train_sub):,} train / {len(val_sub):,} val | "
        f"Test Splits: C2={len(c2_pairs):,}, C3={len(c3_pairs):,}, C4={len(c4_pairs):,}"
    )

    train_loader = DataLoader(
        ARACoFusionDataset(train_sub, embedding_dict, random_swap=True),
        batch_size=args.batch_size,
        shuffle=True,
        num_workers=0,
        drop_last=len(train_sub) > args.batch_size,
    )

    model = ARACoFusion(
        embed_dim=1280,
        num_heads=args.num_heads,
        dropout=args.dropout,
    ).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr)
    focal_ls = FocalLossWithLabelSmoothing(
        gamma=args.gamma, smoothing=args.smoothing
    ).to(device)

    best_val_auprc = -1.0
    best_state = None
    no_improve = 0
    start_total = time.time()

    for epoch in range(args.epochs):
        model.train()
        start_ep = time.time()
        running_loss = 0.0

        for e1, e2, batch_y in train_loader:
            e1 = e1.to(device)
            e2 = e2.to(device)
            batch_y = batch_y.to(device).float()

            # Equation 15: M = 3 Monte Carlo forward passes + variance regularization
            mc_preds = []
            focal_sum = 0.0
            for _ in range(args.mc_passes):
                y_hat = model(e1, e2, temperature=1.0)
                mc_preds.append(y_hat)
                focal_sum = focal_sum + focal_ls(y_hat, batch_y)

            focal_term = focal_sum / args.mc_passes
            stacked_preds = torch.stack(mc_preds, dim=0)  # [M, B]
            var_term = torch.mean(torch.var(stacked_preds, dim=0, unbiased=False))
            loss = focal_term + args.uncertainty_weight * var_term

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            running_loss += loss.item() * len(batch_y)

        ep_loss = running_loss / max(1, len(train_sub))
        val_preds, val_labels = predict_pairs(
            model,
            val_sub,
            embedding_dict,
            device,
            temperature=args.temperature,
            symmetric=True,
        )
        val_auprc = average_precision_score(val_labels, val_preds)
        val_auroc = roc_auc_score(val_labels, val_preds)

        improved = val_auprc > best_val_auprc
        if improved:
            best_val_auprc = val_auprc
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            no_improve = 0
        else:
            no_improve += 1

        print(
            f"  Epoch {epoch + 1:02d}/{args.epochs:02d} | Loss: {ep_loss:.5f} | "
            f"Val AUPRC: {val_auprc:.4f} (Best: {best_val_auprc:.4f}) | "
            f"Val AUROC: {val_auroc:.4f} | Time: {time.time() - start_ep:.1f}s"
        )

        if no_improve >= args.patience:
            print(
                f"Early stopping at epoch {epoch + 1} (no Val AUPRC improvement for {args.patience} epochs)"
            )
            break

    if best_state is not None:
        model.load_state_dict(best_state)

    if args.model_out:
        os.makedirs(os.path.dirname(os.path.abspath(args.model_out)), exist_ok=True)
        torch.save(model.state_dict(), args.model_out)
        print(f"Saved best ARACoFusion checkpoint to: {args.model_out}")

    os.makedirs(args.output_dir, exist_ok=True)
    c1_csv = os.path.join(args.output_dir, "aracofusion_c1_scores.csv")
    c2_csv = os.path.join(args.output_dir, "aracofusion_c2_scores.csv")
    c3_csv = os.path.join(args.output_dir, "aracofusion_c3_scores.csv")
    c4_csv = os.path.join(args.output_dir, "aracofusion_c4_scores.csv")

    print(
        f"\nEvaluating best ARACoFusion checkpoint (Temperature T={args.temperature}, Total time: {time.time() - start_total:.0f}s):"
    )
    for tag, split_pairs, out_csv in [
        ("C1 (Train)", c1_pairs, c1_csv),
        ("C2 (Test)", c2_pairs, c2_csv),
        ("C3 (Test)", c3_pairs, c3_csv),
        ("C4 (Test)", c4_pairs, c4_csv),
    ]:
        pr, roc = save_and_score_split(
            model,
            split_pairs,
            embedding_dict,
            device,
            out_csv,
            temperature=args.temperature,
        )
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
