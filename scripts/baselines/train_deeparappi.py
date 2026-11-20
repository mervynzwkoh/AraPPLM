#!/usr/bin/env python3
"""
Train and Evaluate DeepAraPPI (Zheng et al., The Plant Journal 2023) on AraPPLM (C1, C2, C3, C4).

PyTorch 1:1 port of the official Keras scripts from https://github.com/zjy1125/DeepAraPPI:
  1. `rcnn.py` (Siamese RCNN on 32-d word2vec amino acid embeddings, L = 2000):
     - 5 stacked blocks of:
         Conv1D(in_ch -> 50, kernel_size=3) -> MaxPool1D(3) -> s
         Bidirectional(GRU(50 -> 50)) -> r(s) (100 channels)
         s = concatenate([r(s), s]) (150 channels)
     - Final Conv1D(150 -> 50, kernel_size=3) -> GlobalAveragePooling1D() -> [B, 50]
     - Pair fusion: multiply([s1, s2]) -> [B, 50]
     - Classifier: Dense(50 -> 100) -> LeakyReLU(0.3) -> Dense(100 -> 28) -> LeakyReLU(0.3) -> Dense(28 -> 2, Softmax)
     - Optimizer: Adam(lr=0.001, betas=(0.9, 0.999), eps=1e-8, amsgrad=True), categorical_crossentropy, 50 epochs.
  2. `Domain2vec.py` (Domain2vec 128-d MLP):
     - Pair fusion: concatenate([s1, s2]) -> [B, 256]
     - Classifier: Dense(256 -> 64) -> LeakyReLU(0.1) -> Dense(64 -> 16) -> LeakyReLU(0.1) -> Dense(16 -> 2, Softmax)
     - Optimizer: Adam(lr=0.01, betas=(0.9, 0.999), eps=1e-8, amsgrad=True), categorical_crossentropy, 50 epochs.
  3. `go2vec.py` (GO2vec 128-d MLP):
     - Identical architecture and optimizer to `Domain2vec.py` on 128-d GO2vec embeddings.
  4. Integrated L2 Logistic Regression Meta-Classifier:
     - All 3 sub-models share `KFold(n_splits=5, random_state=63, shuffle=True)` on C1.
     - Concatenated 5-fold out-of-fold validation predictions [S_RCNN, S_Domain2vec, S_GO2vec] on C1
       train an L2-penalized LogisticRegression tuned via 5-fold GridSearchCV.
     - Evaluates all 4 models (`DeepAraPPI_RCNN`, `DeepAraPPI_Domain2vec`, `DeepAraPPI_GO2vec`,
       `DeepAraPPI_Integrated`) across C1, C2, C3, and C4 via `scripts/benchmarking/evaluate_pplm.py`.
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
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GridSearchCV, KFold, cross_val_predict
from torch.utils.data import DataLoader, Dataset


DEFAULT_SPLIT_DIR = (
    "data/dataset_construction/04_benchmark_splits/arabidopsis/benchmark_data"
)
DEFAULT_ASSETS_PKL = "features/baselines/deeparappi_assets.pkl"
DEFAULT_OUTPUT_ROOT = "results/AraPPLM/baselines"
DEFAULT_MODEL_DIR = "models/AraPPLM/baselines/DeepAraPPI"


class DeepAraPPI_VectorMLP(nn.Module):
    """
    Exact PyTorch port of `Domain2vec.py::build_model` and `go2vec.py::build_model`:
      concatenate([s1, s2]) (256-d)
      -> Dense(64) -> LeakyReLU(0.1)
      -> Dense(16) -> LeakyReLU(0.1)
      -> Dense(2, softmax)
    """

    def __init__(self, in_dim=128):
        super().__init__()
        self.fc1 = nn.Linear(2 * in_dim, 64)
        self.act1 = nn.LeakyReLU(negative_slope=0.1)
        self.fc2 = nn.Linear(64, 16)
        self.act2 = nn.LeakyReLU(negative_slope=0.1)
        self.out = nn.Linear(16, 2)

    def forward(self, s1, s2):
        x = torch.cat([s1, s2], dim=-1)
        x = self.act1(self.fc1(x))
        x = self.act2(self.fc2(x))
        return self.out(x)


class DeepAraPPI_SiameseRCNN(nn.Module):
    """
    Exact PyTorch port of `rcnn.py::build_model` (seq_size=2000, dim=32, hidden_dim=50):
      5 stacked [Conv1D(hidden_dim, 3) -> MaxPool1D(3) -> BiGRU(50) + residual concat] blocks
      + final Conv1D(hidden_dim, 3) -> GlobalAveragePooling1D()
      + element-wise multiply([s1, s2])
      + Dense(50 -> 100) -> LeakyReLU(0.3) -> Dense(100 -> 28) -> LeakyReLU(0.3) -> Dense(28 -> 2)
    """

    def __init__(self, dim=32, hidden_dim=50):
        super().__init__()
        self.l1 = nn.Conv1d(dim, hidden_dim, kernel_size=3)
        self.r1 = nn.GRU(
            hidden_dim, hidden_dim, batch_first=True, bidirectional=True
        )

        self.l2 = nn.Conv1d(3 * hidden_dim, hidden_dim, kernel_size=3)
        self.r2 = nn.GRU(
            hidden_dim, hidden_dim, batch_first=True, bidirectional=True
        )

        self.l3 = nn.Conv1d(3 * hidden_dim, hidden_dim, kernel_size=3)
        self.r3 = nn.GRU(
            hidden_dim, hidden_dim, batch_first=True, bidirectional=True
        )

        self.l4 = nn.Conv1d(3 * hidden_dim, hidden_dim, kernel_size=3)
        self.r4 = nn.GRU(
            hidden_dim, hidden_dim, batch_first=True, bidirectional=True
        )

        self.l5 = nn.Conv1d(3 * hidden_dim, hidden_dim, kernel_size=3)
        self.r5 = nn.GRU(
            hidden_dim, hidden_dim, batch_first=True, bidirectional=True
        )

        self.l6 = nn.Conv1d(3 * hidden_dim, hidden_dim, kernel_size=3)

        mid_dim = int((hidden_dim + 7) / 2)  # 28 when hidden_dim=50
        self.fc1 = nn.Linear(hidden_dim, 100)
        self.act1 = nn.LeakyReLU(negative_slope=0.3)
        self.fc2 = nn.Linear(100, mid_dim)
        self.act2 = nn.LeakyReLU(negative_slope=0.3)
        self.out = nn.Linear(mid_dim, 2)

    def _conv_gru_block(self, x, conv, gru):
        # x: [B, C_in, L] -> conv -> maxpool(3) -> [B, 50, L_out]
        s = F.max_pool1d(conv(x), kernel_size=3, stride=3)
        # GRU expects [B, L_out, 50] -> outputs [B, L_out, 100]
        r_out, _ = gru(s.transpose(1, 2))
        # Concatenate [r(s), s] along channel dim -> [B, 150, L_out]
        return torch.cat([r_out.transpose(1, 2), s], dim=1)

    def encode_sequence(self, seq_x):
        # seq_x: [B, 32, 2000]
        s = self._conv_gru_block(seq_x, self.l1, self.r1)
        s = self._conv_gru_block(s, self.l2, self.r2)
        s = self._conv_gru_block(s, self.l3, self.r3)
        s = self._conv_gru_block(s, self.l4, self.r4)
        s = self._conv_gru_block(s, self.l5, self.r5)
        s = self.l6(s)  # [B, 50, 5]
        return s.mean(dim=-1)  # GlobalAveragePooling1D -> [B, 50]

    def forward(self, seq1, seq2):
        s1 = self.encode_sequence(seq1)
        s2 = self.encode_sequence(seq2)
        merge_text = s1 * s2
        x = self.act1(self.fc1(merge_text))
        x = self.act2(self.fc2(x))
        return self.out(x)


class VectorPairDataset(Dataset):
    def __init__(self, pairs, vec_dict):
        self.pairs = pairs
        self.vec_dict = vec_dict

    def __len__(self):
        return len(self.pairs)

    def __getitem__(self, idx):
        p1, p2, lbl = self.pairs[idx]
        return self.vec_dict[p1], self.vec_dict[p2], int(lbl)


class RCNNPairDataset(Dataset):
    """Indexes pre-encoded [32, 2000] protein sequence tensors by pair."""

    def __init__(self, pairs, seq_tensor_dict):
        self.pairs = pairs
        self.seq_tensor_dict = seq_tensor_dict

    def __len__(self):
        return len(self.pairs)

    def __getitem__(self, idx):
        p1, p2, lbl = self.pairs[idx]
        return self.seq_tensor_dict[p1], self.seq_tensor_dict[p2], int(lbl)


def build_rcnn_sequence_tensors(sequences, aa_vec, seq_size=2000, dim=32):
    """
    Pre-encode all unique proteins into channel-first [32, 2000] float32 tensors
    matching `rcnn.py::get_pro_embedding`.
    For 8,426 unique proteins, this takes ~2.16 GB RAM total.
    """
    zero_aa = np.zeros(dim, dtype=np.float32)
    seq_tensors = {}
    for pro, seq in sequences.items():
        mat = np.zeros((seq_size, dim), dtype=np.float32)
        limit = min(len(seq), seq_size)
        for idx in range(limit):
            mat[idx] = aa_vec.get(seq[idx], zero_aa)
        # Store channel-first [32, 2000] for PyTorch Conv1d
        seq_tensors[pro] = torch.from_numpy(mat.T.copy())
    return seq_tensors


def load_pairs_from_tsv(tsv_path):
    pairs = []
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
            pairs.append((parts[0], parts[1], int(parts[2])))
    return pairs


@torch.no_grad()
def predict_loader(model, loader, device):
    model.eval()
    probs = []
    for x1, x2, _ in loader:
        x1 = x1.to(device, non_blocking=True)
        x2 = x2.to(device, non_blocking=True)
        logits = model(x1, x2)
        p_pos = F.softmax(logits, dim=-1)[:, 1]
        probs.extend(p_pos.detach().cpu().numpy().tolist())
    return np.asarray(probs, dtype=np.float64)


def train_5fold_branch(
    branch_name,
    model_fn,
    c1_pairs,
    test_splits_dict,
    dataset_cls,
    feat_lookup,
    kf_splits,
    device,
    lr,
    batch_size,
    eval_batch_size,
    epochs,
    model_dir,
):
    """
    Run 5-fold CV on C1 matching `rcnn.py`, `Domain2vec.py`, and `go2vec.py`:
      - Optimizer: Adam(lr=lr, betas=(0.9, 0.999), eps=1e-8, amsgrad=True)
      - Loss: CrossEntropyLoss (2-class softmax categorical_crossentropy)
      - Collects out-of-fold C1 predictions (`oof_c1`) and 5-fold averaged predictions
        for C1 (full ensemble), C2, C3, and C4.
    """
    print(f"\n>>> Training DeepAraPPI Branch: {branch_name} (5 folds x {epochs} epochs, lr={lr}, batch_size={batch_size}) <<<")
    n_c1 = len(c1_pairs)
    oof_c1 = np.zeros(n_c1, dtype=np.float64)
    ens_c1 = np.zeros(n_c1, dtype=np.float64)
    test_ens = {
        k: np.zeros(len(v), dtype=np.float64) for k, v in test_splits_dict.items()
    }

    full_c1_loader = DataLoader(
        dataset_cls(c1_pairs, feat_lookup),
        batch_size=eval_batch_size,
        shuffle=False,
        num_workers=0,
    )
    test_loaders = {
        k: DataLoader(
            dataset_cls(v, feat_lookup),
            batch_size=eval_batch_size,
            shuffle=False,
            num_workers=0,
        )
        for k, v in test_splits_dict.items()
    }

    os.makedirs(model_dir, exist_ok=True)
    criterion = nn.CrossEntropyLoss().to(device)

    for fold_idx, (train_idx, val_idx) in enumerate(kf_splits, start=1):
        fold_start = time.time()
        train_sub = [c1_pairs[i] for i in train_idx]
        val_sub = [c1_pairs[i] for i in val_idx]
        val_labels = np.array([c1_pairs[i][2] for i in val_idx])

        train_loader = DataLoader(
            dataset_cls(train_sub, feat_lookup),
            batch_size=batch_size,
            shuffle=True,
            num_workers=0,
        )
        val_loader = DataLoader(
            dataset_cls(val_sub, feat_lookup),
            batch_size=eval_batch_size,
            shuffle=False,
            num_workers=0,
        )

        model = model_fn().to(device)
        optimizer = torch.optim.Adam(
            model.parameters(),
            lr=lr,
            betas=(0.9, 0.999),
            eps=1e-8,
            amsgrad=True,
        )

        for epoch in range(epochs):
            model.train()
            running_loss = 0.0
            for x1, x2, y in train_loader:
                x1 = x1.to(device, non_blocking=True)
                x2 = x2.to(device, non_blocking=True)
                y = y.to(device, non_blocking=True).long()

                logits = model(x1, x2)
                loss = criterion(logits, y)
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()
                running_loss += loss.item() * len(y)

            if (epoch + 1) % 10 == 0 or epoch == 0 or (epoch + 1) == epochs:
                val_probs = predict_loader(model, val_loader, device)
                val_pr = average_precision_score(val_labels, val_probs)
                print(
                    f"  [{branch_name} Fold {fold_idx}/5 | Ep {epoch + 1:02d}/{epochs:02d}] "
                    f"Loss: {running_loss / max(1, len(train_sub)):.4f} | "
                    f"Fold-Val AUPRC: {val_pr:.4f} ({time.time() - fold_start:.0f}s)"
                )

        ckpt_path = os.path.join(model_dir, f"{branch_name.lower()}_model_fold{fold_idx}.pt")
        torch.save(model.state_dict(), ckpt_path)

        val_probs = predict_loader(model, val_loader, device)
        oof_c1[val_idx] = val_probs
        ens_c1 += predict_loader(model, full_c1_loader, device) / len(kf_splits)

        for split_k, t_loader in test_loaders.items():
            test_ens[split_k] += predict_loader(model, t_loader, device) / len(kf_splits)

    c1_labels = np.array([lbl for _, _, lbl in c1_pairs])
    print(
        f"  [{branch_name} 5-Fold Complete] "
        f"C1 OOF AUPRC: {average_precision_score(c1_labels, oof_c1):.4f} | "
        f"C1 Ensemble AUPRC: {average_precision_score(c1_labels, ens_c1):.4f}"
    )
    return oof_c1, ens_c1, test_ens


def write_pred_csv(pairs, scores, out_csv):
    os.makedirs(os.path.dirname(os.path.abspath(out_csv)), exist_ok=True)
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Protein1", "Protein2", "pred_score", "true_label"])
        for (p1, p2, lbl), score in zip(pairs, scores):
            writer.writerow([p1, p2, f"{float(score):.6f}", int(lbl)])


def run_evaluate_pplm(c1_csv, c2_csv, c3_csv, c4_csv, output_dir):
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
        output_dir,
    ]
    subprocess.run(cmd, check=True)


def main():
    parser = argparse.ArgumentParser(
        description="Train and evaluate DeepAraPPI (RCNN, Domain2vec, GO2vec, Integrated LR) on AraPPLM"
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
    parser.add_argument("--assets", default=DEFAULT_ASSETS_PKL)
    parser.add_argument("--output_root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--model_dir", default=DEFAULT_MODEL_DIR)
    parser.add_argument("--epochs", type=int, default=50, help="Epochs per fold (default: 50)")
    parser.add_argument(
        "--mlp_batch_size",
        type=int,
        default=32,
        help="Batch size for Domain2vec and GO2vec (Domain2vec.py/go2vec.py: 32)",
    )
    parser.add_argument(
        "--rcnn_batch_size",
        type=int,
        default=32,
        help="Batch size for RCNN (rcnn.py: 32)",
    )
    parser.add_argument("--n_splits", type=int, default=5, help="CV folds (default: 5)")
    parser.add_argument("--seed", type=int, default=63, help="Random state (DeepAraPPI: 63)")
    parser.add_argument("--gpu_id", type=int, default=0)
    args = parser.parse_args()

    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    device = torch.device(
        f"cuda:{args.gpu_id}" if torch.cuda.is_available() else "cpu"
    )
    print(f"=== DeepAraPPI 5-Fold Training & Evaluation (Device: {device}) ===")

    with open(args.assets, "rb") as f:
        assets = pickle.load(f)

    c1_pairs = load_pairs_from_tsv(args.c1)
    c2_pairs = load_pairs_from_tsv(args.c2)
    c3_pairs = load_pairs_from_tsv(args.c3)
    c4_pairs = load_pairs_from_tsv(args.c4)
    test_splits = {"c2": c2_pairs, "c3": c3_pairs, "c4": c4_pairs}

    # Exact 5-fold split from Domain2vec.py / go2vec.py / rcnn.py:
    # kf = KFold(n_splits=5, random_state=63, shuffle=True)
    kf = KFold(n_splits=args.n_splits, random_state=args.seed, shuffle=True)
    kf_splits = list(kf.split(np.zeros(len(c1_pairs))))

    d2v_tensors = {
        k: torch.as_tensor(v, dtype=torch.float32)
        for k, v in assets["domain2vec"].items()
    }
    g2v_tensors = {
        k: torch.as_tensor(v, dtype=torch.float32)
        for k, v in assets["go2vec"].items()
    }

    # 1. Train Domain2vec 5-fold MLP (lr=0.01, amsgrad=True)
    d2v_oof_c1, d2v_ens_c1, d2v_test = train_5fold_branch(
        branch_name="Domain2vec",
        model_fn=lambda: DeepAraPPI_VectorMLP(in_dim=128),
        c1_pairs=c1_pairs,
        test_splits_dict=test_splits,
        dataset_cls=VectorPairDataset,
        feat_lookup=d2v_tensors,
        kf_splits=kf_splits,
        device=device,
        lr=0.01,
        batch_size=args.mlp_batch_size,
        eval_batch_size=2048,
        epochs=args.epochs,
        model_dir=args.model_dir,
    )

    # 2. Train GO2vec 5-fold MLP (lr=0.01, amsgrad=True)
    g2v_oof_c1, g2v_ens_c1, g2v_test = train_5fold_branch(
        branch_name="GO2vec",
        model_fn=lambda: DeepAraPPI_VectorMLP(in_dim=128),
        c1_pairs=c1_pairs,
        test_splits_dict=test_splits,
        dataset_cls=VectorPairDataset,
        feat_lookup=g2v_tensors,
        kf_splits=kf_splits,
        device=device,
        lr=0.01,
        batch_size=args.mlp_batch_size,
        eval_batch_size=2048,
        epochs=args.epochs,
        model_dir=args.model_dir,
    )

    # 3. Train Siamese RCNN 5-fold (lr=0.001, amsgrad=True)
    print("\nEncoding protein sequences into [32, 2000] word2vec tensors...")
    rcnn_seq_tensors = build_rcnn_sequence_tensors(
        assets["sequences"], assets["aa_vec"], seq_size=2000, dim=32
    )
    rcnn_oof_c1, rcnn_ens_c1, rcnn_test = train_5fold_branch(
        branch_name="RCNN",
        model_fn=lambda: DeepAraPPI_SiameseRCNN(dim=32, hidden_dim=50),
        c1_pairs=c1_pairs,
        test_splits_dict=test_splits,
        dataset_cls=RCNNPairDataset,
        feat_lookup=rcnn_seq_tensors,
        kf_splits=kf_splits,
        device=device,
        lr=0.001,
        batch_size=args.rcnn_batch_size,
        eval_batch_size=256,
        epochs=args.epochs,
        model_dir=args.model_dir,
    )

    # 4. Meta-Classifier: L2 Logistic Regression tuned via 5-fold GridSearchCV on C1 OOF predictions
    print("\n>>> Training Integrated L2 Logistic Regression Meta-Classifier (5-fold GridSearchCV) <<<")
    x_oof_c1 = np.column_stack([rcnn_oof_c1, d2v_oof_c1, g2v_oof_c1])
    x_ens_c1 = np.column_stack([rcnn_ens_c1, d2v_ens_c1, g2v_ens_c1])
    y_c1 = np.array([lbl for _, _, lbl in c1_pairs])

    base_lr = LogisticRegression(penalty="l2", solver="liblinear", random_state=args.seed)
    grid = GridSearchCV(
        base_lr,
        param_grid={"C": [0.001, 0.01, 0.1, 1.0, 10.0, 100.0]},
        cv=5,
        scoring="average_precision",
    )
    grid.fit(x_oof_c1, y_c1)
    best_lr = grid.best_estimator_
    print(f"  Best LR C={grid.best_params_['C']} | Coefs (RCNN, Domain2vec, GO2vec): {best_lr.coef_[0]}")

    lr_oof_c1 = cross_val_predict(
        best_lr, x_oof_c1, y_c1, cv=kf, method="predict_proba"
    )[:, 1]
    lr_ens_c1 = best_lr.predict_proba(x_ens_c1)[:, 1]
    lr_test = {}
    for split_k in ("c2", "c3", "c4"):
        x_t = np.column_stack(
            [rcnn_test[split_k], d2v_test[split_k], g2v_test[split_k]]
        )
        lr_test[split_k] = best_lr.predict_proba(x_t)[:, 1]

    with open(os.path.join(args.model_dir, "logistic_regression_meta.pkl"), "wb") as f:
        pickle.dump(best_lr, f)

    # 5. Save predictions & evaluate all 4 variants on C1, C2, C3, C4
    variants = [
        ("DeepAraPPI_RCNN", "deeparappi_rcnn", rcnn_oof_c1, rcnn_ens_c1, rcnn_test),
        (
            "DeepAraPPI_Domain2vec",
            "deeparappi_domain2vec",
            d2v_oof_c1,
            d2v_ens_c1,
            d2v_test,
        ),
        ("DeepAraPPI_GO2vec", "deeparappi_go2vec", g2v_oof_c1, g2v_ens_c1, g2v_test),
        (
            "DeepAraPPI_Integrated",
            "deeparappi_integrated",
            lr_oof_c1,
            lr_ens_c1,
            lr_test,
        ),
    ]

    for dir_name, prefix, c1_oof_scores, c1_ens_scores, test_dict in variants:
        out_dir = os.path.join(args.output_root, dir_name)
        c1_oof_csv = os.path.join(out_dir, f"{prefix}_c1_scores.csv")
        c1_ens_csv = os.path.join(out_dir, f"{prefix}_c1_ensemble_scores.csv")
        c2_csv = os.path.join(out_dir, f"{prefix}_c2_scores.csv")
        c3_csv = os.path.join(out_dir, f"{prefix}_c3_scores.csv")
        c4_csv = os.path.join(out_dir, f"{prefix}_c4_scores.csv")

        write_pred_csv(c1_pairs, c1_oof_scores, c1_oof_csv)
        write_pred_csv(c1_pairs, c1_ens_scores, c1_ens_csv)
        write_pred_csv(c2_pairs, test_dict["c2"], c2_csv)
        write_pred_csv(c3_pairs, test_dict["c3"], c3_csv)
        write_pred_csv(c4_pairs, test_dict["c4"], c4_csv)

        run_evaluate_pplm(c1_oof_csv, c2_csv, c3_csv, c4_csv, out_dir)


if __name__ == "__main__":
    main()
