#!/usr/bin/env python3
"""
PPI Classifier Model Definition for Plant-PPLM Training.

Architecture: PPI_inter_intra_attn_embed_single_pooling
  - 3 linear projection heads (inter-attn, intra-attn, embedding) → 660 dims each
  - 5-way concatenation → 3300 dims
  - 5-layer MLP: 3300→1024→512→256→128→1 (LayerNorm + ReLU, Sigmoid output)

This is architecturally identical to the original PPLM-PPI model
(junliu621/PPLM/pplm_ppi/model.py) to ensure weight compatibility.

Reference: Liu, Chen & Zhang, Nat. Commun. 2026
"""

import torch
from torch import nn
import numpy as np
from sklearn.metrics import (
    roc_auc_score,
    accuracy_score,
    f1_score,
    precision_score,
    recall_score,
    precision_recall_curve,
    auc,
    matthews_corrcoef,
    confusion_matrix,
)


class PPI_inter_intra_attn_embed_single_pooling(nn.Module):
    """
    PPLM-PPI classifier MLP head.

    Takes pre-pooled features from one pooling strategy (mean OR max OR min):
      - inter_attn:  [batch, 660]  inter-protein cross-attention (33 layers × 20 heads)
      - intra_attn_A: [batch, 660]  intra-protein self-attention for protein A
      - intra_attn_B: [batch, 660]  intra-protein self-attention for protein B
      - embed_A:     [batch, 1280] per-residue embedding (mean/max/min pooled) for protein A
      - embed_B:     [batch, 1280] per-residue embedding (mean/max/min pooled) for protein B

    Architecture:
      - linear_inter:  660 → 660  (projects inter-protein attention)
      - linear_intra:  660 → 660  (shared weight for both intra-protein attentions)
      - linear_embed: 1280 → 660  (projects embeddings to match attention dim)
      - Concatenate → 660 × 5 = 3300
      - MLP: 3300 → 1024 → 512 → 256 → 128 → 1 (Sigmoid)
    """

    def __init__(self, input_dim=660, embedding_dim=1280):
        super(PPI_inter_intra_attn_embed_single_pooling, self).__init__()
        prev_dim = input_dim
        self.linear_intra = nn.Linear(prev_dim, prev_dim)
        self.linear_inter = nn.Linear(prev_dim, prev_dim)
        self.linear_embed = nn.Linear(embedding_dim, prev_dim)

        layers = []
        layers.append(nn.Linear(prev_dim * 5, 1024))
        layers.append(nn.LayerNorm(1024))
        layers.append(nn.ReLU())
        layers.append(nn.Linear(1024, 512))
        layers.append(nn.LayerNorm(512))
        layers.append(nn.ReLU())
        layers.append(nn.Linear(512, 256))
        layers.append(nn.LayerNorm(256))
        layers.append(nn.ReLU())
        layers.append(nn.Linear(256, 128))
        layers.append(nn.LayerNorm(128))
        layers.append(nn.ReLU())
        layers.append(nn.Linear(128, 1))
        layers.append(nn.Sigmoid())

        self.mlp = nn.Sequential(*layers)

    def forward(self, inter_attn, intra_attn_A, intra_attn_B, embed_A, embed_B):
        inter_attn = self.linear_inter(inter_attn)
        intra_attn_A = self.linear_intra(intra_attn_A)
        intra_attn_B = self.linear_intra(intra_attn_B)
        embed_A = self.linear_embed(embed_A)
        embed_B = self.linear_embed(embed_B)

        features = torch.cat(
            [inter_attn, intra_attn_A, intra_attn_B, embed_A, embed_B], dim=-1
        )

        preds = self.mlp(features)  # Shape: (batch_size, 1)
        return preds


def evaluate(GT, Pre, thresh=0.5):
    """
    Compute comprehensive evaluation metrics for PPI prediction.

    Returns:
        precision, recall, accuracy, F1_score, specificity, MCC,
        Top10, Top20, Top50, AUC_ROC, AUC_PR
    """
    GT = np.array(GT)
    Pre = np.array(Pre)

    # Top-K precision
    idx_list = np.argsort(-Pre)
    GT_sorted = [GT[idx] for idx in idx_list]
    Top10 = np.sum(GT_sorted[:10]) / 10
    Top20 = np.sum(GT_sorted[:20]) / 20
    Top50 = np.sum(GT_sorted[:50]) / 50

    # Threshold-free metrics
    AUC_ROC = roc_auc_score(GT, Pre)
    precision_list, recall_list, _ = precision_recall_curve(GT, Pre)
    AUC_PR = auc(recall_list, precision_list)

    # Threshold-dependent metrics
    Pre_binary = [1 if item > thresh else 0 for item in Pre]
    accuracy = accuracy_score(GT, Pre_binary)
    recall = recall_score(GT, Pre_binary, zero_division=0)
    precision = precision_score(GT, Pre_binary, zero_division=0)
    F1 = f1_score(GT, Pre_binary, average="binary", zero_division=0)
    MCC = matthews_corrcoef(GT, Pre_binary)

    # Specificity (True Negative Rate)
    tn, fp, fn, tp = confusion_matrix(GT, Pre_binary, labels=[0, 1]).ravel()
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0

    return precision, recall, accuracy, F1, specificity, MCC, Top10, Top20, Top50, AUC_ROC, AUC_PR
