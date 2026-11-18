import torch
from torch import nn
import numpy as np
from sklearn.metrics import roc_auc_score,accuracy_score,f1_score,precision_score,recall_score,precision_recall_curve,auc
import os, sys


#####################################################################################################

class PPI_inter_intra_attn_embed_single_pooling(nn.Module):
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

        # print(inter_attn.shape, intra_attn_A.shape, intra_attn_B.shape, embed_A.shape, embed_B.shape)

        features = torch.cat([inter_attn, intra_attn_A, intra_attn_B, embed_A, embed_B], dim=-1)   #for  batchsize > 1

        preds = self.mlp(features)  # Shape: (batch_size, 1)

        return preds


def evaluate(GT, Pre, thresh=0.5):
    GT = np.array(GT)
    Pre = np.array(Pre)

    # Top 10/20/50
    idx_list = np.argsort(-Pre)
    GT_new = [GT[idx] for idx in idx_list]

    Top10 = np.sum(GT_new[:10]) / 10
    Top20 = np.sum(GT_new[:20]) / 20
    Top50 = np.sum(GT_new[:50]) / 50

    AUC_ROC = roc_auc_score(GT, Pre)
    precision_list, recall_list, _ = precision_recall_curve(GT, Pre)
    AUC_PR = auc(recall_list, precision_list)

    # Pre = [1 if item>0.5 else 0 for item in Pre]
    Pre = [1 if item > thresh else 0 for item in Pre]
    accuracy = accuracy_score(GT, Pre)
    recall = recall_score(GT, Pre)
    precision = precision_score(GT, Pre)
    F1_score = f1_score(GT, Pre, average='binary')

    return precision, recall, accuracy, F1_score, Top10, Top20, Top50, AUC_ROC, AUC_PR
