#!/usr/bin/env python3
"""
PyTorch Dataset for loading pre-extracted PPLM features from .pkl files.

Each .pkl file contains pooled attention and embedding features for one protein pair,
produced by extract_features.py. This dataset is used by train_ppi_head.py for
training the PPI classifier MLP head.

Key design decisions matching the original PPLM training pipeline:
  - 50% random protein-order swapping (A↔B) to enforce interaction symmetry
  - Supports mean/max/min pooling mode selection
  - Resilient to corrupted/missing .pkl files (random re-sample on error)
"""

import torch
from torch.utils.data import Dataset
import numpy as np
import pickle
import random
from os.path import join, isfile
import os


class PlantPPIDataset(Dataset):
    """
    Dataset for pre-extracted PPLM features.

    Args:
        sample_list: List of (pair_id, label) tuples, where pair_id is the
                     filename stem (e.g., "Q9FHY1@Q9LP46")
        feat_dir:    Directory containing {pair_id}.pkl files
        mode:        Pooling mode to load — "mean", "max", or "min"
        shuffle:     Whether to shuffle samples on init (True for training)
        augment:     Whether to apply 50% random A↔B swapping (True for training)
    """

    def __init__(
        self,
        sample_list,
        feat_dir,
        mode="mean",
        shuffle=True,
        augment=True,
    ):
        self.feat_dir = feat_dir
        self.mode = mode
        self.augment = augment

        # Validate that all samples have corresponding .pkl files
        self.samples = []
        missing = 0
        for pair_id, label in sample_list:
            pkl_path = join(self.feat_dir, pair_id + ".pkl")
            if isfile(pkl_path):
                self.samples.append((pair_id, label))
            else:
                missing += 1

        if missing > 0:
            print(f"WARNING: {missing} samples missing .pkl files, skipped")

        if shuffle:
            np.random.shuffle(self.samples)

        print(
            f"PlantPPIDataset | mode={mode} | samples={len(self.samples)} | "
            f"augment={augment}"
        )

    def __len__(self):
        return len(self.samples)

    def __getitem__(self, idx):
        if torch.is_tensor(idx):
            idx = idx.tolist()

        pair_id, label = self.samples[idx]
        success = False
        attempts = 0

        while not success and attempts < 10:
            try:
                pkl_path = join(self.feat_dir, pair_id + ".pkl")
                with open(pkl_path, "rb") as f:
                    feat_data = pickle.load(f)

                # Select the correct pooling mode features
                prefix = self.mode  # "mean", "max", or "min"

                embed_A = feat_data[f"{prefix}_embed_A"]
                embed_B = feat_data[f"{prefix}_embed_B"]
                attn_AA = feat_data[f"{prefix}_attn_AA"]
                attn_BB = feat_data[f"{prefix}_attn_BB"]
                inter_attn = feat_data[f"{prefix}_inter_attn"]

                # 50% random A↔B protein-order swapping (matching PPLM training)
                if self.augment and random.random() < 0.5:
                    data = {
                        "embed_A": embed_B,
                        "embed_B": embed_A,
                        "attn_AA": attn_BB,
                        "attn_BB": attn_AA,
                        "inter_attn": inter_attn,  # symmetric, unchanged
                        "label": label,
                    }
                else:
                    data = {
                        "embed_A": embed_A,
                        "embed_B": embed_B,
                        "attn_AA": attn_AA,
                        "attn_BB": attn_BB,
                        "inter_attn": inter_attn,
                        "label": label,
                    }

                success = True

            except Exception as e:
                attempts += 1
                print(f"WARNING: Error loading {pair_id}: {e}, re-sampling...")
                new_idx = np.random.choice(len(self.samples))
                pair_id, label = self.samples[new_idx]

        return data


def load_pair_list(dataset_path):
    """
    Load a PPI dataset file (TSV format: Protein1\tProtein2\tlabel) and
    return a list of (pair_id, label) tuples.

    The pair_id is formatted as "ProtA@ProtB" to match the .pkl filenames
    produced by extract_features.py.

    Handles both DeepAraPPI format (3 columns with optional header) and
    ESMAraPPI format (4 columns: Prot1, Prot2, label, score).
    """
    pairs = []
    with open(dataset_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("Protein"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            prot_a = parts[0]
            prot_b = parts[1]
            label = int(parts[2])
            pair_id = f"{prot_a}@{prot_b}"
            pairs.append((pair_id, label))

    print(f"Loaded {len(pairs)} pairs from {dataset_path}")
    return pairs
