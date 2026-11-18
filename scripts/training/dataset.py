#!/usr/bin/env python3
"""
PyTorch Dataset and CV Split Utilities for PPLM-PPI Retraining.

Aligned with the original PPLM author's implementation
(ppi_training_code/dataset.py):
  - Dataset_PPI_Interact loads all 15 pooled attention/embedding features
    (mean, max, min) plus "sample" ID and "label" from pre-extracted .pkl files.
  - Supports order_mode ("random", "normal", "reverse", "both") for 50% random
    A<->B protein-order swapping during training and deterministic evaluation
    during testing.
  - Includes generate_cv_split_files() to persist stratified 10-fold
    .trainset / .testset files for per-fold (--cross_round) execution.
  - Retains PlantPPIDataset and load_pair_list() for backward compatibility.
"""

import os
import pickle
import random
from os.path import isfile, join

import numpy as np
import torch
from sklearn.model_selection import StratifiedKFold
from torch.utils.data import Dataset


POOLING_MODES = ("mean", "max", "min")


class Dataset_PPI_Interact(Dataset):
    """
    Dataset for pre-extracted PPLM features, matching ppi_training_code/dataset.py.

    Args:
        target_list:    Path to a split/dataset file (.trainset, .testset, or TSV)
                        OR a Python list of sample IDs / (sample_id, label) tuples.
        esm_dir:        Directory containing {sample}.pkl files.
        esm_flag:       Optional flag for single-chain ESM features (default: False).
        esm_chain_dir:  Directory for single-chain ESM features (if esm_flag=True).
        order_mode:     "random" (50% A<->B swap), "normal" / "both" (A, B order),
                        or "reverse" (B, A order).
        debug:          If True, limit number of samples to n_debug_target.
        n_debug_target: Maximum samples when debug=True.
        crop_size:      Kept for API compatibility with original author's constructor.
        n_layer:        Kept for API compatibility with original author's constructor.
        validation:     Validation flag (disables shuffling).
        testing:        Testing flag (disables shuffling and logs missing files).
        e_45:           Kept for API compatibility with original author's constructor.
    """

    def __init__(
        self,
        target_list,
        esm_dir="/mnt/rna01/zhy/tools/PPLM-PPI_data/D_script/PPLM_feats/",
        esm_flag=False,
        esm_chain_dir="",
        order_mode="random",
        debug=False,
        n_debug_target=100,
        crop_size=512,
        n_layer=12,
        validation=False,
        testing=False,
        e_45=False,
    ):
        self.esm_dir = esm_dir
        self.esm_flag = esm_flag
        self.esm_chain_dir = esm_chain_dir
        self.order_mode = order_mode
        self.testing = testing
        self.validation = validation
        self.n_layer = n_layer
        self.crop_size = crop_size
        self.e_45 = e_45

        print("self.esm_dir:", self.esm_dir)

        raw_samples = self._parse_target_list(target_list)
        self.samples = []
        missing = 0

        for sample in raw_samples:
            esm_path = join(self.esm_dir, sample + ".pkl")
            if isfile(esm_path):
                self.samples.append(sample)
            else:
                missing += 1

            if debug and len(self.samples) >= n_debug_target:
                break

        if missing > 0:
            print(f"WARNING: {missing} samples missing .pkl files in {self.esm_dir}, skipped")

        if not validation and not testing:
            np.random.shuffle(self.samples)

        print(
            f"Dataset_PPI_Interact | number of targets: {len(self.samples)} | "
            f"order_mode={self.order_mode}"
        )

    @staticmethod
    def _parse_target_list(target_list):
        if isinstance(target_list, (list, tuple)):
            parsed = []
            for item in target_list:
                if isinstance(item, (list, tuple)):
                    parsed.append(str(item[0]))
                else:
                    parsed.append(str(item))
            return parsed

        parsed = []
        with open(target_list, "r") as f:
            for line in f:
                line = line.strip()
                if (
                    not line
                    or line.startswith("#")
                    or line.startswith("Protein")
                    or line.startswith("participant_a_id")
                ):
                    continue
                parts = line.split()
                if len(parts) >= 3 and parts[2] == "label":
                    continue
                if len(parts) >= 3 and "@" not in parts[0]:
                    # Raw TSV line: Prot1 Prot2 label [score/partition...]
                    sample = f"{parts[0]}@{parts[1]}"
                else:
                    sample = parts[0]
                parsed.append(sample)
        return parsed

    def __len__(self):
        return len(self.samples)

    def _should_swap(self):
        if self.order_mode == "random":
            return random.random() >= 0.5
        if self.order_mode == "reverse":
            return True
        # "both" and "normal" return canonical (A, B) order;
        # testing loops run the explicit swapped pass when order_mode == "both".
        return False

    def __getitem__(self, idx, pindex=-1):
        if torch.is_tensor(idx):
            idx = idx.tolist()

        sample = self.samples[idx]
        attempts = 0

        while attempts < 10:
            esm_path = join(self.esm_dir, sample + ".pkl")
            try:
                with open(esm_path, "rb") as fr:
                    esm_data = pickle.load(fr)

                swap = self._should_swap()
                data = {
                    "sample": sample,
                    "label": int(esm_data["label"]),
                }

                for mode in POOLING_MODES:
                    embed_a = esm_data[f"{mode}_embed_A"]
                    embed_b = esm_data[f"{mode}_embed_B"]
                    attn_aa = esm_data[f"{mode}_attn_AA"]
                    attn_bb = esm_data[f"{mode}_attn_BB"]
                    inter_attn = esm_data[f"{mode}_inter_attn"]

                    if swap:
                        data[f"{mode}_embed_A"] = embed_b
                        data[f"{mode}_embed_B"] = embed_a
                        data[f"{mode}_attn_AA"] = attn_bb
                        data[f"{mode}_attn_BB"] = attn_aa
                    else:
                        data[f"{mode}_embed_A"] = embed_a
                        data[f"{mode}_embed_B"] = embed_b
                        data[f"{mode}_attn_AA"] = attn_aa
                        data[f"{mode}_attn_BB"] = attn_bb
                    data[f"{mode}_inter_attn"] = inter_attn

                if self.esm_flag and self.esm_chain_dir:
                    chain_a_name, chain_b_name = sample.split("@", 1)
                    with open(join(self.esm_chain_dir, chain_a_name + ".fa.pkl"), "rb") as fr:
                        chain_a_data = pickle.load(fr)
                    with open(join(self.esm_chain_dir, chain_b_name + ".fa.pkl"), "rb") as fr:
                        chain_b_data = pickle.load(fr)

                    for mode in POOLING_MODES:
                        ce_a = chain_a_data[f"{mode}_embed"]
                        ce_b = chain_b_data[f"{mode}_embed"]
                        ca_a = chain_a_data[f"{mode}_attn"]
                        ca_b = chain_b_data[f"{mode}_attn"]
                        if swap:
                            data[f"chain_{mode}_embed_A"] = ce_b
                            data[f"chain_{mode}_embed_B"] = ce_a
                            data[f"chain_{mode}_attn_A"] = ca_b
                            data[f"chain_{mode}_attn_B"] = ca_a
                        else:
                            data[f"chain_{mode}_embed_A"] = ce_a
                            data[f"chain_{mode}_embed_B"] = ce_b
                            data[f"chain_{mode}_attn_A"] = ca_a
                            data[f"chain_{mode}_attn_B"] = ca_b

                return data

            except Exception as e:
                attempts += 1
                print(f"Warning: {sample} ({esm_path}): {e}, re-sampling...")
                new_idx = np.random.choice(len(self.samples))
                sample = self.samples[new_idx]

        raise RuntimeError(f"Failed to load feature file after 10 attempts (last: {esm_path})")


class PlantPPIDataset(Dataset):
    """
    Backward-compatible single-mode dataset wrapper around pre-extracted .pkl files.
    """

    def __init__(self, sample_list, feat_dir, mode="mean", shuffle=True, augment=True):
        order_mode = "random" if augment else "normal"
        self._inner = Dataset_PPI_Interact(
            target_list=sample_list,
            esm_dir=feat_dir,
            order_mode=order_mode,
            validation=not shuffle,
            testing=not shuffle,
        )
        self.mode = mode

    def __len__(self):
        return len(self._inner)

    def __getitem__(self, idx):
        full = self._inner[idx]
        prefix = self.mode
        return {
            "embed_A": full[f"{prefix}_embed_A"],
            "embed_B": full[f"{prefix}_embed_B"],
            "attn_AA": full[f"{prefix}_attn_AA"],
            "attn_BB": full[f"{prefix}_attn_BB"],
            "inter_attn": full[f"{prefix}_inter_attn"],
            "label": full["label"],
        }


def load_pair_list(dataset_path):
    """
    Load a PPI dataset file (TSV format: Protein1\\tProtein2\\tlabel) and
    return a list of (pair_id, label) tuples.
    """
    pairs = []
    with open(dataset_path, "r") as f:
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
            prot_a = parts[0]
            prot_b = parts[1]
            label = int(parts[2])
            pair_id = f"{prot_a}@{prot_b}"
            pairs.append((pair_id, label))

    print(f"Loaded {len(pairs)} pairs from {dataset_path}")
    return pairs


def discover_labeled_samples(feat_dir, pair_list_path=None):
    """
    Discover valid (pair_id, label) pairs that have existing .pkl files in feat_dir.
    Uses pair_list_path if provided to avoid loading every .pkl file during splitting.
    """
    samples = []
    seen = set()

    if pair_list_path and isfile(pair_list_path):
        for pair_id, label in load_pair_list(pair_list_path):
            if pair_id in seen:
                continue
            if isfile(join(feat_dir, pair_id + ".pkl")):
                samples.append((pair_id, int(label)))
                seen.add(pair_id)
        print(f"Matched {len(samples)} existing .pkl samples against {pair_list_path}")
        return samples

    pkl_files = sorted(f for f in os.listdir(feat_dir) if f.endswith(".pkl"))
    print(f"Scanning {len(pkl_files)} .pkl files in {feat_dir} for labels...")
    for fname in pkl_files:
        pair_id = fname[:-4]
        pkl_path = join(feat_dir, fname)
        try:
            with open(pkl_path, "rb") as f:
                data = pickle.load(f)
            samples.append((pair_id, int(data["label"])))
        except Exception as e:
            print(f"WARNING: Could not read {fname}: {e}")

    print(f"Discovered {len(samples)} valid samples in {feat_dir}")
    return samples


def generate_cv_split_files(
    feat_dir,
    split_dir,
    dataset_prefix="plant",
    pair_list_path=None,
    n_folds=10,
    seed=73591024,
    force=False,
):
    """
    Persist stratified K-fold .trainset and .testset split files on disk so each
    fold (--cross_round 0..9) runs as an independent, resumable process matching
    ppi_training_code/run_training.sh.

    Returns:
        List of (trainset_path, testset_path) tuples indexed by fold (0..n_folds-1).
    """
    os.makedirs(split_dir, exist_ok=True)
    split_paths = [
        (
            join(split_dir, f"{dataset_prefix}.cross_validation_{k}.trainset"),
            join(split_dir, f"{dataset_prefix}.cross_validation_{k}.testset"),
        )
        for k in range(n_folds)
    ]

    if not force and all(isfile(tr) and isfile(va) for tr, va in split_paths):
        return split_paths

    samples = discover_labeled_samples(feat_dir, pair_list_path=pair_list_path)
    if len(samples) < n_folds:
        raise ValueError(
            f"Cannot create {n_folds}-fold splits from only {len(samples)} samples in {feat_dir}"
        )

    pair_ids = [s[0] for s in samples]
    labels = np.array([s[1] for s in samples], dtype=np.int64)

    print(
        f"Generating {n_folds}-fold StratifiedKFold split files in {split_dir} "
        f"(seed={seed}, total={len(samples)}, pos={int(np.sum(labels))}, "
        f"neg={len(labels) - int(np.sum(labels))})"
    )

    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=seed)
    for fold_idx, (train_idx, valid_idx) in enumerate(skf.split(pair_ids, labels)):
        train_path, valid_path = split_paths[fold_idx]
        with open(train_path, "w") as f:
            for idx in train_idx:
                f.write(pair_ids[idx] + "\n")
        with open(valid_path, "w") as f:
            for idx in valid_idx:
                f.write(pair_ids[idx] + "\n")

    return split_paths
