#!/usr/bin/env python3
"""
PPLM Feature Extraction for Plant PPI Training.

Extracts pre-pooled PPLM backbone features for all protein pairs in a C1
training dataset and saves them as individual .pkl files. These are consumed
by train_ppi_head.py during Stage 2 of the training pipeline.

For each pair, the frozen PPLM backbone (33-layer Transformer, 650M params)
produces:
  - Per-residue embeddings from layer 33: [L, 1280]
  - Full attention tensor: [33, 20, L, L]

These are immediately pooled into compact vectors:
  - Inter-protein attention: mean/max/min over spatial dims → [660]
  - Intra-protein attention A/B: mean/max/min → [660]
  - Per-protein embeddings A/B: mean/max/min → [1280]

Combined cropping (LA + LB ≤ 1020) and OOM resilience match the existing
benchmarking pipeline (scripts/benchmarking/batch_predict.py).

Usage:
    python scripts/training/extract_features.py \
        --input data/DeepAraPPI/c1_ppi_sample_DeepAraPPI.txt \
        --output_dir features/DeepAraPPI_C1/ \
        --seq_db data/arabidopsis/uniprot_final.pkl \
        --gpu_id 0
"""

import os
import sys
import argparse
import pickle
import torch
import numpy as np
import time

# Path setup
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))
PPLM_DIR = os.path.join(PROJECT_ROOT, "PPLM")
sys.path.insert(0, PPLM_DIR)

from pplm import PPLM, Alphabet


def load_pplm_backbone(device):
    """Load the frozen PPLM backbone model."""
    print("Loading PPLM backbone...")
    alphabet = Alphabet.from_architecture()
    batch_converter = alphabet.get_batch_converter()

    model_path = os.path.join(PPLM_DIR, "weights", "pplm_t33_650M.pt")
    if not os.path.exists(model_path):
        print(f"ERROR: PPLM weights not found at {model_path}")
        print("Run: cd PPLM/weights && bash download_weights.sh")
        sys.exit(1)

    model_data = torch.load(model_path, map_location="cpu", weights_only=False)
    model_param = model_data["param"]
    model_state = model_data["model"]

    model = PPLM(
        num_layers=model_param["encoder_layers"],
        embed_dim=model_param["encoder_embed_dim"],
        attention_heads=model_param["encoder_attention_heads"],
        token_dropout=False,
        alphabet=alphabet,
    )
    model.to(device)
    model.load_state_dict(model_state, strict=False)
    model.eval()

    print(f"PPLM backbone loaded on {device}")
    return model, batch_converter


def extract_pair_features(model, batch_converter, seqA, seqB, device, max_pair_len=1020):
    """
    Extract all pooled features for a single protein pair.

    Returns a dict containing 15 feature tensors (3 pooling × 5 feature types)
    plus the sequence lengths.
    """
    # Combined cropping (matching PPLM pretraining: LA + LB ≤ 1020)
    total_len = len(seqA) + len(seqB)
    if total_len > max_pair_len:
        ratio_A = len(seqA) / total_len
        budget_A = max(50, int(max_pair_len * ratio_A))
        budget_B = max_pair_len - budget_A
        seqA = seqA[:budget_A]
        seqB = seqB[:budget_B]

    # Tokenize
    _, _, seqA_tokens = batch_converter([("seqA", seqA)])
    _, _, seqB_tokens = batch_converter([("seqB", seqB)])
    tokens = torch.cat([seqA_tokens, seqB_tokens], dim=-1).to(device)

    # Inter-chain mask
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

        # Extract per-residue embeddings (layer 33)
        embed_A = out["representations"][33][0, 1 : lenA + 1, :]  # [lenA, 1280]
        embed_B = out["representations"][33][0, -(lenB + 1) : -1, :]  # [lenB, 1280]

        # Extract attention matrices
        attns = out["attentions"].squeeze()  # [33, 20, L, L]
        n_heads = 33 * 20

        attn_AA = attns[:, :, 1 : lenA + 1, 1 : lenA + 1].reshape(n_heads, lenA, lenA)
        attn_AB = attns[:, :, 1 : lenA + 1, -(lenB + 1) : -1].reshape(n_heads, lenA, lenB)
        attn_BA = attns[:, :, -(lenB + 1) : -1, 1 : lenA + 1].reshape(n_heads, lenB, lenA)
        attn_BB = attns[:, :, -(lenB + 1) : -1, -(lenB + 1) : -1].reshape(n_heads, lenB, lenB)

        # Symmetrized inter-protein attention
        inter_attn = (attn_AB + attn_BA.transpose(1, 2)) / 2

        # Pool features across spatial dimensions (3 pooling strategies)
        features = {
            # Mean pooling
            "mean_inter_attn": inter_attn.mean(dim=[1, 2]).cpu(),
            "mean_attn_AA": attn_AA.mean(dim=[1, 2]).cpu(),
            "mean_attn_BB": attn_BB.mean(dim=[1, 2]).cpu(),
            "mean_embed_A": embed_A.mean(dim=0).cpu(),
            "mean_embed_B": embed_B.mean(dim=0).cpu(),
            # Max pooling
            "max_inter_attn": torch.amax(inter_attn, dim=(1, 2)).cpu(),
            "max_attn_AA": torch.amax(attn_AA, dim=(1, 2)).cpu(),
            "max_attn_BB": torch.amax(attn_BB, dim=(1, 2)).cpu(),
            "max_embed_A": torch.amax(embed_A, dim=0).cpu(),
            "max_embed_B": torch.amax(embed_B, dim=0).cpu(),
            # Min pooling
            "min_inter_attn": torch.amin(inter_attn, dim=(1, 2)).cpu(),
            "min_attn_AA": torch.amin(attn_AA, dim=(1, 2)).cpu(),
            "min_attn_BB": torch.amin(attn_BB, dim=(1, 2)).cpu(),
            "min_embed_A": torch.amin(embed_A, dim=0).cpu(),
            "min_embed_B": torch.amin(embed_B, dim=0).cpu(),
            # Metadata
            "lens": (lenA, lenB),
        }

        # Free VRAM immediately
        del out, attns, attn_AA, attn_AB, attn_BA, attn_BB
        del inter_attn, embed_A, embed_B, tokens, inter_chain_mask

    return features


def main():
    parser = argparse.ArgumentParser(
        description="Extract PPLM backbone features for plant PPI training"
    )
    parser.add_argument(
        "--input", required=True, help="Input TSV file (Protein1\\tProtein2\\tlabel)"
    )
    parser.add_argument(
        "--output_dir", required=True, help="Output directory for .pkl feature files"
    )
    parser.add_argument(
        "--seq_db", required=True, help="Path to pickled sequence dictionary (.pkl)"
    )
    parser.add_argument("--gpu_id", type=int, default=0, help="GPU device ID")
    parser.add_argument(
        "--max_pair_len",
        type=int,
        default=1020,
        help="Maximum combined sequence pair length (LA + LB)",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        default=False,
        help="Skip pairs that already have .pkl files (for resumption)",
    )
    args = parser.parse_args()

    # Setup device
    device = torch.device(f"cuda:{args.gpu_id}" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    # Load PPLM backbone
    model, batch_converter = load_pplm_backbone(device)

    # Load sequence database
    print(f"Loading sequence database: {args.seq_db}")
    with open(args.seq_db, "rb") as f:
        seq_db = pickle.load(f)
    print(f"Loaded {len(seq_db)} sequences")

    # Read input pairs
    print(f"Reading input: {args.input}")
    pairs = []
    with open(args.input, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("Protein"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            pairs.append((parts[0], parts[1], int(parts[2])))
    print(f"Total pairs: {len(pairs)}")

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Extract features
    start_time = time.time()
    extracted = 0
    skipped = 0
    errors = 0

    for i, (protA, protB, label) in enumerate(pairs):
        pair_id = f"{protA}@{protB}"
        output_path = os.path.join(args.output_dir, pair_id + ".pkl")

        # Resume support: skip if already extracted
        if args.resume and os.path.isfile(output_path):
            skipped += 1
            continue

        # Lookup sequences
        seqA = seq_db.get(protA)
        seqB = seq_db.get(protB)
        if seqA is None or seqB is None:
            print(f"WARNING: Missing sequence for {protA} or {protB}, skipping")
            errors += 1
            continue

        try:
            features = extract_pair_features(
                model, batch_converter, seqA, seqB, device, args.max_pair_len
            )
            features["label"] = label

            # Save to .pkl
            with open(output_path, "wb") as f:
                pickle.dump(features, f, protocol=pickle.HIGHEST_PROTOCOL)

            extracted += 1

        except torch.cuda.OutOfMemoryError:
            torch.cuda.empty_cache()
            try:
                features = extract_pair_features(
                    model, batch_converter, seqA, seqB, device, max_pair_len=800
                )
                features["label"] = label
                with open(output_path, "wb") as f:
                    pickle.dump(features, f, protocol=pickle.HIGHEST_PROTOCOL)
                extracted += 1
            except Exception as e:
                print(f"WARNING: Failed pair {pair_id} ({e}), skipping")
                errors += 1
                torch.cuda.empty_cache()

        except Exception as e:
            print(f"WARNING: Error for {pair_id} ({e}), skipping")
            errors += 1

        # Periodic VRAM cleanup and progress
        if (i + 1) % 200 == 0:
            if torch.cuda.is_available():
                torch.cuda.empty_cache()

        if (i + 1) % 500 == 0 or (i + 1) == len(pairs):
            elapsed = time.time() - start_time
            rate = (extracted + skipped) / max(0.1, elapsed)
            print(
                f"  Progress: {i + 1}/{len(pairs)} | "
                f"Extracted: {extracted} | Skipped: {skipped} | "
                f"Errors: {errors} | Rate: {rate:.1f} pairs/sec"
            )

    elapsed = time.time() - start_time
    print(f"\nDone! Extracted {extracted} features in {elapsed:.1f}s")
    print(f"Skipped (resume): {skipped} | Errors: {errors}")
    print(f"Output directory: {args.output_dir}")


if __name__ == "__main__":
    main()
