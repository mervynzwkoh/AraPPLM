#!/usr/bin/env python3
"""
PPLM Batch PPI Prediction Script
Predicts protein-protein interactions for a dataset of protein pairs.

Usage:
    python batch_predict.py \
        --input ../data/DeepAraPPI/total_positive_negative_samples_DeepAraPPI.txt \
        --output results/deepara_c1_scores.csv \
        --batch_size 4 \
        --gpu_id 0 \
        --seq_db ../data/uniprot_final.pkl
"""

import os
import sys
import argparse
import csv
import pickle
import torch
import numpy as np
import time
from pathlib import Path

# Get the directory where THIS script is located
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Project root is the parent directory of scripts/
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
# PPLM submodule / package directory
PPLM_DIR = os.path.join(PROJECT_ROOT, "PPLM")
DATA_DIR = os.path.join(PROJECT_ROOT, "data")
sys.path.insert(0, PPLM_DIR)

# Debug: Print paths for verification
print(f"Script directory: {SCRIPT_DIR}")
print(f"Project root:     {PROJECT_ROOT}")
print(f"PPLM directory:   {PPLM_DIR}")
print(f"Data directory:   {DATA_DIR}")
print(f"sys.path[0]:      {sys.path[0]}")

# Verify PPLM can be imported
try:
    from pplm import PPLM, Alphabet
    from pplm_ppi import PPLM_PPI
    print("✅ PPLM imports successful")
except ImportError as e:
    print(f"❌ Import failed: {e}")
    sys.exit(1)

# Verify weights exist
weights_dir = os.path.join(PPLM_DIR, "weights")
if not os.path.exists(weights_dir):
    print(f"❌ Weights directory not found: {weights_dir}")
    sys.exit(1)

pplm_weight = os.path.join(weights_dir, "pplm_t33_650M.pt")
ppi_weight = os.path.join(weights_dir, "ppi_models.pkl")
if not os.path.exists(pplm_weight) or not os.path.exists(ppi_weight):
    print(f"❌ Missing weight files in {weights_dir}")
    sys.exit(1)

def load_models(device):
    """Load PPLM model and PPI classifier weights."""
    print("Loading PPLM model...")
    
    # Load PPLM backbone
    alphabet = Alphabet.from_architecture()
    batch_converter = alphabet.get_batch_converter()
    model_path = os.path.join(PPLM_DIR, "weights", "pplm_t33_650M.pt")
    model_data = torch.load(model_path, map_location="cpu")
    model_param = model_data["param"]
    model_state = model_data["model"]
    
    pplm_model = PPLM(
        num_layers=model_param["encoder_layers"],
        embed_dim=model_param["encoder_embed_dim"],
        attention_heads=model_param["encoder_attention_heads"],
        token_dropout=False,
        alphabet=alphabet,
    )
    pplm_model.to(device)
    pplm_model.load_state_dict(model_state, strict=False)
    pplm_model.eval()
    
    # Load PPI classifier weights
    print("Loading PPI classifier weights...")
    ppi_path = os.path.join(PPLM_DIR, "weights", "ppi_models.pkl")
    ppi_weights = torch.load(ppi_path, map_location=device)
    
    from pplm_ppi import PPLM_PPI
    ppi_model = PPLM_PPI()
    ppi_model.to(device)
    ppi_model.eval()
    
    print(f"Models loaded. Device: {device}")
    return pplm_model, batch_converter, ppi_model, ppi_weights


def get_features(model, batch_converter, seqA, seqB, device):
    """Extract features from PPLM for a protein pair."""
    seqA_labels, seqA_strs, seqA_tokens = batch_converter([("seqA", seqA)])
    seqB_labels, seqB_strs, seqB_tokens = batch_converter([("seqB", seqB)])
    tokens = torch.cat([seqA_tokens, seqB_tokens], dim=-1).to(device)
    
    inter_chain_mask = torch.ones(
        (len(seqA) + 2 + len(seqB) + 2, len(seqA) + 2 + len(seqB) + 2),
        device=device,
    )
    inter_chain_mask[: len(seqA) + 2, : len(seqA) + 2] = 0
    inter_chain_mask[len(seqA) + 2 :, len(seqA) + 2 :] = 0
    
    with torch.no_grad():
        out = model(
            tokens, inter_chain_mask,
            repr_layers=[33],
            need_head_weights=True,
            return_contacts=False,
        )
    
    embed_A = out["representations"][33][0, 1 : len(seqA) + 1, :]
    embed_B = out["representations"][33][0, -(len(seqB) + 1) : -1, :]
    
    attn_AA = out["attentions"].squeeze()[:, :, 1 : len(seqA) + 1, 1 : len(seqA) + 1].reshape(33 * 20, len(seqA), len(seqA))
    attn_AB = out["attentions"].squeeze()[:, :, 1 : len(seqA) + 1, -(len(seqB) + 1) : -1].reshape(33 * 20, len(seqA), len(seqB))
    attn_BA = out["attentions"].squeeze()[:, :, -(len(seqB) + 1) : -1, 1 : len(seqA) + 1].reshape(33 * 20, len(seqB), len(seqA))
    attn_BB = out["attentions"].squeeze()[:, :, -(len(seqB) + 1) : -1, -(len(seqB) + 1) : -1].reshape(33 * 20, len(seqB), len(seqB))
    inter_attn = (attn_AB + attn_BA.transpose(1, 2)) / 2
    
    features = {
        "mean_inter_attn": inter_attn.mean(dim=[1, 2]).unsqueeze(0),
        "mean_attn_AA": attn_AA.mean(dim=[1, 2]).unsqueeze(0),
        "mean_attn_BB": attn_BB.mean(dim=[1, 2]).unsqueeze(0),
        "mean_embed_A": embed_A.mean(dim=[0]).unsqueeze(0),
        "mean_embed_B": embed_B.mean(dim=[0]).unsqueeze(0),
        "max_inter_attn": torch.amax(inter_attn, dim=(1, 2)).unsqueeze(0),
        "max_attn_AA": torch.amax(attn_AA, dim=(1, 2)).unsqueeze(0),
        "max_attn_BB": torch.amax(attn_BB, dim=(1, 2)).unsqueeze(0),
        "max_embed_A": torch.amax(embed_A, dim=0).unsqueeze(0),
        "max_embed_B": torch.amax(embed_B, dim=0).unsqueeze(0),
    }
    
    return features


def predict_with_weights(ppi_model, features, ppi_weights):
    """Predict using the full 10-model ensemble (5 mean + 5 max classifier weights)."""
    with torch.no_grad():
        predictions_list = []
        
        # 1. 5 Mean-pooling classifiers
        for model_weight in ppi_weights['mean']:
            ppi_model.load_state_dict(model_weight)
            pred = ppi_model(
                features["mean_inter_attn"],
                features["mean_attn_AA"],
                features["mean_attn_BB"],
                features["mean_embed_A"],
                features["mean_embed_B"],
            )
            pred_swap = ppi_model(
                features["mean_inter_attn"],
                features["mean_attn_BB"],
                features["mean_attn_AA"],
                features["mean_embed_B"],
                features["mean_embed_A"],
            )
            pred_sym = (pred + pred_swap) / 2
            predictions_list.append(pred_sym)
            
        # 2. 5 Max-pooling classifiers
        for model_weight in ppi_weights['max']:
            ppi_model.load_state_dict(model_weight)
            pred = ppi_model(
                features["max_inter_attn"],
                features["max_attn_AA"],
                features["max_attn_BB"],
                features["max_embed_A"],
                features["max_embed_B"],
            )
            pred_swap = ppi_model(
                features["max_inter_attn"],
                features["max_attn_BB"],
                features["max_attn_AA"],
                features["max_embed_B"],
                features["max_embed_A"],
            )
            pred_sym = (pred + pred_swap) / 2
            predictions_list.append(pred_sym)
            
        # Average over all 10 ensemble models
        predictions = torch.stack(predictions_list)
        final_score = torch.mean(predictions, dim=0).squeeze().detach().cpu().numpy()
        return float(final_score)


def main():
    parser = argparse.ArgumentParser(description="PPLM Batch PPI Prediction")
    parser.add_argument("--input", required=True, help="Input TSV file (Protein1\tProtein2\tlabel)")
    parser.add_argument("--output", required=True, help="Output CSV file")
    parser.add_argument("--batch_size", type=int, default=4, help="Batch size for GPU inference")
    parser.add_argument("--gpu_id", type=int, default=0, help="GPU device ID")
    parser.add_argument("--seq_db", required=True, help="Path to pickled sequence dictionary")
    args = parser.parse_args()
    
    # Setup device
    device = torch.device(f"cuda:{args.gpu_id}" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")
    
    # Load models
    pplm_model, batch_converter, ppi_model, ppi_weights = load_models(device)
    
    # Load sequence database
    print(f"Loading sequence database: {args.seq_db}")
    with open(args.seq_db, "rb") as f:
        seq_db = pickle.load(f)
    print(f"Loaded {len(seq_db)} sequences")
    
    # Read input data
    print(f"Reading input: {args.input}")
    rows = []
    with open(args.input, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("Protein"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            rows.append((parts[0], parts[1], int(parts[2])))
    print(f"Total rows: {len(rows)}")
    
    # Predict
    print("Starting prediction...")
    start_time = time.time()
    results = []
    skipped = 0
    
    for i in range(0, len(rows), args.batch_size):
        batch = rows[i : i + args.batch_size]
        batch_scores = []
        
        for protA, protB, label in batch:
            # Get sequences
            seqA = seq_db.get(protA, None)
            seqB = seq_db.get(protB, None)
            
            if seqA is None or seqB is None:
                print(f"WARNING: Missing sequence for {protA} or {protB}, skipping")
                skipped += 1
                continue
            
            # Get features
            features = get_features(pplm_model, batch_converter, seqA, seqB, device)
            
            # Predict (average over all 10 weight sets)
            score = predict_with_weights(ppi_model, features, ppi_weights)
            batch_scores.append((f"{protA}:{protB}", score, label))
        
        results.extend(batch_scores)
        
        # Progress report
        if (i + len(batch)) % 1000 == 0:
            elapsed = time.time() - start_time
            rate = (i + len(batch)) / elapsed
            print(f"  Progress: {i + len(batch)}/{len(rows)} ({rate:.1f} pairs/sec)")
    
    # Save results
    print(f"Saving results to: {args.output}")
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    
    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["pair_id", "pred_score", "true_label"])
        writer.writerows(results)
    
    elapsed = time.time() - start_time
    print(f"Done! Predicted {len(results)} pairs in {elapsed:.1f} seconds ({len(results)/elapsed:.1f} pairs/sec)")
    print(f"Skipped {skipped} pairs due to missing sequences")


if __name__ == "__main__":
    main()