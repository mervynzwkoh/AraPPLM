#!/usr/bin/env python3
"""
PPLM Paired Embedding Extractor for Arabidopsis Positive PPIs
============================================================
Extracts 2,560-dimensional paired-protein embeddings from the frozen PPLM
backbone (ESM2-650M architecture, 33 Transformer layers with inter-protein
cross-attention) across all deduplicated positive pairs.

Embedding Representation:
  Concatenation of mean-pooled and max-pooled representations from layer 33
  across the full interacting sequence pair [mean(seqA+seqB), max(seqA+seqB)] -> 2,560 dims.

Key Features:
  - Supports streaming from .parquet or .tsv
  - Sequence resolution from uniprot_final.pkl with GN= and isoform fallback
  - Proportional sequence length cropping (LA + LB <= 1020)
  - Memory-mapped checkpointing: resumable after HPC walltime termination
  - Progress reporting and execution audit logging

Usage:
  python scripts/hpc/extract_ara_pplm_embeddings.py \
      --input data/processed/arabidopsis/ara_positive_pairs_deduplicated.parquet \
      --seq_db data/arabidopsis/uniprot_final.pkl \
      --output data/processed/arabidopsis/ara_positive_embeddings_2560d.npy \
      --meta_output data/processed/arabidopsis/ara_positive_embeddings_meta.parquet \
      --batch_size 8 \
      --device cuda
"""

import os
import sys
import time
import argparse
import pickle
import json
import torch
import numpy as np
import pandas as pd

# Path setup to import pplm from submodule
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))
PPLM_DIR = os.path.join(PROJECT_ROOT, "PPLM")
if PPLM_DIR not in sys.path:
    sys.path.insert(0, PPLM_DIR)

from pplm import PPLM, Alphabet


def parse_args():
    parser = argparse.ArgumentParser(description="Extract PPLM paired embeddings for Arabidopsis positive PPIs")
    parser.add_argument("--input", default="data/processed/arabidopsis/ara_positive_pairs_deduplicated.parquet",
                        help="Path to deduplicated pairs parquet or tsv")
    parser.add_argument("--seq_db", default="data/arabidopsis/uniprot_final.pkl",
                        help="Path to indexed UniProt sequence pickle")
    parser.add_argument("--output", default="data/processed/arabidopsis/ara_positive_embeddings_2560d.npy",
                        help="Output path for embeddings (.npy)")
    parser.add_argument("--meta_output", default="data/processed/arabidopsis/ara_positive_embeddings_meta.parquet",
                        help="Output path for aligned pair metadata (.parquet)")
    parser.add_argument("--batch_size", type=int, default=8,
                        help="Batch size for extraction (default: 8)")
    parser.add_argument("--max_pair_len", type=int, default=1020,
                        help="Maximum combined length LA + LB (default: 1020)")
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu",
                        help="Computing device (cuda or cpu)")
    parser.add_argument("--limit", type=int, default=None,
                        help="Optional limit on number of pairs (for testing)")
    return parser.parse_args()


def load_sequence_database(seq_db_path):
    print(f"Loading sequence database from {seq_db_path}...")
    if not os.path.exists(seq_db_path):
        raise FileNotFoundError(f"Sequence database not found: {seq_db_path}")
    with open(seq_db_path, "rb") as f:
        seq_db = pickle.load(f)
    print(f"Loaded {len(seq_db):,} sequence entries.")
    return seq_db


def resolve_sequence(protein_id, seq_db):
    """Resolve protein sequence with case-insensitivity, isoform stripping, and clean formatting."""
    if protein_id in seq_db:
        return seq_db[protein_id]
    
    p_up = protein_id.upper()
    if p_up in seq_db:
        return seq_db[p_up]
    
    # Strip isoform suffix (e.g., Q9LUI9-1 -> Q9LUI9)
    if "-" in protein_id:
        base_id = protein_id.split("-")[0]
        if base_id in seq_db:
            return seq_db[base_id]
        if base_id.upper() in seq_db:
            return seq_db[base_id.upper()]
            
    return None


def load_pplm_backbone(device):
    """Load the frozen PPLM backbone model from weights."""
    print("Initializing PPLM backbone model...")
    alphabet = Alphabet.from_architecture()
    batch_converter = alphabet.get_batch_converter()

    model_path = os.path.join(PPLM_DIR, "weights", "pplm_t33_650M.pt")
    if not os.path.exists(model_path):
        raise FileNotFoundError(
            f"PPLM weights not found at {model_path}. Run: cd PPLM/weights && bash download_weights.sh"
        )

    print(f"Loading model weights from {model_path}...")
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
    print(f"PPLM backbone successfully loaded on {device}.")
    return model, batch_converter


def crop_pair(seqA, seqB, max_len=1020):
    """Crop sequences proportionally if combined length exceeds max_len."""
    total_len = len(seqA) + len(seqB)
    if total_len <= max_len:
        return seqA, seqB
    
    ratio_A = len(seqA) / total_len
    budget_A = max(50, int(max_len * ratio_A))
    budget_B = max_len - budget_A
    return seqA[:budget_A], seqB[:budget_B]


def extract_batch_embeddings(model, batch_converter, pair_batch, device, max_pair_len=1020):
    """
    Extract 2,560-dim paired embeddings for a list of (seqA, seqB) tuples.
    Returns np.ndarray of shape (batch_size, 2560).
    """
    batch_embeddings = []
    
    for seqA, seqB in pair_batch:
        seqA_c, seqB_c = crop_pair(seqA, seqB, max_pair_len)
        lenA, lenB = len(seqA_c), len(seqB_c)
        
        # Tokenize pair
        _, _, seqA_tokens = batch_converter([("seqA", seqA_c)])
        _, _, seqB_tokens = batch_converter([("seqB", seqB_c)])
        tokens = torch.cat([seqA_tokens, seqB_tokens], dim=-1).to(device)
        
        # Inter-chain mask
        total_tok_len = lenA + 2 + lenB + 2
        inter_chain_mask = torch.ones((total_tok_len, total_tok_len), device=device)
        inter_chain_mask[: lenA + 2, : lenA + 2] = 0
        inter_chain_mask[lenA + 2 :, lenA + 2 :] = 0
        
        with torch.no_grad():
            out = model(
                tokens,
                inter_chain_mask,
                repr_layers=[33],
                need_head_weights=False,
                return_contacts=False
            )
            
            # Representation tensor: [1, total_tok_len, 1280]
            # Exclude special tokens (BOS, EOS) from both chains
            rep_A = out["representations"][33][0, 1:lenA + 1, :]
            rep_B = out["representations"][33][0, -(lenB + 1):-1, :]
            rep_pair = torch.cat([rep_A, rep_B], dim=0) # [lenA + lenB, 1280]
            
            # Mean and Max pooling across all paired amino acids
            mean_embed = rep_pair.mean(dim=0)
            max_embed = torch.amax(rep_pair, dim=0)
            
            paired_vec = torch.cat([mean_embed, max_embed], dim=-1).cpu().numpy().astype(np.float32)
            batch_embeddings.append(paired_vec)
            
            del out, rep_A, rep_B, rep_pair, tokens, inter_chain_mask
            
    return np.stack(batch_embeddings, axis=0)


def main():
    args = parse_args()
    print("=" * 70)
    print("PPLM PAIRED EMBEDDING EXTRACTION FOR ARABIDOPSIS POSITIVE PPIs")
    print(f"Device: {args.device} | Batch Size: {args.batch_size} | Max Length: {args.max_pair_len}")
    print("=" * 70)

    # 1. Load input dataset
    input_path = os.path.join(PROJECT_ROOT, args.input) if not os.path.isabs(args.input) else args.input
    print(f"Loading positive pairs from: {input_path}")
    if input_path.endswith(".parquet"):
        df = pd.read_parquet(input_path)
    else:
        df = pd.read_csv(input_path, sep="\t")
        
    if args.limit:
        df = df.head(args.limit).copy()
        print(f"Limiting to first {len(df):,} pairs for testing.")
        
    total_pairs = len(df)
    print(f"Total positive pairs to process: {total_pairs:,}")

    # 2. Load sequence database
    seq_db_path = os.path.join(PROJECT_ROOT, args.seq_db) if not os.path.isabs(args.seq_db) else args.seq_db
    seq_db = load_sequence_database(seq_db_path)

    # 3. Filter pairs with resolvable sequences
    print("\nResolving protein sequences for all pairs...")
    valid_indices = []
    pair_sequences = []
    missing_count = 0

    for idx, row in df.iterrows():
        p_a = str(row["participant_a_id"])
        p_b = str(row["participant_b_id"])
        
        seq_a = resolve_sequence(p_a, seq_db)
        seq_b = resolve_sequence(p_b, seq_db)
        
        if seq_a is not None and seq_b is not None:
            valid_indices.append(idx)
            pair_sequences.append((seq_a, seq_b))
        else:
            missing_count += 1

    df_valid = df.loc[valid_indices].reset_index(drop=True)
    num_valid = len(df_valid)
    print(f"Sequence resolution complete: {num_valid:,} valid pairs ({num_valid/total_pairs*100:.2f}%)")
    if missing_count > 0:
        print(f"Warning: {missing_count:,} pairs could not be resolved in the sequence database.")

    # 4. Setup output directory & memmap array for checkpointing
    output_path = os.path.join(PROJECT_ROOT, args.output) if not os.path.isabs(args.output) else args.output
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Progress file to track resumes
    progress_file = output_path + ".progress.json"
    start_idx = 0
    
    if os.path.exists(output_path) and os.path.exists(progress_file):
        with open(progress_file, "r") as f:
            prog = json.load(f)
            start_idx = prog.get("completed_pairs", 0)
            if start_idx < num_valid:
                print(f"Resuming from previous checkpoint at index {start_idx:,} / {num_valid:,}")
                embeddings_mmap = np.memmap(output_path, dtype="float32", mode="r+", shape=(num_valid, 2560))
            else:
                print(f"Extraction already completed ({start_idx:,} pairs). Exiting.")
                return
    else:
        print(f"Creating memory-mapped array for {num_valid:,} pairs: {output_path}")
        embeddings_mmap = np.memmap(output_path, dtype="float32", mode="w+", shape=(num_valid, 2560))

    # 5. Load model
    device = torch.device(args.device)
    model, batch_converter = load_pplm_backbone(device)

    # 6. Extraction Loop
    print("\nStarting batch extraction...")
    start_time = time.time()
    checkpoint_interval = 1000

    for i in range(start_idx, num_valid, args.batch_size):
        batch_end = min(i + args.batch_size, num_valid)
        curr_batch = pair_sequences[i:batch_end]
        
        try:
            vecs = extract_batch_embeddings(model, batch_converter, curr_batch, device, args.max_pair_len)
            embeddings_mmap[i:batch_end] = vecs
        except torch.cuda.OutOfMemoryError:
            print(f"\n[Warning] CUDA OOM at batch {i}:{batch_end}. Clearing cache and retrying individually with aggressive cropping...")
            torch.cuda.empty_cache()
            for sub_i, (sA, sB) in enumerate(curr_batch):
                single_vec = extract_batch_embeddings(model, batch_converter, [(sA, sB)], device, max_pair_len=800)
                embeddings_mmap[i + sub_i] = single_vec[0]

        # Checkpoint flush
        if batch_end % checkpoint_interval == 0 or batch_end == num_valid:
            embeddings_mmap.flush()
            with open(progress_file, "w") as f:
                json.dump({"completed_pairs": batch_end, "total_pairs": num_valid}, f)
                
            elapsed = time.time() - start_time
            rate = (batch_end - start_idx) / max(elapsed, 1.0)
            remaining = (num_valid - batch_end) / max(rate, 0.001)
            print(f"[{batch_end:,}/{num_valid:,} | {batch_end/num_valid*100:5.1f}%] - {rate:.1f} pairs/sec - ETA: {remaining/60:.1f} mins")

    embeddings_mmap.flush()
    print(f"\nExtraction complete in {(time.time() - start_time)/60:.2f} minutes.")

    # 7. Export aligned metadata
    meta_path = os.path.join(PROJECT_ROOT, args.meta_output) if not os.path.isabs(args.meta_output) else args.meta_output
    try:
        df_valid.to_parquet(meta_path, compression="snappy", index=False)
    except Exception as e:
        tsv_fallback = meta_path.replace(".parquet", ".tsv")
        print(f"Parquet export failed ({e}). Falling back to TSV: {tsv_fallback}")
        df_valid.to_csv(tsv_fallback, sep="\t", index=False)

    if os.path.exists(progress_file):
        os.remove(progress_file)
        
    print(f"Done! Embeddings saved to: {output_path}")
    print(f"Aligned Metadata saved to: {meta_path}")


if __name__ == "__main__":
    main()
