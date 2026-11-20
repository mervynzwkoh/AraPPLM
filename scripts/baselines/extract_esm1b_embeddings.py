#!/usr/bin/env python3
"""
Extract Single-Chain ESM-1b (esm1b_t33_650M_UR50S) Embeddings for AraPPLM Baselines.

Shared feature extraction script for both ESMAraPPI (Zhou et al., 2023) and
ARACoFusion (Sarkar & Sarkar, 2026).

Following ESMAraPPI (MLP.py & Plant Methods 2023) and ARACoFusion (Eq. 1):
  - Loads frozen ESM-1b (esm1b_t33_650M_UR50S, 33 layers, 650M params).
  - Truncates sequences > 1022 residues to 1022 (ESM-1b max positional length = 1024 tokens
    including <cls> and <eos>).
  - Extracts final layer (layer 33) residue embeddings [1 : L+1] and averages across L:
      z = (1 / L) * sum_{i=1}^L x_i  in R^{1280}
  - Saves {protein_id: torch.Tensor(1280, dtype=torch.float32)} to a single .pkl file
    directly compatible with ESMAraPPI's MLP.py (`embeddingDict[prot_id]`).
"""

import argparse
import glob
import os
import pickle
import sys
import time

import torch


DEFAULT_SPLIT_DIR = (
    "data/dataset_construction/04_benchmark_splits/arabidopsis/benchmark_data"
)
DEFAULT_SEQ_DB = "data/arabidopsis/uniprot_final.pkl"
DEFAULT_AUX_FASTA = (
    "data/dataset_construction/01_positive_pairs/arabidopsis/"
    "ara_positive_pairs_high_only_hetero.fasta"
)
DEFAULT_OUTPUT = "features/baselines/esm1b_t33_650M_UR50S_ara.pkl"


def load_fasta_dict(fasta_path):
    """Parse a FASTA file into {record_id: sequence}."""
    seqs = {}
    if not fasta_path or not os.path.exists(fasta_path):
        return seqs
    cur_id, cur_seq = None, []
    with open(fasta_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id:
                    seqs[cur_id] = "".join(cur_seq)
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line)
        if cur_id:
            seqs[cur_id] = "".join(cur_seq)
    return seqs


def collect_unique_proteins(tsv_paths):
    """Collect all unique protein IDs from the benchmark TSV files."""
    proteins = set()
    for path in tsv_paths:
        with open(path, "r", encoding="utf-8") as f:
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
                proteins.add(parts[0])
                proteins.add(parts[1])
    return sorted(proteins)


def load_esm1b_model(device, weights_path=None):
    """Load frozen esm1b_t33_650M_UR50S model and alphabet."""
    try:
        import esm
    except ImportError:
        print(
            "ERROR: fair-esm (`import esm`) is not installed. "
            "Install with `pip install fair-esm`."
        )
        sys.exit(1)

    print("Loading ESM-1b (esm1b_t33_650M_UR50S)...")
    if weights_path and os.path.exists(weights_path):
        model, alphabet = esm.pretrained.load_model_and_alphabet_local(weights_path)
    else:
        model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()

    model.to(device)
    model.eval()
    batch_converter = alphabet.get_batch_converter()
    return model, batch_converter


def main():
    parser = argparse.ArgumentParser(
        description="Extract ESM-1b 1280-d mean-pooled embeddings for AraPPLM proteins"
    )
    parser.add_argument(
        "--datasets",
        nargs="+",
        default=None,
        help="List of benchmark TSV files (defaults to all C1-C4 TSVs in benchmark_data)",
    )
    parser.add_argument(
        "--seq_db",
        default=DEFAULT_SEQ_DB,
        help="Path to pickled sequence dictionary (default: data/arabidopsis/uniprot_final.pkl)",
    )
    parser.add_argument(
        "--aux_fasta",
        default=DEFAULT_AUX_FASTA,
        help="Auxiliary FASTA file for curated positive pairs",
    )
    parser.add_argument(
        "--output",
        default=DEFAULT_OUTPUT,
        help="Output .pkl path for {protein_id: torch.Tensor(1280)}",
    )
    parser.add_argument(
        "--esm_weights",
        default=None,
        help="Optional local path to esm1b_t33_650M_UR50S.pt",
    )
    parser.add_argument(
        "--max_seq_len",
        type=int,
        default=1022,
        help="Maximum sequence length before truncation (ESM-1b limit: 1022)",
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        default=16,
        help="Number of sequences per forward batch",
    )
    parser.add_argument("--gpu_id", type=int, default=0, help="GPU device ID")
    parser.add_argument(
        "--resume",
        action="store_true",
        default=False,
        help="Resume from existing output .pkl if present",
    )
    args = parser.parse_args()

    tsv_paths = args.datasets or sorted(
        glob.glob(os.path.join(DEFAULT_SPLIT_DIR, "*.tsv"))
    )
    if not tsv_paths:
        print(f"ERROR: No TSV files found in {DEFAULT_SPLIT_DIR}")
        sys.exit(1)

    print("Scanning benchmark datasets:")
    for p in tsv_paths:
        print(f"  - {p}")
    unique_proteins = collect_unique_proteins(tsv_paths)
    print(f"Total unique proteins required: {len(unique_proteins):,}")

    # Load sequences from seq_db + aux_fasta
    seq_db = {}
    if args.seq_db and os.path.exists(args.seq_db):
        print(f"Loading sequence DB: {args.seq_db}")
        with open(args.seq_db, "rb") as f:
            seq_db.update(pickle.load(f))
        print(f"  Loaded {len(seq_db):,} sequences from seq_db")

    if args.aux_fasta and os.path.exists(args.aux_fasta):
        aux_seqs = load_fasta_dict(args.aux_fasta)
        seq_db.update(aux_seqs)
        print(
            f"  Merged {len(aux_seqs):,} sequences from {args.aux_fasta} "
            f"(total DB size: {len(seq_db):,})"
        )

    missing_prots = [p for p in unique_proteins if p not in seq_db]
    if missing_prots:
        print(
            f"WARNING: {len(missing_prots)} proteins missing from sequence DB: "
            f"{missing_prots[:10]}"
        )
    else:
        print("100.00% sequence coverage verified for all unique proteins!")

    embedding_dict = {}
    if args.resume and os.path.exists(args.output):
        with open(args.output, "rb") as f:
            embedding_dict = pickle.load(f)
        print(f"Resumed {len(embedding_dict):,} existing embeddings from {args.output}")

    to_extract = [
        (p, seq_db[p][: args.max_seq_len])
        for p in unique_proteins
        if p in seq_db and p not in embedding_dict
    ]
    print(f"Proteins remaining to extract: {len(to_extract):,}")

    if to_extract:
        device = torch.device(
            f"cuda:{args.gpu_id}" if torch.cuda.is_available() else "cpu"
        )
        print(f"Using device: {device}")
        model, batch_converter = load_esm1b_model(device, args.esm_weights)

        # Sort by length for efficient padding within batches
        to_extract.sort(key=lambda x: len(x[1]))

        start_time = time.time()
        processed = 0
        for start_idx in range(0, len(to_extract), args.batch_size):
            batch = to_extract[start_idx : start_idx + args.batch_size]
            labels, strs, toks = batch_converter(batch)
            toks = toks.to(device)

            with torch.no_grad():
                out = model(toks, repr_layers=[33], return_contacts=False)
                reps = out["representations"][33]

                for i, (prot_id, seq_str) in enumerate(batch):
                    seq_len = len(seq_str)
                    mean_emb = reps[i, 1 : seq_len + 1].mean(dim=0).detach().cpu().float()
                    embedding_dict[prot_id] = mean_emb

            processed += len(batch)
            if processed % 500 == 0 or processed == len(to_extract):
                elapsed = time.time() - start_time
                rate = processed / max(0.1, elapsed)
                print(
                    f"  Extracted {processed:,}/{len(to_extract):,} proteins "
                    f"({rate:.1f} prots/sec)"
                )

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "wb") as f:
        pickle.dump(embedding_dict, f, protocol=pickle.HIGHEST_PROTOCOL)

    size_mb = os.path.getsize(args.output) / (1024 * 1024)
    print(
        f"Saved {len(embedding_dict):,} ESM-1b embeddings to {args.output} ({size_mb:.2f} MB)"
    )


if __name__ == "__main__":
    main()
