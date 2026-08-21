#!/usr/bin/env python3
"""
Verify Protein Sequence Database Coverage
Checks a protein interaction dataset against a pickled sequence database to ensure
all interacting pairs have sequences available before launching GPU benchmarks.

Usage:
    python scripts/data_preparation/verify_coverage.py \
        --dataset data/DeepAraPPI/total_positive_negative_samples_DeepAraPPI.txt \
        --seq_db data/uniprot_final.pkl
"""

import os
import sys
import pickle
import argparse
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description="Verify Sequence Database Coverage for a PPI Dataset")
    parser.add_argument("--dataset", required=True, help="Path to PPI dataset file (TSV with Protein1 <tab> Protein2 <tab> [label])")
    parser.add_argument("--seq_db", required=True, help="Path to pickled sequence database (.pkl)")
    return parser.parse_args()

def main():
    args = parse_args()
    
    if not os.path.exists(args.seq_db):
        print(f"❌ Error: Sequence database '{args.seq_db}' does not exist.")
        sys.exit(1)
        
    if not os.path.exists(args.dataset):
        print(f"❌ Error: Dataset file '{args.dataset}' does not exist.")
        sys.exit(1)
        
    print(f"Loading sequence database: {args.seq_db}")
    with open(args.seq_db, "rb") as f:
        seq_db = pickle.load(f)
    print(f"Loaded {len(seq_db):,} sequence/accession entries.")
    
    print(f"\nChecking dataset: {args.dataset}")
    total_pairs = 0
    found_pairs = 0
    missing_pairs = []
    missing_ids = set()
    
    with open(args.dataset, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("Protein"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            protA, protB = parts[0], parts[1]
            total_pairs += 1
            
            in_a = protA in seq_db
            in_b = protB in seq_db
            
            if in_a and in_b:
                found_pairs += 1
            else:
                missing_pairs.append((protA, protB, in_a, in_b))
                if not in_a:
                    missing_ids.add(protA)
                if not in_b:
                    missing_ids.add(protB)
                    
    pct = (found_pairs / total_pairs * 100) if total_pairs > 0 else 0
    print("=" * 60)
    print(f"Total pairs in dataset:    {total_pairs:,}")
    print(f"Pairs with BOTH sequences: {found_pairs:,} ({pct:.2f}%)")
    print(f"Pairs with missing seqs:   {len(missing_pairs):,}")
    print(f"Unique missing protein IDs: {len(missing_ids):,}")
    print("=" * 60)
    
    if missing_ids:
        print("\n⚠️ Missing Protein IDs (first 10):")
        for mid in sorted(list(missing_ids))[:10]:
            print(f"  - {mid}")
        sys.exit(1)
    else:
        print("\n✅ Verification SUCCESS: 100% of protein pairs are present in the sequence database!")

if __name__ == "__main__":
    main()
