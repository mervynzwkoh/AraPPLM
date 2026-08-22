#!/usr/bin/env python3
"""
Merge Missing Rice UniParc Sequences into the Rice Sequence Database

Takes the missing UniParc sequences in data/rice/missing_rice.fasta,
maps them to the corresponding DeepAraPPI accessions in data/rice/missing_rice_ids.txt,
and merges them into data/rice/uniprot_rice_final.pkl.

Usage:
    python scripts/data_preparation/merge_rice_uniparc.py \
        --db data/rice/uniprot_rice_final.pkl \
        --missing_ids data/rice/missing_rice_ids.txt \
        --fasta data/rice/missing_rice.fasta \
        --verify data/DeepAraPPI/all_rice_positive_negative_DeepAraPPI.txt
"""

import os
import sys
import pickle
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Merge missing Rice UniParc sequences into PKL database")
    parser.add_argument("--db", default="data/rice/uniprot_rice_final.pkl", help="Path to rice PKL sequence database")
    parser.add_argument("--missing_ids", default="data/rice/missing_rice_ids.txt", help="Path to missing IDs text file")
    parser.add_argument("--fasta", default="data/rice/missing_rice.fasta", help="Path to missing UniParc FASTA file")
    parser.add_argument("--output", default="data/rice/uniprot_rice_final.pkl", help="Output path for final merged PKL database")
    parser.add_argument("--verify", default="data/DeepAraPPI/all_rice_positive_negative_DeepAraPPI.txt", help="Dataset to verify coverage against")
    return parser.parse_args()

def parse_fasta(fasta_path):
    """Parse FASTA into a list of (id, sequence) tuples."""
    entries = []
    current_id = None
    current_seq = []
    with open(fasta_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    entries.append((current_id, "".join(current_seq)))
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            entries.append((current_id, "".join(current_seq)))
    return entries

def main():
    args = parse_args()
    
    print("=" * 60)
    print("MERGING MISSING RICE SEQUENCES INTO DATABASE")
    print("=" * 60)
    
    # 1. Load missing IDs
    print(f"\n1. Loading missing ID list: {args.missing_ids}")
    with open(args.missing_ids, "r", encoding="utf-8") as f:
        missing_ids = [line.strip() for line in f if line.strip()]
    print(f"   Loaded {len(missing_ids):,} missing accession IDs")
    
    # 2. Load missing FASTA
    print(f"\n2. Loading UniParc FASTA: {args.fasta}")
    fasta_records = parse_fasta(args.fasta)
    print(f"   Loaded {len(fasta_records):,} UniParc sequences")
    
    # 3. Load existing database
    print(f"\n3. Loading existing sequence database: {args.db}")
    with open(args.db, "rb") as f:
        seq_db = pickle.load(f)
    print(f"   Original entries in database: {len(seq_db):,}")
    
    # 4. Merge sequences
    print("\n4. Merging sequences...")
    merged_count = 0
    for i, orig_id in enumerate(missing_ids):
        if i < len(fasta_records):
            upi_id, seq = fasta_records[i]
            seq_db[orig_id] = seq
            seq_db[upi_id] = seq
            merged_count += 1
            
    # Also add all UniParc IDs directly
    for upi_id, seq in fasta_records:
        seq_db[upi_id] = seq
        
    print(f"   Successfully merged {merged_count:,} missing accessions into database")
    print(f"   Total entries in updated database: {len(seq_db):,}")
    
    # 5. Save updated database
    print(f"\n5. Saving updated database to: {args.output}")
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "wb") as f:
        pickle.dump(seq_db, f)
    file_size_mb = os.path.getsize(args.output) / (1024 * 1024)
    print(f"   Saved successfully ({file_size_mb:.2f} MB)")
    
    # 6. Verify coverage
    if args.verify and os.path.exists(args.verify):
        print(f"\n6. Verifying coverage against dataset: {args.verify}")
        total = 0
        found = 0
        missing = []
        with open(args.verify, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("Protein"):
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                p1, p2 = parts[0].strip(), parts[1].strip()
                total += 1
                in_1 = p1 in seq_db
                in_2 = p2 in seq_db
                if in_1 and in_2:
                    found += 1
                else:
                    missing.append((p1, p2, in_1, in_2))
                    
        pct = (found / total * 100) if total > 0 else 0
        print("=" * 60)
        print(f"Total pairs checked:       {total:,}")
        print(f"Pairs with BOTH sequences: {found:,} ({pct:.2f}%)")
        print(f"Pairs with missing seqs:   {len(missing):,}")
        print("=" * 60)
        
        if missing:
            print(f"[!] Warning: {len(missing)} pairs are still missing sequences.")
        else:
            print("[SUCCESS] 100% Sequence Coverage Confirmed for Rice Benchmark!")

if __name__ == "__main__":
    main()
