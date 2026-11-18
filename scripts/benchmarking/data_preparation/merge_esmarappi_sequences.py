#!/usr/bin/env python3
"""
Merge Missing ESMAraPPI FASTA Sequences into the Arabidopsis Sequence Database

Usage:
    python scripts/data_preparation/merge_esmarappi_sequences.py \
        --db data/arabidopsis/uniprot_final.pkl \
        --fasta data/ESMAraPPI/missing_esmarappi.fasta \
        --output data/arabidopsis/uniprot_final.pkl
"""

import os
import sys
import pickle
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Merge missing ESMAraPPI FASTA sequences into Arabidopsis PKL database")
    parser.add_argument("--db", default="data/arabidopsis/uniprot_final.pkl", help="Path to existing Arabidopsis PKL database")
    parser.add_argument("--fasta", default="data/ESMAraPPI/missing_esmarappi.fasta", help="Path to downloaded FASTA file containing missing sequences")
    parser.add_argument("--output", default="data/arabidopsis/uniprot_final.pkl", help="Output path for updated PKL database")
    return parser.parse_args()

def parse_fasta(fasta_path):
    """Parse FASTA into a dictionary of {accession: sequence}."""
    entries = {}
    cur_id = None
    cur_seq = []
    with open(fasta_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if cur_id:
                    entries[cur_id] = "".join(cur_seq)
                header = line[1:].strip()
                # Parse bare accession from compound header >sp|P25069|... or >P25069
                if "|" in header:
                    parts = header.split("|")
                    cur_id = parts[1] if len(parts) >= 2 else header.split()[0]
                else:
                    cur_id = header.split()[0]
                cur_seq = []
            else:
                cur_seq.append(line)
        if cur_id:
            entries[cur_id] = "".join(cur_seq)
    return entries

def main():
    args = parse_args()
    
    print("=" * 60)
    print("MERGING MISSING ESMARAPPI SEQUENCES")
    print("=" * 60)
    
    if not os.path.exists(args.fasta):
        print(f"[!] Error: FASTA file not found at: {args.fasta}")
        print("    Please save your downloaded FASTA file to that path.")
        sys.exit(1)
        
    print(f"\n1. Loading FASTA sequences: {args.fasta}")
    fasta_seqs = parse_fasta(args.fasta)
    print(f"   Loaded {len(fasta_seqs)} sequences from FASTA")
    for acc, seq in fasta_seqs.items():
        print(f"     - {acc}: {len(seq)} amino acids")
        
    print(f"\n2. Loading existing sequence database: {args.db}")
    with open(args.db, "rb") as f:
        seq_db = pickle.load(f)
    print(f"   Original entries in database: {len(seq_db):,}")
    
    print("\n3. Merging sequences into database...")
    for acc, seq in fasta_seqs.items():
        seq_db[acc] = seq
        
    print(f"   Updated total entries: {len(seq_db):,}")
    
    print(f"\n4. Saving updated database to: {args.output}")
    with open(args.output, "wb") as f:
        pickle.dump(seq_db, f)
    print("[SUCCESS] Database updated successfully!")

if __name__ == "__main__":
    main()
