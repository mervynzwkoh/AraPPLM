#!/usr/bin/env python3
"""
Build Protein Sequence Database for PPLM Inference
Extracts Swiss-Prot / TrEMBL sequences from FASTA, applies UniParc mappings (if provided),
and exports a single pickled sequence database dictionary with 100% ID coverage.

Usage:
    # Build complete Arabidopsis database from project root:
    python scripts/data_preparation/build_sequence_db.py \
        --fasta data/uniprot_arabidopsis.fasta \
        --uniparc_fasta data/uniparc_348.fasta \
        --mapping data/uniparc_id_mapping.json \
        --output data/uniprot_final.pkl \
        --verify data/DeepAraPPI/total_positive_negative_samples_DeepAraPPI.txt
"""

import os
import re
import json
import pickle
import argparse
from pathlib import Path
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description="Build Sequence PKL Database for PPLM")
    parser.add_argument("--fasta", required=True, help="Path to primary UniProt FASTA file")
    parser.add_argument("--uniparc_fasta", default=None, help="Optional path to UniParc auxiliary FASTA")
    parser.add_argument("--mapping", default=None, help="Optional JSON file mapping original IDs to UniParc IDs")
    parser.add_argument("--output", required=True, help="Output PKL database file path")
    parser.add_argument("--verify", default=None, help="Optional dataset TSV file to verify pair coverage against")
    return parser.parse_args()

def build_database(fasta_path, uniparc_fasta=None, mapping_path=None, output_path=None):
    print("=" * 60)
    print("BUILDING PPLM SEQUENCE DATABASE")
    print("=" * 60)
    
    seq_db = {}
    
    # 1. Parse primary UniProt FASTA
    print(f"\n1. Parsing primary FASTA: {fasta_path}")
    count = 0
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq)
        full_id = record.id
        
        # Store full ID (e.g. sp|Q9LUI9|PIAL1_ARATH)
        seq_db[full_id] = seq
        
        # Extract simple accession ID (e.g. Q9LUI9 or A0A178VKB9)
        match = re.match(r"(?:sp|tr)\|([A-Z0-9]+)\|", full_id)
        if match:
            simple_id = match.group(1)
            seq_db[simple_id] = seq
        else:
            # If no pipe format, use record.id directly
            seq_db[full_id.split()[0]] = seq
            
        count += 1

    print(f"   Parsed {count:,} records -> {len(seq_db):,} ID mappings")
    print(f"   Unique sequences: {len(set(seq_db.values())):,}")

    # 2. Merge UniParc auxiliary sequences if provided
    if uniparc_fasta and mapping_path and os.path.exists(uniparc_fasta) and os.path.exists(mapping_path):
        print(f"\n2. Merging UniParc sequences from: {uniparc_fasta}")
        with open(mapping_path, "r") as f:
            id_mapping = json.load(f)
            
        uniparc_seqs = {}
        for record in SeqIO.parse(uniparc_fasta, "fasta"):
            uniparc_seqs[record.id] = str(record.seq)
            
        merged_count = 0
        for orig_id, uniparc_id in id_mapping.items():
            if uniparc_id in uniparc_seqs:
                seq = uniparc_seqs[uniparc_id]
                seq_db[orig_id] = seq
                seq_db[uniparc_id] = seq
                merged_count += 1

        print(f"   Merged {merged_count:,} UniParc entries into database")

    # 3. Save database
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    print(f"\n3. Saving final database to: {output_path}")
    with open(output_path, "wb") as f:
        pickle.dump(seq_db, f)

    file_size_mb = os.path.getsize(output_path) / (1024 * 1024)
    print(f"   Database saved successfully ({file_size_mb:.2f} MB, {len(seq_db):,} total mappings)")

    return seq_db

def verify_coverage(seq_db, dataset_path):
    print(f"\n4. Verifying coverage against dataset: {dataset_path}")
    total_pairs = 0
    found_pairs = 0
    missing_pairs = []
    
    with open(dataset_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("Protein"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            protA, protB = parts[0], parts[1]
            total_pairs += 1
            
            if protA in seq_db and protB in seq_db:
                found_pairs += 1
            else:
                missing_pairs.append((protA, protB, protA in seq_db, protB in seq_db))

    pct = (found_pairs / total_pairs * 100) if total_pairs > 0 else 0
    print(f"   Total protein pairs:   {total_pairs:,}")
    print(f"   Pairs with BOTH seqs:  {found_pairs:,} ({pct:.2f}%)")
    print(f"   Pairs with MISSING:    {len(missing_pairs):,}")
    
    if missing_pairs:
        print("\n   Sample missing pairs:")
        for pA, pB, fA, fB in missing_pairs[:5]:
            print(f"     {pA} (found={fA}) + {pB} (found={fB})")
    else:
        print("   ✅ 100% Sequence Coverage Confirmed!")

def main():
    args = parse_args()
    seq_db = build_database(
        fasta_path=args.fasta,
        uniparc_fasta=args.uniparc_fasta,
        mapping_path=args.mapping,
        output_path=args.output
    )
    
    if args.verify and os.path.exists(args.verify):
        verify_coverage(seq_db, args.verify)

    print("\n" + "=" * 60)
    print("DONE")
    print("=" * 60)

if __name__ == "__main__":
    main()
