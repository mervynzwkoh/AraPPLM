#!/usr/bin/env python3
"""
Build Canonical FASTA File for Arabidopsis High-Quality Heteromeric Positive Pairs
==================================================================================
Reads:
  data/dataset_construction/01_positive_pairs/arabidopsis/ara_positive_pairs_high_only_hetero.tsv
Resolves sequences for all 8,996 unique proteins from:
  1. Primary UniProt reference proteome stream (TaxID 3702)
  2. Gene locus mappings (TAIR OLN/symbol -> UniProt)
  3. External interactors (fungal, viral, human, yeast, mammalian orthologs/preys)
  4. UniSave historical entry for Q8RVQ5
  5. UniParc & ESMAraPPI auxiliary datasets
Outputs:
  data/dataset_construction/01_positive_pairs/arabidopsis/ara_positive_pairs_high_only_hetero.fasta
"""

import os
import sys
import csv
import gzip
import json
import re
import urllib.request

WORKSPACE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
TSV_PATH = os.path.join(WORKSPACE_ROOT, "data", "dataset_construction", "01_positive_pairs", "arabidopsis", "ara_positive_pairs_high_only_hetero.tsv")
STREAM_TSV_GZ = os.path.join(WORKSPACE_ROOT, "data", "arabidopsis", "uniprot_arabidopsis_stream.tsv.gz")
OUTPUT_FASTA = os.path.join(WORKSPACE_ROOT, "data", "dataset_construction", "01_positive_pairs", "arabidopsis", "ara_positive_pairs_high_only_hetero.fasta")
UNIPARC_FASTA = os.path.join(WORKSPACE_ROOT, "data", "arabidopsis", "uniparc_348.fasta")
UNIPARC_MAP = os.path.join(WORKSPACE_ROOT, "data", "arabidopsis", "uniparc_id_mapping.json")
ESM_FASTA = os.path.join(WORKSPACE_ROOT, "data", "ESMAraPPI", "missing_esmarappi.fasta")


def get_tsv_proteins(tsv_path):
    unique_prots = set()
    total_pairs = 0
    with open(tsv_path, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        for row in reader:
            total_pairs += 1
            unique_prots.add(row[0].strip())
            unique_prots.add(row[1].strip())
    print(f"Dataset TSV: {tsv_path}")
    print(f"Total pairs: {total_pairs:,} | Unique proteins: {len(unique_prots):,}")
    return unique_prots


def load_auxiliary_sources():
    aux = {}
    if os.path.exists(UNIPARC_FASTA) and os.path.exists(UNIPARC_MAP):
        with open(UNIPARC_MAP, "r", encoding="utf-8") as f:
            mapping = json.load(f)
        uniparc_seqs = {}
        cur_id = None
        cur_seq = []
        with open(UNIPARC_FASTA, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if cur_id:
                        uniparc_seqs[cur_id] = "".join(cur_seq)
                    cur_id = line[1:].split()[0]
                    cur_seq = []
                else:
                    cur_seq.append(line)
            if cur_id:
                uniparc_seqs[cur_id] = "".join(cur_seq)
        for orig_id, upi in mapping.items():
            if upi in uniparc_seqs:
                aux[orig_id] = uniparc_seqs[upi]
                aux[orig_id.upper()] = uniparc_seqs[upi]

    if os.path.exists(ESM_FASTA):
        cur_id = None
        cur_seq = []
        with open(ESM_FASTA, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if cur_id:
                        aux[cur_id] = "".join(cur_seq)
                    header = line[1:]
                    acc = header.split()[0]
                    aux[acc] = ""
                    cur_id = acc
                    cur_seq = []
                    if "secondary_acc=" in header:
                        sec = header.split("secondary_acc=")[1].split()[0]
                        aux[sec] = ""
                        cur_id = sec
                else:
                    cur_seq.append(line)
            if cur_id:
                aux[cur_id] = "".join(cur_seq)
    return aux


def parse_uniprot_stream(stream_gz_path):
    print("Loading UniProt reference stream...")
    acc_to_seq = {}
    oln_to_seq = {}
    with gzip.open(stream_gz_path, "rt", encoding="utf-8", errors="replace") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        idx_acc = header.index("Entry")
        idx_oln = header.index("Gene Names (ordered locus)") if "Gene Names (ordered locus)" in header else -1
        idx_gene = header.index("Gene Names") if "Gene Names" in header else -1
        idx_seq = header.index("Sequence")
        
        for row in reader:
            acc = row[idx_acc].strip()
            seq = row[idx_seq].strip()
            acc_to_seq[acc] = seq
            acc_to_seq[acc.upper()] = seq
            
            if idx_oln >= 0 and row[idx_oln].strip():
                for oln in row[idx_oln].strip().split():
                    oln_upper = oln.upper()
                    if oln_upper not in oln_to_seq or len(seq) > len(oln_to_seq[oln_upper]):
                        oln_to_seq[oln_upper] = seq
                        
            if idx_gene >= 0 and row[idx_gene].strip():
                for g in row[idx_gene].strip().split():
                    g_upper = g.upper()
                    if g_upper.startswith("AT") and "G" in g_upper:
                        if g_upper not in oln_to_seq or len(seq) > len(oln_to_seq[g_upper]):
                            oln_to_seq[g_upper] = seq
    return acc_to_seq, oln_to_seq


def fetch_direct_uniprot(accession):
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=15) as resp:
        content = resp.read().decode("utf-8")
        lines = content.strip().splitlines()
        seq = "".join(l.strip() for l in lines[1:])
        return seq


def fetch_unisave(accession):
    url = f"https://rest.uniprot.org/unisave/{accession}?format=fasta"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=15) as resp:
        content = resp.read().decode("utf-8")
        lines = content.strip().splitlines()
        seq_lines = []
        first_record = False
        for line in lines:
            if line.startswith(">"):
                if first_record:
                    break
                first_record = True
            elif first_record:
                seq_lines.append(line.strip())
        seq = "".join(seq_lines)
        return seq


def build_fasta():
    target_prots = get_tsv_proteins(TSV_PATH)
    aux = load_auxiliary_sources()
    acc_to_seq, oln_to_seq = parse_uniprot_stream(STREAM_TSV_GZ)
    
    # Specific known mappings
    specific_mappings = {
        "AT1G56070": "Q9ASR1",     # LOS1 / EF2
        "AT1G67370": "F4HRV8",     # ASY1
        "EBI-15192079": "Q9FL05",  # AT5G24470
    }
    
    resolved = {}
    missing = []
    
    for p in target_prots:
        seq = None
        p_up = p.upper()
        
        # 1. Direct accession in stream
        if p in acc_to_seq:
            seq = acc_to_seq[p]
        elif p_up in acc_to_seq:
            seq = acc_to_seq[p_up]
        # 2. Ordered locus name in stream
        elif p_up in oln_to_seq:
            seq = oln_to_seq[p_up]
        # 3. Specific mapping
        elif p in specific_mappings:
            mapped = specific_mappings[p]
            if mapped in acc_to_seq:
                seq = acc_to_seq[mapped]
            else:
                seq = fetch_direct_uniprot(mapped)
        # 4. Aux (UniParc / ESMAraPPI)
        elif p in aux:
            seq = aux[p]
        elif p_up in aux:
            seq = aux[p_up]
        # 5. Isoform fallback (e.g. Q9LUI9-1 -> Q9LUI9)
        elif "-" in p:
            base = p.split("-")[0]
            if base in acc_to_seq:
                seq = acc_to_seq[base]
            elif base in aux:
                seq = aux[base]
        
        if seq:
            resolved[p] = seq
        else:
            missing.append(p)
            
    print(f"Stream + Aux resolution: {len(resolved):,} / {len(target_prots):,} (Missing: {len(missing)})")
    
    # Fetch any remaining from UniProt REST / UniSave
    if missing:
        print(f"Fetching remaining {len(missing)} entries from UniProt/UniSave...")
        for m in missing:
            try:
                if m == "Q8RVQ5":
                    seq = fetch_unisave(m)
                else:
                    seq = fetch_direct_uniprot(m)
                if seq:
                    resolved[m] = seq
                    print(f"  Fetched {m}: {len(seq)} aa")
            except Exception as e:
                print(f"  [ERROR] Failed to fetch {m}: {e}")
                
    # Final assertion
    unresolved = [p for p in target_prots if p not in resolved]
    if unresolved:
        print(f"[FATAL] Could not resolve {len(unresolved)} proteins: {unresolved}")
        sys.exit(1)
        
    print(f"\n[SUCCESS] 100.00% Coverage achieved: All {len(resolved):,} proteins resolved!")
    
    # Write to FASTA
    os.makedirs(os.path.dirname(OUTPUT_FASTA), exist_ok=True)
    print(f"\nWriting fasta to: {OUTPUT_FASTA}")
    
    # Keep deterministic order matching unique proteins sorted
    with open(OUTPUT_FASTA, "w", encoding="utf-8") as out:
        for p in sorted(target_prots):
            seq = resolved[p]
            out.write(f">{p}\n")
            # Wrap at 60 characters
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")
                
    file_size_mb = os.path.getsize(OUTPUT_FASTA) / (1024 * 1024)
    print(f"FASTA written successfully ({file_size_mb:.2f} MB)")
    
    # Verification check
    print("\n" + "=" * 60)
    print("RUNNING INDEPENDENT VERIFICATION CHECK")
    print("=" * 60)
    fasta_proteins = set()
    fasta_records = 0
    with open(OUTPUT_FASTA, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                fasta_records += 1
                fasta_proteins.add(line[1:].strip())
                
    print(f"1. Total pairs in dataset TSV: {len(target_prots):,} unique proteins represented")
    print(f"2. Total FASTA records: {fasta_records:,}")
    print(f"3. Unique proteins in FASTA: {len(fasta_proteins):,}")
    assert len(fasta_proteins) == len(target_prots), f"Count mismatch: {len(fasta_proteins)} vs {len(target_prots)}"
    assert fasta_proteins == target_prots, "Set mismatch between TSV and FASTA!"
    print("4. Exact 1-to-1 set equality check: PASSED (100% matched)")
    print("=" * 60)


if __name__ == "__main__":
    build_fasta()
