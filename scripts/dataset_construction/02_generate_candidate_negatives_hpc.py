#!/usr/bin/env python3
"""
02_generate_candidate_negatives_hpc.py
======================================
Generates the Candidate Negative Sample Dataset for Arabidopsis thaliana on the HPC (NSCC ASPIRE 2A).
Execution Steps:
1. Loads canonical reference proteome (Option A: 1 protein per locus, N ~ 28,609).
2. Maps TAIR Cellular Component GO annotations from tair.gaf.gz to 10 major compartments as a 10-bit uint16 mask.
3. Loads the exhaustive positive blacklist (ara_all_databases_positive_blacklist.parquet) into a 64-bit integer hash set.
4. Performs vectorized 2D block scanning across all N*(N-1)/2 pairs:
   - Rejects unannotated / ambiguous proteins.
   - Rejects pairs sharing any subcellular compartment: (mask_a & mask_b) != 0.
   - Rejects pairs recorded in the positive blacklist.
5. Streams surviving candidate negative pairs directly into chunked Parquet files (1M pairs per chunk).
6. Outputs an audit summary JSON file.

Outputs:
  data/dataset_construction/02_negative_sampling/candidates/ara_candidate_negatives_chunk_*.parquet
  data/dataset_construction/02_negative_sampling/ara_candidate_negatives_audit.json
"""

import os
import sys
import re
import csv
import gzip
import json
import time
import argparse
from collections import defaultdict, Counter
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

WORKSPACE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

DEFAULT_STREAM = os.path.join(WORKSPACE_ROOT, "data", "arabidopsis", "uniprot_arabidopsis_stream.tsv.gz")
DEFAULT_TAIR_GAF = os.path.join(WORKSPACE_ROOT, "data", "raw", "arabidopsis", "tair", "tair.gaf.gz")
DEFAULT_BLACKLIST = os.path.join(WORKSPACE_ROOT, "data", "dataset_construction", "02_negative_sampling", "ara_all_databases_positive_blacklist.parquet")
DEFAULT_OUT_DIR = os.path.join(WORKSPACE_ROOT, "data", "dataset_construction", "02_negative_sampling", "candidates")
DEFAULT_AUDIT_PATH = os.path.join(WORKSPACE_ROOT, "data", "dataset_construction", "02_negative_sampling", "ara_candidate_negatives_audit.json")

# 10 Major Cellular Compartments in Arabidopsis (Bits 0-9)
COMPARTMENTS = [
    "nucleus",                 # Bit 0
    "cytosol_cytoplasm",       # Bit 1
    "chloroplast_plastid",     # Bit 2
    "mitochondrion",           # Bit 3
    "plasma_membrane",         # Bit 4
    "endoplasmic_reticulum",   # Bit 5
    "golgi_apparatus",         # Bit 6
    "vacuole",                 # Bit 7
    "peroxisome",              # Bit 8
    "extracellular_apoplast",  # Bit 9
]

COMP_TO_BIT = {comp: (1 << i) for i, comp in enumerate(COMPARTMENTS)}

# Mapping of Cellular Component GO terms to 10 major compartments
GO_TO_COMPARTMENT = {
    # 1. Nucleus
    "GO:0005634": "nucleus", "GO:0005730": "nucleus", "GO:0005654": "nucleus", "GO:0005635": "nucleus",
    "GO:0005694": "nucleus", "GO:0016607": "nucleus", "GO:0005638": "nucleus", "GO:0001650": "nucleus",
    "GO:0005667": "nucleus", "GO:0000785": "nucleus", "GO:0000790": "nucleus", "GO:0000792": "nucleus",
    "GO:0005681": "nucleus", "GO:0005686": "nucleus", "GO:0005687": "nucleus", "GO:0005689": "nucleus",
    "GO:0005690": "nucleus", "GO:0016604": "nucleus", "GO:0031965": "nucleus", "GO:0000228": "nucleus",
    "GO:0005643": "nucleus", "GO:0034399": "nucleus", "GO:0000786": "nucleus", "GO:0016592": "nucleus",
    "GO:0000781": "nucleus", "GO:0000784": "nucleus", "GO:0005669": "nucleus", "GO:0000118": "nucleus",

    # 2. Cytoplasm / Cytosol
    "GO:0005737": "cytosol_cytoplasm", "GO:0005829": "cytosol_cytoplasm", "GO:0005856": "cytosol_cytoplasm",
    "GO:0005840": "cytosol_cytoplasm", "GO:0022625": "cytosol_cytoplasm", "GO:0022626": "cytosol_cytoplasm",
    "GO:0022627": "cytosol_cytoplasm", "GO:0005884": "cytosol_cytoplasm", "GO:0005874": "cytosol_cytoplasm",
    "GO:0005875": "cytosol_cytoplasm", "GO:0005813": "cytosol_cytoplasm", "GO:0005815": "cytosol_cytoplasm",
    "GO:0005819": "cytosol_cytoplasm", "GO:0005732": "cytosol_cytoplasm", "GO:0030529": "cytosol_cytoplasm",
    "GO:0015934": "cytosol_cytoplasm", "GO:0015935": "cytosol_cytoplasm", "GO:0043687": "cytosol_cytoplasm",
    "GO:0005844": "cytosol_cytoplasm", "GO:0010494": "cytosol_cytoplasm", "GO:0000932": "cytosol_cytoplasm",
    "GO:0019005": "cytosol_cytoplasm", "GO:0005871": "cytosol_cytoplasm", "GO:0005811": "cytosol_cytoplasm",
    "GO:0000502": "cytosol_cytoplasm", "GO:0044444": "cytosol_cytoplasm", "GO:0044445": "cytosol_cytoplasm",

    # 3. Chloroplast / Plastid
    "GO:0009507": "chloroplast_plastid", "GO:0009536": "chloroplast_plastid", "GO:0009570": "chloroplast_plastid",
    "GO:0009941": "chloroplast_plastid", "GO:0009535": "chloroplast_plastid", "GO:0009579": "chloroplast_plastid",
    "GO:0009534": "chloroplast_plastid", "GO:0009705": "chloroplast_plastid", "GO:0009706": "chloroplast_plastid",
    "GO:0009524": "chloroplast_plastid", "GO:0045271": "chloroplast_plastid", "GO:0042651": "chloroplast_plastid",
    "GO:0009532": "chloroplast_plastid", "GO:0009521": "chloroplast_plastid", "GO:0009522": "chloroplast_plastid",
    "GO:0009523": "chloroplast_plastid", "GO:0009525": "chloroplast_plastid", "GO:0009526": "chloroplast_plastid",
    "GO:0009527": "chloroplast_plastid", "GO:0009528": "chloroplast_plastid", "GO:0009529": "chloroplast_plastid",
    "GO:0009530": "chloroplast_plastid", "GO:0009531": "chloroplast_plastid", "GO:0009533": "chloroplast_plastid",
    "GO:0009538": "chloroplast_plastid", "GO:0009539": "chloroplast_plastid", "GO:0009540": "chloroplast_plastid",
    "GO:0009541": "chloroplast_plastid", "GO:0009542": "chloroplast_plastid", "GO:0009543": "chloroplast_plastid",
    "GO:0009707": "chloroplast_plastid", "GO:0010287": "chloroplast_plastid", "GO:0042644": "chloroplast_plastid",
    "GO:0009544": "chloroplast_plastid", "GO:0009545": "chloroplast_plastid", "GO:0009546": "chloroplast_plastid",

    # 4. Mitochondrion
    "GO:0005739": "mitochondrion", "GO:0005743": "mitochondrion", "GO:0005759": "mitochondrion",
    "GO:0005740": "mitochondrion", "GO:0005741": "mitochondrion", "GO:0005742": "mitochondrion",
    "GO:0005744": "mitochondrion", "GO:0005745": "mitochondrion", "GO:0005746": "mitochondrion",
    "GO:0005747": "mitochondrion", "GO:0005748": "mitochondrion", "GO:0005749": "mitochondrion",
    "GO:0005750": "mitochondrion", "GO:0005751": "mitochondrion", "GO:0005752": "mitochondrion",
    "GO:0005753": "mitochondrion", "GO:0005754": "mitochondrion", "GO:0005755": "mitochondrion",
    "GO:0005756": "mitochondrion", "GO:0005757": "mitochondrion", "GO:0005758": "mitochondrion",
    "GO:0005760": "mitochondrion", "GO:0005761": "mitochondrion", "GO:0005762": "mitochondrion",
    "GO:0005763": "mitochondrion", "GO:0031966": "mitochondrion",

    # 5. Plasma Membrane
    "GO:0005886": "plasma_membrane", "GO:0005887": "plasma_membrane", "GO:0009506": "plasma_membrane",

    # 6. Endoplasmic Reticulum
    "GO:0005783": "endoplasmic_reticulum", "GO:0005789": "endoplasmic_reticulum", "GO:0005790": "endoplasmic_reticulum",
    "GO:0005788": "endoplasmic_reticulum", "GO:0005791": "endoplasmic_reticulum", "GO:0005792": "endoplasmic_reticulum",
    "GO:0005793": "endoplasmic_reticulum", "GO:0031969": "endoplasmic_reticulum",

    # 7. Golgi Apparatus
    "GO:0005794": "golgi_apparatus", "GO:0000139": "golgi_apparatus", "GO:0005802": "golgi_apparatus",
    "GO:0005801": "golgi_apparatus",

    # 8. Vacuole
    "GO:0000325": "vacuole", "GO:0005773": "vacuole", "GO:0005774": "vacuole", "GO:0000323": "vacuole",
    "GO:0000324": "vacuole", "GO:0000326": "vacuole", "GO:0000327": "vacuole",

    # 9. Peroxisome
    "GO:0005777": "peroxisome", "GO:0005778": "peroxisome", "GO:0005779": "peroxisome",
    "GO:0005780": "peroxisome", "GO:0005781": "peroxisome", "GO:0005782": "peroxisome",

    # 10. Extracellular / Apoplast
    "GO:0005576": "extracellular_apoplast", "GO:0009505": "extracellular_apoplast",
    "GO:0048046": "extracellular_apoplast", "GO:0005618": "extracellular_apoplast",
    "GO:0005615": "extracellular_apoplast", "GO:0009514": "extracellular_apoplast"
}


def load_canonical_proteome(stream_path):
    print(f"[{time.strftime('%X')}] Loading canonical reference proteome from: {stream_path}")
    if not os.path.exists(stream_path):
        raise FileNotFoundError(f"Stream file missing: {stream_path}")

    locus_candidates = defaultdict(list)
    with gzip.open(stream_path, "rt", encoding="utf-8", errors="replace") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        idx_acc = header.index("Entry")
        idx_oln = header.index("Gene Names (ordered locus)") if "Gene Names (ordered locus)" in header else -1
        idx_gene = header.index("Gene Names") if "Gene Names" in header else -1
        idx_seq = header.index("Sequence") if "Sequence" in header else -1

        for row in reader:
            acc = row[idx_acc].strip()
            seq = row[idx_seq].strip() if idx_seq >= 0 else ""
            seq_len = len(seq)
            is_trembl = acc.startswith("A0A")

            loci = []
            if idx_oln >= 0 and row[idx_oln].strip():
                for g in row[idx_oln].strip().split():
                    g_up = g.upper()
                    if g_up.startswith("AT") and ("G" in g_up or "M" in g_up or "C" in g_up):
                        loci.append(g_up)
            if not loci and idx_gene >= 0 and row[idx_gene].strip():
                for g in row[idx_gene].strip().split():
                    g_up = g.upper()
                    if g_up.startswith("AT") and ("G" in g_up or "M" in g_up or "C" in g_up):
                        loci.append(g_up)

            if loci:
                primary_locus = loci[0]
                locus_candidates[primary_locus].append((acc, is_trembl, seq_len))

    # Canonical selection per locus
    locus_to_acc = {}
    for locus, cands in locus_candidates.items():
        sorted_cands = sorted(cands, key=lambda x: (x[1], -x[2], x[0]))
        canonical_acc = sorted_cands[0][0]
        locus_to_acc[locus] = canonical_acc

    # Canonical proteins list sorted
    sorted_loci = sorted(locus_to_acc.keys())
    canonical_proteins = [locus_to_acc[loc] for loc in sorted_loci]

    prot_to_idx = {p: np.uint32(i) for i, p in enumerate(canonical_proteins)}
    locus_to_idx = {loc: np.uint32(i) for i, loc in enumerate(sorted_loci)}
    idx_to_prot = {np.uint32(i): p for i, p in enumerate(canonical_proteins)}
    idx_to_locus = {np.uint32(i): loc for i, loc in enumerate(sorted_loci)}

    print(f"[{time.strftime('%X')}] Canonical proteome: {len(canonical_proteins):,} proteins.")
    return {
        "canonical_proteins": canonical_proteins,
        "sorted_loci": sorted_loci,
        "prot_to_idx": prot_to_idx,
        "locus_to_idx": locus_to_idx,
        "idx_to_prot": idx_to_prot,
        "idx_to_locus": idx_to_locus,
    }


def load_subcellular_masks(tair_gaf_path, proteome_info):
    print(f"\n[{time.strftime('%X')}] Assigning subcellular bitmasks from TAIR: {tair_gaf_path}")
    N = len(proteome_info["canonical_proteins"])
    masks = np.zeros(N, dtype=np.uint16)
    locus_to_idx = proteome_info["locus_to_idx"]

    annotated_loci = set()
    total_annotations = 0

    with gzip.open(tair_gaf_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!"):
                continue
            parts = line.strip().split("\t")
            if len(parts) > 8 and parts[8] == "C":
                qualifier = parts[3].strip()
                if qualifier.startswith("NOT"):
                    continue

                locus = parts[1].strip().upper()
                go_id = parts[4].strip()

                if locus in locus_to_idx and go_id in GO_TO_COMPARTMENT:
                    comp = GO_TO_COMPARTMENT[go_id]
                    bit = COMP_TO_BIT[comp]
                    idx = locus_to_idx[locus]
                    masks[idx] |= np.uint16(bit)
                    annotated_loci.add(locus)
                    total_annotations += 1

    annotated_count = np.count_nonzero(masks > 0)
    print(f"[{time.strftime('%X')}] Subcellular annotations parsed: {total_annotations:,} mapped.")
    print(f"[{time.strftime('%X')}] Proteins with confident compartment mask (>0): {annotated_count:,} / {N:,} ({annotated_count / N * 100:.2f}%)")
    return masks


def mask_to_str(mask: int) -> str:
    assigned = []
    for comp in COMPARTMENTS:
        if mask & COMP_TO_BIT[comp]:
            assigned.append(comp)
    return "|".join(assigned) if assigned else "unannotated"


def load_positive_blacklist(blacklist_path, proteome_info):
    print(f"\n[{time.strftime('%X')}] Loading positive blacklist from: {blacklist_path}")
    if not os.path.exists(blacklist_path):
        raise FileNotFoundError(f"Blacklist file not found: {blacklist_path}")

    table = pq.read_table(blacklist_path, columns=["participant_a_id", "participant_b_id"])
    col_a = table["participant_a_id"].to_pylist()
    col_b = table["participant_b_id"].to_pylist()

    prot_to_idx = proteome_info["prot_to_idx"]
    blacklist_set = set()

    mapped_count = 0
    for a, b in zip(col_a, col_b):
        if a in prot_to_idx and b in prot_to_idx:
            u = prot_to_idx[a]
            v = prot_to_idx[b]
            if u != v:
                u_min = min(u, v)
                u_max = max(u, v)
                packed = (np.uint64(u_min) << np.uint64(32)) | np.uint64(u_max)
                blacklist_set.add(packed)
                mapped_count += 1

    print(f"[{time.strftime('%X')}] Loaded {len(blacklist_set):,} canonical positive pairs into integer blacklist hash set.")
    return blacklist_set


def generate_candidate_negatives(args):
    start_time = time.time()
    os.makedirs(args.out_dir, exist_ok=True)

    proteome_info = load_canonical_proteome(args.stream_path)
    masks = load_subcellular_masks(args.tair_gaf_path, proteome_info)
    blacklist_set = load_positive_blacklist(args.blacklist_path, proteome_info)

    N = len(proteome_info["canonical_proteins"])
    total_pairs_evaluated = (N * (N - 1)) // 2
    block_size = args.block_size
    num_blocks = (N + block_size - 1) // block_size

    print(f"\n[{time.strftime('%X')}] Initiating 2D Block Streamer: N={N:,} ({total_pairs_evaluated:,} non-self combinations)")
    print(f"[{time.strftime('%X')}] Block size: {block_size} ({num_blocks}x{num_blocks} blocks) | Chunk size: {args.chunk_size:,} pairs")

    idx_to_prot = proteome_info["idx_to_prot"]
    idx_to_locus = proteome_info["idx_to_locus"]

    # Pre-compute string representations for masks to accelerate output formatting
    mask_str_cache = {int(m): mask_to_str(int(m)) for m in np.unique(masks)}

    # Buffer for streaming chunks
    buf_a = []
    buf_b = []
    buf_loc_a = []
    buf_loc_b = []
    buf_comp_a = []
    buf_comp_b = []

    chunk_idx = 1
    total_valid_negatives = 0
    total_rejected_unannotated = 0
    total_rejected_colocalized = 0
    total_rejected_blacklisted = 0

    def flush_chunk():
        nonlocal chunk_idx, buf_a, buf_b, buf_loc_a, buf_loc_b, buf_comp_a, buf_comp_b
        if not buf_a:
            return

        out_chunk_path = os.path.join(args.out_dir, f"ara_candidate_negatives_chunk_{chunk_idx:04d}.parquet")
        chunk_table = pa.Table.from_arrays(
            [
                pa.array(buf_a, type=pa.string()),
                pa.array(buf_b, type=pa.string()),
                pa.array(buf_loc_a, type=pa.string()),
                pa.array(buf_loc_b, type=pa.string()),
                pa.array(buf_comp_a, type=pa.string()),
                pa.array(buf_comp_b, type=pa.string()),
            ],
            names=["participant_a_id", "participant_b_id", "locus_a", "locus_b", "compartments_a", "compartments_b"],
        )
        pq.write_table(chunk_table, out_chunk_path, compression="snappy")
        print(f"[{time.strftime('%X')}] Written chunk {chunk_idx:04d}: {len(buf_a):,} pairs -> {out_chunk_path}")
        chunk_idx += 1
        buf_a = []
        buf_b = []
        buf_loc_a = []
        buf_loc_b = []
        buf_comp_a = []
        buf_comp_b = []

    # Iterate over 2D triangular blocks
    block_eval_count = 0
    total_blocks_to_eval = (num_blocks * (num_blocks + 1)) // 2

    for bi in range(num_blocks):
        r0 = bi * block_size
        r1 = min(r0 + block_size, N)
        row_m = masks[r0:r1]

        for bj in range(bi, num_blocks):
            c0 = bj * block_size
            c1 = min(c0 + block_size, N)
            col_m = masks[c0:c1]
            block_eval_count += 1

            # 1. Annotation check
            annotated_2d = (row_m[:, None] > 0) & (col_m[None, :] > 0)

            # 2. Disjoint subcellular check
            disjoint_2d = (row_m[:, None] & col_m[None, :]) == 0

            valid_2d = annotated_2d & disjoint_2d

            if bi == bj:
                # Diagonal block: exclude lower triangle and diagonal
                R = r1 - r0
                C = c1 - c0
                triu_mask = np.triu(np.ones((R, C), dtype=bool), k=1)
                valid_2d &= triu_mask

            ii, jj = np.where(valid_2d)
            if len(ii) == 0:
                continue

            u_global = (r0 + ii).astype(np.uint32)
            v_global = (c0 + jj).astype(np.uint32)

            # 3. Positive blacklist check
            packed_pairs = (u_global.astype(np.uint64) << np.uint64(32)) | v_global.astype(np.uint64)

            for idx in range(len(packed_pairs)):
                pk = packed_pairs[idx]
                if pk in blacklist_set:
                    total_rejected_blacklisted += 1
                else:
                    u = u_global[idx]
                    v = v_global[idx]
                    p_a = idx_to_prot[u]
                    p_b = idx_to_prot[v]
                    loc_a = idx_to_locus[u]
                    loc_b = idx_to_locus[v]
                    m_a = masks[u]
                    m_b = masks[v]

                    buf_a.append(p_a)
                    buf_b.append(p_b)
                    buf_loc_a.append(loc_a)
                    buf_loc_b.append(loc_b)
                    buf_comp_a.append(mask_str_cache[int(m_a)])
                    buf_comp_b.append(mask_str_cache[int(m_b)])
                    total_valid_negatives += 1

                    if len(buf_a) >= args.chunk_size:
                        flush_chunk()

        if (bi + 1) % 2 == 0 or (bi + 1) == num_blocks:
            pct = (block_eval_count / total_blocks_to_eval) * 100
            elapsed = time.time() - start_time
            print(f"[{time.strftime('%X')}] Progress: block {bi + 1}/{num_blocks} ({pct:.1f}%) | Valid Negatives: {total_valid_negatives:,} | Elapsed: {elapsed:.1f}s")

    # Flush remaining records
    flush_chunk()

    total_chunks = chunk_idx - 1
    total_time = time.time() - start_time
    print(f"\n[{time.strftime('%X')}] Scanning Complete!")
    print(f"[{time.strftime('%X')}] Total candidate negatives generated: {total_valid_negatives:,} in {total_chunks} chunks.")
    print(f"[{time.strftime('%X')}] Blacklist exclusions: {total_rejected_blacklisted:,}")
    print(f"[{time.strftime('%X')}] Total execution time: {total_time:.2f} seconds.")

    # Audit JSON
    audit_data = {
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "proteome_scope": "Option A (Canonical UP000006548)",
        "total_canonical_proteins": N,
        "total_pairwise_combinations": total_pairs_evaluated,
        "total_valid_candidate_negatives": total_valid_negatives,
        "blacklist_exclusions_count": total_rejected_blacklisted,
        "total_chunks_written": total_chunks,
        "chunk_size": args.chunk_size,
        "output_directory": args.out_dir,
        "execution_time_seconds": round(total_time, 2),
    }

    with open(args.audit_path, "w", encoding="utf-8") as f:
        json.dump(audit_data, f, indent=2)
    print(f"[{time.strftime('%X')}] Audit report saved: {args.audit_path}")


def main():
    parser = argparse.ArgumentParser(description="Generate Arabidopsis Candidate Negatives on HPC")
    parser.add_argument("--stream_path", default=DEFAULT_STREAM, help="Path to uniprot_arabidopsis_stream.tsv.gz")
    parser.add_argument("--tair_gaf_path", default=DEFAULT_TAIR_GAF, help="Path to tair.gaf.gz")
    parser.add_argument("--blacklist_path", default=DEFAULT_BLACKLIST, help="Path to ara_all_databases_positive_blacklist.parquet")
    parser.add_argument("--out_dir", default=DEFAULT_OUT_DIR, help="Output directory for chunks")
    parser.add_argument("--audit_path", default=DEFAULT_AUDIT_PATH, help="Output path for audit JSON")
    parser.add_argument("--chunk_size", type=int, default=1000000, help="Number of pairs per parquet chunk")
    parser.add_argument("--block_size", type=int, default=2000, help="Block size for 2D vectorized evaluation")

    args = parser.parse_args()
    generate_candidate_negatives(args)


if __name__ == "__main__":
    main()
