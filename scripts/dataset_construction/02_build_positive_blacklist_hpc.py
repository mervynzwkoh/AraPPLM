#!/usr/bin/env python3
"""
02_build_positive_blacklist_hpc.py
==================================
Builds an exhaustive, zero-filter Positive Interaction Blacklist for Arabidopsis thaliana
on the HPC (NSCC ASPIRE 2A) by aggregating all positive interaction evidence across:
1. BioGRID (v5.0.261): All records (physical, genetic, proximity-MS, colocalization).
2. IntAct (Release 252 / IMEx): All records (includes MINT and DIP curation).
3. PRIDE / pLink XL-MS (2026): All cross-linked peptide pairs without FDR cutoff.
4. STRING (v12.0 Physical Full): All physical links (experiments, databases, textmining, transfers).
5. STRING (v12.0 Links Detailed): All links with experimental > 0 or database > 0.
6. Legacy Benchmarks: DeepAraPPI and ESMAraPPI positive interactions.

All interactors are resolved to Canonical UniProtKB accessions (Option A: 1 protein per locus).
Outputs:
  data/dataset_construction/02_negative_sampling/ara_all_databases_positive_blacklist.parquet
  data/dataset_construction/02_negative_sampling/ara_blacklist_audit.json
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
DEFAULT_BIOGRID = os.path.join(WORKSPACE_ROOT, "data", "raw", "arabidopsis", "biogrid", "BIOGRID-ORGANISM-Arabidopsis_thaliana_Columbia-5.0.261.tab3.txt")
DEFAULT_INTACT = os.path.join(WORKSPACE_ROOT, "data", "raw", "arabidopsis", "intact", "intact_arabidopsis_taxid3702_release252.mitab27.txt")
DEFAULT_XLMS = os.path.join(WORKSPACE_ROOT, "data", "raw", "arabidopsis", "xlms_2026", "Total_XL_plink3-2_v3.csv")
DEFAULT_STRING_PHYS = os.path.join(WORKSPACE_ROOT, "data", "raw", "arabidopsis", "string", "3702.protein.physical.links.full.v12.0.txt.gz")
DEFAULT_STRING_DET = os.path.join(WORKSPACE_ROOT, "data", "raw", "arabidopsis", "string", "3702.protein.links.detailed.v12.0.txt.gz")
DEFAULT_DEEPARAPPI = os.path.join(WORKSPACE_ROOT, "data", "raw", "legacy_benchmarks", "deeparappi", "total_positive_samples_DeepAraPPI.txt")
DEFAULT_ESMARAPPI_DIR = os.path.join(WORKSPACE_ROOT, "data", "raw", "legacy_benchmarks", "esmarappi")

DEFAULT_OUT_DIR = os.path.join(WORKSPACE_ROOT, "data", "dataset_construction", "02_negative_sampling")


def parse_uniprot_mappings(stream_path):
    """
    Parses uniprot_arabidopsis_stream.tsv.gz and establishes:
    - 1-to-1 mapping from AGI locus code to Canonical UniProt Accession.
    - Mapping from any UniProt accession to Canonical Accession for that locus.
    - Index mappings for ultra-compact 64-bit pair hashing.
    """
    print(f"[{time.strftime('%X')}] Loading UniProt reference stream: {stream_path}")
    if not os.path.exists(stream_path):
        raise FileNotFoundError(f"UniProt stream not found: {stream_path}")

    # locus -> list of (acc, is_trembl, seq_len)
    locus_candidates = defaultdict(list)
    acc_to_locus = {}

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
                acc_to_locus[acc] = primary_locus
                acc_to_locus[acc.upper()] = primary_locus

    # Choose canonical accession for each locus:
    # Rule: Swiss-Prot (is_trembl == False) over TrEMBL, then longest sequence, then lexicographic
    locus_to_canonical = {}
    for locus, cands in locus_candidates.items():
        sorted_cands = sorted(cands, key=lambda x: (x[1], -x[2], x[0]))
        canonical_acc = sorted_cands[0][0]
        locus_to_canonical[locus] = canonical_acc

    # Also map all known accessions to the canonical accession of their locus
    acc_to_canonical = {}
    for acc, locus in acc_to_locus.items():
        if locus in locus_to_canonical:
            acc_to_canonical[acc] = locus_to_canonical[locus]

    # Sorted list of all canonical proteins
    canonical_proteins = sorted(list(set(locus_to_canonical.values())))
    prot_to_idx = {p: np.uint32(i) for i, p in enumerate(canonical_proteins)}
    idx_to_prot = {np.uint32(i): p for i, p in enumerate(canonical_proteins)}
    canonical_to_locus = {v: k for k, v in locus_to_canonical.items()}

    print(f"[{time.strftime('%X')}] Loaded {len(canonical_proteins):,} canonical proteins across {len(locus_to_canonical):,} AGI loci.")
    return {
        "locus_to_canonical": locus_to_canonical,
        "acc_to_canonical": acc_to_canonical,
        "canonical_to_locus": canonical_to_locus,
        "canonical_proteins": canonical_proteins,
        "prot_to_idx": prot_to_idx,
        "idx_to_prot": idx_to_prot,
    }


def make_normalizer(mappings):
    """
    Returns a fast normalizer that resolves raw strings to canonical integer indices.
    """
    locus_map = mappings["locus_to_canonical"]
    acc_map = mappings["acc_to_canonical"]
    prot_to_idx = mappings["prot_to_idx"]

    # Pre-compile regex cleaners
    re_clean = re.compile(r"^(uniprotkb:|uniprot:|intact:|refseq:|ncbi:)", re.IGNORECASE)
    re_strip_string = re.compile(r"^3702\.")

    def resolve(raw_id):
        if not raw_id:
            return None
        s = raw_id.strip()
        s = re_clean.sub("", s)
        s = re_strip_string.sub("", s)

        # Check if locus code (AT...G... or AT...M... or AT...C...)
        s_up = s.upper()
        base_locus = s_up.split(".")[0].split("-")[0]
        if base_locus in locus_map:
            can_acc = locus_map[base_locus]
            return prot_to_idx.get(can_acc)

        # Check if direct accession (or isoform)
        base_acc = s_up.split("-")[0].split(".")[0]
        if base_acc in acc_map:
            can_acc = acc_map[base_acc]
            return prot_to_idx.get(can_acc)
        if s in acc_map:
            can_acc = acc_map[s]
            return prot_to_idx.get(can_acc)

        return None

    return resolve


def pack_pair(u: np.uint32, v: np.uint32) -> np.uint64:
    """Packs two 32-bit indices into a single 64-bit integer where u < v."""
    if u < v:
        return (np.uint64(u) << np.uint64(32)) | np.uint64(v)
    elif v < u:
        return (np.uint64(v) << np.uint64(32)) | np.uint64(u)
    return np.uint64(0)  # Self-interaction


def parse_biogrid(biogrid_path, resolve_fn):
    print(f"\n[{time.strftime('%X')}] Parsing BioGRID: {biogrid_path}")
    pairs = set()
    total_rows = 0
    if not os.path.exists(biogrid_path):
        print(f"WARNING: BioGRID file missing at {biogrid_path}, skipping.")
        return pairs

    with open(biogrid_path, "r", encoding="utf-8", errors="replace") as f:
        header = f.readline()
        for line in f:
            total_rows += 1
            cols = line.strip().split("\t")
            if len(cols) > 16:
                tax_a = cols[15].strip()
                tax_b = cols[16].strip()
                if tax_a == "3702" and tax_b == "3702":
                    sp_a = cols[23].strip() if len(cols) > 23 else ""
                    sp_b = cols[26].strip() if len(cols) > 26 else ""
                    sys_a = cols[5].strip() if len(cols) > 5 else ""
                    sys_b = cols[6].strip() if len(cols) > 6 else ""

                    raw_a = sp_a.split("|")[0] if (sp_a and sp_a != "-") else sys_a.split("|")[0]
                    raw_b = sp_b.split("|")[0] if (sp_b and sp_b != "-") else sys_b.split("|")[0]

                    idx_a = resolve_fn(raw_a)
                    idx_b = resolve_fn(raw_b)
                    if idx_a is not None and idx_b is not None:
                        packed = pack_pair(idx_a, idx_b)
                        if packed != 0:
                            pairs.add(packed)

    print(f"[{time.strftime('%X')}] BioGRID: {total_rows:,} rows -> {len(pairs):,} unique canonical pairs.")
    return pairs


def parse_intact(intact_path, resolve_fn):
    print(f"\n[{time.strftime('%X')}] Parsing IntAct / IMEx: {intact_path}")
    pairs = set()
    total_rows = 0
    if not os.path.exists(intact_path):
        print(f"WARNING: IntAct file missing at {intact_path}, skipping.")
        return pairs

    with open(intact_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            total_rows += 1
            cols = line.strip().split("\t")
            if len(cols) > 10:
                tax_a = cols[9]
                tax_b = cols[10]
                if "3702" in tax_a and "3702" in tax_b:
                    raw_a = cols[0].split(":")[-1].strip()
                    raw_b = cols[1].split(":")[-1].strip()
                    idx_a = resolve_fn(raw_a)
                    idx_b = resolve_fn(raw_b)
                    if idx_a is not None and idx_b is not None:
                        packed = pack_pair(idx_a, idx_b)
                        if packed != 0:
                            pairs.add(packed)

    print(f"[{time.strftime('%X')}] IntAct: {total_rows:,} rows -> {len(pairs):,} unique canonical pairs.")
    return pairs


def parse_xlms(xlms_path, resolve_fn):
    print(f"\n[{time.strftime('%X')}] Parsing PRIDE / pLink XL-MS: {xlms_path}")
    pairs = set()
    total_rows = 0
    if not os.path.exists(xlms_path):
        print(f"WARNING: XL-MS file missing at {xlms_path}, skipping.")
        return pairs

    with open(xlms_path, "r", encoding="utf-8", errors="replace") as f:
        reader = csv.reader(f)
        header = next(reader)
        p1_idx = 28
        p2_idx = 30
        for i, col in enumerate(header):
            if "first_valid_protein_pair" in col.lower():
                p1_idx = i
            elif col.lower() == "proteinb":
                p2_idx = i

        for row in reader:
            total_rows += 1
            if len(row) > max(p1_idx, p2_idx):
                raw_a = row[p1_idx].strip()
                raw_b = row[p2_idx].strip()
                idx_a = resolve_fn(raw_a)
                idx_b = resolve_fn(raw_b)
                if idx_a is not None and idx_b is not None:
                    packed = pack_pair(idx_a, idx_b)
                    if packed != 0:
                        pairs.add(packed)

    print(f"[{time.strftime('%X')}] XL-MS: {total_rows:,} rows -> {len(pairs):,} unique canonical pairs.")
    return pairs


def parse_string_physical(string_phys_path, resolve_fn):
    print(f"\n[{time.strftime('%X')}] Parsing STRING physical links full: {string_phys_path}")
    pairs = set()
    total_rows = 0
    if not os.path.exists(string_phys_path):
        print(f"WARNING: STRING physical file missing at {string_phys_path}, skipping.")
        return pairs

    with gzip.open(string_phys_path, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline()
        for line in f:
            total_rows += 1
            parts = line.strip().split()
            if len(parts) >= 2:
                idx_a = resolve_fn(parts[0])
                idx_b = resolve_fn(parts[1])
                if idx_a is not None and idx_b is not None:
                    packed = pack_pair(idx_a, idx_b)
                    if packed != 0:
                        pairs.add(packed)

    print(f"[{time.strftime('%X')}] STRING Physical: {total_rows:,} rows -> {len(pairs):,} unique canonical pairs.")
    return pairs


def parse_string_detailed(string_det_path, resolve_fn):
    print(f"\n[{time.strftime('%X')}] Parsing STRING detailed links (experimental > 0 or database > 0): {string_det_path}")
    pairs = set()
    total_rows = 0
    matched_rows = 0
    if not os.path.exists(string_det_path):
        print(f"WARNING: STRING detailed file missing at {string_det_path}, skipping.")
        return pairs

    with gzip.open(string_det_path, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline().strip().split()
        idx_exp = header.index("experimental") if "experimental" in header else -1
        idx_db = header.index("database") if "database" in header else -1

        for line in f:
            total_rows += 1
            parts = line.strip().split()
            if len(parts) > max(idx_exp, idx_db):
                exp_score = int(parts[idx_exp]) if idx_exp >= 0 else 0
                db_score = int(parts[idx_db]) if idx_db >= 0 else 0
                if exp_score > 0 or db_score > 0:
                    matched_rows += 1
                    idx_a = resolve_fn(parts[0])
                    idx_b = resolve_fn(parts[1])
                    if idx_a is not None and idx_b is not None:
                        packed = pack_pair(idx_a, idx_b)
                        if packed != 0:
                            pairs.add(packed)
            if total_rows % 2500000 == 0:
                print(f"  Processed {total_rows:,} / 12,557,494 rows... ({len(pairs):,} unique pairs)")

    print(f"[{time.strftime('%X')}] STRING Detailed: {total_rows:,} total rows ({matched_rows:,} exp/db > 0) -> {len(pairs):,} unique canonical pairs.")
    return pairs


def parse_legacy_benchmarks(deeparappi_path, esmarappi_dir, resolve_fn):
    print(f"\n[{time.strftime('%X')}] Parsing Legacy Benchmarks (DeepAraPPI & ESMAraPPI)...")
    pairs = set()

    # DeepAraPPI
    if os.path.exists(deeparappi_path):
        with open(deeparappi_path, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("#") or line.startswith("Protein"):
                    continue
                parts = line.strip().split()
                if len(parts) >= 2:
                    idx_a = resolve_fn(parts[0])
                    idx_b = resolve_fn(parts[1])
                    if idx_a is not None and idx_b is not None:
                        packed = pack_pair(idx_a, idx_b)
                        if packed != 0:
                            pairs.add(packed)
        print(f"[{time.strftime('%X')}] DeepAraPPI loaded.")

    # ESMAraPPI
    if os.path.exists(esmarappi_dir):
        for fname in ["c1Train.txt", "c2Pred.txt", "c3Pred.txt"]:
            fpath = os.path.join(esmarappi_dir, fname)
            if os.path.exists(fpath):
                with open(fpath, "r", encoding="utf-8") as f:
                    for line in f:
                        parts = line.strip().split()
                        if len(parts) >= 3 and parts[2] == "1":
                            idx_a = resolve_fn(parts[0])
                            idx_b = resolve_fn(parts[1])
                            if idx_a is not None and idx_b is not None:
                                packed = pack_pair(idx_a, idx_b)
                                if packed != 0:
                                    pairs.add(packed)
        print(f"[{time.strftime('%X')}] ESMAraPPI loaded.")

    print(f"[{time.strftime('%X')}] Total Legacy Benchmarks: {len(pairs):,} unique canonical pairs.")
    return pairs


def build_unified_blacklist(args):
    start_time = time.time()
    os.makedirs(args.out_dir, exist_ok=True)

    mappings = parse_uniprot_mappings(args.stream_path)
    resolve_fn = make_normalizer(mappings)

    # Harvest from each source
    source_pairs = {}
    source_pairs["BioGRID"] = parse_biogrid(args.biogrid_path, resolve_fn)
    source_pairs["IntAct"] = parse_intact(args.intact_path, resolve_fn)
    source_pairs["XL-MS"] = parse_xlms(args.xlms_path, resolve_fn)
    source_pairs["STRING_physical"] = parse_string_physical(args.string_phys_path, resolve_fn)
    source_pairs["STRING_detailed_exp_db"] = parse_string_detailed(args.string_det_path, resolve_fn)
    source_pairs["Legacy_benchmarks"] = parse_legacy_benchmarks(args.deeparappi_path, args.esmarappi_dir, resolve_fn)

    print(f"\n[{time.strftime('%X')}] Consolidating and annotating unified blacklist...")
    # pair_packed -> list of source names
    pair_sources = defaultdict(list)
    for src_name, p_set in source_pairs.items():
        for p in p_set:
            pair_sources[p].append(src_name)

    total_unique_pairs = len(pair_sources)
    print(f"[{time.strftime('%X')}] Total Unique Positive Interactions Across All Databases: {total_unique_pairs:,}")

    # Build Parquet output
    print(f"[{time.strftime('%X')}] Unpacking pairs to DataFrame...")
    idx_to_prot = mappings["idx_to_prot"]
    can_to_locus = mappings["canonical_to_locus"]

    packed_list = list(pair_sources.keys())
    part_a_list = []
    part_b_list = []
    locus_a_list = []
    locus_b_list = []
    src_list = []

    for packed in packed_list:
        u = np.uint32(packed >> np.uint64(32))
        v = np.uint32(packed & np.uint64(0xFFFFFFFF))
        p_a = idx_to_prot[u]
        p_b = idx_to_prot[v]
        loc_a = can_to_locus.get(p_a, "NA")
        loc_b = can_to_locus.get(p_b, "NA")
        src_str = "|".join(pair_sources[packed])

        part_a_list.append(p_a)
        part_b_list.append(p_b)
        locus_a_list.append(loc_a)
        locus_b_list.append(loc_b)
        src_list.append(src_str)

    table_blacklist = pa.Table.from_arrays(
        [
            pa.array(part_a_list, type=pa.string()),
            pa.array(part_b_list, type=pa.string()),
            pa.array(locus_a_list, type=pa.string()),
            pa.array(locus_b_list, type=pa.string()),
            pa.array(src_list, type=pa.string()),
        ],
        names=["participant_a_id", "participant_b_id", "locus_a", "locus_b", "sources"],
    )

    out_parquet = os.path.join(args.out_dir, "ara_all_databases_positive_blacklist.parquet")
    print(f"[{time.strftime('%X')}] Saving Parquet: {out_parquet}")
    pq.write_table(table_blacklist, out_parquet, compression="snappy")

    # Generate Audit Summary JSON
    source_counts = {k: len(v) for k, v in source_pairs.items()}
    overlap_counts = Counter(len(srcs) for srcs in pair_sources.values())

    audit_data = {
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "total_canonical_proteins": len(mappings["canonical_proteins"]),
        "total_unique_blacklisted_pairs": total_unique_pairs,
        "source_breakdown": source_counts,
        "support_count_distribution": {f"{k}_sources": v for k, v in sorted(overlap_counts.items())},
        "output_file": out_parquet,
        "execution_time_seconds": round(time.time() - start_time, 2),
    }

    out_audit = os.path.join(args.out_dir, "ara_blacklist_audit.json")
    with open(out_audit, "w", encoding="utf-8") as f:
        json.dump(audit_data, f, indent=2)
    print(f"[{time.strftime('%X')}] Audit report saved: {out_audit}")
    print(f"[{time.strftime('%X')}] Blacklist compilation complete in {time.time() - start_time:.2f} seconds.")


def main():
    parser = argparse.ArgumentParser(description="Build Comprehensive Zero-Filter Arabidopsis Positive PPI Blacklist on HPC")
    parser.add_argument("--stream_path", default=DEFAULT_STREAM, help="Path to uniprot_arabidopsis_stream.tsv.gz")
    parser.add_argument("--biogrid_path", default=DEFAULT_BIOGRID, help="Path to BioGRID tab3 file")
    parser.add_argument("--intact_path", default=DEFAULT_INTACT, help="Path to IntAct mitab file")
    parser.add_argument("--xlms_path", default=DEFAULT_XLMS, help="Path to XL-MS CSV file")
    parser.add_argument("--string_phys_path", default=DEFAULT_STRING_PHYS, help="Path to STRING physical links full")
    parser.add_argument("--string_det_path", default=DEFAULT_STRING_DET, help="Path to STRING links detailed")
    parser.add_argument("--deeparappi_path", default=DEFAULT_DEEPARAPPI, help="Path to DeepAraPPI positive file")
    parser.add_argument("--esmarappi_dir", default=DEFAULT_ESMARAPPI_DIR, help="Directory containing ESMAraPPI splits")
    parser.add_argument("--out_dir", default=DEFAULT_OUT_DIR, help="Output directory")

    args = parser.parse_args()
    build_unified_blacklist(args)


if __name__ == "__main__":
    main()
