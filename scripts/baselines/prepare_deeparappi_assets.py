#!/usr/bin/env python3
"""
Prepare DeepAraPPI Assets (word2vec 32-d AA embeddings, Domain2vec 128-d, GO2vec 128-d)
for Training and Evaluating DeepAraPPI on the Curated AraPPLM Dataset.

Reads from the cloned DeepAraPPI repository (`external_repos/DeepAraPPI` or `--repo_dir`):
  - `embeddings/aa_embeddings.txt` (20 standard amino acids -> 32-d float vectors)
  - `embeddings/ara_domain2vec_embeddings.txt` (FASTA-like 2-line format -> 128-d float vectors)
  - `embeddings/ara_go2vec_embeddings.txt` (FASTA-like 2-line format -> 128-d float vectors)

Resolves all 8,426 unique AraPPLM proteins by checking:
  1. Exact participant ID (e.g. UniProt accession `P43273`)
  2. Isoform-stripped base UniProt accession (e.g. `Q94CK9-3` -> `Q94CK9`)
  3. TAIR AGI locus ID from the `locus_a` / `locus_b` column (e.g. `AT5G06950`)
  4. Zero-vector fallback (`np.zeros(128, dtype=np.float32)`) for out-of-corpus proteins,
     while logging exact corpus coverage statistics per split (C1, C2, C3, C4).
"""

import argparse
import glob
import os
import pickle
import sys
import numpy as np


DEFAULT_SPLIT_DIR = (
    "data/dataset_construction/04_benchmark_splits/arabidopsis/benchmark_data"
)
DEFAULT_SEQ_DB = "data/arabidopsis/uniprot_final.pkl"
DEFAULT_AUX_FASTA = (
    "data/dataset_construction/01_positive_pairs/arabidopsis/"
    "ara_positive_pairs_high_only_hetero.fasta"
)
DEFAULT_REPO_DIR = "external_repos/DeepAraPPI"
DEFAULT_OUTPUT = "features/baselines/deeparappi_assets.pkl"

STANDARD_AAS = "ACDEFGHIKLMNPQRSTVWY"


def load_fasta_dict(fasta_path):
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


def parse_aa_embeddings(aa_emb_path):
    """
    Parse `embeddings/aa_embeddings.txt` matching `rcnn.py::get_pro_embedding`:
        line = line.strip().split('\\t')
        aa = line[0]
        vec = np.array([float(x) for x in line[1].split()])
    """
    aa_vec = {}
    with open(aa_emb_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            aa = parts[0].strip()
            if len(parts) == 2:
                vals = [float(x) for x in parts[1].split()]
            else:
                vals = [float(x) for x in parts[1:]]
            aa_vec[aa] = np.asarray(vals, dtype=np.float32)
    return aa_vec


def build_fallback_word2vec(sequences, dim=32, window=3):
    """Fallback CBOW word2vec (dim=32, window=3) if aa_embeddings.txt is not present."""
    try:
        from gensim.models import Word2Vec

        print("Training fallback Gensim CBOW word2vec (dim=32, window=3)...")
        sentences = [list(seq) for seq in sequences.values() if seq]
        w2v = Word2Vec(
            sentences=sentences,
            vector_size=dim,
            window=window,
            min_count=1,
            sg=0,
            seed=63,
            workers=4,
        )
        return {
            aa: np.asarray(w2v.wv[aa], dtype=np.float32)
            for aa in w2v.wv.key_to_index
        }
    except ImportError:
        print(
            "WARNING: gensim not installed and aa_embeddings.txt missing; "
            "using deterministic orthogonal 32-d AA basis."
        )
        rng = np.random.RandomState(63)
        return {
            aa: rng.normal(0, 1, size=(dim,)).astype(np.float32)
            for aa in STANDARD_AAS
        }


def parse_node2vec_embeddings(embedding_file, expected_dim=128):
    """
    Parse `ara_domain2vec_embeddings.txt` or `ara_go2vec_embeddings.txt`
    matching `Domain2vec.py::get_embeddings` and `go2vec.py::get_embeddings`.
    """
    all_pro_vec = {}
    if not embedding_file or not os.path.exists(embedding_file):
        return all_pro_vec

    pro = None
    with open(embedding_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                pro = line[1:].strip().split()[0]
            elif pro is not None:
                vec = [float(x) for x in line.split()]
                if len(vec) == expected_dim:
                    all_pro_vec[pro] = np.asarray(vec, dtype=np.float32)
                pro = None
    return all_pro_vec


def collect_proteins_and_loci(tsv_paths):
    """Collect unique proteins, their TAIR locus mappings, and per-split pair lists."""
    prot_to_locus = {}
    split_pairs = {}
    for path in tsv_paths:
        split_name = os.path.basename(path)
        pairs = []
        with open(path, "r", encoding="utf-8") as f:
            for i, line in enumerate(f):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 3:
                    parts = line.split()
                if len(parts) < 3:
                    continue
                if i == 0 and (
                    parts[0] in ("Protein1", "participant_a_id", "protein1")
                    or parts[2] == "label"
                ):
                    continue
                p_a, p_b = parts[0], parts[1]
                lbl = int(parts[2])
                pairs.append((p_a, p_b, lbl))
                prot_to_locus.setdefault(p_a, None)
                prot_to_locus.setdefault(p_b, None)
                if len(parts) >= 6:
                    loc_a, loc_b = parts[4], parts[5]
                    if loc_a and loc_a != "None":
                        prot_to_locus[p_a] = loc_a
                    if loc_b and loc_b != "None":
                        prot_to_locus[p_b] = loc_b
        split_pairs[split_name] = pairs
    return prot_to_locus, split_pairs


def resolve_vector(prot_id, locus_id, raw_dict, dim=128):
    """Lookup protein vector by exact ID, isoform-stripped ID, or TAIR locus ID."""
    if prot_id in raw_dict:
        return raw_dict[prot_id], True
    base_id = prot_id.split("-")[0]
    if base_id in raw_dict:
        return raw_dict[base_id], True
    if locus_id and locus_id in raw_dict:
        return raw_dict[locus_id], True
    return np.zeros(dim, dtype=np.float32), False


def find_file_in_repo(repo_dir, filename):
    """Search for filename under repo_dir."""
    direct = os.path.join(repo_dir, "embeddings", filename)
    if os.path.exists(direct):
        return direct
    matches = glob.glob(os.path.join(repo_dir, "**", filename), recursive=True)
    return matches[0] if matches else None


def main():
    parser = argparse.ArgumentParser(
        description="Prepare DeepAraPPI embeddings and sequence lookup for AraPPLM"
    )
    parser.add_argument(
        "--repo_dir",
        default=DEFAULT_REPO_DIR,
        help="Path to cloned zjy1125/DeepAraPPI repository",
    )
    parser.add_argument(
        "--datasets",
        nargs="+",
        default=None,
        help="Benchmark TSV files (defaults to all C1-C4 TSVs in benchmark_data)",
    )
    parser.add_argument(
        "--seq_db",
        default=DEFAULT_SEQ_DB,
        help="Path to pickled sequence dictionary",
    )
    parser.add_argument(
        "--aux_fasta",
        default=DEFAULT_AUX_FASTA,
        help="Path to ara_positive_pairs_high_only_hetero.fasta",
    )
    parser.add_argument(
        "--output",
        default=DEFAULT_OUTPUT,
        help="Output .pkl bundle path",
    )
    args = parser.parse_args()

    tsv_paths = args.datasets or sorted(
        glob.glob(os.path.join(DEFAULT_SPLIT_DIR, "*.tsv"))
    )
    if not tsv_paths:
        print(f"ERROR: No TSV files found in {DEFAULT_SPLIT_DIR}")
        sys.exit(1)

    prot_to_locus, split_pairs = collect_proteins_and_loci(tsv_paths)
    unique_prots = sorted(prot_to_locus.keys())
    print(f"Unique proteins across splits: {len(unique_prots):,}")

    # Load sequences
    seq_db = {}
    if args.seq_db and os.path.exists(args.seq_db):
        with open(args.seq_db, "rb") as f:
            seq_db.update(pickle.load(f))
    if args.aux_fasta and os.path.exists(args.aux_fasta):
        seq_db.update(load_fasta_dict(args.aux_fasta))

    sequences = {}
    missing_seq = 0
    for p in unique_prots:
        if p in seq_db:
            sequences[p] = seq_db[p]
        else:
            sequences[p] = ""
            missing_seq += 1
    print(
        f"Sequence coverage: {len(unique_prots) - missing_seq:,}/{len(unique_prots):,} "
        f"({(len(unique_prots) - missing_seq) / max(1, len(unique_prots)) * 100:.2f}%)"
    )

    # 1. Load aa_embeddings.txt (32-d)
    aa_path = find_file_in_repo(args.repo_dir, "aa_embeddings.txt")
    if aa_path and os.path.exists(aa_path):
        print(f"Loading word2vec AA embeddings from: {aa_path}")
        aa_vec = parse_aa_embeddings(aa_path)
        print(f"  Loaded {len(aa_vec)} amino acid vectors (dim={len(next(iter(aa_vec.values())))})")
    else:
        print(f"WARNING: aa_embeddings.txt not found under {args.repo_dir}")
        aa_vec = build_fallback_word2vec(sequences, dim=32, window=3)

    # 2. Load ara_domain2vec_embeddings.txt (128-d)
    d2v_path = find_file_in_repo(args.repo_dir, "ara_domain2vec_embeddings.txt")
    raw_d2v = parse_node2vec_embeddings(d2v_path, expected_dim=128)
    print(f"Loaded raw Domain2vec corpus ({d2v_path}): {len(raw_d2v):,} entries")

    # 3. Load ara_go2vec_embeddings.txt (128-d)
    g2v_path = find_file_in_repo(args.repo_dir, "ara_go2vec_embeddings.txt")
    raw_g2v = parse_node2vec_embeddings(g2v_path, expected_dim=128)
    print(f"Loaded raw GO2vec corpus ({g2v_path}): {len(raw_g2v):,} entries")

    domain2vec = {}
    go2vec = {}
    d2v_found = set()
    g2v_found = set()

    for p in unique_prots:
        loc = prot_to_locus.get(p)
        d_vec, d_ok = resolve_vector(p, loc, raw_d2v, dim=128)
        g_vec, g_ok = resolve_vector(p, loc, raw_g2v, dim=128)
        domain2vec[p] = d_vec
        go2vec[p] = g_vec
        if d_ok:
            d2v_found.add(p)
        if g_ok:
            g2v_found.add(p)

    print(
        f"Unique Protein Corpus Coverage:\n"
        f"  - RCNN (word2vec sequence): {len(unique_prots) - missing_seq:,}/{len(unique_prots):,} "
        f"({(len(unique_prots) - missing_seq) / max(1, len(unique_prots)) * 100:.2f}%)\n"
        f"  - Domain2vec (128-d):       {len(d2v_found):,}/{len(unique_prots):,} "
        f"({len(d2v_found) / max(1, len(unique_prots)) * 100:.2f}%)\n"
        f"  - GO2vec (128-d):           {len(g2v_found):,}/{len(unique_prots):,} "
        f"({len(g2v_found) / max(1, len(unique_prots)) * 100:.2f}%)"
    )

    coverage_stats = {
        "unique_proteins": len(unique_prots),
        "rcnn_covered_proteins": len(unique_prots) - missing_seq,
        "domain2vec_covered_proteins": len(d2v_found),
        "go2vec_covered_proteins": len(g2v_found),
        "splits": {},
    }

    print("\nPair-Level Corpus Coverage by Split (both interactors in corpus):")
    for split_name, pairs in split_pairs.items():
        tot = len(pairs)
        d_both = sum(1 for a, b, _ in pairs if a in d2v_found and b in d2v_found)
        g_both = sum(1 for a, b, _ in pairs if a in g2v_found and b in g2v_found)
        coverage_stats["splits"][split_name] = {
            "total_pairs": tot,
            "domain2vec_both_covered": d_both,
            "go2vec_both_covered": g_both,
        }
        print(
            f"  {split_name}: total={tot:,} | "
            f"Domain2vec both={d_both:,} ({d_both / max(1, tot) * 100:.2f}%) | "
            f"GO2vec both={g_both:,} ({g_both / max(1, tot) * 100:.2f}%)"
        )

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    bundle = {
        "aa_vec": aa_vec,
        "sequences": sequences,
        "domain2vec": domain2vec,
        "go2vec": go2vec,
        "coverage_stats": coverage_stats,
    }
    with open(args.output, "wb") as f:
        pickle.dump(bundle, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"\nSaved DeepAraPPI assets bundle to: {args.output}")


if __name__ == "__main__":
    main()
