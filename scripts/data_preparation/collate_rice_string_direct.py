"""
Collate Direct Rice Interactions from STRING v12
=================================================
Identifies genuine experimentally supported physical interactions in Rice (Oryza sativa, TaxID 39947)
from the STRING v12 database, strictly separating directly observed experimental evidence
(experiments > 0) from orthology-transferred interologs (experiments_transferred > 0).

Exports clean provenance datasets:
1. string_direct_rice_evidence.tsv / .parquet (directed records with full raw channel scores)
2. string_direct_rice_pairs_deduplicated.tsv / .parquet (canonical undirected pairs)
3. string_direct_rice_benchmark_2col.tsv (canonical 2-column format: participant_a, participant_b)
4. string_direct_rice_audit.json (audit metadata and statistics)
"""

import os
import sys
import gzip
import json
from collections import defaultdict
from datetime import datetime, timezone
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

WORKSPACE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

def collate_string_direct_rice(output_dir=None):
    if output_dir is None:
        output_dir = os.path.join(WORKSPACE_ROOT, "data", "processed", "rice")
    os.makedirs(output_dir, exist_ok=True)

    aliases_path = os.path.join(WORKSPACE_ROOT, "data", "raw", "rice", "string", "39947.protein.aliases.v12.0.txt.gz")
    phys_links_path = os.path.join(WORKSPACE_ROOT, "data", "raw", "rice", "string", "39947.protein.physical.links.full.v12.0.txt.gz")

    if not os.path.exists(phys_links_path):
        raise FileNotFoundError(f"Missing required file: {phys_links_path}")

    print("Step 1: Loading protein aliases...")
    aliases = defaultdict(lambda: {"gene_name": None, "locus_tag": None, "description": None})
    if os.path.exists(aliases_path):
        with gzip.open(aliases_path, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    sid, alias, src = parts[0], parts[1], parts[2]
                    if src == "UniProt_GN_Name" and not aliases[sid]["gene_name"]:
                        aliases[sid]["gene_name"] = alias
                    elif src == "UniProt_GN_OrderedLocusNames" and not aliases[sid]["locus_tag"]:
                        aliases[sid]["locus_tag"] = alias
                    elif src == "UniProt_DE_RecName_Full" and not aliases[sid]["description"]:
                        aliases[sid]["description"] = alias

    print("Step 2: Parsing physical links and isolating direct experimental records...")
    evidence_records = []
    undirected_pairs_dict = {}

    total_records = 0
    exp_gt_0_records = 0
    exp_trans_only_records = 0
    neither_exp_records = 0

    with gzip.open(phys_links_path, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline().strip().split()
        for line in f:
            total_records += 1
            parts = line.strip().split()
            p1_full, p2_full = parts[0], parts[1]
            exp = int(parts[3])
            exp_tr = int(parts[4])
            db = int(parts[5])
            db_tr = int(parts[6])
            tm = int(parts[7])
            tm_tr = int(parts[8])
            comb = int(parts[9])

            if exp > 0:
                exp_gt_0_records += 1
                p1 = p1_full.replace("39947.", "")
                p2 = p2_full.replace("39947.", "")

                gene1 = aliases[p1_full]["gene_name"] or aliases[p1_full]["locus_tag"] or p1
                gene2 = aliases[p2_full]["gene_name"] or aliases[p2_full]["locus_tag"] or p2

                evidence_id = f"STRING_RICE_DIR_{exp_gt_0_records:04d}"
                rec = {
                    "evidence_id": evidence_id,
                    "source_database": "STRING v12.0",
                    "source_file": "39947.protein.physical.links.full.v12.0.txt.gz",
                    "ncbi_taxid": 39947,
                    "species_name": "Oryza sativa Japonica Group",
                    "participant_a_id": p1,
                    "participant_b_id": p2,
                    "participant_a_name": gene1,
                    "participant_b_name": gene2,
                    "participant_a_locus": aliases[p1_full]["locus_tag"] or "NA",
                    "participant_b_locus": aliases[p2_full]["locus_tag"] or "NA",
                    "interaction_type": "physical association",
                    "experiments_score": exp,
                    "experiments_transferred_score": exp_tr,
                    "database_score": db,
                    "database_transferred_score": db_tr,
                    "textmining_score": tm,
                    "textmining_transferred_score": tm_tr,
                    "combined_score": comb,
                    "is_direct_experiment": True,
                    "is_interolog_transfer": exp_tr > 0,
                    "source_pair_key": f"{p1_full}-{p2_full}"
                }
                evidence_records.append(rec)

                # Canonical undirected pair
                pair_key = tuple(sorted([p1, p2]))
                if pair_key not in undirected_pairs_dict:
                    undirected_pairs_dict[pair_key] = {
                        "participant_a_id": pair_key[0],
                        "participant_b_id": pair_key[1],
                        "participant_a_name": gene1 if p1 == pair_key[0] else gene2,
                        "participant_b_name": gene2 if p1 == pair_key[0] else gene1,
                        "participant_a_locus": (aliases[p1_full]["locus_tag"] if p1 == pair_key[0] else aliases[p2_full]["locus_tag"]) or "NA",
                        "participant_b_locus": (aliases[p2_full]["locus_tag"] if p1 == pair_key[0] else aliases[p1_full]["locus_tag"]) or "NA",
                        "experiments_score": exp,
                        "experiments_transferred_score": exp_tr,
                        "database_score": db,
                        "database_transferred_score": db_tr,
                        "textmining_score": tm,
                        "textmining_transferred_score": tm_tr,
                        "combined_score": comb,
                        "is_homomeric": pair_key[0] == pair_key[1],
                        "evidence_record_count": 1
                    }
                else:
                    undirected_pairs_dict[pair_key]["evidence_record_count"] += 1
            else:
                if exp_tr > 0:
                    exp_trans_only_records += 1
                else:
                    neither_exp_records += 1

    print(f"Total physical links scanned: {total_records}")
    print(f"Direct experiments > 0 records: {exp_gt_0_records}")
    print(f"Transferred experimental only records: {exp_trans_only_records}")
    print(f"Neither experimental records: {neither_exp_records}")
    print(f"Unique undirected direct pairs: {len(undirected_pairs_dict)}")

    # Step 3: Check overlap with existing IntAct/BioGRID curated dataset
    proc_pairs_file = os.path.join(output_dir, "rice_positive_pairs_deduplicated.tsv")
    curated_pairs = set()
    if os.path.exists(proc_pairs_file):
        df_curated = pd.read_csv(proc_pairs_file, sep="\t")
        for _, row in df_curated.iterrows():
            p_a = str(row["participant_a_id"]).strip()
            p_b = str(row["participant_b_id"]).strip()
            curated_pairs.add(tuple(sorted([p_a, p_b])))

    # Build deduplicated pairs dataframe
    dedup_rows = []
    overlap_count = 0
    novel_count = 0
    unique_proteins = set()

    for idx, (pair_key, pdata) in enumerate(sorted(undirected_pairs_dict.items()), 1):
        pair_id = f"STRING_RICE_POS_{idx:04d}"
        in_curated = pair_key in curated_pairs
        if in_curated:
            overlap_count += 1
        else:
            novel_count += 1

        unique_proteins.add(pair_key[0])
        unique_proteins.add(pair_key[1])

        dedup_rows.append({
            "pair_id": pair_id,
            "canonical_pair_key": f"{pair_key[0]}|{pair_key[1]}",
            "participant_a_id": pdata["participant_a_id"],
            "participant_b_id": pdata["participant_b_id"],
            "participant_a_name": pdata["participant_a_name"],
            "participant_b_name": pdata["participant_b_name"],
            "participant_a_locus": pdata["participant_a_locus"],
            "participant_b_locus": pdata["participant_b_locus"],
            "is_homomeric": pdata["is_homomeric"],
            "source_database": "STRING v12.0",
            "interaction_type": "physical association",
            "experiments_score": pdata["experiments_score"],
            "experiments_transferred_score": pdata["experiments_transferred_score"],
            "database_score": pdata["database_score"],
            "database_transferred_score": pdata["database_transferred_score"],
            "textmining_score": pdata["textmining_score"],
            "textmining_transferred_score": pdata["textmining_transferred_score"],
            "combined_score": pdata["combined_score"],
            "in_curated_intact_or_biogrid": in_curated,
            "evidence_record_count": pdata["evidence_record_count"]
        })

    df_evidence = pd.DataFrame(evidence_records)
    df_dedup = pd.DataFrame(dedup_rows)

    # 4. Save evidence dataset
    ev_tsv = os.path.join(output_dir, "string_direct_rice_evidence.tsv")
    ev_parquet = os.path.join(output_dir, "string_direct_rice_evidence.parquet")
    df_evidence.to_csv(ev_tsv, sep="\t", index=False)
    df_evidence.to_parquet(ev_parquet, index=False)
    print(f"Saved: {ev_tsv} ({len(df_evidence)} records)")

    # 5. Save deduplicated pairs
    dedup_tsv = os.path.join(output_dir, "string_direct_rice_pairs_deduplicated.tsv")
    dedup_parquet = os.path.join(output_dir, "string_direct_rice_pairs_deduplicated.parquet")
    df_dedup.to_csv(dedup_tsv, sep="\t", index=False)
    df_dedup.to_parquet(dedup_parquet, index=False)
    print(f"Saved: {dedup_tsv} ({len(df_dedup)} unique pairs)")

    # 6. Save benchmark 2-column pairs (participant_a, participant_b)
    bench_tsv = os.path.join(output_dir, "string_direct_rice_benchmark_2col.tsv")
    df_dedup[["participant_a_id", "participant_b_id"]].to_csv(bench_tsv, sep="\t", index=False)
    print(f"Saved: {bench_tsv}")

    # 7. Generate Audit JSON
    audit_data = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "species": "Oryza sativa Japonica Group (TaxID 39947)",
        "source": "STRING v12.0 Physical Links Full (39947.protein.physical.links.full.v12.0.txt.gz)",
        "full_network_context": {
            "all_functional_links_total": 67438350,
            "all_functional_links_experiments_gt_0": 862,
            "all_functional_links_experiments_transferred_gt_0": 61272824
        },
        "physical_links_breakdown": {
            "total_physical_records": total_records,
            "total_physical_undirected_pairs": total_records // 2,
            "direct_experiments_records": exp_gt_0_records,
            "direct_experiments_undirected_pairs": len(df_dedup),
            "direct_experiments_unique_proteins": len(unique_proteins),
            "transferred_experiments_only_records": exp_trans_only_records,
            "transferred_experiments_only_undirected_pairs": exp_trans_only_records // 2,
            "neither_experimental_records": neither_exp_records,
            "neither_experimental_undirected_pairs": neither_exp_records // 2
        },
        "overlap_with_curated_intact_biogrid": {
            "curated_total_unique_pairs": len(curated_pairs),
            "string_direct_shared_pairs": overlap_count,
            "string_direct_novel_pairs": novel_count,
            "shared_percentage": round(overlap_count / len(df_dedup) * 100, 2)
        },
        "score_distribution": {
            "experiments_score_counts": df_dedup["experiments_score"].value_counts().to_dict(),
            "combined_score_summary": {
                "min": int(df_dedup["combined_score"].min()),
                "max": int(df_dedup["combined_score"].max()),
                "mean": round(float(df_dedup["combined_score"].mean()), 2),
                "ge_700_high_confidence_pairs": int((df_dedup["combined_score"] >= 700).sum())
            }
        },
        "exported_files": [
            "string_direct_rice_evidence.tsv",
            "string_direct_rice_evidence.parquet",
            "string_direct_rice_pairs_deduplicated.tsv",
            "string_direct_rice_pairs_deduplicated.parquet",
            "string_direct_rice_benchmark_2col.tsv"
        ]
    }

    audit_json = os.path.join(output_dir, "string_direct_rice_audit.json")
    with open(audit_json, "w", encoding="utf-8") as f:
        json.dump(audit_data, f, indent=2)
    print(f"Saved: {audit_json}")

    return audit_data

if __name__ == "__main__":
    collate_string_direct_rice()
