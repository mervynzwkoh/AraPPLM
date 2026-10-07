"""
Collate Positive Sample Dataset for Arabidopsis thaliana
=========================================================
Aggregates positive PPI evidence records across:
1. IntAct (Release 252) - physical association interactions (excluding proximity)
2. BioGRID (Release 5.0.261) - physical interactions (excluding proximity & genetic)
3. PRIDE PhoX XL-MS (Trinh et al., 2026) - high-confidence cross-linking mass spectrometry

Performs canonical pair deduplication (participant_a <= participant_b),
summarizes pure database provenance (sources, detection methods, PMIDs, authors, MIscore),
and exports:
- Master Evidence Tables (.tsv and .parquet)
- Deduplicated Positive Pairs (.tsv and .parquet)
- Benchmark Datasets (3-column and 2-column, protein-only, and inter-protein)
- Comprehensive Audit JSON
"""

import os
import sys
import json
import argparse
import collections
from datetime import datetime, timezone
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

WORKSPACE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if WORKSPACE_ROOT not in sys.path:
    sys.path.insert(0, WORKSPACE_ROOT)

CANONICAL_EVIDENCE_COLUMNS = [
    "source_resource", "source_release", "source_file", "source_row_number", "source_record_id",
    "canonical_pair_key", "participant_a_id", "participant_a_namespace", "participant_a_original_id",
    "participant_a_name", "participant_a_aliases", "participant_a_taxid", "participant_a_species",
    "participant_a_biological_role", "participant_a_experimental_role", "participant_a_molecule_type",
    "participant_b_id", "participant_b_namespace", "participant_b_original_id", "participant_b_name",
    "participant_b_aliases", "participant_b_taxid", "participant_b_species", "participant_b_biological_role",
    "participant_b_experimental_role", "participant_b_molecule_type", "is_protein_protein",
    "same_species_pair", "interaction_type_raw", "interaction_type_psi_mi", "interaction_type_name",
    "detection_method_raw", "detection_method_psi_mi", "detection_method_name", "first_author",
    "publication_id_raw", "pmid", "publication_year", "throughput_raw", "source_confidence_raw",
    "intact_miscore", "xlms_protein_type", "xlms_peptide_a", "xlms_peptide_b", "xlms_linker",
    "xlms_score", "xlms_evalue"
]

CANONICAL_DEDUP_COLUMNS = [
    "pair_id", "canonical_pair_key", "participant_a_id", "participant_b_id",
    "participant_a_name", "participant_b_name", "participant_a_type", "participant_b_type",
    "is_protein_protein_pair", "is_homomeric", "pair_type", "evidence_count",
    "sources", "interaction_types", "detection_methods", "first_authors",
    "pmid_count", "pmids", "publication_years", "throughput_classes",
    "max_intact_miscore", "source_record_ids"
]


def collate_arabidopsis_positives(intermediate_dir=None, output_dir=None, mirror_dir=None):
    if intermediate_dir is None:
        intermediate_dir = os.path.join(WORKSPACE_ROOT, "data", "processed", "arabidopsis", "intermediate")
    if output_dir is None:
        output_dir = os.path.join(WORKSPACE_ROOT, "data", "processed", "arabidopsis")
    if mirror_dir is None:
        mirror_dir = os.path.join(WORKSPACE_ROOT, "data", "arabidopsis")

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(mirror_dir, exist_ok=True)

    print("=" * 70)
    print("COLLATING ARABIDOPSIS POSITIVE PPI DATASET")
    print(f"Intermediate Directory: {intermediate_dir}")
    print(f"Target Output Directory: {output_dir}")
    print("=" * 70)

    # -------------------------------------------------------------
    # 1. Load Intermediate Datasets
    # -------------------------------------------------------------
    intact_pq = os.path.join(intermediate_dir, "intact_positive_evidence.parquet")
    biogrid_pq = os.path.join(intermediate_dir, "biogrid_positive_evidence.parquet")
    xlms_pq = os.path.join(intermediate_dir, "xlms_positive_evidence.parquet")

    for fpath, name in [(intact_pq, "IntAct"), (biogrid_pq, "BioGRID"), (xlms_pq, "XL-MS")]:
        if not os.path.exists(fpath):
            raise FileNotFoundError(f"Missing intermediate file for {name}: {fpath}. Run process script first.")

    print(f"Loading IntAct intermediate parquet: {intact_pq}")
    df_intact = pq.read_table(intact_pq).to_pandas()
    print(f"  -> Loaded {len(df_intact):,} IntAct evidence rows.")

    print(f"Loading BioGRID intermediate parquet: {biogrid_pq}")
    df_biogrid = pq.read_table(biogrid_pq).to_pandas()
    print(f"  -> Loaded {len(df_biogrid):,} BioGRID evidence rows.")

    print(f"Loading XL-MS intermediate parquet: {xlms_pq}")
    df_xlms = pq.read_table(xlms_pq).to_pandas()
    print(f"  -> Loaded {len(df_xlms):,} XL-MS evidence rows.")

    # -------------------------------------------------------------
    # 2. Concatenate Master Evidence Table
    # -------------------------------------------------------------
    df_evidence = pd.concat([df_intact, df_biogrid, df_xlms], ignore_index=True)
    df_evidence = df_evidence[CANONICAL_EVIDENCE_COLUMNS]

    total_evidence_count = len(df_evidence)
    print(f"\nTotal Master Evidence Rows: {total_evidence_count:,}")

    evidence_tsv_path = os.path.join(output_dir, "ara_positive_evidence.tsv")
    evidence_parquet_path = os.path.join(output_dir, "ara_positive_evidence.parquet")

    print(f"Writing Master Evidence TSV: {evidence_tsv_path}")
    df_evidence.to_csv(evidence_tsv_path, sep="\t", index=False)

    print(f"Writing Master Evidence Parquet: {evidence_parquet_path}")
    arrow_schema = pa.schema([(c, pa.string()) for c in CANONICAL_EVIDENCE_COLUMNS])
    table_evidence = pa.Table.from_pandas(df_evidence.astype(str), schema=arrow_schema)
    pq.write_table(table_evidence, evidence_parquet_path, compression="snappy")

    # -------------------------------------------------------------
    # 3. Deduplicate into Unique Positive Interacting Pairs
    # -------------------------------------------------------------
    print("\nDeduplicating evidence into canonical pairs...")
    pair_evidence_map = collections.defaultdict(list)

    # Convert to lightweight dict records for aggregation speed
    for r in df_evidence.to_dict(orient="records"):
        pair_evidence_map[r["canonical_pair_key"]].append(r)

    sorted_pair_keys = sorted(pair_evidence_map.keys())
    total_unique_pairs = len(sorted_pair_keys)
    print(f"Total Unique Canonical Pairs: {total_unique_pairs:,}")

    dedup_rows = []
    for idx, pair_key in enumerate(sorted_pair_keys, start=1):
        ev_list = pair_evidence_map[pair_key]
        parts = pair_key.split("|")
        p_a = parts[0].strip().strip('"').strip("'")
        p_b = parts[1].strip().strip('"').strip("'") if len(parts) > 1 else p_a
        is_homo = (p_a == p_b)

        sources = sorted(list(set(str(r["source_resource"]) for r in ev_list)))
        int_types = sorted(list(set(str(r["interaction_type_name"]) for r in ev_list if str(r["interaction_type_name"]) not in ("NA", "None", ""))))
        det_methods = sorted(list(set(str(r["detection_method_name"]) for r in ev_list if str(r["detection_method_name"]) not in ("NA", "None", ""))))
        authors = sorted(list(set(str(r["first_author"]) for r in ev_list if str(r["first_author"]) not in ("NA", "None", ""))))

        pmids = sorted(list(set(
            str(r["pmid"]) for r in ev_list
            if str(r["pmid"]) not in ("NA", "None", "unassigned", "")
        )))

        years = sorted(list(set(
            str(r["publication_year"]) for r in ev_list
            if str(r["publication_year"]) not in ("NA", "None", "")
        )))

        tp = sorted(list(set(str(r["throughput_raw"]) for r in ev_list if str(r["throughput_raw"]) not in ("NA", "None", ""))))
        src_rec_ids = sorted(list(set(str(r["source_record_id"]) for r in ev_list if str(r["source_record_id"]) not in ("NA", "None", ""))))

        # Participant names & types
        names_a = set(r["participant_a_name"] for r in ev_list if r["participant_a_id"] == p_a)
        names_a.update(r["participant_b_name"] for r in ev_list if r["participant_b_id"] == p_a)
        names_a.discard("NA")

        names_b = set(r["participant_a_name"] for r in ev_list if r["participant_a_id"] == p_b)
        names_b.update(r["participant_b_name"] for r in ev_list if r["participant_b_id"] == p_b)
        names_b.discard("NA")

        types_a = set(r["participant_a_molecule_type"] for r in ev_list if r["participant_a_id"] == p_a)
        types_a.update(r["participant_b_molecule_type"] for r in ev_list if r["participant_b_id"] == p_a)
        types_b = set(r["participant_a_molecule_type"] for r in ev_list if r["participant_a_id"] == p_b)
        types_b.update(r["participant_b_molecule_type"] for r in ev_list if r["participant_b_id"] == p_b)

        type_a_str = "|".join(sorted(types_a)) if types_a else "protein"
        type_b_str = "|".join(sorted(types_b)) if types_b else "protein"

        # Max IntAct MIscore
        miscores = []
        for r in ev_list:
            v = r.get("intact_miscore")
            if v not in ("NA", "None", None, ""):
                try:
                    miscores.append(float(v))
                except ValueError:
                    pass
        max_miscore = f"{max(miscores):.4f}" if miscores else "NA"

        dedup_rows.append({
            "pair_id": f"ARA_POS_{idx:06d}",
            "canonical_pair_key": pair_key,
            "participant_a_id": p_a,
            "participant_b_id": p_b,
            "participant_a_name": ";".join(sorted(names_a)) if names_a else p_a,
            "participant_b_name": ";".join(sorted(names_b)) if names_b else p_b,
            "participant_a_type": type_a_str,
            "participant_b_type": type_b_str,
            "is_protein_protein_pair": "TRUE" if (type_a_str == "protein" and type_b_str == "protein") else "FALSE",
            "is_homomeric": "TRUE" if is_homo else "FALSE",
            "pair_type": "homomer" if is_homo else "heteromer",
            "evidence_count": len(ev_list),
            "sources": "|".join(sources),
            "interaction_types": "|".join(int_types) if int_types else "physical association",
            "detection_methods": "|".join(det_methods),
            "first_authors": "|".join(authors) if authors else "NA",
            "pmid_count": len(pmids),
            "pmids": "|".join(pmids) if pmids else "NA",
            "publication_years": "|".join(years) if years else "NA",
            "throughput_classes": "|".join(tp) if tp else "unassigned",
            "max_intact_miscore": max_miscore,
            "source_record_ids": "|".join(src_rec_ids)
        })

    df_dedup = pd.DataFrame(dedup_rows, columns=CANONICAL_DEDUP_COLUMNS)

    dedup_tsv_path = os.path.join(output_dir, "ara_positive_pairs_deduplicated.tsv")
    dedup_parquet_path = os.path.join(output_dir, "ara_positive_pairs_deduplicated.parquet")

    print(f"Writing Deduplicated Pairs TSV: {dedup_tsv_path}")
    df_dedup.to_csv(dedup_tsv_path, sep="\t", index=False)

    print(f"Writing Deduplicated Pairs Parquet: {dedup_parquet_path}")
    arrow_schema_dedup = pa.schema([(c, pa.string()) for c in CANONICAL_DEDUP_COLUMNS])
    table_dedup = pa.Table.from_pandas(df_dedup.astype(str), schema=arrow_schema_dedup)
    pq.write_table(table_dedup, dedup_parquet_path, compression="snappy")

    # -------------------------------------------------------------
    # 4. Benchmark Dataset Export (Deferred)
    # -------------------------------------------------------------
    # Benchmark dataset files are deferred pending criteria finalization.
    # Master evidence and deduplicated pairs tables are maintained.
    df_prot = df_dedup[df_dedup["is_protein_protein_pair"] == "TRUE"]
    df_hetero = df_dedup[df_dedup["pair_type"] == "heteromer"]

    # -------------------------------------------------------------
    # 5. Audit & Provenance Verification
    # -------------------------------------------------------------
    homo_pairs = df_dedup[df_dedup["pair_type"] == "homomer"]
    unique_proteins = sorted(list(set(df_dedup["participant_a_id"]).union(set(df_dedup["participant_b_id"]))))
    unique_pure_proteins = sorted(list(set(df_prot["participant_a_id"]).union(set(df_prot["participant_b_id"]))))

    all_pmids = set(
        p for p_str in df_dedup["pmids"] if p_str != "NA"
        for p in p_str.split("|")
    )

    evidence_dist = collections.Counter(df_dedup["evidence_count"])

    # Source overlap breakdown
    source_counts = collections.Counter(df_dedup["sources"])

    audit_summary = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "species": "Arabidopsis thaliana (taxid 3702)",
        "provenance_standard": "Database-Only Provenance (Internal Guidelines / Holdout Flags Stripped)",
        "filtering_rules": {
            "IntAct": "Retain interaction_semantics == 'physical_association' (exact Rice match; exclude proximity & direct_binary; no MIscore filter)",
            "BioGRID": "Exclude proximity and genetic interactions (retain direct_binary and physical_association)",
            "XL_MS": "pLink 3.2 < 1% PSM FDR; physical cross-links (<35A constraint); unique pairs retained (total PSMs preserved)"
        },
        "evidence_table_columns_count": len(df_evidence.columns) + 1,  # including record_hippie_score
        "dedup_table_columns_count": len(df_dedup.columns) + 2,        # including hippie_score & level
        "filtering_reconciliation": {
            "intact_positive_retained": len(df_intact),
            "biogrid_positive_retained": len(df_biogrid),
            "xlms_positive_retained": len(df_xlms),
            "total_positive_evidence_records": total_evidence_count
        },
        "deduplicated_pairs_summary": {
            "total_unique_pairs": len(df_dedup),
            "heteromeric_pairs": len(df_hetero),
            "homomeric_pairs": len(homo_pairs),
            "pure_protein_protein_pairs": len(df_prot),
            "non_protein_associated_pairs": len(df_dedup) - len(df_prot),
            "unique_proteins_all": len(unique_proteins),
            "unique_proteins_pure_protein_pairs": len(unique_pure_proteins),
            "source_combination_breakdown": dict(source_counts.most_common(20))
        },
        "literature_provenance": {
            "total_unique_pmids": len(all_pmids),
            "intact_pmids_count": len(set(p for p in df_intact["pmid"] if p not in ("NA", "None", ""))),
            "biogrid_pmids_count": len(set(p for p in df_biogrid["pmid"] if p not in ("NA", "None", ""))),
            "xlms_pmids_count": len(set(p for p in df_xlms["pmid"] if p not in ("NA", "None", ""))),
            "unique_authors_count": len(set(df_evidence["first_author"]))
        },
        "output_files": {
            "evidence_tsv": evidence_tsv_path,
            "evidence_parquet": evidence_parquet_path,
            "dedup_pairs_tsv": dedup_tsv_path,
            "dedup_pairs_parquet": dedup_parquet_path,
            "pair_hippie_provenance_tsv": os.path.join(output_dir, "ara_positive_hippie_provenance.tsv"),
            "pair_hippie_provenance_parquet": os.path.join(output_dir, "ara_positive_hippie_provenance.parquet"),
            "evidence_hippie_provenance_tsv": os.path.join(output_dir, "ara_positive_evidence_hippie_provenance.tsv"),
            "evidence_hippie_provenance_parquet": os.path.join(output_dir, "ara_positive_evidence_hippie_provenance.parquet")
        }
    }

    # -------------------------------------------------------------
    # 6. Compute HIPPIE Quality Scores & Intermediate Provenance
    # -------------------------------------------------------------
    from scripts.data_preparation.compute_hippie_scores import run_hippie_pipeline
    hippie_stats = run_hippie_pipeline(data_dir=output_dir, species="ara")
    audit_summary["hippie_scoring_summary"] = hippie_stats

    audit_path = os.path.join(output_dir, "ara_positive_audit.json")
    with open(audit_path, "w", encoding="utf-8") as f:
        json.dump(audit_summary, f, indent=2)

    print("\n" + "=" * 70)
    print("ARABIDOPSIS POSITIVE DATASET COLLATION COMPLETE:")
    print(f"  - Master Evidence Rows:       {total_evidence_count:,}")
    print(f"  - Deduplicated Unique Pairs:  {len(df_dedup):,}")
    print(f"    * Heteromeric Pairs:        {len(df_hetero):,}")
    print(f"    * Homomeric Pairs:          {len(homo_pairs):,}")
    print(f"    * Pure Protein Pairs:       {len(df_prot):,}")
    print(f"    * High Confidence HIPPIE:   {hippie_stats['high_confidence_count']:,} ({hippie_stats['high_confidence_count']/len(df_dedup)*100:.1f}%)")
    print(f"  - Unique Contributing PMIDs:  {len(all_pmids):,}")
    print("=" * 70)

    return audit_summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Collate Arabidopsis positive PPI dataset")
    parser.add_argument("--intermediate-dir", default=None, help="Directory with intermediate parquets")
    parser.add_argument("--output-dir", default=None, help="Directory to save final processed files")
    parser.add_argument("--mirror-dir", default=None, help="Directory to mirror benchmark files")
    args = parser.parse_args()

    collate_arabidopsis_positives(
        intermediate_dir=args.intermediate_dir,
        output_dir=args.output_dir,
        mirror_dir=args.mirror_dir
    )
