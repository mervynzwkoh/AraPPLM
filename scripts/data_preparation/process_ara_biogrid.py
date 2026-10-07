"""
Process Arabidopsis BioGRID Positive PPI Dataset
================================================
Parses raw interaction data for Arabidopsis thaliana (taxid 3702) from:
BioGRID Release 5.0.261 (data/raw/arabidopsis/biogrid/BIOGRID-ORGANISM-Arabidopsis_thaliana_Columbia-5.0.261.tab3.txt)

Applies interaction type filtering (identical to Rice filter):
- Excludes proximity interaction types (PCA/BiFC, Proximity Label-MS)
- Excludes genetic interaction types (synthetic lethality, suppression)
- Retains physical interactions (direct_binary and physical_association)

Retains ONLY genuine database provenance fields (stripping internal guideline flags,
synthetic hashes, and placeholder columns).
Exports clean intermediate Parquet and TSV files.
"""

import os
import sys
import json
import argparse
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

WORKSPACE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if WORKSPACE_ROOT not in sys.path:
    sys.path.insert(0, WORKSPACE_ROOT)

from scripts.phase2.harmonize_assays import AssayHarmonizer
from scripts.phase2.parse_biogrid import BIOGRID_SOURCES, parse_biogrid_source

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


def extract_first_author(r):
    """Extract primary author citation from raw database fields."""
    try:
        raw_dict = json.loads(r.get("raw_fields_json", "{}"))
        return raw_dict.get("first_author", raw_dict.get("author", "NA"))
    except Exception:
        return "NA"


def process_biogrid(output_dir=None, fast_parquet_fallback=True):
    if output_dir is None:
        output_dir = os.path.join(WORKSPACE_ROOT, "data", "processed", "arabidopsis", "intermediate")
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 60)
    print("INGESTING ARABIDOPSIS BIOGRID DATASET")
    print("=" * 60)

    # Check for fast pre-parsed parquet in AraPPLM_Data
    parquet_path = "C:/Users/User/AraPPLM_Data/phase2/evidence_dataset/species_code=ARA/source_resource=BIOGRID"
    loaded_from_parquet = False

    raw_records = []
    if fast_parquet_fallback and os.path.exists(parquet_path):
        import glob
        parquet_files = sorted(glob.glob(os.path.join(parquet_path, "*.parquet")))
        if parquet_files:
            print(f"Found pre-parsed Parquet cache at: {parquet_path} ({len(parquet_files)} parts)")
            df_pq = pd.concat([pq.read_table(f).to_pandas() for f in parquet_files], ignore_index=True)
            raw_records = df_pq.to_dict(orient="records")
            loaded_from_parquet = True
            print(f"Loaded {len(raw_records):,} records from Parquet cache.")

    if not loaded_from_parquet:
        print("Parsing raw BioGRID TAB 3.0 source file...")
        harmonizer = AssayHarmonizer()
        ara_cfg = [s for s in BIOGRID_SOURCES if s["species_code"] == "ARA"][0]
        raw_records, rejections = parse_biogrid_source(ara_cfg, harmonizer)
        print(f"Parsed {len(raw_records):,} raw records ({len(rejections):,} rejected).")

    # Filtering BioGRID records: exclude proximity and genetic
    total_raw = len(raw_records)
    direct_bin = [r for r in raw_records if r["interaction_semantics"] == "direct_binary"]
    phys_assoc = [r for r in raw_records if r["interaction_semantics"] == "physical_association"]
    proximity = [r for r in raw_records if r["interaction_semantics"] == "proximity"]
    genetic = [r for r in raw_records if r["interaction_semantics"] == "genetic"]
    other_sem = [r for r in raw_records if r["interaction_semantics"] not in ("direct_binary", "physical_association", "proximity", "genetic")]

    print(f"BioGRID Raw Breakdown:")
    print(f"  - direct_binary:        {len(direct_bin):,} ({len(direct_bin)/total_raw*100:.2f}%)")
    print(f"  - physical_association: {len(phys_assoc):,} ({len(phys_assoc)/total_raw*100:.2f}%)")
    print(f"  - proximity (EXCLUDED): {len(proximity):,} ({len(proximity)/total_raw*100:.2f}%)")
    print(f"  - genetic (EXCLUDED):   {len(genetic):,} ({len(genetic)/total_raw*100:.2f}%)")
    if other_sem:
        print(f"  - other:                {len(other_sem):,}")

    retained = [r for r in raw_records if r["interaction_semantics"] not in ("proximity", "genetic")]
    filter_rule = "interaction_semantics not in ('proximity', 'genetic') (retain physical interactions)"

    print(f"Retained Positive Records: {len(retained):,} ({len(retained)/total_raw*100:.2f}%)")

    # Build clean database-only evidence records
    clean_records = []
    for r in retained:
        p_a = str(r["participant_a_id_value"]).strip().strip('"').strip("'")
        p_b = str(r["participant_b_id_value"]).strip().strip('"').strip("'")
        pair_key = f"{min(p_a, p_b)}|{max(p_a, p_b)}"

        author = extract_first_author(r)

        rec = {
            "source_resource": "BioGRID",
            "source_release": str(r.get("source_release", "Release 5.0.261")),
            "source_file": str(r.get("source_file", "BIOGRID-ORGANISM-Arabidopsis_thaliana_Columbia-5.0.261.tab3.txt")),
            "source_row_number": str(r.get("source_row_number", "NA")),
            "source_record_id": str(r.get("source_record_id", "NA")),
            "canonical_pair_key": pair_key,
            "participant_a_id": p_a,
            "participant_a_namespace": str(r.get("participant_a_id_namespace", "uniprotkb")),
            "participant_a_original_id": str(r.get("participant_a_original_id", f"uniprotkb:{p_a}")),
            "participant_a_name": str(r.get("participant_a_original_name", p_a)),
            "participant_a_aliases": str(r.get("participant_a_aliases_raw", "NA")),
            "participant_a_taxid": str(r.get("participant_a_taxid", 3702)),
            "participant_a_species": str(r.get("participant_a_species", "Arabidopsis thaliana (Columbia)")),
            "participant_a_biological_role": str(r.get("participant_a_biological_role", "unspecified")),
            "participant_a_experimental_role": str(r.get("participant_a_experimental_role", "unspecified")),
            "participant_a_molecule_type": "protein",
            "participant_b_id": p_b,
            "participant_b_namespace": str(r.get("participant_b_id_namespace", "uniprotkb")),
            "participant_b_original_id": str(r.get("participant_b_original_id", f"uniprotkb:{p_b}")),
            "participant_b_name": str(r.get("participant_b_original_name", p_b)),
            "participant_b_aliases": str(r.get("participant_b_aliases_raw", "NA")),
            "participant_b_taxid": str(r.get("participant_b_taxid", 3702)),
            "participant_b_species": str(r.get("participant_b_species", "Arabidopsis thaliana (Columbia)")),
            "participant_b_biological_role": str(r.get("participant_b_biological_role", "unspecified")),
            "participant_b_experimental_role": str(r.get("participant_b_experimental_role", "unspecified")),
            "participant_b_molecule_type": "protein",
            "is_protein_protein": "TRUE",
            "same_species_pair": "TRUE" if r.get("same_species_pair", True) else "FALSE",
            "interaction_type_raw": str(r.get("interaction_type_raw", "physical")),
            "interaction_type_psi_mi": str(r.get("interaction_type_psi_mi", "MI:0218")),
            "interaction_type_name": str(r.get("interaction_type_psi_mi_name", "physical")),
            "detection_method_raw": str(r.get("detection_method_raw", "NA")),
            "detection_method_psi_mi": str(r.get("detection_method_psi_mi", "NA")),
            "detection_method_name": str(r.get("detection_method_psi_mi_name", "NA")),
            "first_author": author,
            "publication_id_raw": str(r.get("publication_id_raw", "NA")),
            "pmid": str(r.get("pmid", "NA")),
            "publication_year": str(r.get("publication_year", "NA")),
            "throughput_raw": str(r.get("throughput_raw", "NA")),
            "source_confidence_raw": str(r.get("source_confidence_raw", "NA")),
            "intact_miscore": "NA",
            "xlms_protein_type": "NA",
            "xlms_peptide_a": "NA",
            "xlms_peptide_b": "NA",
            "xlms_linker": "NA",
            "xlms_score": "NA",
            "xlms_evalue": "NA"
        }
        clean_records.append(rec)

    df = pd.DataFrame(clean_records, columns=CANONICAL_EVIDENCE_COLUMNS)

    # Export intermediate files
    tsv_path = os.path.join(output_dir, "biogrid_positive_evidence.tsv")
    parquet_path_out = os.path.join(output_dir, "biogrid_positive_evidence.parquet")
    summary_path = os.path.join(output_dir, "biogrid_summary.json")

    print(f"Writing TSV to: {tsv_path}")
    df.to_csv(tsv_path, sep="\t", index=False)

    print(f"Writing Parquet to: {parquet_path_out}")
    arrow_schema = pa.schema([(c, pa.string()) for c in CANONICAL_EVIDENCE_COLUMNS])
    table = pa.Table.from_pandas(df.astype(str), schema=arrow_schema)
    pq.write_table(table, parquet_path_out, compression="snappy")

    # Compute unique pairs
    unique_pairs = set(df["canonical_pair_key"])
    hetero_pairs = set(r["canonical_pair_key"] for r in clean_records if r["participant_a_id"] != r["participant_b_id"])
    homo_pairs = unique_pairs - hetero_pairs

    summary = {
        "dataset": "BioGRID",
        "species": "Arabidopsis thaliana (taxid 3702)",
        "filter_rule": filter_rule,
        "raw_record_count": total_raw,
        "retained_evidence_count": len(df),
        "excluded_record_count": total_raw - len(df),
        "filtering_breakdown": {
            "direct_binary_retained": len(direct_bin),
            "physical_association_retained": len(phys_assoc),
            "proximity_excluded": len(proximity),
            "genetic_excluded": len(genetic)
        },
        "unique_pairs_count": len(unique_pairs),
        "heteromeric_pairs_count": len(hetero_pairs),
        "homomeric_pairs_count": len(homo_pairs),
        "unique_pmids_count": len(set(p for p in df["pmid"] if p not in ("NA", "None", ""))),
        "columns_count": len(df.columns),
        "output_files": {
            "tsv": tsv_path,
            "parquet": parquet_path_out
        }
    }

    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(f"BioGRID processing complete: {len(df):,} evidence rows, {len(unique_pairs):,} unique pairs.")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process Arabidopsis BioGRID PPI data")
    parser.add_argument("--output-dir", default=None, help="Directory to save intermediate files")
    args = parser.parse_args()

    process_biogrid(output_dir=args.output_dir)
