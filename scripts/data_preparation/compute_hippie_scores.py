"""
Compute HIPPIE Scores for Rice Positive PPI Interactome
======================================================
Implements the HIPPIE (Human Integrated Protein-Protein Interaction rEference)
scoring framework (Schaefer et al., 2012, PLoS ONE) for the Rice (Oryza sativa)
interactome.

Mathematical Formulation:
-------------------------
S = w_s * s_s(n_s) + w_t * s_t(n_t) + w_o * s_o(n_o)

where subscores use the saturating sigmoid function:
s_i(n) = (2 / (1 + exp(-a_i * n))) - 1

Optimized Parameters:
- Study weight w_s = 0.6, steepness a_s = 2.3
- Technique weight w_t = 0.3, steepness a_t = 0.2
- Orthology weight w_o = 0.1, steepness a_o = 1.6

Inputs:
- n_s: Number of independent studies (distinct PMIDs)
- n_t: Sum of technique quality weights (0 to 10 scale)
- n_o: Cross-species orthology reproducibility count

Outputs:
- Updates rice_positive_evidence.tsv/.parquet with `record_hippie_score`
- Updates rice_positive_pairs_deduplicated.tsv/.parquet with `hippie_score` and `hippie_confidence_level`
- Generates dedicated provenance table: rice_positive_hippie_provenance.tsv/.parquet (logging all intermediate variables)
- Generates record-level provenance table: rice_positive_evidence_hippie_provenance.tsv/.parquet
- Exports high-confidence benchmark pairs (HIPPIE >= 0.72)
"""

import os
import sys
import re
import math
import json
from collections import defaultdict
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

WORKSPACE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

# -------------------------------------------------------------
# HIPPIE Experimental Technique Weights (0 - 10 scale)
# Based on Schaefer et al. (2012) PSI-MI mapping
# -------------------------------------------------------------
TECHNIQUE_WEIGHT_MAP = {
    # Highest Confidence (Direct Structural / Biophysical Reconstitution = 10.0)
    "reconstituted complex": 10.0,
    "protein-peptide": 10.0,
    "co-crystal structure": 10.0,
    "x-ray crystallography": 10.0,
    "surface plasmon resonance": 10.0,
    "nuclear magnetic resonance": 10.0,
    
    # Near-Physiological / Direct Biophysical / Cross-linking (9.0 - 8.0)
    "far western": 9.0,
    "cross-linking": 8.0,
    "cross-link": 8.0,
    "xl-ms": 8.0,
    "electron microscopy": 8.0,
    
    # Biochemical / Enzymatic (7.5)
    "biochemical activity": 7.5,
    "biochemical": 7.5,
    "enzymatic study": 7.5,
    "protein kinase assay": 7.5,
    "in-gel kinase assay": 7.5,
    "protease assay": 7.5,
    "gtpase assay": 7.5,
    "phosphatase assay": 7.5,
    "phosphotransferase assay": 7.5,
    
    # Biophysical Proximity / Energy Transfer (6.0)
    "fret": 6.0,
    "bret": 6.0,
    "classical fluorescence spectroscopy": 6.0,
    "competition binding": 6.0,
    "enzyme linked immunosorbent assay": 6.0,
    
    # Standard Binary & Affinity Platforms (5.0)
    "two hybrid": 5.0,
    "two-hybrid": 5.0,
    "affinity capture-ms": 5.0,
    "affinity chromatography technology": 5.0,
    "affinity technology": 5.0,
    "anti bait coimmunoprecipitation": 5.0,
    "anti tag coimmunoprecipitation": 5.0,
    "affinity capture-western": 5.0,
    "co-purification": 5.0,
    "affinity capture-luminescence": 5.0,
    "solid phase assay": 5.0,
    "tandem affinity purification": 5.0,
    "phage display": 5.0,
    "ubiquitin reconstruction": 5.0,
    "three hybrid": 5.0,
    
    # Low Confidence / Proximity / Complementation (2.5 - 2.0)
    "pull down": 2.5,
    "bimolecular fluorescence complementation": 2.5,
    "split firefly luciferase": 2.5,
    "split renilla luciferase": 2.5,
    "pca": 2.0,
    
    # Co-fractionation / Molecular sieving / Hydrodynamic (1.0)
    "co-fractionation": 1.0,
    "molecular sieving": 1.0,
    "blue native page": 1.0,
    "comigration": 1.0,
    "cosedimentation": 1.0,
}

# HIPPIE Hyperparameters
W_S = 0.6  # Weight for studies
W_T = 0.3  # Weight for experimental techniques
W_O = 0.1  # Weight for orthology

A_S = 2.3  # Saturation constant for studies
A_T = 0.2  # Saturation constant for techniques
A_O = 1.6  # Saturation constant for orthology

HIGH_CONFIDENCE_THRESHOLD = 0.72


def map_technique_weight(detection_method_name, detection_method_raw, source_resource):
    """Maps an assay detection method to HIPPIE technique quality score (0 to 10)."""
    method_lower = str(detection_method_name).strip().lower()
    raw_lower = str(detection_method_raw).strip().lower()

    if source_resource == "STRING":
        m = re.search(r"score\s+(\d+)", raw_lower)
        if m:
            sc = int(m.group(1))
            if sc >= 900:
                return 10.0  # Crystallographic / structural complex (PDB)
            elif sc >= 800:
                return 8.0   # High-confidence verified binary
            elif sc >= 400:
                return 6.0   # Multi-assay or intermediate validation
            else:
                return 5.0   # Standard screen (AP-MS 292, Y2H 230)
        return 5.0

    for term, weight in TECHNIQUE_WEIGHT_MAP.items():
        if term in method_lower or term in raw_lower:
            return weight

    return 5.0  # Default standard experimental assay


def calc_saturating_subscore(n, a):
    """Sigmoidal saturating subscore: s_i(n) = (2 / (1 + exp(-a * n))) - 1"""
    if n <= 0:
        return 0.0
    return (2.0 / (1.0 + math.exp(-a * float(n)))) - 1.0


def compute_hippie_score(n_s, n_t, n_o):
    """Computes overall HIPPIE score and returns (score, s_s, s_t, s_o)."""
    s_s = calc_saturating_subscore(n_s, A_S)
    s_t = calc_saturating_subscore(n_t, A_T)
    s_o = calc_saturating_subscore(n_o, A_O)
    score = (W_S * s_s) + (W_T * s_t) + (W_O * s_o)
    return round(score, 4), round(s_s, 4), round(s_t, 4), round(s_o, 4)


def run_hippie_pipeline(data_dir=None, species="rice"):
    prefix = "ara" if str(species).lower() in ("ara", "arabidopsis") else "rice"
    if data_dir is None:
        sub_folder = "arabidopsis" if prefix == "ara" else "rice"
        data_dir = os.path.join(WORKSPACE_ROOT, "data", "processed", sub_folder)

    evidence_tsv = os.path.join(data_dir, f"{prefix}_positive_evidence.tsv")
    evidence_parquet = os.path.join(data_dir, f"{prefix}_positive_evidence.parquet")
    dedup_tsv = os.path.join(data_dir, f"{prefix}_positive_pairs_deduplicated.tsv")
    dedup_parquet = os.path.join(data_dir, f"{prefix}_positive_pairs_deduplicated.parquet")

    if not os.path.exists(evidence_tsv) or not os.path.exists(dedup_tsv):
        raise FileNotFoundError(f"Missing master files in {data_dir}")

    print("Step 1: Loading master evidence and deduplicated datasets...")
    df_evidence = pd.read_csv(evidence_tsv, sep="\t", dtype=str)
    df_dedup = pd.read_csv(dedup_tsv, sep="\t", dtype=str)

    # ---------------------------------------------------------
    # 2. Record-Level HIPPIE Computation
    # ---------------------------------------------------------
    print("Step 2: Computing record-level HIPPIE scores and intermediate variables...")
    evidence_provenance_rows = []
    record_hippie_scores = []
    evidence_records = df_evidence.to_dict(orient="records")

    for r in evidence_records:
        source_res = r["source_resource"]
        method_name = r["detection_method_name"]
        method_raw = r["detection_method_raw"]
        pmid = str(r["pmid"]).strip()

        tech_weight = map_technique_weight(method_name, method_raw, source_res)
        has_pmid = 1 if pmid not in ("NA", "None", "", "unassigned") else 0
        
        # Check cross-species orthology transfer
        conf_raw = str(r["source_confidence_raw"])
        has_orthology = 0
        if "exp_transferred:" in conf_raw:
            m = re.search(r"exp_transferred:(\d+)", conf_raw)
            if m and int(m.group(1)) > 0:
                has_orthology = 1

        rec_score, rec_s_s, rec_s_t, rec_s_o = compute_hippie_score(
            n_s=has_pmid,
            n_t=tech_weight,
            n_o=has_orthology
        )
        record_hippie_scores.append(f"{rec_score:.4f}")

        evidence_provenance_rows.append({
            "source_record_id": r["source_record_id"],
            "canonical_pair_key": r["canonical_pair_key"],
            "source_resource": source_res,
            "detection_method_name": method_name,
            "technique_weight": tech_weight,
            "n_studies_record": has_pmid,
            "s_studies_record": rec_s_s,
            "n_techniques_record": tech_weight,
            "s_techniques_record": rec_s_t,
            "n_orthology_record": has_orthology,
            "s_orthology_record": rec_s_o,
            "record_hippie_score": rec_score
        })

    # Update df_evidence
    df_evidence["record_hippie_score"] = record_hippie_scores

    # ---------------------------------------------------------
    # 3. Pair-Level Cumulative HIPPIE Computation
    # ---------------------------------------------------------
    print("Step 3: Computing pair-level cumulative HIPPIE scores and intermediate variables...")
    # Group evidence records by canonical_pair_key
    pair_evidence_groups = defaultdict(list)
    for r in evidence_records:
        pair_evidence_groups[r["canonical_pair_key"]].append(r)

    pair_provenance_rows = []
    pair_hippie_scores = []
    pair_confidence_levels = []
    dedup_records = df_dedup.to_dict(orient="records")

    for r in dedup_records:
        pair_key = r["canonical_pair_key"]
        ev_list = pair_evidence_groups[pair_key]

        # 1. Distinct PMIDs
        pmids = set()
        methods_list = []
        weights_list = []
        has_orthology = 0

        for ev in ev_list:
            pmid = str(ev["pmid"]).strip()
            if pmid not in ("NA", "None", "", "unassigned"):
                pmids.add(pmid)

            w = map_technique_weight(ev["detection_method_name"], ev["detection_method_raw"], ev["source_resource"])
            methods_list.append(str(ev["detection_method_name"]))
            weights_list.append(str(w))

            conf_raw = str(ev["source_confidence_raw"])
            if "exp_transferred:" in conf_raw:
                m = re.search(r"exp_transferred:(\d+)", conf_raw)
                if m and int(m.group(1)) > 0:
                    has_orthology = 1

        n_s = len(pmids) if len(pmids) > 0 else 1
        n_t = sum(float(w) for w in weights_list)
        n_o = has_orthology

        pair_score, s_s, s_t, s_o = compute_hippie_score(n_s, n_t, n_o)
        
        conf_level = "high" if pair_score >= HIGH_CONFIDENCE_THRESHOLD else "medium"
        is_high = "TRUE" if pair_score >= HIGH_CONFIDENCE_THRESHOLD else "FALSE"

        pair_hippie_scores.append(f"{pair_score:.4f}")
        pair_confidence_levels.append(conf_level)

        pair_provenance_rows.append({
            "pair_id": r["pair_id"],
            "canonical_pair_key": pair_key,
            "participant_a_id": r["participant_a_id"],
            "participant_b_id": r["participant_b_id"],
            "participant_a_name": r["participant_a_name"],
            "participant_b_name": r["participant_b_name"],
            "sources": r["sources"],
            "evidence_count": r["evidence_count"],
            "n_studies": n_s,
            "s_studies": s_s,
            "detection_methods": "|".join(methods_list),
            "technique_weights": "|".join(weights_list),
            "n_techniques": round(n_t, 2),
            "s_techniques": s_t,
            "n_orthology": n_o,
            "s_orthology": s_o,
            "weight_studies": W_S,
            "weight_techniques": W_T,
            "weight_orthology": W_O,
            "hippie_score": pair_score,
            "hippie_confidence_level": conf_level,
            "is_high_confidence_hippie": is_high
        })

    # Update df_dedup
    df_dedup["hippie_score"] = pair_hippie_scores
    df_dedup["hippie_confidence_level"] = pair_confidence_levels

    # ---------------------------------------------------------
    # 4. Save Updated Master Tables & Dedicated Provenance Tables
    # ---------------------------------------------------------
    print("Step 4: Saving updated master tables and dedicated HIPPIE provenance tables...")
    # Master Evidence Table
    df_evidence.to_csv(evidence_tsv, sep="\t", index=False)
    arrow_schema_ev = pa.schema([(c, pa.string()) for c in df_evidence.columns])
    pq.write_table(pa.Table.from_pandas(df_evidence.astype(str), schema=arrow_schema_ev), evidence_parquet, compression="snappy")
    print(f"Updated master evidence table: {evidence_tsv} ({len(df_evidence)} records, {len(df_evidence.columns)} columns)")

    # Master Deduplicated Table
    df_dedup.to_csv(dedup_tsv, sep="\t", index=False)
    arrow_schema_dedup = pa.schema([(c, pa.string()) for c in df_dedup.columns])
    pq.write_table(pa.Table.from_pandas(df_dedup.astype(str), schema=arrow_schema_dedup), dedup_parquet, compression="snappy")
    print(f"Updated master deduplicated table: {dedup_tsv} ({len(df_dedup)} pairs, {len(df_dedup.columns)} columns)")

    # Dedicated Pair-Level HIPPIE Provenance Table
    df_pair_prov = pd.DataFrame(pair_provenance_rows)
    pair_prov_tsv = os.path.join(data_dir, f"{prefix}_positive_hippie_provenance.tsv")
    pair_prov_parquet = os.path.join(data_dir, f"{prefix}_positive_hippie_provenance.parquet")
    df_pair_prov.to_csv(pair_prov_tsv, sep="\t", index=False)
    df_pair_prov.to_parquet(pair_prov_parquet, index=False)
    print(f"Saved pair HIPPIE provenance table: {pair_prov_tsv} ({len(df_pair_prov)} rows, {len(df_pair_prov.columns)} columns)")

    # Dedicated Evidence-Level HIPPIE Provenance Table
    df_ev_prov = pd.DataFrame(evidence_provenance_rows)
    ev_prov_tsv = os.path.join(data_dir, f"{prefix}_positive_evidence_hippie_provenance.tsv")
    ev_prov_parquet = os.path.join(data_dir, f"{prefix}_positive_evidence_hippie_provenance.parquet")
    df_ev_prov.to_csv(ev_prov_tsv, sep="\t", index=False)
    df_ev_prov.to_parquet(ev_prov_parquet, index=False)
    print(f"Saved evidence HIPPIE provenance table: {ev_prov_tsv} ({len(df_ev_prov)} rows, {len(df_ev_prov.columns)} columns)")

    # ---------------------------------------------------------
    # 5. Benchmark Dataset Export (Deferred)
    # ---------------------------------------------------------
    # Benchmark dataset files are deferred pending criteria finalization.
    # Master evidence, deduplicated pairs, and HIPPIE provenance tables are maintained.

    # ---------------------------------------------------------
    # 6. Summary Statistics
    # ---------------------------------------------------------
    scores = df_pair_prov["hippie_score"]
    stats = {
        "species": prefix,
        "total_pairs": len(scores),
        "high_confidence_count": int((scores >= HIGH_CONFIDENCE_THRESHOLD).sum()),
        "medium_confidence_count": int(((scores >= 0.45) & (scores < HIGH_CONFIDENCE_THRESHOLD)).sum()),
        "low_confidence_count": int((scores < 0.45).sum()),
        "mean_score": round(float(scores.mean()), 4),
        "median_score": round(float(scores.median()), 4),
        "min_score": round(float(scores.min()), 4),
        "max_score": round(float(scores.max()), 4),
        "iqr": round(float(scores.quantile(0.75) - scores.quantile(0.25)), 4)
    }

    print(f"\nHIPPIE Scoring Summary ({prefix.upper()}):")
    print(f"  - Total Evaluated Pairs:       {stats['total_pairs']:,}")
    print(f"  - High Confidence (S >= {HIGH_CONFIDENCE_THRESHOLD}):  {stats['high_confidence_count']:,} ({stats['high_confidence_count']/stats['total_pairs']*100:.1f}%)")
    print(f"  - Medium Confidence:           {stats['medium_confidence_count']:,} ({stats['medium_confidence_count']/stats['total_pairs']*100:.1f}%)")
    print(f"  - Low Confidence:              {stats['low_confidence_count']:,} ({stats['low_confidence_count']/stats['total_pairs']*100:.1f}%)")
    print(f"  - Mean Score:                  {stats['mean_score']:.4f}")
    print(f"  - Score Range:                 [{stats['min_score']:.4f}, {stats['max_score']:.4f}]")

    return stats


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Compute HIPPIE Scores for PPI interactome datasets")
    parser.add_argument("--species", choices=["rice", "ara", "arabidopsis"], default="rice", help="Species code (default: rice)")
    parser.add_argument("--data-dir", default=None, help="Directory containing master files")
    args = parser.parse_args()

    run_hippie_pipeline(data_dir=args.data_dir, species=args.species)
