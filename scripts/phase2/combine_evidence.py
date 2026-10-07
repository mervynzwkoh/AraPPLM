"""
Phase 2: Streaming, Resumable ETL Evidence Normalization Pipeline.

Key architectural properties:
1. Canonical Warehouse: Hive-style partitioned Parquet dataset at local non-OneDrive root:
     <ARAPPLM_DATA_ROOT>/evidence_dataset/species_code=<SPECIES>/source_resource=<RESOURCE>/part_<CHUNK>.parquet
2. No 20 GB master TSV or duplicate training-candidate datasets.
3. Training eligibility is preserved strictly as the boolean column 'eligible_for_training'.
4. Full provenance preserved: source_id, source_file, source_row_number, source_record_id.
5. raw_fields_json omitted for STRING (reconstruction possible from preserved columns + checksummed raw).
6. Resumable at source/chunk level with atomic writes (.parquet.tmp -> .parquet) and progress manifest.
7. Strictly memory-bounded: ceiling of 1.5 GB, actively monitored via psutil.
8. Generates small human-readable TSV audit tables and a 5,000-record sample in the repository.
"""
import os
import sys
import gc
import json
import csv
import collections
from datetime import datetime, timezone
from pathlib import Path
import psutil
import pyarrow as pa
import pyarrow.parquet as pq

WORKSPACE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if WORKSPACE_ROOT not in sys.path:
    sys.path.insert(0, WORKSPACE_ROOT)

from scripts.phase2.harmonize_assays import AssayHarmonizer
from scripts.phase2.parse_intact import INTACT_SOURCES, parse_intact_source
from scripts.phase2.parse_biogrid import BIOGRID_SOURCES, parse_biogrid_source
from scripts.phase2.parse_xlms import parse_xlms_records, parse_xlms_stream
from scripts.phase2.parse_string import STRING_SOURCES, parse_string_stream
from scripts.phase2.parse_legacy import parse_legacy_benchmarks
from scripts.phase2.parse_bipseq import audit_bipseq

# Default to non-OneDrive local directory
DEFAULT_DATA_ROOT = os.path.join(os.path.expanduser("~"), "AraPPLM_Data", "phase2")
DATA_ROOT = os.environ.get("ARAPPLM_DATA_ROOT", DEFAULT_DATA_ROOT)
EVIDENCE_DATASET_DIR = os.path.join(DATA_ROOT, "evidence_dataset")

MEMORY_CEILING_MB = 1536.0  # 1.5 GB ceiling
CHUNK_SIZE_LARGE = 100000   # For STRING (100k rows/chunk)
CHUNK_SIZE_MED = 50000      # For XL-MS / IntAct / BioGRID

SCHEMA_COLUMNS = [
    "evidence_id",
    "source_id",
    "source_resource",
    "source_release",
    "source_file",
    "source_row_number",
    "source_record_id",
    "retrieval_date",
    "participant_a_original_id",
    "participant_a_id_namespace",
    "participant_a_id_value",
    "participant_a_original_name",
    "participant_a_aliases_raw",
    "participant_a_taxid",
    "participant_a_species",
    "participant_a_biological_role",
    "participant_a_experimental_role",
    "participant_b_original_id",
    "participant_b_id_namespace",
    "participant_b_id_value",
    "participant_b_original_name",
    "participant_b_aliases_raw",
    "participant_b_taxid",
    "participant_b_species",
    "participant_b_biological_role",
    "participant_b_experimental_role",
    "same_species_pair",
    "plant_plant_pair",
    "host_pathogen_pair",
    "taxon_ambiguous",
    "unordered_pair_key_provisional",
    "interaction_semantics",
    "interaction_type_raw",
    "interaction_type_psi_mi",
    "interaction_type_psi_mi_name",
    "detection_method_raw",
    "detection_method_psi_mi",
    "detection_method_psi_mi_name",
    "assay_family",
    "assay_supports_direct_binding",
    "assay_supports_physical_association",
    "assay_supports_proximity_only",
    "assay_is_genetic",
    "assay_is_computational",
    "native_in_planta",
    "publication_id_raw",
    "pmid",
    "doi",
    "publication_year",
    "publication_source",
    "throughput_raw",
    "throughput_class",
    "screen_id",
    "study_id",
    "intact_miscore",
    "biogrid_score",
    "string_combined_score",
    "string_physical_combined_score",
    "string_experimental_score",
    "string_database_score",
    "string_textmining_score",
    "native_experimental_status",
    "source_confidence_raw",
    "source_confidence_name",
    "xlms_interprotein_flag",
    "xlms_intraprotein_flag",
    "xlms_peptide_a",
    "xlms_peptide_b",
    "xlms_linker",
    "xlms_score",
    "xlms_evalue",
    "xlms_fdr",
    "is_external_holdout",
    "is_temporal_holdout",
    "eligible_for_training",
    "evidence_fingerprint",
    "possible_cross_database_duplicate",
    "needs_manual_review",
    "review_reason",
    "raw_fields_json"
]

ARROW_SCHEMA = pa.schema([(col, pa.string()) for col in SCHEMA_COLUMNS])

def format_val_str(v):
    if v is None:
        return "NA"
    if isinstance(v, bool):
        return "TRUE" if v else "FALSE"
    s = str(v).replace("\t", " ").replace("\n", " ").replace("\r", " ")
    return s if s != "" else "NA"

class MemoryMonitor:
    def __init__(self, ceiling_mb=MEMORY_CEILING_MB):
        self.ceiling_mb = ceiling_mb
        self.process = psutil.Process()
        self.peak_rss_mb = 0.0

    def check(self, context=""):
        rss_mb = self.process.memory_info().rss / (1024 * 1024)
        if rss_mb > self.peak_rss_mb:
            self.peak_rss_mb = rss_mb

        if rss_mb > self.ceiling_mb * 0.8:
            gc.collect()
            rss_mb = self.process.memory_info().rss / (1024 * 1024)

        if rss_mb > self.ceiling_mb:
            raise MemoryError(
                f"[FATAL] Peak RSS {rss_mb:.1f} MB exceeded ceiling of {self.ceiling_mb} MB at '{context}'"
            )
        return rss_mb

class ProgressTracker:
    def __init__(self, manifest_path):
        self.manifest_path = manifest_path
        self.state = {
            "dataset_root": EVIDENCE_DATASET_DIR,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "updated_utc": datetime.now(timezone.utc).isoformat(),
            "peak_rss_mb": 0.0,
            "sources": {}
        }
        self.load()

    def load(self):
        if os.path.exists(self.manifest_path):
            try:
                with open(self.manifest_path, "r", encoding="utf-8") as f:
                    self.state = json.load(f)
            except Exception:
                pass

    def save(self, peak_rss=0.0):
        self.state["updated_utc"] = datetime.now(timezone.utc).isoformat()
        if peak_rss > self.state.get("peak_rss_mb", 0.0):
            self.state["peak_rss_mb"] = round(peak_rss, 2)
        os.makedirs(os.path.dirname(self.manifest_path), exist_ok=True)
        tmp_path = self.manifest_path + ".tmp"
        with open(tmp_path, "w", encoding="utf-8") as f:
            json.dump(self.state, f, indent=2)
        os.replace(tmp_path, self.manifest_path)

    def is_source_completed(self, source_id):
        s = self.state.get("sources", {}).get(source_id)
        return s is not None and s.get("status") == "COMPLETED"

    def record_source_completed(self, source_id, stats):
        self.state.setdefault("sources", {})[source_id] = {
            "status": "COMPLETED",
            "source_resource": stats.get("source_resource"),
            "file": stats.get("file"),
            "n_input": stats.get("n_input"),
            "n_normalized": stats.get("n_normalized"),
            "n_rejected": stats.get("n_rejected"),
            "rejection_rate_pct": stats.get("rejection_rate_pct"),
            "chunks_written": stats.get("chunks_written"),
            "completed_at": datetime.now(timezone.utc).isoformat()
        }
        self.save()

def write_parquet_chunk_atomic(records_chunk, out_dir, chunk_index):
    """
    Writes a chunk of records to Hive partition directory using atomic rename.
    """
    os.makedirs(out_dir, exist_ok=True)
    chunk_filename = f"part_{chunk_index:05d}.parquet"
    final_path = os.path.join(out_dir, chunk_filename)
    tmp_path = os.path.join(out_dir, f".part_{chunk_index:05d}.parquet.tmp")

    # Build columnar data
    col_dict = {col: [] for col in SCHEMA_COLUMNS}
    for r in records_chunk:
        for col in SCHEMA_COLUMNS:
            col_dict[col].append(format_val_str(r.get(col)))

    tbl = pa.Table.from_pydict(col_dict, schema=ARROW_SCHEMA)
    pq.write_table(tbl, tmp_path, compression="snappy")
    os.replace(tmp_path, final_path)
    return final_path

def run_pipeline():
    print("================================================================================")
    print("           AraPPLM PHASE 2: STREAMING & RESUMABLE ETL PIPELINE                 ")
    print("================================================================================")
    print(f"Non-OneDrive Data Root: {DATA_ROOT}")
    print(f"Canonical Evidence Dataset: {EVIDENCE_DATASET_DIR}")
    print(f"Memory Ceiling: {MEMORY_CEILING_MB} MB")
    print("================================================================================\n")

    mem = MemoryMonitor(ceiling_mb=MEMORY_CEILING_MB)
    progress_file = os.path.join(WORKSPACE_ROOT, "manifests", "phase2_progress.json")
    tracker = ProgressTracker(progress_file)

    harmonizer = AssayHarmonizer()

    # Local repo directories for audits, sample, schemas, rejected
    repo_p2_dir = os.path.join(WORKSPACE_ROOT, "data", "interim", "phase2")
    audits_dir = os.path.join(repo_p2_dir, "audits")
    rejected_dir = os.path.join(repo_p2_dir, "rejected")
    sample_dir = os.path.join(repo_p2_dir, "sample")

    os.makedirs(audits_dir, exist_ok=True)
    os.makedirs(rejected_dir, exist_ok=True)
    os.makedirs(sample_dir, exist_ok=True)
    os.makedirs(EVIDENCE_DATASET_DIR, exist_ok=True)

    # Pointer in repository
    pointer_path = os.path.join(repo_p2_dir, "DATA_LOCATION.txt")
    with open(pointer_path, "w", encoding="utf-8") as f:
        f.write("AraPPLM Phase 2 Canonical Evidence Dataset\n")
        f.write("==========================================\n\n")
        f.write(f"Dataset Root: {EVIDENCE_DATASET_DIR}\n")
        f.write("Partitioning Format: Hive-style partitioned Parquet\n")
        f.write("Partition Layout: species_code=<SPECIES>/source_resource=<RESOURCE>/part_<CHUNK>.parquet\n")
        f.write("Training Eligibility: Preserved as boolean column 'eligible_for_training'\n")
        f.write("Last Updated UTC: " + datetime.now(timezone.utc).isoformat() + "\n")

    # Global tracking accumulators
    source_stats = {}
    all_rejected = []
    taxon_counter = collections.Counter()
    species_counter = collections.Counter()
    assay_counter = collections.Counter()
    semantics_counter = collections.Counter()
    sample_records = []
    total_normalized = 0
    total_training_candidates = 0

    def sample_accumulator(recs, target_sample=350):
        # Accumulate a stratified slice for sample TSV
        if len(sample_records) < 5000 and recs:
            step = max(1, len(recs) // target_sample)
            for i in range(0, len(recs), step):
                if len(sample_records) >= 5000:
                    break
                sample_records.append(recs[i])

    def accumulate_stats(recs):
        nonlocal total_normalized, total_training_candidates
        for r in recs:
            total_normalized += 1
            if r.get("eligible_for_training"):
                total_training_candidates += 1
            taxon_counter[str(r.get("participant_a_taxid"))] += 1
            taxon_counter[str(r.get("participant_b_taxid"))] += 1
            species_counter[str(r.get("participant_a_species"))] += 1
            species_counter[str(r.get("participant_b_species"))] += 1
            assay_counter[str(r.get("assay_family"))] += 1
            semantics_counter[str(r.get("interaction_semantics"))] += 1

    # =========================================================================
    # 1. INTACT SOURCES
    # =========================================================================
    print("--- [1/5] Processing IntAct Sources ---")
    for src in INTACT_SOURCES:
        src_id = src["source_id"]
        species_code = src["species_code"].upper()
        if species_code.startswith("RICE"):
            species_code = "RICE"
        part_dir = os.path.join(EVIDENCE_DATASET_DIR, f"species_code={species_code}", "source_resource=INTACT")

        if tracker.is_source_completed(src_id):
            print(f"  [SKIPPED] {src_id} already completed (resumed).")
            # Load stored stats
            st = tracker.state["sources"][src_id]
            source_stats[src_id] = st
            continue

        print(f"  Parsing {src_id} ({src['file_path']})...")
        recs, rej = parse_intact_source(src, harmonizer)
        all_rejected.extend(rej)

        # Write in chunks of 50k, offsetting by existing parts in directory
        existing_parts = glob.glob(os.path.join(part_dir, "part_*.parquet"))
        chunks_written = len(existing_parts)
        for i in range(0, len(recs), CHUNK_SIZE_MED):
            chunk = recs[i:i + CHUNK_SIZE_MED]
            write_parquet_chunk_atomic(chunk, part_dir, chunks_written)
            chunks_written += 1

        accumulate_stats(recs)
        sample_accumulator(recs, target_sample=350)

        n_in = len(recs) + len(rej)
        st = {
            "source_resource": "IntAct",
            "file": src["file_path"],
            "n_input": n_in,
            "n_normalized": len(recs),
            "n_rejected": len(rej),
            "rejection_rate_pct": (len(rej) / n_in * 100) if n_in > 0 else 0.0,
            "chunks_written": chunks_written
        }
        source_stats[src_id] = st
        tracker.record_source_completed(src_id, st)
        mem.check(src_id)
        print(f"    Normalized: {len(recs):,}, Rejected: {len(rej)}, Chunks: {chunks_written}")
        del recs, rej
        gc.collect()

    # =========================================================================
    # 2. BIOGRID SOURCES
    # =========================================================================
    print("\n--- [2/5] Processing BioGRID Sources ---")
    for src in BIOGRID_SOURCES:
        src_id = src["source_id"]
        species_code = src["species_code"].upper()
        part_dir = os.path.join(EVIDENCE_DATASET_DIR, f"species_code={species_code}", "source_resource=BIOGRID")

        if tracker.is_source_completed(src_id):
            print(f"  [SKIPPED] {src_id} already completed (resumed).")
            st = tracker.state["sources"][src_id]
            source_stats[src_id] = st
            continue

        print(f"  Parsing {src_id} ({src['file_path']})...")
        recs, rej = parse_biogrid_source(src, harmonizer)
        all_rejected.extend(rej)

        chunks_written = 0
        for i in range(0, len(recs), CHUNK_SIZE_MED):
            chunk = recs[i:i + CHUNK_SIZE_MED]
            write_parquet_chunk_atomic(chunk, part_dir, chunks_written)
            chunks_written += 1

        accumulate_stats(recs)
        sample_accumulator(recs, target_sample=350)

        n_in = len(recs) + len(rej)
        st = {
            "source_resource": "BioGRID",
            "file": src["file_path"],
            "n_input": n_in,
            "n_normalized": len(recs),
            "n_rejected": len(rej),
            "rejection_rate_pct": (len(rej) / n_in * 100) if n_in > 0 else 0.0,
            "chunks_written": chunks_written
        }
        source_stats[src_id] = st
        tracker.record_source_completed(src_id, st)
        mem.check(src_id)
        print(f"    Normalized: {len(recs):,}, Rejected: {len(rej)}, Chunks: {chunks_written}")
        del recs, rej
        gc.collect()

    # =========================================================================
    # 3. 2026 PhoX XL-MS (Trinh et al., 2026)
    # =========================================================================
    print("\n--- [3/5] Processing 2026 PhoX XL-MS Source ---")
    xlms_id = "SRC_A4_XLMS_2026_PLINK_SEARCH"
    part_dir_xlms = os.path.join(EVIDENCE_DATASET_DIR, "species_code=ARA", "source_resource=XLMS")

    if tracker.is_source_completed(xlms_id):
        print(f"  [SKIPPED] {xlms_id} already completed (resumed).")
        st = tracker.state["sources"][xlms_id]
        source_stats[xlms_id] = st
    else:
        print(f"  Streaming {xlms_id}...")
        norm_count = 0
        rej_count = 0
        chunks_written = 0

        for batch_recs, batch_rej in parse_xlms_stream(batch_size=CHUNK_SIZE_MED):
            all_rejected.extend(batch_rej)
            rej_count += len(batch_rej)
            if batch_recs:
                write_parquet_chunk_atomic(batch_recs, part_dir_xlms, chunks_written)
                chunks_written += 1
                norm_count += len(batch_recs)
                accumulate_stats(batch_recs)
                sample_accumulator(batch_recs, target_sample=50)

            cur_rss = mem.check(f"{xlms_id}_chunk_{chunks_written}")
            tracker.save(peak_rss=cur_rss)

        n_in = norm_count + rej_count
        st = {
            "source_resource": "PRIDE / pLink 3.2",
            "file": "data/raw/arabidopsis/xlms_2026/Total_XL_plink3-2_v3.csv",
            "n_input": n_in,
            "n_normalized": norm_count,
            "n_rejected": rej_count,
            "rejection_rate_pct": (rej_count / n_in * 100) if n_in > 0 else 0.0,
            "chunks_written": chunks_written
        }
        source_stats[xlms_id] = st
        tracker.record_source_completed(xlms_id, st)
        print(f"    Normalized: {norm_count:,}, Rejected: {rej_count}, Chunks: {chunks_written}, Peak RSS: {mem.peak_rss_mb:.1f} MB")
        gc.collect()

    # =========================================================================
    # 4. STRING PHYSICAL LINKS SOURCES (Streaming in bounded chunks)
    # =========================================================================
    print("\n--- [4/5] Processing STRING Physical Links (Streaming 100k Bounded Chunks) ---")
    for src in STRING_SOURCES:
        src_id = src["source_id"]
        species_code = src["species_code"].upper()
        part_dir = os.path.join(EVIDENCE_DATASET_DIR, f"species_code={species_code}", "source_resource=STRING")

        if tracker.is_source_completed(src_id):
            print(f"  [SKIPPED] {src_id} already completed (resumed).")
            st = tracker.state["sources"][src_id]
            source_stats[src_id] = st
            continue

        print(f"  Streaming {src_id} ({src['file_path']})...")
        norm_count = 0
        rej_count = 0
        chunks_written = 0

        for batch_recs, batch_rej in parse_string_stream(src, batch_size=CHUNK_SIZE_LARGE):
            all_rejected.extend(batch_rej)
            rej_count += len(batch_rej)

            if batch_recs:
                write_parquet_chunk_atomic(batch_recs, part_dir, chunks_written)
                chunks_written += 1
                norm_count += len(batch_recs)
                accumulate_stats(batch_recs)
                sample_accumulator(batch_recs, target_sample=100)

            # Monitor memory and force collection per chunk
            cur_rss = mem.check(f"{src_id}_chunk_{chunks_written}")
            tracker.save(peak_rss=cur_rss)

        n_in = norm_count + rej_count
        st = {
            "source_resource": "STRING",
            "file": src["file_path"],
            "n_input": n_in,
            "n_normalized": norm_count,
            "n_rejected": rej_count,
            "rejection_rate_pct": (rej_count / n_in * 100) if n_in > 0 else 0.0,
            "chunks_written": chunks_written
        }
        source_stats[src_id] = st
        tracker.record_source_completed(src_id, st)
        print(f"    Normalized: {norm_count:,}, Rejected: {rej_count}, Chunks: {chunks_written}, Peak RSS: {mem.peak_rss_mb:.1f} MB")
        gc.collect()

    # =========================================================================
    # 5. LEGACY BENCHMARKS & BIP-SEQ AUDITS
    # =========================================================================
    print("\n--- [5/5] Normalizing Legacy Benchmarks & Auditing BIP-seq ---")
    legacy_stats = parse_legacy_benchmarks()
    bipseq_status = audit_bipseq()

    # =========================================================================
    # 6. WRITE COMPACT HUMAN-READABLE SAMPLE TSV
    # =========================================================================
    sample_tsv_path = os.path.join(sample_dir, "sample_normalized_evidence.tsv")
    print(f"\nWriting sample evidence TSV ({len(sample_records):,} rows) to: {sample_tsv_path}...")
    with open(sample_tsv_path, "w", encoding="utf-8", newline="") as f_samp:
        writer = csv.writer(f_samp, delimiter="\t")
        writer.writerow(SCHEMA_COLUMNS)
        for r in sample_records:
            writer.writerow([format_val_str(r.get(col)) for col in SCHEMA_COLUMNS])

    # =========================================================================
    # 7. WRITE REJECTED RECORDS LOG
    # =========================================================================
    rej_path = os.path.join(rejected_dir, "rejected_records.tsv")
    with open(rej_path, "w", encoding="utf-8", newline="") as f:
        f.write("source_id\tsource_file\tsource_row\treason\traw_record\n")
        for r in all_rejected:
            clean_raw = str(r.get("raw_record", "")).replace("\t", " ").replace("\n", " ").replace("\r", " ")
            f.write(f"{r.get('source_id')}\t{r.get('source_file')}\t{r.get('source_row')}\t{r.get('reason')}\t{clean_raw}\n")

    # =========================================================================
    # 8. WRITE AUDIT SUMMARY TABLES
    # =========================================================================
    # 8.1 source_record_counts.tsv
    with open(os.path.join(audits_dir, "source_record_counts.tsv"), "w", encoding="utf-8", newline="") as f:
        f.write("source_id\tsource_resource\tfile_path\tn_input\tn_normalized\tn_rejected\trejection_rate_pct\n")
        for sid, st in source_stats.items():
            f.write(f"{sid}\t{st['source_resource']}\t{st['file']}\t{st['n_input']}\t{st['n_normalized']}\t{st['n_rejected']}\t{st['rejection_rate_pct']:.4f}\n")

    # 8.2 taxon_summary.tsv
    with open(os.path.join(audits_dir, "taxon_summary.tsv"), "w", encoding="utf-8", newline="") as f:
        f.write("taxid\tparticipant_count\tnotes\n")
        for taxid, cnt in taxon_counter.most_common():
            f.write(f"{taxid}\t{cnt}\tAudited participant taxid\n")

    # 8.3 assay_summary.tsv
    with open(os.path.join(audits_dir, "assay_summary.tsv"), "w", encoding="utf-8", newline="") as f:
        f.write("assay_family\trecord_count\n")
        for fam, cnt in assay_counter.most_common():
            f.write(f"{fam}\t{cnt}\n")

    # 8.4 interaction_semantics_summary.tsv
    with open(os.path.join(audits_dir, "interaction_semantics_summary.tsv"), "w", encoding="utf-8", newline="") as f:
        f.write("interaction_semantics\trecord_count\n")
        for sem, cnt in semantics_counter.most_common():
            f.write(f"{sem}\t{cnt}\n")

    # 8.5 unmapped_psi_mi_terms.tsv
    with open(os.path.join(audits_dir, "unmapped_psi_mi_terms.tsv"), "w", encoding="utf-8", newline="") as f:
        f.write("unmapped_term\tnotes\n")
        for term in sorted(harmonizer.unmapped_terms):
            f.write(f"{term}\tRequires manual review\n")

    tracker.save(peak_rss=mem.peak_rss_mb)

    print("\n================================================================================")
    print("                    PHASE 2 ETL EXECUTION SUMMARY                               ")
    print("================================================================================")
    print(f"Canonical Evidence Warehouse: {EVIDENCE_DATASET_DIR}")
    print(f"Total Normalized Evidence Records: {total_normalized:,}")
    print(f"Total Training Candidates (eligible_for_training=TRUE): {total_training_candidates:,}")
    print(f"Total Holdout Records (eligible_for_training=FALSE): {total_normalized - total_training_candidates:,}")
    print(f"Total Rejected Records: {len(all_rejected)}")
    print(f"Unmapped Assay Terms: {len(harmonizer.unmapped_terms)}")
    print(f"Peak Physical RSS Memory: {mem.peak_rss_mb:.1f} MB (Ceiling: {MEMORY_CEILING_MB} MB)")
    print(f"Progress Manifest: {progress_file}")
    print(f"Sample Evidence Table: {sample_tsv_path}")
    print("================================================================================\n")

    return {
        "total_normalized": total_normalized,
        "total_training_candidates": total_training_candidates,
        "total_rejected": len(all_rejected),
        "peak_rss_mb": mem.peak_rss_mb,
        "source_stats": source_stats
    }

if __name__ == "__main__":
    run_pipeline()
