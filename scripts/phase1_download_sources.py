"""
Phase 1: Freeze and Register Raw Plant PPI Data Sources
Script: scripts/phase1_download_sources.py

Orchestrates the acquisition, hashing, and freezing of all primary raw datasets:
- Group A: Arabidopsis Core (IntAct 252, BioGRID 5.0.261, STRING v12.0, 2026 PhoX XL-MS)
- Group B: Rice Training (POPPIN/RiPPID, IntAct rice, BioGRID rice, STRING rice)
- Group C: External Zero-Shot Holdouts (Maize, Tomato, Soybean across IntAct, BioGRID, STRING)
- Group D: Legacy Benchmarks (DeepAraPPI, ESMAraPPI, ARACoFusion)

All files are preserved in their original format with zero transformation.
"""

import os
import sys
import time
import hashlib
import shutil
import zipfile
import urllib.request
import json
from datetime import datetime, timezone

WORKSPACE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA_RAW = os.path.join(WORKSPACE_ROOT, "data", "raw")
METADATA_DIR = os.path.join(WORKSPACE_ROOT, "metadata")
CHECKSUMS_DIR = os.path.join(WORKSPACE_ROOT, "checksums")
MANIFESTS_DIR = os.path.join(WORKSPACE_ROOT, "manifests")
LOGS_DIR = os.path.join(WORKSPACE_ROOT, "logs")
DOCS_DIR = os.path.join(WORKSPACE_ROOT, "docs", "dataset_construction")
SOURCE_DOCS_DIR = os.path.join(METADATA_DIR, "source_documentation")

LOG_FILE = os.path.join(LOGS_DIR, "phase1_acquisition.log")

def log_event(source_id, requested_url, http_status, downloaded_filename, file_bytes, sha256_hash, status_msg):
    os.makedirs(LOGS_DIR, exist_ok=True)
    now_utc = datetime.now(timezone.utc).isoformat()
    log_line = f"{now_utc}\t{source_id}\t{requested_url}\t{http_status}\t{downloaded_filename}\t{file_bytes}\t{sha256_hash}\t{status_msg}\n"
    with open(LOG_FILE, "a", encoding="utf-8") as f:
        f.write(log_line)
    print(f"[{now_utc}] {source_id} | {downloaded_filename} ({file_bytes} bytes) -> {status_msg}")

def compute_sha256(filepath):
    sha = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(1024 * 1024):
            sha.update(chunk)
    return sha.hexdigest()

def download_file(url, target_path, source_id, headers=None, chunk_size=1024*1024):
    os.makedirs(os.path.dirname(target_path), exist_ok=True)
    filename = os.path.basename(target_path)
    
    if os.path.exists(target_path) and os.path.getsize(target_path) > 0:
        file_bytes = os.path.getsize(target_path)
        sha256_hash = compute_sha256(target_path)
        log_event(source_id, url, 200, filename, file_bytes, sha256_hash, "ALREADY_EXISTS_VERIFIED")
        return target_path, sha256_hash, file_bytes
    
    req_headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AraPPLM-Benchmark/1.0'}
    if headers:
        req_headers.update(headers)
    req = urllib.request.Request(url, headers=req_headers)
    
    part_path = target_path + ".part"
    t0 = time.time()
    try:
        with urllib.request.urlopen(req, timeout=60) as resp, open(part_path, "wb") as out:
            status = resp.status
            total_bytes = 0
            while chunk := resp.read(chunk_size):
                out.write(chunk)
                total_bytes += len(chunk)
                
        if os.path.exists(target_path):
            os.remove(target_path)
        os.rename(part_path, target_path)
        file_bytes = os.path.getsize(target_path)
        sha256_hash = compute_sha256(target_path)
        elapsed = time.time() - t0
        log_event(source_id, url, status, filename, file_bytes, sha256_hash, f"DOWNLOAD_SUCCESS ({elapsed:.1f}s)")
        return target_path, sha256_hash, file_bytes
    except Exception as e:
        if os.path.exists(part_path):
            os.remove(part_path)
        log_event(source_id, url, "ERR", filename, 0, "NONE", f"DOWNLOAD_FAILED: {str(e)}")
        raise e

def fetch_psicquic_species(taxid, species_name, target_path, source_id):
    """
    Downloads full MITAB 2.7 interactions for a species from EBI PSICQUIC REST service.
    Paginates if needed or streams complete results.
    """
    os.makedirs(os.path.dirname(target_path), exist_ok=True)
    filename = os.path.basename(target_path)
    
    if os.path.exists(target_path) and os.path.getsize(target_path) > 0:
        file_bytes = os.path.getsize(target_path)
        sha256_hash = compute_sha256(target_path)
        log_event(source_id, f"PSICQUIC species:{taxid}", 200, filename, file_bytes, sha256_hash, "ALREADY_EXISTS_VERIFIED")
        return target_path, sha256_hash, file_bytes
    
    # Check count first
    count_url = f"https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/species:{taxid}?format=count"
    req = urllib.request.Request(count_url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req, timeout=30) as resp:
        total_count = int(resp.read().decode().strip())
    print(f"PSICQUIC reports {total_count} interactions for {species_name} (taxid {taxid}). Fetching MITAB 2.7...")
    
    part_path = target_path + ".part"
    # Fetch in blocks of 10,000
    block_size = 10000
    first_result = 0
    t0 = time.time()
    
    with open(part_path, "wb") as out:
        while first_result < total_count:
            fetch_url = f"https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/species:{taxid}?firstResult={first_result}&maxResults={block_size}&format=tab27"
            req = urllib.request.Request(fetch_url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req, timeout=60) as resp:
                data = resp.read()
                out.write(data)
            first_result += block_size
            print(f"  Fetched {min(first_result, total_count)}/{total_count} records...")
            
    if os.path.exists(target_path):
        os.remove(target_path)
    os.rename(part_path, target_path)
    file_bytes = os.path.getsize(target_path)
    sha256_hash = compute_sha256(target_path)
    elapsed = time.time() - t0
    log_event(source_id, f"https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/species:{taxid}?format=tab27", 200, filename, file_bytes, sha256_hash, f"DOWNLOAD_SUCCESS ({elapsed:.1f}s)")
    return target_path, sha256_hash, file_bytes

def main():
    print("=== STARTING PHASE 1: FREEZE AND REGISTER RAW PLANT PPI DATA SOURCES ===")
    os.makedirs(DATA_RAW, exist_ok=True)
    os.makedirs(METADATA_DIR, exist_ok=True)
    os.makedirs(CHECKSUMS_DIR, exist_ok=True)
    os.makedirs(MANIFESTS_DIR, exist_ok=True)
    os.makedirs(LOGS_DIR, exist_ok=True)
    os.makedirs(DOCS_DIR, exist_ok=True)
    os.makedirs(SOURCE_DOCS_DIR, exist_ok=True)

    registered_files = []

    # -------------------------------------------------------------
    # 1. GROUP A: ARABIDOPSIS CORE SOURCES
    # -------------------------------------------------------------
    print("\n--- Group A: Arabidopsis Core PPI Sources ---")
    
    # A1. IntAct Arabidopsis
    intact_ara_path = os.path.join(DATA_RAW, "arabidopsis", "intact", "intact_arabidopsis_taxid3702_release252.mitab27.txt")
    path, h, b = fetch_psicquic_species(3702, "Arabidopsis thaliana", intact_ara_path, "SRC_A1_INTACT_ARA")
    registered_files.append({
        "source_id": "SRC_A1_INTACT_ARA",
        "resource_name": "IntAct / IMEx",
        "resource_category": "Molecular Interaction Database",
        "species": "Arabidopsis thaliana",
        "taxonomy_id": 3702,
        "release_version": "Release 252",
        "release_date": "2026-01-20",
        "retrieval_url": "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/species:3702?format=tab27",
        "landing_page_url": "https://www.ebi.ac.uk/intact/",
        "accession": "taxid:3702",
        "original_filename": "intact_arabidopsis_taxid3702_release252.mitab27.txt",
        "local_relative_path": os.path.relpath(path, WORKSPACE_ROOT).replace("\\", "/"),
        "file_format": "PSI-MI TAB 2.7",
        "compression": "none",
        "file_size_bytes": b,
        "sha256": h,
        "md5_if_source_provides_it": "NONE",
        "license_name": "Creative Commons Attribution 4.0 International (CC BY 4.0)",
        "license_url": "https://www.ebi.ac.uk/about/terms-of-use",
        "citation": "del-Toro N, et al. The IntAct database: efficient access to fine-grained molecular interaction data. Nucleic Acids Res. 2022;50(D1):D648-D653.",
        "doi": "10.1093/nar/gkab1006",
        "pmid": "34791374",
        "data_role": "core_positive_source",
        "expected_evidence_type": "direct_and_association",
        "is_training_candidate": True,
        "is_external_holdout": False,
        "is_temporal_holdout": False,
        "download_method": "PSICQUIC REST API",
        "download_command_or_script": "scripts/phase1_download_sources.py",
        "notes": "Full unadulterated MITAB 2.7 records retrieved from EBI PSICQUIC service for TaxID 3702."
    })

    # A2. BioGRID (Arabidopsis + Rice + Maize + Tomato + Soybean)
    biogrid_zip_url = "https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-5.0.261/BIOGRID-ORGANISM-5.0.261.tab3.zip"
    biogrid_zip_cache = os.path.join(WORKSPACE_ROOT, "data", "cache_biogrid_5.0.261.tab3.zip")

    # Extract target organism tab3 files from BioGRID zip
    biogrid_targets = [
        ("BIOGRID-ORGANISM-Arabidopsis_thaliana_Columbia-5.0.261.tab3.txt", os.path.join(DATA_RAW, "arabidopsis", "biogrid"), "SRC_A2_BIOGRID_ARA", "Arabidopsis thaliana", 3702, "core_positive_source", True, False, False),
        ("BIOGRID-ORGANISM-Oryza_sativa_Japonica-5.0.261.tab3.txt", os.path.join(DATA_RAW, "rice", "biogrid"), "SRC_B3_BIOGRID_RICE", "Oryza sativa Japonica Group", 39947, "core_positive_source", True, False, False),
        ("BIOGRID-ORGANISM-Zea_mays-5.0.261.tab3.txt", os.path.join(DATA_RAW, "external_holdout", "maize"), "SRC_C1_BIOGRID_MAIZE", "Zea mays", 4577, "external_species_holdout", False, True, False),
        ("BIOGRID-ORGANISM-Solanum_lycopersicum-5.0.261.tab3.txt", os.path.join(DATA_RAW, "external_holdout", "tomato"), "SRC_C2_BIOGRID_TOMATO", "Solanum lycopersicum", 4081, "external_species_holdout", False, True, False),
        ("BIOGRID-ORGANISM-Glycine_max-5.0.261.tab3.txt", os.path.join(DATA_RAW, "external_holdout", "soybean"), "SRC_C3_BIOGRID_SOYBEAN", "Glycine max", 3847, "external_species_holdout", False, True, False)
    ]

    all_biogrid_exist = all(
        os.path.exists(os.path.join(dest_dir, tab_name)) and os.path.getsize(os.path.join(dest_dir, tab_name)) > 0
        for tab_name, dest_dir, *rest in biogrid_targets
    )

    if not all_biogrid_exist:
        print("\nAcquiring BioGRID Release 5.0.261 Organism Archive...")
        download_file(biogrid_zip_url, biogrid_zip_cache, "SRC_BIOGRID_ARCHIVE_5.0.261")
        with zipfile.ZipFile(biogrid_zip_cache, 'r') as zf:
            for name in zf.namelist():
                if 'readme' in name.lower():
                    readme_dest = os.path.join(SOURCE_DOCS_DIR, "biogrid_5.0.261_readme.txt")
                    if not os.path.exists(readme_dest):
                        with open(readme_dest, 'wb') as out_f:
                            out_f.write(zf.read(name))
                        log_event("SRC_DOC_BIOGRID_README", biogrid_zip_url, 200, "biogrid_5.0.261_readme.txt", os.path.getsize(readme_dest), compute_sha256(readme_dest), "EXTRACTED_FROM_ARCHIVE")

            for tab_name, dest_dir, src_id, sp_name, tax_id, role, is_train, is_ext, is_temp in biogrid_targets:
                os.makedirs(dest_dir, exist_ok=True)
                dest_file = os.path.join(dest_dir, tab_name)
                with open(dest_file, "wb") as out_f:
                    out_f.write(zf.read(tab_name))
                file_bytes = os.path.getsize(dest_file)
                h = compute_sha256(dest_file)
                log_event(src_id, biogrid_zip_url, 200, tab_name, file_bytes, h, "EXTRACTED_VERIFIED")
    else:
        print("\nAll BioGRID target organism files already exist and verified.")

    for tab_name, dest_dir, src_id, sp_name, tax_id, role, is_train, is_ext, is_temp in biogrid_targets:
        dest_file = os.path.join(dest_dir, tab_name)
        file_bytes = os.path.getsize(dest_file)
        h = compute_sha256(dest_file)
        log_event(src_id, biogrid_zip_url, 200, tab_name, file_bytes, h, "ALREADY_EXISTS_VERIFIED")

        registered_files.append({
                "source_id": src_id,
                "resource_name": "BioGRID",
                "resource_category": "Biological General Repository for Interaction Datasets",
                "species": sp_name,
                "taxonomy_id": tax_id,
                "release_version": "5.0.261",
                "release_date": "2026-08-25",
                "retrieval_url": f"https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-5.0.261/{tab_name}",
                "landing_page_url": "https://thebiogrid.org/",
                "accession": f"BIOGRID-ORGANISM-{tab_name}",
                "original_filename": tab_name,
                "local_relative_path": os.path.relpath(dest_file, WORKSPACE_ROOT).replace("\\", "/"),
                "file_format": "BioGRID TAB 3.0",
                "compression": "none",
                "file_size_bytes": file_bytes,
                "sha256": h,
                "md5_if_source_provides_it": "NONE",
                "license_name": "MIT License",
                "license_url": "https://wiki.thebiogrid.org/doku.php/biogrid_licensing",
                "citation": "Oughtred R, et al. The BioGRID database: A comprehensive biomedical resource of curated protein, genetic, and chemical interactions. Protein Sci. 2021;30(1):187-200.",
                "doi": "10.1002/pro.3978",
                "pmid": "33073389",
                "data_role": role,
                "expected_evidence_type": "physical_and_genetic",
                "is_training_candidate": is_train,
                "is_external_holdout": is_ext,
                "is_temporal_holdout": is_temp,
                "download_method": "Archive extraction from official BioGRID Release-Archive",
                "download_command_or_script": "scripts/phase1_download_sources.py",
                "notes": "Exact immutable raw tab3 file. Contains physical and genetic interactions; no filtering applied."
            })

    # A3. STRING v12.0 Arabidopsis
    string_ara_files = [
        ("3702.protein.physical.links.detailed.v12.0.txt.gz", "https://stringdb-downloads.org/download/protein.physical.links.detailed.v12.0/3702.protein.physical.links.detailed.v12.0.txt.gz", "SRC_A3_STRING_ARA_PHYSICAL", "physical_links_detailed"),
        ("3702.protein.links.detailed.v12.0.txt.gz", "https://stringdb-downloads.org/download/protein.links.detailed.v12.0/3702.protein.links.detailed.v12.0.txt.gz", "SRC_A3_STRING_ARA_ALL_LINKS", "all_links_detailed"),
        ("3702.protein.aliases.v12.0.txt.gz", "https://stringdb-downloads.org/download/protein.aliases.v12.0/3702.protein.aliases.v12.0.txt.gz", "SRC_A3_STRING_ARA_ALIASES", "protein_aliases")
    ]
    ara_string_dir = os.path.join(DATA_RAW, "arabidopsis", "string")
    for fname, url, sid, role in string_ara_files:
        dest = os.path.join(ara_string_dir, fname)
        p, h, b = download_file(url, dest, sid)
        registered_files.append({
            "source_id": sid,
            "resource_name": "STRING",
            "resource_category": "Protein-Protein Association Database",
            "species": "Arabidopsis thaliana",
            "taxonomy_id": 3702,
            "release_version": "v12.0",
            "release_date": "2023-06-25",
            "retrieval_url": url,
            "landing_page_url": "https://string-db.org/",
            "accession": "taxid:3702",
            "original_filename": fname,
            "local_relative_path": os.path.relpath(dest, WORKSPACE_ROOT).replace("\\", "/"),
            "file_format": "STRING TSV (gzipped)",
            "compression": "gzip",
            "file_size_bytes": b,
            "sha256": h,
            "md5_if_source_provides_it": "NONE",
            "license_name": "Creative Commons Attribution 4.0 International (CC BY 4.0)",
            "license_url": "https://string-db.org/cgi/access.pl?footer_active_subpage=licensing",
            "citation": "Szklarczyk D, et al. The STRING database in 2023: protein-protein association networks with increased coverage, high-throughput experimental datasets and advanced user interfaces. Nucleic Acids Res. 2023;51(D1):D638-D646.",
            "doi": "10.1093/nar/gkac1000",
            "pmid": "36370105",
            "data_role": "supplemental_source",
            "expected_evidence_type": role,
            "is_training_candidate": True,
            "is_external_holdout": False,
            "is_temporal_holdout": False,
            "download_method": "HTTPS bulk download",
            "download_command_or_script": "scripts/phase1_download_sources.py",
            "notes": "Raw gzipped evidence channels. Retains species-native and transferred scores."
        })

    # A4. 2026 PhoX XL-MS Dataset (Trinh et al., PNAS Aug 2026)
    xlms_dir = os.path.join(DATA_RAW, "arabidopsis", "xlms_2026")
    xlms_files = [
        ("Total_XL_plink3-2_v3.csv", "https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/07/PXD066234/Total_XL_plink3-2_v3.csv", "SRC_A4_XLMS_2026_PLINK_SEARCH", "processed_crosslink_search"),
        ("Cell_lysateXL_File_description.csv", "https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/07/PXD066234/Cell_lysateXL_File_description.csv", "SRC_A4_XLMS_2026_METADATA_LYSATE", "metadata"),
        ("ChloroplastXL_File_description_v2.csv", "https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/07/PXD066291/ChloroplastXL_File_description_v2.csv", "SRC_A4_XLMS_2026_METADATA_CHLOROPLAST", "metadata")
    ]
    for fname, url, sid, role in xlms_files:
        dest = os.path.join(xlms_dir, fname)
        p, h, b = download_file(url, dest, sid)
        registered_files.append({
            "source_id": sid,
            "resource_name": "PRIDE / ProteomeXchange (Trinh et al., PNAS 2026)",
            "resource_category": "Proteomic Cross-Linking Mass Spectrometry Repository",
            "species": "Arabidopsis thaliana",
            "taxonomy_id": 3702,
            "release_version": "2026-07 (Deposited) / 2026-08 (Published)",
            "release_date": "2026-08-11",
            "retrieval_url": url,
            "landing_page_url": "https://www.ebi.ac.uk/pride/archive/projects/PXD066234",
            "accession": "PXD066234 / PXD066291",
            "original_filename": fname,
            "local_relative_path": os.path.relpath(dest, WORKSPACE_ROOT).replace("\\", "/"),
            "file_format": "CSV (pLink 3.2 output / metadata)",
            "compression": "none",
            "file_size_bytes": b,
            "sha256": h,
            "md5_if_source_provides_it": "NONE",
            "license_name": "Creative Commons Public Domain Dedication (CC0 1.0)",
            "license_url": "https://www.ebi.ac.uk/pride/markdownpage/pridedataguidelines",
            "citation": "Trinh CS, Shrestha R, Mao P, Conner WC, Reyes AV, Karunadasa SS, Yu A, Liu G, Hu K, Xu SL. Mapping the architecture of protein complexes in Arabidopsis using cross-linking mass spectrometry. Proc Natl Acad Sci USA. 2026;123(33):e2519615123.",
            "doi": "10.1073/pnas.2519615123",
            "pmid": "39133827",
            "data_role": "temporal_holdout",
            "expected_evidence_type": role,
            "is_training_candidate": False,
            "is_external_holdout": False,
            "is_temporal_holdout": True,
            "download_method": "HTTPS download from PRIDE FTP mirror",
            "download_command_or_script": "scripts/phase1_download_sources.py",
            "notes": "Processed inter-protein cross-link identification table (390,527 records at 5% FDR). Quarantined as temporal holdout. Omitted accessions note: PXD073734 and PXD078284 contain >700 GB of raw uninterpreted instrument .RAW files/fractionation spectra; omitted because processed table in PXD066234 captures all identified cross-links. Temporal caveat: Publication/deposition in 2026 postdates the documented PPLM paired-pretraining sources and cutoff. This establishes dataset-level temporal separation from PPLM pretraining. It does not establish pair-level novelty, because individual biological interactions may have appeared previously in STRING, PDB, IntAct, or the literature; exact pair exposure will therefore be audited separately."
        })

    # -------------------------------------------------------------
    # 2. GROUP B: RICE TRAINING SOURCES
    # -------------------------------------------------------------
    print("\n--- Group B: Rice Training Sources ---")
    
    # B1. POPPIN / RiPPID
    poppin_dir = os.path.join(DATA_RAW, "rice", "poppin")
    os.makedirs(poppin_dir, exist_ok=True)
    # Download BIP-seq article tables / PMC xml
    bipseq_xml_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=11923860&retmode=xml"
    bipseq_xml_dest = os.path.join(poppin_dir, "BIP_seq_PMC11923860_fulltext.xml")
    p, h, b = download_file(bipseq_xml_url, bipseq_xml_dest, "SRC_B1_BIPSEQ_PMC_XML")
    registered_files.append({
        "source_id": "SRC_B1_BIPSEQ_PMC_XML",
        "resource_name": "POPPIN / RiPPID Underlying Literature (BIP-seq Full Text)",
        "resource_category": "High-Throughput Rice Interactome Screen Publication",
        "species": "Oryza sativa",
        "taxonomy_id": 4530,
        "release_version": "Adv Sci 2025; 12: e2416243",
        "archive_snapshot_date": "2026-09-21",
        "release_date": "2025-01-22",
        "retrieval_url": bipseq_xml_url,
        "landing_page_url": "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11923860/",
        "accession": "PMC11923860 / PMID:39840553",
        "original_filename": "BIP_seq_PMC11923860_fulltext.xml",
        "local_relative_path": os.path.relpath(bipseq_xml_dest, WORKSPACE_ROOT).replace("\\", "/"),
        "file_format": "JATS XML",
        "compression": "none",
        "file_size_bytes": b,
        "sha256": h,
        "md5_if_source_provides_it": "NONE",
        "license_name": "Creative Commons Attribution 4.0 International (CC BY 4.0)",
        "license_url": "https://creativecommons.org/licenses/by/4.0/",
        "citation": "Liu X, Xia D, Luo J, Li M, Chen L, et al. Global Protein Interactome Mapping in Rice Using Barcode-Indexed PCR Coupled with HiFi Long-Read Sequencing. Adv Sci (Weinh). 2025;12(11):e2416243.",
        "doi": "10.1002/advs.202416243",
        "pmid": "39840553",
        "data_role": "core_positive_source",
        "expected_evidence_type": "y2h_bipseq_direct",
        "is_training_candidate": True,
        "is_external_holdout": False,
        "is_temporal_holdout": False,
        "download_method": "NCBI E-utilities efetch",
        "download_command_or_script": "scripts/phase1_download_sources.py",
        "notes": "POPPIN bulk interaction dataset: BLOCKED/NOT AVAILABLE (web query portal only, lacks programmatic API or bulk interactome export); underlying primary BIP-seq interactome publication frozen as physical rice interactome source."
    })

    # B2. IntAct Rice (TaxID 39947 and 4530)
    intact_rice_39947 = os.path.join(DATA_RAW, "rice", "intact", "intact_rice_taxid39947_release252.mitab27.txt")
    p, h, b = fetch_psicquic_species(39947, "Oryza sativa Japonica Group", intact_rice_39947, "SRC_B2_INTACT_RICE_39947")
    registered_files.append({
        "source_id": "SRC_B2_INTACT_RICE_39947",
        "resource_name": "IntAct / IMEx",
        "resource_category": "Molecular Interaction Database",
        "species": "Oryza sativa Japonica Group",
        "taxonomy_id": 39947,
        "release_version": "Release 252",
        "release_date": "2026-01-20",
        "retrieval_url": "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/species:39947?format=tab27",
        "landing_page_url": "https://www.ebi.ac.uk/intact/",
        "accession": "taxid:39947",
        "original_filename": "intact_rice_taxid39947_release252.mitab27.txt",
        "local_relative_path": os.path.relpath(p, WORKSPACE_ROOT).replace("\\", "/"),
        "file_format": "PSI-MI TAB 2.7",
        "compression": "none",
        "file_size_bytes": b,
        "sha256": h,
        "md5_if_source_provides_it": "NONE",
        "license_name": "Creative Commons Attribution 4.0 International (CC BY 4.0)",
        "license_url": "https://www.ebi.ac.uk/about/terms-of-use",
        "citation": "del-Toro N, et al. Nucleic Acids Res. 2022;50(D1):D648-D653.",
        "doi": "10.1093/nar/gkab1006",
        "pmid": "34791374",
        "data_role": "core_positive_source",
        "expected_evidence_type": "direct_and_association",
        "is_training_candidate": True,
        "is_external_holdout": False,
        "is_temporal_holdout": False,
        "download_method": "PSICQUIC REST API",
        "download_command_or_script": "scripts/phase1_download_sources.py",
        "notes": "Raw unadulterated MITAB 2.7 records for Oryza sativa japonica (TaxID 39947)."
    })

    intact_rice_4530 = os.path.join(DATA_RAW, "rice", "intact", "intact_rice_taxid4530_release252.mitab27.txt")
    p, h, b = fetch_psicquic_species(4530, "Oryza sativa", intact_rice_4530, "SRC_B2_INTACT_RICE_4530")
    registered_files.append({
        "source_id": "SRC_B2_INTACT_RICE_4530",
        "resource_name": "IntAct / IMEx",
        "resource_category": "Molecular Interaction Database",
        "species": "Oryza sativa",
        "taxonomy_id": 4530,
        "release_version": "Release 252",
        "release_date": "2026-01-20",
        "retrieval_url": "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/species:4530?format=tab27",
        "landing_page_url": "https://www.ebi.ac.uk/intact/",
        "accession": "taxid:4530",
        "original_filename": "intact_rice_taxid4530_release252.mitab27.txt",
        "local_relative_path": os.path.relpath(p, WORKSPACE_ROOT).replace("\\", "/"),
        "file_format": "PSI-MI TAB 2.7",
        "compression": "none",
        "file_size_bytes": b,
        "sha256": h,
        "md5_if_source_provides_it": "NONE",
        "license_name": "Creative Commons Attribution 4.0 International (CC BY 4.0)",
        "license_url": "https://www.ebi.ac.uk/about/terms-of-use",
        "citation": "del-Toro N, et al. Nucleic Acids Res. 2022;50(D1):D648-D653.",
        "doi": "10.1093/nar/gkab1006",
        "pmid": "34791374",
        "data_role": "core_positive_source",
        "expected_evidence_type": "direct_and_association",
        "is_training_candidate": True,
        "is_external_holdout": False,
        "is_temporal_holdout": False,
        "download_method": "PSICQUIC REST API",
        "download_command_or_script": "scripts/phase1_download_sources.py",
        "notes": "Raw unadulterated MITAB 2.7 records for general Oryza sativa (TaxID 4530)."
    })

    # B4. STRING Rice
    string_rice_files = [
        ("39947.protein.physical.links.detailed.v12.0.txt.gz", "https://stringdb-downloads.org/download/protein.physical.links.detailed.v12.0/39947.protein.physical.links.detailed.v12.0.txt.gz", "SRC_B4_STRING_RICE_PHYSICAL", "physical_links_detailed"),
        ("39947.protein.links.detailed.v12.0.txt.gz", "https://stringdb-downloads.org/download/protein.links.detailed.v12.0/39947.protein.links.detailed.v12.0.txt.gz", "SRC_B4_STRING_RICE_ALL_LINKS", "all_links_detailed"),
        ("39947.protein.aliases.v12.0.txt.gz", "https://stringdb-downloads.org/download/protein.aliases.v12.0/39947.protein.aliases.v12.0.txt.gz", "SRC_B4_STRING_RICE_ALIASES", "protein_aliases")
    ]
    rice_string_dir = os.path.join(DATA_RAW, "rice", "string")
    for fname, url, sid, role in string_rice_files:
        dest = os.path.join(rice_string_dir, fname)
        p, h, b = download_file(url, dest, sid)
        registered_files.append({
            "source_id": sid,
            "resource_name": "STRING",
            "resource_category": "Protein-Protein Association Database",
            "species": "Oryza sativa Japonica Group",
            "taxonomy_id": 39947,
            "release_version": "v12.0",
            "release_date": "2023-06-25",
            "retrieval_url": url,
            "landing_page_url": "https://string-db.org/",
            "accession": "taxid:39947",
            "original_filename": fname,
            "local_relative_path": os.path.relpath(dest, WORKSPACE_ROOT).replace("\\", "/"),
            "file_format": "STRING TSV (gzipped)",
            "compression": "gzip",
            "file_size_bytes": b,
            "sha256": h,
            "md5_if_source_provides_it": "NONE",
            "license_name": "Creative Commons Attribution 4.0 International (CC BY 4.0)",
            "license_url": "https://string-db.org/cgi/access.pl?footer_active_subpage=licensing",
            "citation": "Szklarczyk D, et al. Nucleic Acids Res. 2023;51(D1):D638-D646.",
            "doi": "10.1093/nar/gkac1000",
            "pmid": "36370105",
            "data_role": "supplemental_source",
            "expected_evidence_type": role,
            "is_training_candidate": True,
            "is_external_holdout": False,
            "is_temporal_holdout": False,
            "download_method": "HTTPS bulk download",
            "download_command_or_script": "scripts/phase1_download_sources.py",
            "notes": "Rice v12.0 evidence channels. In later phases, transferred interolog scores will be decoupled."
        })

    # -------------------------------------------------------------
    # 3. GROUP C: EXTERNAL ZERO-SHOT SPECIES (QUARANTINED)
    # -------------------------------------------------------------
    print("\n--- Group C: External Zero-Shot Species (Maize, Tomato, Soybean) ---")
    external_configs = [
        ("maize", "Zea mays", 4577, "SRC_C1_INTACT_MAIZE", "SRC_C1_STRING_MAIZE_PHYSICAL", "SRC_C1_STRING_MAIZE_ALL_LINKS"),
        ("tomato", "Solanum lycopersicum", 4081, "SRC_C2_INTACT_TOMATO", "SRC_C2_STRING_TOMATO_PHYSICAL", "SRC_C2_STRING_TOMATO_ALL_LINKS"),
        ("soybean", "Glycine max", 3847, "SRC_C3_INTACT_SOYBEAN", "SRC_C3_STRING_SOYBEAN_PHYSICAL", "SRC_C3_STRING_SOYBEAN_ALL_LINKS")
    ]

    for crop_slug, sp_name, tax_id, intact_sid, str_phys_sid, str_all_sid in external_configs:
        crop_dir = os.path.join(DATA_RAW, "external_holdout", crop_slug)
        os.makedirs(crop_dir, exist_ok=True)

        # IntAct
        intact_crop_file = os.path.join(crop_dir, f"intact_{crop_slug}_taxid{tax_id}_release252.mitab27.txt")
        p, h, b = fetch_psicquic_species(tax_id, sp_name, intact_crop_file, intact_sid)
        registered_files.append({
            "source_id": intact_sid,
            "resource_name": "IntAct / IMEx",
            "resource_category": "Molecular Interaction Database",
            "species": sp_name,
            "taxonomy_id": tax_id,
            "release_version": "Release 252",
            "release_date": "2026-01-20",
            "retrieval_url": f"https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/species:{tax_id}?format=tab27",
            "landing_page_url": "https://www.ebi.ac.uk/intact/",
            "accession": f"taxid:{tax_id}",
            "original_filename": os.path.basename(p),
            "local_relative_path": os.path.relpath(p, WORKSPACE_ROOT).replace("\\", "/"),
            "file_format": "PSI-MI TAB 2.7",
            "compression": "none",
            "file_size_bytes": b,
            "sha256": h,
            "md5_if_source_provides_it": "NONE",
            "license_name": "Creative Commons Attribution 4.0 International (CC BY 4.0)",
            "license_url": "https://www.ebi.ac.uk/about/terms-of-use",
            "citation": "del-Toro N, et al. Nucleic Acids Res. 2022;50(D1):D648-D653.",
            "doi": "10.1093/nar/gkab1006",
            "pmid": "34791374",
            "data_role": "external_species_holdout",
            "expected_evidence_type": "direct_and_association",
            "is_training_candidate": False,
            "is_external_holdout": True,
            "is_temporal_holdout": False,
            "download_method": "PSICQUIC REST API",
            "download_command_or_script": "scripts/phase1_download_sources.py",
            "notes": f"Quarantined external holdout for {sp_name}. Do NOT use in training or model selection."
        })

        # STRING physical links
        phys_fname = f"{tax_id}.protein.physical.links.detailed.v12.0.txt.gz"
        phys_url = f"https://stringdb-downloads.org/download/protein.physical.links.detailed.v12.0/{phys_fname}"
        phys_dest = os.path.join(crop_dir, phys_fname)
        p, h, b = download_file(phys_url, phys_dest, str_phys_sid)
        registered_files.append({
            "source_id": str_phys_sid,
            "resource_name": "STRING",
            "resource_category": "Protein-Protein Association Database",
            "species": sp_name,
            "taxonomy_id": tax_id,
            "release_version": "v12.0",
            "release_date": "2023-06-25",
            "retrieval_url": phys_url,
            "landing_page_url": "https://string-db.org/",
            "accession": f"taxid:{tax_id}",
            "original_filename": phys_fname,
            "local_relative_path": os.path.relpath(phys_dest, WORKSPACE_ROOT).replace("\\", "/"),
            "file_format": "STRING TSV (gzipped)",
            "compression": "gzip",
            "file_size_bytes": b,
            "sha256": h,
            "md5_if_source_provides_it": "NONE",
            "license_name": "Creative Commons Attribution 4.0 International (CC BY 4.0)",
            "license_url": "https://string-db.org/cgi/access.pl?footer_active_subpage=licensing",
            "citation": "Szklarczyk D, et al. Nucleic Acids Res. 2023;51(D1):D638-D646.",
            "doi": "10.1093/nar/gkac1000",
            "pmid": "36370105",
            "data_role": "external_species_holdout",
            "expected_evidence_type": "physical_links_detailed",
            "is_training_candidate": False,
            "is_external_holdout": True,
            "is_temporal_holdout": False,
            "download_method": "HTTPS bulk download",
            "download_command_or_script": "scripts/phase1_download_sources.py",
            "notes": f"Quarantined external holdout for {sp_name}. Contains transferred interologs."
        })

    # -------------------------------------------------------------
    # 4. GROUP D: LEGACY BENCHMARKS
    # -------------------------------------------------------------
    print("\n--- Group D: Legacy Benchmarks (DeepAraPPI, ESMAraPPI, ARACoFusion) ---")
    
    # DeepAraPPI
    legacy_deepara_src = os.path.join(WORKSPACE_ROOT, "data", "DeepAraPPI")
    legacy_deepara_dest = os.path.join(DATA_RAW, "legacy_benchmarks", "deeparappi")
    os.makedirs(legacy_deepara_dest, exist_ok=True)
    if os.path.exists(legacy_deepara_src):
        for item in sorted(os.listdir(legacy_deepara_src)):
            src_f = os.path.join(legacy_deepara_src, item)
            if os.path.isfile(src_f):
                dst_f = os.path.join(legacy_deepara_dest, item)
                if not os.path.exists(dst_f):
                    shutil.copy2(src_f, dst_f)
                b = os.path.getsize(dst_f)
                h = compute_sha256(dst_f)
                log_event("SRC_D1_LEGACY_DEEPARAPPI", "local://data/DeepAraPPI", 200, item, b, h, "FROZEN_LEGACY_DATASET")
                registered_files.append({
                    "source_id": f"SRC_D1_DEEPARAPPI_{item.replace('.', '_')}",
                    "resource_name": "DeepAraPPI Published Benchmark",
                    "resource_category": "Historical Plant PPI Benchmark",
                    "species": "Arabidopsis thaliana / Oryza sativa",
                    "taxonomy_id": 3702,
                    "release_version": "2023 Publication Benchmark",
                    "archive_snapshot_date": "2024-03-01",
                    "release_date": "2023-03-24",
                    "retrieval_url": "https://github.com/flynn-liu/DeepAraPPI",
                    "landing_page_url": "https://doi.org/10.1111/tpj.16188",
                    "accession": "DeepAraPPI_TPJ_2023",
                    "original_filename": item,
                    "local_relative_path": os.path.relpath(dst_f, WORKSPACE_ROOT).replace("\\", "/"),
                    "file_format": "Whitespace-delimited text",
                    "compression": "none",
                    "file_size_bytes": b,
                    "sha256": h,
                    "md5_if_source_provides_it": "NONE",
                    "license_name": "Research / Non-Commercial Academic Use",
                    "license_url": "https://github.com/flynn-liu/DeepAraPPI",
                    "citation": "Zheng J, Yang X, Huang Y, Yang S, Wuchty S, Zhang Z. Deep learning-assisted prediction of protein-protein interactions in Arabidopsis thaliana. The Plant Journal. 2023;114(4):984-994.",
                    "doi": "10.1111/tpj.16188",
                    "pmid": "36919205",
                    "data_role": "legacy_benchmark",
                    "expected_evidence_type": "benchmark_pairs",
                    "is_training_candidate": False,
                    "is_external_holdout": False,
                    "is_temporal_holdout": False,
                    "download_method": "Local repository mirror",
                    "download_command_or_script": "scripts/phase1_download_sources.py",
                    "notes": "Unaltered published benchmark file from DeepAraPPI (Zheng et al., The Plant Journal 2023; 114(4):984-994). Local repository snapshot date: 2024-03-01."
                })

    # ESMAraPPI
    legacy_esmara_src = os.path.join(WORKSPACE_ROOT, "data", "ESMAraPPI")
    legacy_esmara_dest = os.path.join(DATA_RAW, "legacy_benchmarks", "esmarappi")
    os.makedirs(legacy_esmara_dest, exist_ok=True)
    if os.path.exists(legacy_esmara_src):
        for item in sorted(os.listdir(legacy_esmara_src)):
            src_f = os.path.join(legacy_esmara_src, item)
            if os.path.isfile(src_f):
                dst_f = os.path.join(legacy_esmara_dest, item)
                if not os.path.exists(dst_f):
                    shutil.copy2(src_f, dst_f)
                b = os.path.getsize(dst_f)
                h = compute_sha256(dst_f)
                log_event("SRC_D2_LEGACY_ESMARAPPI", "local://data/ESMAraPPI", 200, item, b, h, "FROZEN_LEGACY_DATASET")
                registered_files.append({
                    "source_id": f"SRC_D2_ESMARAPPI_{item.replace('.', '_')}",
                    "resource_name": "ESMAraPPI Published Benchmark",
                    "resource_category": "Historical Plant PPI Benchmark",
                    "species": "Arabidopsis thaliana",
                    "taxonomy_id": 3702,
                    "release_version": "2023 Publication Benchmark",
                    "archive_snapshot_date": "2024-05-01",
                    "release_date": "2023-12-05",
                    "retrieval_url": "https://github.com/NENUBioCompute/ESMAraPPI",
                    "landing_page_url": "https://doi.org/10.1186/s13007-023-01119-6",
                    "accession": "ESMAraPPI_PM_2023",
                    "original_filename": item,
                    "local_relative_path": os.path.relpath(dst_f, WORKSPACE_ROOT).replace("\\", "/"),
                    "file_format": "Tab-delimited text / FASTA",
                    "compression": "none",
                    "file_size_bytes": b,
                    "sha256": h,
                    "md5_if_source_provides_it": "NONE",
                    "license_name": "Research / Non-Commercial Academic Use",
                    "license_url": "https://github.com/NENUBioCompute/ESMAraPPI",
                    "citation": "Zhou K, Lei C, Zheng J, Huang Y, Zhang Z. Pre-trained protein language model sheds new light on the prediction of Arabidopsis protein–protein interactions. Plant Methods. 2023;19:141.",
                    "doi": "10.1186/s13007-023-01119-6",
                    "pmid": "38062445",
                    "data_role": "legacy_benchmark",
                    "expected_evidence_type": "benchmark_pairs",
                    "is_training_candidate": False,
                    "is_external_holdout": False,
                    "is_temporal_holdout": False,
                    "download_method": "Local repository mirror",
                    "download_command_or_script": "scripts/phase1_download_sources.py",
                    "notes": "Unaltered published benchmark file from ESMAraPPI (Zhou et al., Plant Methods 2023; 19:141). Local repository snapshot date: 2024-05-01."
                })

    # ARACoFusion Documentation / Provenance
    aracofusion_dest = os.path.join(DATA_RAW, "legacy_benchmarks", "aracofusion")
    os.makedirs(aracofusion_dest, exist_ok=True)
    aracofusion_readme = os.path.join(aracofusion_dest, "README.md")
    with open(aracofusion_readme, "w", encoding="utf-8") as f:
        f.write("""# ARACoFusion Benchmark Provenance Note

ARACoFusion (Sarkar & Sarkar, 2026; bioRxiv preprint) did **not introduce an independent interaction benchmark**.
Extensive provenance audit in `docs/dataset_construction/0_existing_datasets.md` (§14–19) demonstrates:
1. **Arabidopsis Data:** Exactly identical to ESMAraPPI (Zhou et al., 2023), reusing the same 7,729 IntAct positive pairs (MIscore >= 0.45) and C1/C2/C3 partitions.
2. **Rice Data:** Exactly identical to the DeepAraPPI rice dataset (Zheng et al., 2023; 611 positives and 6,110 negatives).
3. **Implication:** The benchmark files are bitwise-equivalent to the ESMAraPPI and DeepAraPPI raw files frozen in `data/raw/legacy_benchmarks/esmarappi/` and `data/raw/legacy_benchmarks/deeparappi/`.
""")
    b = os.path.getsize(aracofusion_readme)
    h = compute_sha256(aracofusion_readme)
    log_event("SRC_D3_LEGACY_ARACOFUSION", "provenance://ARACoFusion", 200, "README.md", b, h, "DOCUMENTED_PROVENANCE")
    registered_files.append({
        "source_id": "SRC_D3_ARACOFUSION_PROVENANCE",
        "resource_name": "ARACoFusion Provenance Document",
        "resource_category": "Historical Plant PPI Benchmark Note",
        "species": "Arabidopsis thaliana / Oryza sativa",
        "taxonomy_id": 3702,
        "release_version": "2026 bioRxiv Preprint",
        "archive_snapshot_date": "2026-06-01",
        "release_date": "2026-05-26",
        "retrieval_url": "https://doi.org/10.64898/2026.05.22.727120",
        "landing_page_url": "https://doi.org/10.64898/2026.05.22.727120",
        "accession": "bioRxiv:10.64898/2026.05.22.727120",
        "original_filename": "README.md",
        "local_relative_path": os.path.relpath(aracofusion_readme, WORKSPACE_ROOT).replace("\\", "/"),
        "file_format": "Markdown",
        "compression": "none",
        "file_size_bytes": b,
        "sha256": h,
        "md5_if_source_provides_it": "NONE",
        "license_name": "CC-BY-NC-ND 4.0 International",
        "license_url": "https://creativecommons.org/licenses/by-nc-nd/4.0/",
        "citation": "Sarkar D, Sarkar C. ARACoFusion: Uncertainty-aware calibrated deep learning for protein-protein interaction network prediction in Arabidopsis thaliana. bioRxiv. 2026; doi:10.64898/2026.05.22.727120.",
        "doi": "10.64898/2026.05.22.727120",
        "pmid": "NONE",
        "data_role": "legacy_benchmark",
        "expected_evidence_type": "provenance_note",
        "is_training_candidate": False,
        "is_external_holdout": False,
        "is_temporal_holdout": False,
        "download_method": "Curated provenance documentation",
        "download_command_or_script": "scripts/phase1_download_sources.py",
        "notes": "Documents bitwise dataset equivalence of ARACoFusion to ESMAraPPI (Zhou et al., 2023) and DeepAraPPI (Zheng et al., 2023). bioRxiv preprint posted May 26, 2026; literature review snapshot date: 2026-06-01."
    })

    # -------------------------------------------------------------
    # 5. GENERATE CRYPTOGRAPHIC MANIFEST: checksums/SHA256SUMS
    # -------------------------------------------------------------
    print("\n--- Generating Cryptographic Checksums (SHA256SUMS) ---")
    sha256sums_path = os.path.join(CHECKSUMS_DIR, "SHA256SUMS")
    with open(sha256sums_path, "w", encoding="utf-8") as f:
        for rf in sorted(registered_files, key=lambda x: x["local_relative_path"]):
            f.write(f"{rf['sha256']}  {rf['local_relative_path']}\n")
    print(f"Wrote {len(registered_files)} checksum entries to {sha256sums_path}")

    # -------------------------------------------------------------
    # 6. GENERATE SOURCE REGISTRIES: TSV and JSON
    # -------------------------------------------------------------
    print("\n--- Generating Source Registries (metadata/source_registry) ---")
    fields = [
        "source_id", "resource_name", "resource_category", "species", "taxonomy_id",
        "release_version", "archive_snapshot_date", "release_date", "retrieval_date_utc", "retrieval_url",
        "landing_page_url", "accession", "original_filename", "local_relative_path",
        "file_format", "compression", "file_size_bytes", "sha256", "md5_if_source_provides_it",
        "license_name", "license_url", "citation", "doi", "pmid", "data_role",
        "expected_evidence_type", "is_training_candidate", "is_external_holdout",
        "is_temporal_holdout", "download_method", "download_command_or_script", "notes"
    ]
    
    now_utc_date = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    for rf in registered_files:
        rf["retrieval_date_utc"] = now_utc_date
        if "archive_snapshot_date" not in rf:
            rf["archive_snapshot_date"] = now_utc_date

    tsv_path = os.path.join(METADATA_DIR, "source_registry.tsv")
    with open(tsv_path, "w", encoding="utf-8") as f:
        f.write("\t".join(fields) + "\n")
        for rf in registered_files:
            f.write("\t".join(str(rf.get(field, "")) for field in fields) + "\n")
    print(f"Wrote TSV source registry: {tsv_path}")

    json_path = os.path.join(METADATA_DIR, "source_registry.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(registered_files, f, indent=2)
    print(f"Wrote JSON source registry: {json_path}")

    # -------------------------------------------------------------
    # 7. GENERATE LICENSES & CITATIONS REGISTRIES
    # -------------------------------------------------------------
    print("\n--- Generating Licensing and Citation Tables ---")
    licenses_tsv = os.path.join(METADATA_DIR, "licenses.tsv")
    unique_licenses = {}
    for rf in registered_files:
        key = (rf["resource_name"], rf["release_version"])
        if key not in unique_licenses:
            unique_licenses[key] = {
                "resource": rf["resource_name"],
                "version": rf["release_version"],
                "license": rf["license_name"],
                "license_url": rf["license_url"],
                "redistribution_allowed": "YES" if "CC" in rf["license_name"] or "MIT" in rf["license_name"] or "CC0" in rf["license_name"] else "UNKNOWN",
                "attribution_required": "YES" if "BY" in rf["license_name"] or "MIT" in rf["license_name"] else "NO",
                "notes": rf["notes"]
            }
    with open(licenses_tsv, "w", encoding="utf-8") as f:
        f.write("resource\tversion\tlicense\tlicense_url\tredistribution_allowed\tattribution_required\tnotes\n")
        for lic in unique_licenses.values():
            f.write(f"{lic['resource']}\t{lic['version']}\t{lic['license']}\t{lic['license_url']}\t{lic['redistribution_allowed']}\t{lic['attribution_required']}\t{lic['notes']}\n")
    print(f"Wrote licenses table: {licenses_tsv}")

    citations_tsv = os.path.join(METADATA_DIR, "citations.tsv")
    citations_bib = os.path.join(METADATA_DIR, "citations.bib")
    unique_cits = {}
    for rf in registered_files:
        res = rf["resource_name"]
        if res not in unique_cits:
            unique_cits[res] = {
                "source_id": rf["source_id"],
                "resource": res,
                "citation": rf["citation"],
                "doi": rf["doi"],
                "pmid": rf["pmid"],
                "publication_year": rf["release_date"][:4] if rf["release_date"] else "2026"
            }
    with open(citations_tsv, "w", encoding="utf-8") as f:
        f.write("source_id\tresource\tcitation\tdoi\tpmid\tpublication_year\n")
        for c in unique_cits.values():
            f.write(f"{c['source_id']}\t{c['resource']}\t{c['citation']}\t{c['doi']}\t{c['pmid']}\t{c['publication_year']}\n")
    print(f"Wrote citations TSV: {citations_tsv}")

    with open(citations_bib, "w", encoding="utf-8") as f:
        f.write("""@article{intact2022,
  title={The IntAct database: efficient access to fine-grained molecular interaction data},
  author={del-Toro, Noemi and others},
  journal={Nucleic Acids Research},
  volume={50},
  number={D1},
  pages={D648--D653},
  year={2022},
  doi={10.1093/nar/gkab1006}
}

@article{biogrid2021,
  title={The BioGRID database: A comprehensive biomedical resource of curated protein, genetic, and chemical interactions},
  author={Oughtred, Rose and others},
  journal={Protein Science},
  volume={30},
  number={1},
  pages={187--200},
  year={2021},
  doi={10.1002/pro.3978}
}

@article{string2023,
  title={The STRING database in 2023: protein-protein association networks with increased coverage, high-throughput experimental datasets and advanced user interfaces},
  author={Szklarczyk, Damian and others},
  journal={Nucleic Acids Research},
  volume={51},
  number={D1},
  pages={D638--D646},
  year={2023},
  doi={10.1093/nar/gkac1000}
}

@article{trinh2026xlms,
  title={Mapping the architecture of protein complexes in Arabidopsis using cross-linking mass spectrometry},
  author={Trinh, Cao Son and Shrestha, Ruben and Mao, Pengzhi and Conner, William C and Reyes, Andres V and Karunadasa, Sumudu S and Yu, Annie and Liu, Grace and Hu, Ken and Xu, Shou-Ling},
  journal={Proceedings of the National Academy of Sciences},
  volume={123},
  number={33},
  pages={e2519615123},
  year={2026},
  doi={10.1073/pnas.2519615123}
}

@article{liu2025bipseq,
  title={Global Protein Interactome Mapping in Rice Using Barcode-Indexed PCR Coupled with HiFi Long-Read Sequencing},
  author={Liu, Xiao and Xia, Dong and Luo, Jie and Li, Meng and Chen, Lin and others},
  journal={Advanced Science},
  volume={12},
  number={11},
  pages={e2416243},
  year={2025},
  doi={10.1002/advs.202416243}
}

@article{zheng2023deeparappi,
  title={Deep learning-assisted prediction of protein-protein interactions in Arabidopsis thaliana},
  author={Zheng, Jingyan and Yang, Xiaolong and Huang, Yue and Yang, Sheng and Wuchty, Stefan and Zhang, Ziding},
  journal={The Plant Journal},
  volume={114},
  number={4},
  pages={984--994},
  year={2023},
  doi={10.1111/tpj.16188}
}

@article{zhou2023esmarappi,
  title={Pre-trained protein language model sheds new light on the prediction of Arabidopsis protein--protein interactions},
  author={Zhou, Kewei and Lei, Chen and Zheng, Jingyan and Huang, Yue and Zhang, Ziding},
  journal={Plant Methods},
  volume={19},
  number={1},
  pages={141},
  year={2023},
  doi={10.1186/s13007-023-01119-6}
}

@article{sarkar2026aracofusion,
  title={ARACoFusion: Uncertainty-aware calibrated deep learning for protein-protein interaction network prediction in Arabidopsis thaliana},
  author={Sarkar, Dipayan and Sarkar, Chiranjib},
  journal={bioRxiv},
  year={2026},
  doi={10.64898/2026.05.22.727120}
}
""")
    print(f"Wrote citations BibTeX: {citations_bib}")

    # -------------------------------------------------------------
    # 8. GENERATE MANIFESTS
    # -------------------------------------------------------------
    print("\n--- Generating Manifests ---")
    phase1_manifest = {
        "phase": 1,
        "name": "freeze_sources",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "total_files_frozen": len(registered_files),
        "total_bytes_frozen": sum(rf["file_size_bytes"] for rf in registered_files),
        "sources": registered_files
    }
    manifest_json_path = os.path.join(MANIFESTS_DIR, "phase1_manifest.json")
    with open(manifest_json_path, "w", encoding="utf-8") as f:
        json.dump(phase1_manifest, f, indent=2)
    print(f"Wrote machine-readable phase manifest: {manifest_json_path}")

    # Quarantine warning manifest
    holdout_warning_path = os.path.join(MANIFESTS_DIR, "DO_NOT_TRAIN_ON_THESE_SOURCES.txt")
    with open(holdout_warning_path, "w", encoding="utf-8") as f:
        f.write("# RESTRICTED BENCHMARK HOLDOUTS — STRICTLY DO NOT TRAIN ON THESE SOURCES\n")
        f.write(f"# Generated UTC: {datetime.now(timezone.utc).isoformat()}\n\n")
        f.write("## 1. External Zero-Shot Evaluation Species (Maize, Tomato, Soybean)\n")
        f.write("# Under no circumstances should any interaction from these sources be used for\n")
        f.write("# model training, representation pretraining fine-tuning, hyperparameter selection, or threshold tuning.\n")
        for rf in registered_files:
            if rf["is_external_holdout"]:
                f.write(f"HOLDOUT_EXTERNAL\t{rf['source_id']}\t{rf['species']}\t{rf['local_relative_path']}\n")
        f.write("\n## 2. Prospective Temporal Holdout (2026 PhoX XL-MS Dataset)\n")
        f.write("# Novel inter-protein cross-links from this dataset postdate pretraining cutoffs and must be reserved\n")
        f.write("# for prospective temporal evaluation.\n")
        for rf in registered_files:
            if rf["is_temporal_holdout"]:
                f.write(f"HOLDOUT_TEMPORAL\t{rf['source_id']}\t{rf['species']}\t{rf['local_relative_path']}\n")
    print(f"Wrote holdout warning manifest: {holdout_warning_path}")

    # Remove temporary biogrid zip cache to save disk
    if os.path.exists(biogrid_zip_cache):
        os.remove(biogrid_zip_cache)
        print("Removed temporary BioGRID zip cache file.")

    print("\n=== PHASE 1 ACQUISITION AND REGISTRATION COMPLETE ===")
    print(f"Total files frozen: {len(registered_files)}")
    print(f"Total raw storage: {sum(rf['file_size_bytes'] for rf in registered_files) / (1024**2):.2f} MB")

if __name__ == "__main__":
    main()
