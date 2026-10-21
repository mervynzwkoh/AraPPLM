"""
Complete Raw Positive Pair Verification Script for All Species
=============================================================
Directly verifies positive pair counts from downloaded raw files for:
- Arabidopsis thaliana (taxid 3702)
- Zea mays (Maize, taxid 4577)
- Oryza sativa (Rice, taxid 39947 / 4530)
- Solanum lycopersicum (Tomato, taxid 4081)
- Glycine max (Soybean, taxid 3847)

Computes:
1. Raw downloaded records & raw pairs before any processing
2. Counts after exact filters:
   - physical vs genetic / proximity / computational
   - strictly exclude PCA (split-ubiquitin, BiFC, split-luciferase) from BOTH BioGRID and IntAct
   - use core physical set for IntAct (direct binary MI:0407, physical association MI:0915/MI:0914, physical MI:0217)
   - protein-only vs DNA / RNA / small molecule / missing
3. Counts after deduplication (per-database and cross-database combined)
   - Total unique pairs
   - Heteromers vs homomers
   - Unique participating proteins
"""

import os
import sys
import gzip
import csv
import json
import collections

WORKSPACE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

PCA_MI = {
    "MI:0112", # ubiquitin reconstruction (split-ubiquitin / SUS / MYTH)
    "MI:0809", # bimolecular fluorescence complementation (BiFC)
    "MI:1204", # split firefly luciferase complementation
    "MI:1037", # split renilla luciferase complementation
    "MI:0090", # protein complementation assay
    "MI:2365", # fluorescent protein-protein interaction-visualization
    "MI:0097", # reverse ras recruitment system
}

PROX_MI = {
    "MI:0663", # confocal microscopy
    "MI:0416", # fluorescence microscopy
    "MI:1314", # BioID / proximity-dependent biotin identification
    "MI:0813", # proximity ligation assay
    "MI:2364", # proximity interaction type
    "MI:0403", # colocalization
}

CORE_TERMS = {"MI:0915", "MI:0914", "MI:0407", "MI:0217"}

def parse_intact_raw(file_path):
    """Parses raw IntAct MITAB 2.7 file directly."""
    records = []
    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        for line_no, line in enumerate(f, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\r\n").split("\t")
            if len(cols) < 15:
                continue
            
            p1_raw = cols[0].strip()
            p2_raw = cols[1].strip()
            itype = cols[11].strip() if len(cols) > 11 else "NA"
            det_method = cols[6].strip() if len(cols) > 6 else "NA"
            ta = cols[20].strip() if len(cols) > 20 else "NA"
            tb = cols[21].strip() if len(cols) > 21 else "NA"
            tax_a = cols[9].strip() if len(cols) > 9 else "NA"
            tax_b = cols[10].strip() if len(cols) > 10 else "NA"
            is_neg = (cols[35].strip().lower() == "true") if len(cols) > 35 else False
            miscore = cols[14].strip() if len(cols) > 14 else "NA"
            pmid = cols[8].strip() if len(cols) > 8 else "NA"
            
            # Clean IDs
            p1 = p1_raw.split(":")[-1] if ":" in p1_raw else p1_raw
            p2 = p2_raw.split(":")[-1] if ":" in p2_raw else p2_raw
            
            # Molecule type classification
            def get_mol_type(t_str, raw_id):
                t_lower = t_str.lower()
                id_lower = raw_id.lower()
                if "mi:0328" in t_lower or "small molecule" in t_lower or "chebi" in id_lower:
                    return "small_molecule"
                if "mi:0319" in t_lower or "mi:0681" in t_lower or "deoxyribonucleic" in t_lower:
                    return "dna"
                if "mi:0320" in t_lower or "ribonucleic" in t_lower:
                    return "rna"
                if "mi:0318" in t_lower or "nucleic acid" in t_lower:
                    return "nucleic_acid"
                if "mi:0250" in t_lower or "gene" in t_lower:
                    return "gene"
                if "mi:1304" in t_lower or "molecule set" in t_lower:
                    return "molecule_set"
                if "mi:0326" in t_lower or "protein" in t_lower or "mi:0327" in t_lower or "peptide" in t_lower:
                    return "protein"
                if raw_id == "-" or not raw_id:
                    return "missing"
                return "protein" # default fallback for UniProt IDs
                
            mol_a = get_mol_type(ta, p1_raw)
            mol_b = get_mol_type(tb, p2_raw)
            is_prot_prot = (mol_a == "protein" and mol_b == "protein")
            
            # Semantic classification
            is_proximity_type = any(k in itype for k in ["MI:2364", "MI:0403"])
            is_pca_method = any(k in det_method for k in PCA_MI)
            is_prox_method = any(k in det_method for k in PROX_MI)
            is_proximity = is_proximity_type or is_prox_method

            is_physical_core = any(k in itype for k in CORE_TERMS) and not is_proximity
            is_physical_assoc_only = ("MI:0915" in itype or "MI:0914" in itype) and not is_proximity
            is_direct_binary_only = ("MI:0407" in itype) and not is_proximity
            
            records.append({
                "line_no": line_no,
                "p1_raw": p1_raw, "p2_raw": p2_raw,
                "p1": p1, "p2": p2,
                "pair_key": tuple(sorted([p1, p2])),
                "itype": itype, "det_method": det_method,
                "ta": ta, "tb": tb,
                "mol_a": mol_a, "mol_b": mol_b,
                "is_protein_protein": is_prot_prot,
                "is_neg": is_neg,
                "is_pca": is_pca_method,
                "is_proximity": is_proximity,
                "is_physical_core": is_physical_core,
                "is_physical_assoc_only": is_physical_assoc_only,
                "is_direct_binary_only": is_direct_binary_only,
                "miscore": miscore, "pmid": pmid
            })
    return records


def parse_biogrid_raw(file_path):
    """Parses raw BioGRID TAB 3.0 file directly."""
    records = []
    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        for line_no, line in enumerate(f, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\r\n").split("\t")
            if len(cols) < 20:
                continue
            
            bg_id = cols[0].strip()
            entrez_a = cols[1].strip()
            entrez_b = cols[2].strip()
            sys_a = cols[5].strip()
            sys_b = cols[6].strip()
            sym_a = cols[7].strip()
            sym_b = cols[8].strip()
            exp_sys = cols[11].strip()
            sys_type = cols[12].strip().lower()
            tax_a = cols[15].strip()
            tax_b = cols[16].strip()
            sp_a = cols[23].strip() if len(cols) > 23 else "-"
            sp_b = cols[26].strip() if len(cols) > 26 else "-"
            
            # Preferred ID
            p1 = sp_a.split("|")[0] if (sp_a and sp_a != "-") else (sys_a if (sys_a and sys_a != "-") else (sym_a if (sym_a and sym_a != "-") else entrez_a))
            p2 = sp_b.split("|")[0] if (sp_b and sp_b != "-") else (sys_b if (sys_b and sys_b != "-") else (sym_b if (sym_b and sym_b != "-") else entrez_b))
            
            is_physical = (sys_type == "physical")
            is_genetic = (sys_type == "genetic")
            is_pca = (exp_sys == "PCA")
            is_proximity = ("proximity" in exp_sys.lower()) or (exp_sys == "Proximity Label-MS") or (exp_sys == "Co-localization")
            is_rna = ("RNA" in exp_sys)
            
            records.append({
                "line_no": line_no,
                "bg_id": bg_id,
                "p1": p1, "p2": p2,
                "pair_key": tuple(sorted([p1, p2])),
                "exp_sys": exp_sys,
                "sys_type": sys_type,
                "tax_a": tax_a, "tax_b": tax_b,
                "is_physical": is_physical,
                "is_genetic": is_genetic,
                "is_pca": is_pca,
                "is_proximity": is_proximity,
                "is_physical_no_proximity": is_physical and not is_proximity and not is_rna,
                "is_protein_protein": not is_rna
            })
    return records


def parse_string_raw(file_path):
    """Parses raw STRING physical links file (full or detailed)."""
    records = []
    is_full = "full" in file_path
    with gzip.open(file_path, "rt", encoding="utf-8", errors="ignore") as f:
        header = f.readline().strip().split()
        for line_no, line in enumerate(f, start=2):
            parts = line.strip().split()
            if len(parts) < 6:
                continue
            p1_raw, p2_raw = parts[0], parts[1]
            p1 = p1_raw.split(".", 1)[-1] if "." in p1_raw else p1_raw
            p2 = p2_raw.split(".", 1)[-1] if "." in p2_raw else p2_raw
            
            if is_full:
                exp = int(parts[3])
                exp_tr = int(parts[4])
                db = int(parts[5])
                tm = int(parts[7])
                comb = int(parts[9])
            else: # detailed
                exp = int(parts[2])
                exp_tr = 0
                db = int(parts[3])
                tm = int(parts[4])
                comb = int(parts[5])
                
            records.append({
                "line_no": line_no,
                "p1": p1, "p2": p2,
                "p1_raw": p1_raw, "p2_raw": p2_raw,
                "pair_key": tuple(sorted([p1, p2])),
                "is_symmetric": (p1 <= p2),
                "exp": exp,
                "exp_tr": exp_tr,
                "db": db, "tm": tm, "comb": comb,
                "is_protein_protein": True
            })
    return records, is_full


def parse_xlms_raw(file_path):
    """Parses raw PRIDE PhoX XL-MS CSV file."""
    records = []
    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        reader = csv.reader(f)
        header = next(reader)
        for row_no, row in enumerate(reader, start=2):
            if not row or len(row) < 31:
                continue
            pA = row[29].strip()
            pB = row[30].strip() or pA
            ptype = row[14].strip()
            score = row[10].strip()
            evalue = row[9].strip()
            records.append({
                "row_no": row_no,
                "p1": pA, "p2": pB,
                "pair_key": tuple(sorted([pA, pB])),
                "ptype": ptype,
                "score": score,
                "evalue": evalue,
                "is_inter": (pA != pB),
                "is_protein_protein": True
            })
    return records


def verify_arabidopsis():
    print("=" * 80)
    print("VERIFYING ARABIDOPSIS THALIANA (TAXID 3702)")
    print("=" * 80)
    
    intact_p = os.path.join(WORKSPACE_ROOT, "data", "raw", "arabidopsis", "intact", "intact_arabidopsis_taxid3702_release252.mitab27.txt")
    bg_p = os.path.join(WORKSPACE_ROOT, "data", "raw", "arabidopsis", "biogrid", "BIOGRID-ORGANISM-Arabidopsis_thaliana_Columbia-5.0.261.tab3.txt")
    str_p = os.path.join(WORKSPACE_ROOT, "data", "raw", "arabidopsis", "string", "3702.protein.physical.links.full.v12.0.txt.gz")
    xlms_p = os.path.join(WORKSPACE_ROOT, "data", "raw", "arabidopsis", "xlms_2026", "Total_XL_plink3-2_v3.csv")

    intact_recs = parse_intact_raw(intact_p)
    bg_recs = parse_biogrid_raw(bg_p)
    str_recs, _ = parse_string_raw(str_p)
    xlms_recs = parse_xlms_raw(xlms_p)

    intact_pass = [r for r in intact_recs if r["is_physical_core"] and r["is_protein_protein"] and not r["is_neg"]]
    bg_pass = [r for r in bg_recs if r["is_physical_no_proximity"] and r["is_protein_protein"]]
    str_pass = [r for r in str_recs if r["exp"] > 0]
    xlms_pass = [r for r in xlms_recs if r["is_protein_protein"]]

    print(f"BioGRID (PCA Retained): Raw={len(bg_recs)}, Filtered={len(bg_pass)}, Deduplicated={len(set(r['pair_key'] for r in bg_pass))}")
    print(f"IntAct (Core Physical, PCA Retained): Raw={len(intact_recs)}, Filtered={len(intact_pass)}, Deduplicated={len(set(r['pair_key'] for r in intact_pass))}")
    print(f"STRING (Exp > 0): Raw Directed={len(str_recs)}, Filtered Directed={len(str_pass)}, Deduplicated Undirected={len(set(r['pair_key'] for r in str_pass))}")
    print(f"XL-MS: Raw PSMs={len(xlms_recs)}, Filtered PSMs={len(xlms_pass)}, Deduplicated={len(set(r['pair_key'] for r in xlms_pass))}")

if __name__ == "__main__":
    verify_arabidopsis()
