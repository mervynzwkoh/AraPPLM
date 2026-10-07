"""
Phase 2: Harmonize Experimental Assays and Interaction Semantics.
Uses data/interim/phase2/schemas/assay_mapping.tsv as the authoritative mapping table.
"""
import os
import re

WORKSPACE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
ASSAY_MAPPING_PATH = os.path.join(WORKSPACE_ROOT, "data", "interim", "phase2", "schemas", "assay_mapping.tsv")

class AssayHarmonizer:
    def __init__(self, mapping_path=ASSAY_MAPPING_PATH):
        self.mapping_path = mapping_path
        self.psi_mi_map = {}   # keyed by MI:XXXX
        self.term_map = {}     # keyed by normalized source_term string
        self.unmapped_terms = set()
        self._load_mapping()

    def _load_mapping(self):
        if not os.path.exists(self.mapping_path):
            raise FileNotFoundError(f"Assay mapping table not found: {self.mapping_path}")

        with open(self.mapping_path, "r", encoding="utf-8") as f:
            header = f.readline().strip().split("\t")
            for line in f:
                if not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 10:
                    continue
                source_term = parts[0]
                psi_mi_acc = parts[1]
                psi_mi_name = parts[2]
                assay_fam = parts[3]
                direct = parts[4].upper() == "TRUE"
                phys = parts[5].upper() == "TRUE"
                prox = parts[6].upper() == "TRUE"
                native = parts[7]
                rationale = parts[8]

                entry = {
                    "psi_mi_accession": psi_mi_acc if psi_mi_acc != "NONE" else "NA",
                    "psi_mi_name": psi_mi_name if psi_mi_name != "NONE" else "NA",
                    "assay_family": assay_fam,
                    "supports_direct_binding": direct,
                    "supports_physical_association": phys,
                    "supports_proximity_only": prox,
                    "assay_is_genetic": (assay_fam == "genetic"),
                    "assay_is_computational": (assay_fam == "computational"),
                    "native_in_planta": native,
                    "mapping_rationale": rationale
                }

                if psi_mi_acc.startswith("MI:"):
                    self.psi_mi_map[psi_mi_acc] = entry
                self.term_map[source_term.strip().lower()] = entry
                self.term_map[psi_mi_name.strip().lower()] = entry

    def harmonize_psi_mi(self, raw_method_str):
        """
        Parses PSI-MI method string e.g. psi-mi:"MI:0018"(two hybrid)
        """
        if not raw_method_str or raw_method_str == "-":
            return self._unknown_entry(raw_method_str, "empty_or_missing_method")

        # Extract MI:XXXX
        m = re.search(r'MI:\d{4}', raw_method_str)
        if m:
            mi_acc = m.group(0)
            if mi_acc in self.psi_mi_map:
                res = dict(self.psi_mi_map[mi_acc])
                res["detection_method_raw"] = raw_method_str
                res["detection_method_psi_mi"] = mi_acc
                res["detection_method_psi_mi_name"] = res["psi_mi_name"]
                res["needs_manual_review"] = False
                res["review_reason"] = "NA"
                return res

        # Try term lookup
        norm = raw_method_str.strip().lower()
        if norm in self.term_map:
            res = dict(self.term_map[norm])
            res["detection_method_raw"] = raw_method_str
            res["detection_method_psi_mi"] = res["psi_mi_accession"]
            res["detection_method_psi_mi_name"] = res["psi_mi_name"]
            res["needs_manual_review"] = False
            res["review_reason"] = "NA"
            return res

        self.unmapped_terms.add(raw_method_str)
        return self._unknown_entry(raw_method_str, f"unmapped_psi_mi_method:{raw_method_str}")

    def harmonize_biogrid_system(self, system_name, system_type):
        """
        Maps BioGRID Experimental System and Type
        """
        norm = system_name.strip().lower()
        if norm in self.term_map:
            res = dict(self.term_map[norm])
            res["detection_method_raw"] = system_name
            res["detection_method_psi_mi"] = res["psi_mi_accession"]
            res["detection_method_psi_mi_name"] = res["psi_mi_name"]
            res["needs_manual_review"] = False
            res["review_reason"] = "NA"
            return res

        self.unmapped_terms.add(f"{system_name} ({system_type})")
        return self._unknown_entry(system_name, f"unmapped_biogrid_system:{system_name}")

    def _unknown_entry(self, raw_str, reason):
        return {
            "detection_method_raw": raw_str if raw_str else "NA",
            "detection_method_psi_mi": "NA",
            "detection_method_psi_mi_name": "NA",
            "assay_family": "unknown",
            "supports_direct_binding": False,
            "supports_physical_association": False,
            "supports_proximity_only": False,
            "assay_is_genetic": False,
            "assay_is_computational": False,
            "native_in_planta": "unresolved",
            "mapping_rationale": "No controlled mapping match found",
            "needs_manual_review": True,
            "review_reason": reason
        }

    def infer_interaction_semantics(self, assay_info, raw_type_str="", is_genetic=False):
        """
        Provisional interaction semantics: direct_binary, physical_association,
        co_complex, proximity, genetic, computational, or unknown.
        """
        if is_genetic or assay_info.get("assay_is_genetic", False):
            return "genetic"
        if assay_info.get("assay_is_computational", False):
            return "computational"

        # Check PSI-MI interaction type if available
        if "MI:0407" in raw_type_str: # direct interaction
            return "direct_binary"
        if "MI:0914" in raw_type_str: # association
            return "physical_association"
        if "MI:0915" in raw_type_str: # physical association
            return "physical_association"
        if "MI:0462" in raw_type_str: # enzymatic reaction
            return "functional"

        # Fallback to assay support
        if assay_info.get("supports_direct_binding", False):
            return "direct_binary"
        if assay_info.get("supports_proximity_only", False):
            return "proximity"
        if assay_info.get("supports_physical_association", False):
            return "physical_association"

        return "unknown"
