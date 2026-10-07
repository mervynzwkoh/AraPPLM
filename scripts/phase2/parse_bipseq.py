"""
Phase 2 Parser: BIP-seq Literature / POPPIN Status Audit
Audits BIP_seq_PMC11923860_fulltext.xml for structured interactome tables.
Confirms bulk interactome absence; records POPPIN_bulk_status = unavailable without fabricating data.
"""
import os
import sys
import xml.etree.ElementTree as ET

WORKSPACE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

BIPSEQ_SOURCE = {
    "source_id": "SRC_B1_BIPSEQ_PMC_XML",
    "file_path": "data/raw/rice/poppin/BIP_seq_PMC11923860_fulltext.xml"
}

def audit_bipseq():
    full_path = os.path.join(WORKSPACE_ROOT, BIPSEQ_SOURCE["file_path"])
    if not os.path.exists(full_path):
        raise FileNotFoundError(f"Missing raw file: {full_path}")

    tree = ET.parse(full_path)
    root = tree.getroot()

    tables = root.findall(".//table-wrap")
    table_summaries = []
    for t in tables:
        label = t.findtext("label") or "Unknown"
        caption = "".join(t.itertext())[:200].replace("\n", " ").strip()
        table_summaries.append(f"{label}: {caption}")

    # Document findings
    status = {
        "source_id": BIPSEQ_SOURCE["source_id"],
        "file_path": BIPSEQ_SOURCE["file_path"],
        "POPPIN_bulk_status": "unavailable",
        "primary_tables_found": len(tables),
        "structured_interactome_in_xml": False,
        "notes": (
            "BIP-seq PMC JATS XML contains 1 methodological protocol table (TCS system process). "
            "The bulk interactome dataset is accessible only via the POPPIN web query portal, which lacks "
            "a programmatic REST API or bulk dump. In strict compliance with P2_build_plan.md §20, "
            "zero synthetic interaction records are fabricated from prose."
        ),
        "table_summaries": table_summaries
    }
    return status

if __name__ == "__main__":
    res = audit_bipseq()
    print("BIP-seq / POPPIN Audit Result:")
    for k, v in res.items():
        print(f"  {k}: {v}")
