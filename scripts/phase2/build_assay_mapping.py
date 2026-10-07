"""
Build Phase 2 Controlled Assay Mapping Table.
File: data/interim/phase2/schemas/assay_mapping.tsv
"""
import os

MAPPINGS = [
    # IntAct PSI-MI terms: (psi_mi_accession, psi_mi_name, assay_family, direct, physical, proximity, native_in_planta, rationale)
    ("MI:0018", "two hybrid", "Y2H", True, True, False, "FALSE", "Classic Gal4/LexA yeast two-hybrid binary interaction assay; heterologous host"),
    ("MI:0397", "two hybrid array", "Y2H", True, True, False, "FALSE", "High-throughput matrix/array Y2H binary screen"),
    ("MI:1112", "two hybrid prey pooling approach", "Y2H", True, True, False, "FALSE", "Pooled prey Y2H binary screen"),
    ("MI:1356", "validated two hybrid", "Y2H", True, True, False, "FALSE", "Retested/independently validated Y2H assay"),
    ("MI:2277", "Cr-two hybrid", "Y2H", True, True, False, "FALSE", "Chlamydomonas or chloroplast-specific two-hybrid"),
    ("MI:0726", "reverse two hybrid", "Y2H", True, True, False, "FALSE", "Reverse two-hybrid for interaction disruption/validation"),
    ("MI:0112", "ubiquitin reconstruction", "split_ubiquitin", True, True, False, "FALSE", "Split-ubiquitin membrane yeast two-hybrid (SUS/MYTH)"),
    ("MI:0006", "anti bait coimmunoprecipitation", "coIP", False, True, False, "TRUE", "Immunoprecipitation of endogenous/tagged bait to detect associated prey"),
    ("MI:0007", "anti tag coimmunoprecipitation", "coIP", False, True, False, "TRUE", "Tag-directed co-immunoprecipitation of complex components"),
    ("MI:0019", "coimmunoprecipitation", "coIP", False, True, False, "TRUE", "Co-immunoprecipitation from cell/tissue lysate"),
    ("MI:0096", "pull down", "pull_down", True, True, False, "FALSE", "In vitro or cell-free affinity pull-down assay (GST, His, MBP, etc.)"),
    ("MI:0676", "tandem affinity purification", "AP_MS", False, True, False, "TRUE", "TAP-tag tandem affinity purification followed by MS identification"),
    ("MI:0400", "affinity technology", "AP_MS", False, True, False, "unresolved", "General affinity capture / chromatography technology"),
    ("MI:0004", "affinity chromatography technology", "AP_MS", False, True, False, "unresolved", "Affinity chromatography purification of protein complexes"),
    ("MI:0809", "bimolecular fluorescence complementation", "BiFC", False, True, True, "TRUE", "Split-FP complementation in planta; irreversibility can stabilize weak/proximity associations"),
    ("MI:2365", "fluorescent protein-protein interaction-visualization", "BiFC", False, True, True, "TRUE", "Fluorophore complementation or direct visual interaction assay"),
    ("MI:0090", "protein complementation assay", "BiFC", False, True, True, "unresolved", "General protein fragment complementation assay (PCA)"),
    ("MI:1037", "Split renilla luciferase complementation", "split_luciferase", False, True, True, "TRUE", "Split Renilla luciferase complementation assay"),
    ("MI:1204", "split firefly luciferase complementation", "split_luciferase", False, True, True, "TRUE", "Split Firefly luciferase complementation assay"),
    ("MI:0055", "fluorescent resonance energy transfer", "biophysical_binding", True, True, True, "TRUE", "FRET / FLIM-FRET distance-dependent energy transfer assay"),
    ("MI:0051", "fluorescence technology", "biophysical_binding", True, True, False, "unresolved", "Biophysical fluorescence spectroscopy or polarization"),
    ("MI:0017", "classical fluorescence spectroscopy", "biophysical_binding", True, True, False, "FALSE", "In vitro fluorescence emission/quenching assay"),
    ("MI:0107", "surface plasmon resonance", "biophysical_binding", True, True, False, "FALSE", "In vitro biophysical surface plasmon resonance (Biacore) direct kinetic binding"),
    ("MI:0065", "isothermal titration calorimetry", "biophysical_binding", True, True, False, "FALSE", "Direct thermodynamic measurement of binding stoichiometry and affinity"),
    ("MI:0016", "circular dichroism", "biophysical_binding", True, True, False, "FALSE", "Secondary structure alteration upon binding"),
    ("MI:0067", "light scattering", "biophysical_binding", True, True, False, "FALSE", "Static or dynamic light scattering for complex formation"),
    ("MI:0077", "nuclear magnetic resonance", "biophysical_binding", True, True, False, "FALSE", "NMR chemical shift perturbation or structural interaction"),
    ("MI:0114", "x-ray crystallography", "biophysical_binding", True, True, False, "FALSE", "Atomic co-crystal structure determination"),
    ("MI:0040", "electron microscopy", "biophysical_binding", False, True, False, "FALSE", "Cryo-EM or negative-stain electron microscopy of complex"),
    ("MI:0826", "x ray scattering", "biophysical_binding", False, True, False, "FALSE", "Small-angle X-ray scattering (SAXS) of complex"),
    ("MI:0030", "cross-linking study", "XL_MS", True, True, False, "unresolved", "Covalent cross-linking of interacting proteins"),
    ("MI:0031", "protein cross-linking with a bifunctional reagent", "XL_MS", True, True, False, "unresolved", "Bifunctional chemical cross-linking (DSS, BS3, etc.)"),
    ("MI:1314", "proximity-dependent biotin identification", "proximity_labeling", False, False, True, "TRUE", "BioID / TurboID / APEX proximity-dependent biotinylation"),
    ("MI:0813", "proximity ligation assay", "proximity_labeling", False, False, True, "TRUE", "In situ proximity ligation assay (PLA, <40 nm distance)"),
    ("MI:0089", "protein array", "protein_microarray", True, True, False, "FALSE", "Recombinant protein microarray / chip binding screen"),
    ("MI:0084", "phage display", "biophysical_binding", True, True, False, "FALSE", "Phage display selection of interacting polypeptides"),
    ("MI:0066", "lambda phage display", "biophysical_binding", True, True, False, "FALSE", "Lambda phage display library screening"),
    ("MI:0047", "far western blotting", "biophysical_binding", True, True, False, "FALSE", "In vitro Far-Western blotting on immobilized target protein"),
    ("MI:0049", "filter binding", "biophysical_binding", True, True, False, "FALSE", "Nitrocellulose or membrane filter retention assay"),
    ("MI:0097", "reverse ras recruitment system", "other", True, True, False, "FALSE", "Yeast cytoplasmic Ras recruitment system (RRS)"),
    ("MI:0276", "blue native page", "biophysical_binding", False, True, False, "unresolved", "Blue Native polyacrylamide gel electrophoresis of intact complexes"),
    ("MI:0404", "comigration in non denaturing gel electrophoresis", "biophysical_binding", False, True, False, "unresolved", "Native PAGE co-migration of binding partners"),
    ("MI:0807", "comigration in gel electrophoresis", "biophysical_binding", False, True, False, "unresolved", "Electrophoretic co-migration of complexes"),
    ("MI:0028", "cosedimentation in solution", "biophysical_binding", False, True, False, "unresolved", "Analytical ultracentrifugation / sedimentation in solution"),
    ("MI:0029", "cosedimentation through density gradient", "biophysical_binding", False, True, False, "unresolved", "Sucrose or glycerol gradient density co-sedimentation"),
    ("MI:0071", "molecular sieving", "biophysical_binding", False, True, False, "unresolved", "Size exclusion chromatography / gel filtration"),
    ("MI:0401", "biochemical", "biophysical_binding", False, True, False, "unresolved", "In vitro enzymatic or biochemical binding assay"),
    ("MI:0013", "biophysical", "biophysical_binding", True, True, False, "unresolved", "General biophysical binding assay"),
    ("MI:0405", "competition binding", "biophysical_binding", True, True, False, "FALSE", "In vitro competitive ligand/protein binding displacement"),
    ("MI:0411", "enzyme linked immunosorbent assay", "biophysical_binding", True, True, False, "FALSE", "ELISA-based direct interaction binding assay"),
    ("MI:0412", "electrophoretic mobility supershift assay", "biophysical_binding", False, True, False, "FALSE", "EMSA antibody supershift proving complex involvement"),
    ("MI:0413", "electrophoretic mobility shift assay", "biophysical_binding", False, True, False, "FALSE", "Gel shift assay demonstrating association"),
    ("MI:0415", "enzymatic study", "biophysical_binding", False, True, False, "unresolved", "Enzymatic substrate or regulator binding assay"),
    ("MI:0416", "fluorescence microscopy", "proximity_labeling", False, False, True, "TRUE", "Subcellular co-localization by fluorescence microscopy (proximity only)"),
    ("MI:0663", "confocal microscopy", "proximity_labeling", False, False, True, "TRUE", "Subcellular co-localization by confocal microscopy (proximity only)"),
    ("MI:0419", "gtpase assay", "biophysical_binding", False, True, False, "FALSE", "GTPase activating or nucleotide exchange protein interaction"),
    ("MI:0423", "in-gel kinase assay", "biophysical_binding", False, True, False, "FALSE", "In-gel phosphorylation assay"),
    ("MI:0424", "protein kinase assay", "biophysical_binding", False, True, False, "FALSE", "In vitro kinase-substrate phosphorylation interaction"),
    ("MI:0434", "phosphatase assay", "biophysical_binding", False, True, False, "FALSE", "Phosphatase-substrate dephosphorylation interaction"),
    ("MI:0435", "protease assay", "biophysical_binding", False, True, False, "FALSE", "Protease-substrate cleavage interaction"),
    ("MI:0515", "methyltransferase assay", "biophysical_binding", False, True, False, "FALSE", "Methyltransferase enzyme-substrate interaction"),
    ("MI:0841", "phosphotransferase assay", "biophysical_binding", False, True, False, "FALSE", "Phosphotransferase reaction monitoring interaction"),
    ("MI:0880", "atpase assay", "biophysical_binding", False, True, False, "FALSE", "ATPase stimulation/inhibition by interacting factor"),
    ("MI:0440", "saturation binding", "biophysical_binding", True, True, False, "FALSE", "Equilibrium saturation binding measurement"),
    ("MI:0892", "solid phase assay", "biophysical_binding", True, True, False, "FALSE", "Solid phase immobilized binding assay"),
    ("MI:0437", "protein three hybrid", "other", False, True, False, "FALSE", "Yeast three-hybrid mediated by ternary bridge"),
    ("MI:0588", "three hybrid", "other", False, True, False, "FALSE", "Three-hybrid system for ternary complex"),
    ("MI:0432", "one hybrid", "other", False, False, False, "FALSE", "Yeast one-hybrid (protein-DNA, indirect PPI)"),
    ("MI:0606", "DNase I footprinting", "other", False, False, False, "FALSE", "DNA protection footprinting assay"),
    ("MI:0402", "chromatin immunoprecipitation assay", "other", False, True, False, "TRUE", "ChIP assay showing chromatin co-occupancy"),
    ("MI:1017", "rna immunoprecipitation", "other", False, True, False, "TRUE", "RIP-seq/RIP identifying RNA-associated protein complex"),
    ("MI:0001", "interaction detection method", "unknown", False, False, False, "unresolved", "Unspecified root interaction detection method term"),

    # BioGRID Experimental Systems: (source_term, system_type, assay_family, direct, physical, proximity, native_in_planta, rationale)
    ("Two-hybrid", "physical", "Y2H", True, True, False, "FALSE", "BioGRID Two-hybrid physical interaction screen"),
    ("Affinity Capture-Western", "physical", "coIP", False, True, False, "TRUE", "BioGRID Affinity Capture with Western blot detection (co-IP/pull-down)"),
    ("Affinity Capture-MS", "physical", "AP_MS", False, True, False, "TRUE", "BioGRID Affinity Capture with mass spectrometry identification (AP-MS)"),
    ("Affinity Capture-Luminescence", "physical", "coIP", False, True, False, "TRUE", "BioGRID Affinity Capture detected by luminescence/reporter"),
    ("Affinity Capture-RNA", "physical", "other", False, True, False, "TRUE", "BioGRID Affinity Capture of ribonucleoprotein complex"),
    ("Co-fractionation", "physical", "AP_MS", False, True, False, "unresolved", "BioGRID biochemical co-fractionation across chromatographic fractions"),
    ("Co-localization", "physical", "proximity_labeling", False, False, True, "TRUE", "BioGRID subcellular co-localization by microscopy (proximity only)"),
    ("Co-purification", "physical", "AP_MS", False, True, False, "unresolved", "BioGRID multi-step co-purification of native complex"),
    ("Co-crystal Structure", "physical", "biophysical_binding", True, True, False, "FALSE", "BioGRID atomic resolution co-crystal structure"),
    ("Cross-Linking-MS (XL-MS)", "physical", "XL_MS", True, True, False, "unresolved", "BioGRID chemical cross-linking coupled with mass spectrometry"),
    ("Proximity Label-MS", "physical", "proximity_labeling", False, False, True, "TRUE", "BioGRID proximity-dependent labeling (BioID, APEX, TurboID)"),
    ("Reconstituted Complex", "physical", "biophysical_binding", True, True, False, "FALSE", "BioGRID in vitro reconstitution of purified recombinant components"),
    ("Biochemical Activity", "physical", "biophysical_binding", False, True, False, "unresolved", "BioGRID in vitro enzymatic activity alteration or reaction"),
    ("Far Western", "physical", "biophysical_binding", True, True, False, "FALSE", "BioGRID in vitro Far-Western overlay binding assay"),
    ("FRET", "physical", "biophysical_binding", True, True, True, "TRUE", "BioGRID Förster resonance energy transfer in living cells"),
    ("PCA", "physical", "BiFC", False, True, True, "unresolved", "BioGRID protein complementation assay (BiFC, split-luc, DHFR PCA)"),
    ("Protein-peptide", "physical", "biophysical_binding", True, True, False, "FALSE", "BioGRID direct binding of protein to synthetic peptide"),
    ("Protein-RNA", "physical", "other", False, True, False, "unresolved", "BioGRID protein-RNA mediated interaction"),
    ("Dosage Growth Defect", "genetic", "genetic", False, False, False, "TRUE", "BioGRID dosage lethality/growth defect genetic interaction"),
    ("Dosage Rescue", "genetic", "genetic", False, False, False, "TRUE", "BioGRID high-copy suppressor dosage rescue genetic interaction"),
    ("Phenotypic Enhancement", "genetic", "genetic", False, False, False, "TRUE", "BioGRID double mutant phenotypic enhancement genetic interaction"),
    ("Phenotypic Suppression", "genetic", "genetic", False, False, False, "TRUE", "BioGRID suppressor mutant phenotypic suppression genetic interaction"),
    ("Synthetic Growth Defect", "genetic", "genetic", False, False, False, "TRUE", "BioGRID synthetic sickness double mutant interaction"),
    ("Synthetic Lethality", "genetic", "genetic", False, False, False, "TRUE", "BioGRID synthetic lethality genetic interaction"),
    ("Synthetic Rescue", "genetic", "genetic", False, False, False, "TRUE", "BioGRID synthetic rescue genetic interaction"),

    # STRING, XL-MS, BIP-seq
    ("string_experimental", "computational", "computational", False, False, False, "unresolved", "STRING experimental channel (mix of native and transferred evidence)"),
    ("string_database", "computational", "computational", False, False, False, "unresolved", "STRING curated database channel"),
    ("string_textmining", "computational", "computational", False, False, False, "unresolved", "STRING automated literature text-mining channel"),
    ("plink_search_table", "XL_MS", "XL_MS", True, True, False, "TRUE", "Processed pLink 3.2 identified cross-link peptide spectra (5% FDR)"),
    ("bip_seq_article", "other", "other", False, False, False, "TRUE", "Primary literature XML for BIP-seq interactome screen")
]

def main():
    out_dir = os.path.join("data", "interim", "phase2", "schemas")
    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, "assay_mapping.tsv")

    header = [
        "source_term",
        "psi_mi_accession",
        "psi_mi_name",
        "normalized_assay_family",
        "supports_direct_binding",
        "supports_physical_association",
        "supports_proximity_only",
        "native_in_planta",
        "mapping_rationale",
        "mapping_version"
    ]

    with open(out_file, "w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for item in MAPPINGS:
            if item[0].startswith("MI:"):
                acc, name, fam, direct, phys, prox, nat, rat = item
                src = f'psi-mi:"{acc}"({name})'
                row = [src, acc, name, fam, str(direct).upper(), str(phys).upper(), str(prox).upper(), nat, rat, "v1.0"]
            else:
                src, typ, fam, direct, phys, prox, nat, rat = item
                acc = "NONE"
                name = src
                row = [src, acc, name, fam, str(direct).upper(), str(phys).upper(), str(prox).upper(), nat, rat, "v1.0"]
            f.write("\t".join(row) + "\n")

    print(f"[OK] Wrote {out_file} with {len(MAPPINGS)} controlled assay mappings.")

if __name__ == "__main__":
    main()
