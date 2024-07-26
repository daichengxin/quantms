#!/usr/bin/env python

import pandas as pd
from lxml import etree
import pyopenms as oms
import re
import os
import sys

mass_proton = 1.007276


def idXML2tsv(idXML):
    prot_ids = []
    pep_ids = []
    oms.IdXMLFile().load(idXML, prot_ids, pep_ids)
    # basename = os.path.splitext(os.path.basename(mzML))[0]
    oms.PepXMLFile().store("test.out.pepxml", prot_ids, pep_ids)

    # for peptide_id in pep_ids:
    #     retention_time = peptide_id.getRT()
    #     exp_mass_to_charge = peptide_id.getMZ()
    #     scan_number = int(re.findall(r"(spectrum|scan)=(\d+)", peptide_id.getMetaValue("spectrum_reference"))[0][1])
    #     for hit in peptide_id.getHits():
    #         charge = hit.getCharge()
    #         spectrum = ".".join(["test", str(scan_number), str(scan_number), str(charge)])
    #         peptide = hit.getSequence().toUnmodifiedString()
    #         pass


def pepXML2tsv(pepXML, idXML):
    tsv_df = pd.DataFrame()
    prot_ids = []
    pep_ids = []
    oms.IdXMLFile().load(idXML, prot_ids, pep_ids)
    # basename = os.path.splitext(os.path.basename(mzML))[0]
    # oms.PepXMLFile().store("test.out.pepxml", prot_ids, pep_ids)

    tree = etree.parse(pepXML, parser=etree.XMLParser())
    root_el = tree.getroot()
    idx = 0
    # ed = oms.EnzymaticDigestion()
    # enzyme = root_el.xpath(".//*[local-name()='sample_enzyme']")[0]
    # enzyme = enzyme.xpath("@name")[0]
    # print(enzyme)
    # dig = oms.ProteaseDigestion()
    # print(dig.getEnzymeName())  # Trypsin'
    #
    # de = oms.DigestionEnzyme()
    # de.setRegEx()
    # de.setName(enzyme)
    # ed.setEnzyme(de)
    print(len(root_el.xpath(".//*[local-name()='spectrum_query']")))
    records = []
    for spectrum_query in root_el.xpath(".//*[local-name()='spectrum_query']"):
        spectrum = spectrum_query.xpath("@start_scan")[0]
        # print(spectrum)
        charge = spectrum_query.xpath("@assumed_charge")[0]
        RetentionTime = spectrum_query.xpath("@retention_time_sec")[0]
        ObservedMZ = pep_ids[idx].getMZ()
        precursor_neutral_mass = pep_ids[idx].getMZ() * int(charge) - int(charge) * mass_proton
        # print(pep_ids[idx].getHits()[0].getSequence().getAverageWeight())
        spectrum_query.set("precursor_neutral_mass", str(precursor_neutral_mass))
        Hyperscore = pep_ids[idx].getHits()[0].getMetaValue("hyperscore")
        Probability = 1 - pep_ids[idx].getHits()[0].getMetaValue("MS:1001493")
        Peptide = pep_ids[idx].getHits()[0].getSequence().toUnmodifiedString()
        PeptideLength = len(Peptide)
        hit = spectrum_query.xpath(".//*[local-name()='search_hit']")[0]
        # for hit in spectrum_query.xpath(".//*[local-name()='search_hit']"):
        calc_neutral_pep_mass = float(hit.xpath("@calc_neutral_pep_mass")[0])
        CalculatedMZ = pep_ids[idx].getHits()[0].getSequence().getMZ(int(charge))
        massdiff = precursor_neutral_mass - calc_neutral_pep_mass
        hit.set("massdiff", str(massdiff))
        num_matched_ions = pep_ids[idx].getHits()[0].getMetaValue("SAGE:matched_peaks")
        hit.set("num_matched_ions", num_matched_ions)
        modification_infos = hit.xpath(".//*[local-name()='modification_info']")
        Modifications = []
        if len(modification_infos) > 0:
            modified_peptide = modification_infos[0].xpath("@modified_peptide")[0]
            ifnterm = modification_infos[0].xpath("@mod_nterm_mass")
            if len(ifnterm) > 0:
                Modifications.append("N-term(" + str(ifnterm) + ")")
            for mod in modification_infos[0].xpath(".//*[local-name()='mod_aminoacid_mass']"):
                AssignedModifications = mod.xpath("@position")[0] + Peptide[
                    int(mod.xpath("@position")[0]) - 1] + "(" + mod.xpath("@mass")[0] + ")"
                Modifications.append(AssignedModifications)
        else:
            modified_peptide = ""
        PrevAA = hit.xpath("@peptide_prev_aa")[0]
        NextAA = hit.xpath("@peptide_next_aa")[0]
        NumberofEnzymaticTermini = hit.xpath("@num_tol_term")[0]  # 这两个要重新计算
        NumberofMissedCleavages = 0  # ed.countInternalCleavageSites(Peptide)  # 这两个要重新计算
        ProteinStart = pep_ids[idx].getHits()[0].getMetaValue("start")
        ProteinEnd = pep_ids[idx].getHits()[0].getMetaValue("end")
        Protein = hit.xpath("@protein")[0]
        ProteinID = Protein.split("|")[1]
        EntryName = Protein.split("|")[2]

        Spectrum = spectrum_query.xpath("@spectrum")[0]
        SpectrumFile = pepXML
        AssignedModifications = ", ".join(Modifications) if len(Modifications) > 0 else None
        IsUnique = "TRUE" if "unique" == pep_ids[idx].getHits()[0].getMetaValue("protein_references") else "FALSE"

        records.append({"Spectrum": Spectrum, "Spectrum File": SpectrumFile, "Peptide": Peptide,
                        "Modified Peptide": modified_peptide, "Prev AA": PrevAA, "Next AA": NextAA,
                        "Peptide Length": PeptideLength, "Charge": charge,
                        "Retention": RetentionTime,
                        "Observed Mass": precursor_neutral_mass, "Observed M/Z": ObservedMZ,
                        "Calculated Peptide Mass": calc_neutral_pep_mass,
                        "Calculated M/Z": CalculatedMZ,
                        "Delta Mass": massdiff, "Hyperscore": Hyperscore,
                        "Probability": Probability,
                        "Number of Enzymatic Termini": 2,
                        "Number of Missed Cleavages": NumberofMissedCleavages,
                        "Protein Start": ProteinStart, "Protein End": ProteinEnd, "Intensity": 0,
                        "Assigned Modifications": AssignedModifications, "Is Unique": IsUnique,
                        "Protein": Protein, "Protein ID": ProteinID, "Entry Name": EntryName})

        idx += 1

    tsv_df = pd.concat(
        [tsv_df, pd.DataFrame.from_records(records)],
        ignore_index=True)
    tsv_df.to_csv("psm.tsv", index=False, sep="\t")
    tree.write("corrected.pepXML")


def corrected_pepXML(pepXML, idXML):
    prot_ids = []
    pep_ids = []
    oms.IdXMLFile().load(idXML, prot_ids, pep_ids)
    # basename = os.path.splitext(os.path.basename(mzML))[0]
    # oms.PepXMLFile().store("test.out.pepxml", prot_ids, pep_ids)

    tree = etree.parse(pepXML, parser=etree.XMLParser())
    root_el = tree.getroot()
    idx = 0
    for spectrum_query in root_el.xpath(".//*[local-name()='spectrum_query']"):
        spectrum = spectrum_query.xpath("@start_scan")[0]
        if spectrum != "52037":
            idx += 1
            continue
        charge = spectrum_query.xpath("@assumed_charge")[0]
        precursor_neutral_mass = pep_ids[idx].getMZ() * int(charge) - int(charge) * mass_proton
        print(pep_ids[idx].getHits()[0].getSequence().getAverageWeight())
        spectrum_query.set("precursor_neutral_mass", str(precursor_neutral_mass))
        for hit in spectrum_query.xpath(".//*[local-name()='search_hit']"):
            massdiff = precursor_neutral_mass - float(hit.xpath("@calc_neutral_pep_mass")[0])
            hit.set("massdiff", str(massdiff))
            num_matched_ions = pep_ids[idx].getHits()[0].getMetaValue("SAGE:matched_peaks")
            hit.set("num_matched_ions", num_matched_ions)

        idx += 1

    # tree.write("D:/bigdata/sage_opensearch/test.pepXML")


def idXML2PepXML(idXML):
    prot_ids = []
    pep_ids = []
    oms.IdXMLFile().load(idXML, prot_ids, pep_ids)
    for peptide_id in pep_ids:
        retention_time = peptide_id.getRT()
        exp_mass_to_charge = peptide_id.getMZ()
        scan_number = int(re.findall(r"(spectrum|scan)=(\d+)", peptide_id.getMetaValue("spectrum_reference"))[0][1])
        for hit in peptide_id.getHits():
            charge = hit.getCharge()
            spectrum = ".".join(["test", str(scan_number), str(scan_number), str(charge)])
            peptide = hit.getSequence().toUnmodifiedString()
            pass


def main():
    pepxml_path = sys.argv[1]
    idxml_path = sys.argv[2]
    pepXML2tsv(pepxml_path, idxml_path)


if __name__ == "__main__":
    sys.exit(main())

# if __name__ == "__main__":
#     # corrected_pepXML("D:/bigdata/sage_opensearch/b1929_293T_proteinID_09A_QE3_122212_sage_perc_pep_filter.pepXML",
#     #                  "D:/bigdata/sage_opensearch/b1929_293T_proteinID_09A_QE3_122212_sage_perc_pep_filter.idXML")
#     pepXML2tsv("D:/bigdata/sage_opensearch/b1929_293T_proteinID_09A_QE3_122212_sage_perc_pep_filter.pepXML",
#                "D:/bigdata/sage_opensearch/b1929_293T_proteinID_09A_QE3_122212_sage_perc_pep_filter.idXML")
#     # idXML2PepXML("D:/bigdata/sage_opensearch/b1929_293T_proteinID_09A_QE3_122212_sage_perc_pep_filter.idXML")
