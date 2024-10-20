#!/usr/bin/env python

import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyopenms as oms

_parquet_field = [
    "sequence",
    "protein_accessions",
    "protein_start_positions",
    "protein_end_positions",
    "modifications",
    "mass_offset_proforma",
    "retention_time",
    "charge",
    "exp_mass_to_charge",
    "cal_mass_to_charge",
    "collision_energy",
    "reference_file_name",
    "scan_number",
    "USI",
    "peptidoform",
    "posterior_error_probability",
    "global_qvalue",
    "is_decoy",
    "consensus_support",
    "mz_array",
    "intensity_array",
    "num_peaks",
    "search_engines",
    "id_scores",
    "hit_rank",
]


def mods_position(sequence, peptide, mass_offset_dict):
    mass_offset_proforma = ""
    if ".(Glu->pyro-Glu)" in peptide:
        peptide = peptide.replace(".(Glu->pyro-Glu)", "")[0] + "(Glu->pyro-Glu)" + peptide.replace(".(Glu->pyro-Glu)", "")[1:]
    if ".(Gln->pyro-Glu)" in peptide:
        peptide = peptide.replace(".(Gln->pyro-Glu)", "")[0] + "(Gln->pyro-Glu)" + peptide.replace(".(Gln->pyro-Glu)", "")[1:]
    sub_peptide = peptide

    if peptide.startswith("."):
        sub_peptide = peptide[1:]
    pattern = re.compile(r"\((.*?)\)")
    original_mods = pattern.findall(sub_peptide)
    sub_peptide = re.sub(r"\(.*?\)", ".", sub_peptide)
    position = [i.start() for i in re.finditer(r"\.", sub_peptide)]
    for j in range(1, len(position)):
        position[j] -= j

    for k in range(0, len(original_mods)):
        if position[k] == 0:
            original_mods[k] = str(position[k]) + "-" + mass_offset_dict[original_mods[k]][1]
        else:
            original_mods[k] = str(position[k]) + "-" + original_mods[k] + " (" + sequence[position[k] - 1] + ")"

    original_mods = (
        [str(i) for i in original_mods] if len(original_mods) > 0 else np.nan
    )

    for name, offset in mass_offset_dict.items():
        if "-" in offset[0]:
            mass_offset_proforma = peptide.replace(name, offset[0])
        else:
            mass_offset_proforma = peptide.replace(name, "+" + offset[0])
        peptide = mass_offset_proforma
    mass_offset_proforma = mass_offset_proforma.replace("(", "[").replace(")", "]")
    if peptide.startswith("."):
        mass_offset_proforma = (mass_offset_proforma[:mass_offset_proforma.index(sequence[0])] + "-" +
                                mass_offset_proforma[mass_offset_proforma.index(sequence[0]):])
        mass_offset_proforma = mass_offset_proforma[1:]

    return original_mods, mass_offset_proforma


def convert_psm(idxml, spectra_file, exp_design, export_flr, export_decoy_psm):
    prot_ids = []
    pep_ids = []
    parquet_data = []
    consensus_support = np.nan
    mz_array = []
    intensity_array = []
    num_peaks = np.nan
    id_scores = []
    search_engines = []
    sdrf = pd.read_csv(exp_design, sep="\t")
    PXD = sdrf["source name"].tolist()[0].split("-")[0]
    enable_timstof = False

    oms.IdXMLFile().load(idxml, prot_ids, pep_ids)
    if "ConsensusID" in prot_ids[0].getSearchEngine():
        if prot_ids[0].getSearchParameters().metaValueExists("SE:MS-GF+"):
            search_engines = ["MS-GF+"]
        if prot_ids[0].getSearchParameters().metaValueExists("SE:Comet"):
            search_engines.append("Comet")
        if prot_ids[0].getSearchParameters().metaValueExists("SE:Sage"):
            search_engines.append("Sage")
    else:
        search_engines = [prot_ids[0].getSearchEngine()]

    reference_file_name = os.path.splitext(
        prot_ids[0].getMetaValue("spectra_data")[0].decode("UTF-8")
    )[0]
    spectra_df = pd.read_parquet(spectra_file) if spectra_file else None
    fixed_mods = prot_ids[0].getSearchParameters().fixed_modifications
    var_mods = prot_ids[0].getSearchParameters().variable_modifications
    mdb = oms.ModificationsDB()
    mass_offset_dict = {}
    for m in fixed_mods + var_mods:
        r = mdb.getModification(m)
        mass_offset_dict[m.decode('utf-8').split(" ")[0]] = ['%.4f' % r.getDiffMonoMass(), m.decode('utf-8')]

    for peptide_id in pep_ids:
        retention_time = peptide_id.getRT()
        exp_mass_to_charge = peptide_id.getMZ()
        scan_number = re.findall(
            r"(spectrum|scan)=(\d+)", peptide_id.getMetaValue("spectrum_reference")
        )[0][1]

        if isinstance(spectra_df, pd.DataFrame):
            spectra = spectra_df[spectra_df["scan"] == scan_number]
            mz_array = spectra["mz"].values[0].tolist()
            intensity_array = spectra["intensity"].values[0].tolist()
            num_peaks = len(mz_array)
            collision_energy = spectra["collision energy"].values[0]
            if "ccs" in spectra:
                ccs = spectra["collision energy"].values[0]
                enable_timstof = True

        for hit in peptide_id.getHits():
            # if remove decoy when mapped to target+decoy?
            is_decoy = 0 if hit.getMetaValue("target_decoy") == "target" else 1
            if export_decoy_psm == "false" and is_decoy:
                continue
            global_qvalue = np.nan
            if len(search_engines) > 1:
                if "q-value" in peptide_id.getScoreType():
                    global_qvalue = hit.getScore()
                elif hit.metaValueExists("q-value"):
                    global_qvalue = hit.getScore()
                consensus_support = hit.getMetaValue("consensus_support")
            elif search_engines == "Comet":
                id_scores = ["Comet:Expectation value: " + str(hit.getScore())]
            elif search_engines == "MS-GF+":
                id_scores = ["MS-GF:SpecEValue: " + str(hit.getScore())]
            elif search_engines == "Sage":
                id_scores = ["Sage:hyperscore: " + str(hit.getScore())]

            if hit.metaValueExists("MS:1001491"):
                global_qvalue = hit.getMetaValue("MS:1001491")

            charge = hit.getCharge()
            peptidoform = hit.getSequence().toString()
            cal_mass_to_charge = hit.getSequence().getMZ(charge)
            spectrum_name = os.path.basename(idxml).replace("_consensus_fdr_filter.idXML", "")
            sequence = hit.getSequence().toUnmodifiedString()
            modifications, mass_offset_proforma = mods_position(sequence, peptidoform, mass_offset_dict)
            if ".(Glu->pyro-Glu)" in peptidoform or ".(Gln->pyro-Glu)" in peptidoform:
                if ".(Glu->pyro-Glu)" in peptidoform:
                    peptidoform = peptidoform.replace(".(Glu->pyro-Glu)", "")[0] + "(Glu->pyro-Glu)" + peptidoform.replace(
                        ".(Glu->pyro-Glu)", "")[1:]
                if ".(Gln->pyro-Glu)" in peptidoform:
                    peptidoform = peptidoform.replace(".(Gln->pyro-Glu)", "")[0] + "(Gln->pyro-Glu)" + peptidoform.replace(".(Gln->pyro-Glu)", "")[1:]
                usi_peptide = peptidoform.replace("(", "[").replace(")", "]")
            else:
                usi_peptide = peptidoform.replace("(", "[").replace(")", "]")
            usi = "mzspec:{0}:{1}:scan:{2}:{3}/{4}".format(PXD, spectrum_name, str(scan_number), usi_peptide,
                                                           str(charge))

            protein_accessions = [
                ev.getProteinAccession() for ev in hit.getPeptideEvidences()
            ]
            posterior_error_probability = hit.getMetaValue(
                "Posterior Error Probability_score"
            )
            protein_start_positions = [
                ev.getStart() for ev in hit.getPeptideEvidences()
            ]
            protein_end_positions = [ev.getEnd() for ev in hit.getPeptideEvidences()]
            hit_rank = hit.getRank()

            if export_flr == "true":
                Luciphor_global_flr = hit.getMetaValue("Luciphor_global_flr")
                Luciphor_local_flr = hit.getMetaValue("Luciphor_local_flr")
                parquet_data.append(
                    [
                        sequence,
                        protein_accessions,
                        protein_start_positions,
                        protein_end_positions,
                        modifications,
                        mass_offset_proforma,
                        retention_time,
                        charge,
                        exp_mass_to_charge,
                        cal_mass_to_charge,
                        collision_energy,
                        reference_file_name,
                        scan_number,
                        usi,
                        peptidoform,
                        posterior_error_probability,
                        global_qvalue,
                        is_decoy,
                        consensus_support,
                        mz_array,
                        intensity_array,
                        num_peaks,
                        search_engines,
                        id_scores,
                        hit_rank,
                        Luciphor_global_flr,
                        Luciphor_local_flr
                    ]
                )
            else:
                if enable_timstof:
                    parquet_data.append(
                        [
                            sequence,
                            protein_accessions,
                            protein_start_positions,
                            protein_end_positions,
                            modifications,
                            mass_offset_proforma,
                            retention_time,
                            charge,
                            exp_mass_to_charge,
                            cal_mass_to_charge,
                            collision_energy,
                            ccs,
                            reference_file_name,
                            scan_number,
                            usi,
                            peptidoform,
                            posterior_error_probability,
                            global_qvalue,
                            is_decoy,
                            consensus_support,
                            mz_array,
                            intensity_array,
                            num_peaks,
                            search_engines,
                            id_scores,
                            hit_rank,
                        ]
                    )
                else:
                    parquet_data.append(
                        [
                            sequence,
                            protein_accessions,
                            protein_start_positions,
                            protein_end_positions,
                            modifications,
                            mass_offset_proforma,
                            retention_time,
                            charge,
                            exp_mass_to_charge,
                            cal_mass_to_charge,
                            collision_energy,
                            reference_file_name,
                            scan_number,
                            usi,
                            peptidoform,
                            posterior_error_probability,
                            global_qvalue,
                            is_decoy,
                            consensus_support,
                            mz_array,
                            intensity_array,
                            num_peaks,
                            search_engines,
                            id_scores,
                            hit_rank,
                        ]
                    )
    if export_flr == "true":
        _parquet_field.extend(["Luciphor_global_flr", "Luciphor_local_flr"])
        pd.DataFrame(parquet_data, columns=_parquet_field).to_csv(
            f"{Path(idxml).stem}_psm.csv", index=False, sep="\t"
        )
    else:
        if enable_timstof:
            _parquet_field.extend(["ccs"])
            pd.DataFrame(parquet_data, columns=_parquet_field).to_csv(
                f"{Path(idxml).stem}_psm.csv", index=False, sep="\t"
            )
        else:
            pd.DataFrame(parquet_data, columns=_parquet_field).to_csv(
                f"{Path(idxml).stem}_psm.csv", index=False, sep="\t"
            )


def main():
    idxml_path = sys.argv[1]
    spectra_file = sys.argv[2]
    exp_design = sys.argv[3]
    export_flr = sys.argv[4]
    export_decoy_psm = sys.argv[5]
    convert_psm(idxml_path, spectra_file, exp_design, export_flr, export_decoy_psm)


if __name__ == "__main__":
    sys.exit(main())
