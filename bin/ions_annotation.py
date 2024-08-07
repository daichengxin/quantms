#!/usr/bin/env python

import argparse
import errno
import os
import sys
from pathlib import Path
import pandas as pd
import spectrum_utils.spectrum as sus


def ions_annotation(csv_file):
    df = pd.read_csv(csv_file, sep="\t")

    def ion_annotation(row):
        t1 = eval(row["mz_array"])
        mz_array = [i for i in t1 if i != ""]
        t2 = eval(row["intensity_array"])
        intensity_array = [i for i in t2 if i != ""]
        usi_spectrum = sus.MsmsSpectrum(identifier=row["USI"], precursor_mz=row["exp_mass_to_charge"],
                                        precursor_charge=row["charge"],
                                        mz=mz_array,
                                        intensity=intensity_array,
                                        retention_time=row["retention_time"] / 60)
        peptide = row["peptidoform"].replace("(", "[").replace(")", "]").replace("UniMod", "UNIMOD").replace(".[Acetyl]", "[Acetyl]-")

        usi_spectrum.annotate_proforma(
            peptide,
            fragment_tol_mass=20,
            fragment_tol_mode="ppm",
            ion_types="by",
            neutral_losses={"NH3": -17.026549, "H2O": -18.010565},
        )
        peak_annotate = []
        for idx in range(len(usi_spectrum.mz)):
            if str(usi_spectrum.annotation[idx]) != "?":
                peak_annotate.append(f"({usi_spectrum.mz[idx]:.4f},{usi_spectrum.intensity[idx]:.2f},{str(usi_spectrum.annotation[idx])})")
        annotations = ";".join(peak_annotate)
        annotations = peak_annotate
        return annotations

    df.drop(columns=["protein_accessions", "protein_start_positions", "protein_end_positions",
                     "is_decoy", "id_scores", "hit_rank"], inplace=True)
    df["ions_matched"] = df.apply(lambda row: ion_annotation(row), axis=1)
    df.to_csv(f"{Path(csv_file).stem}_ion_psm.csv", index=False)


def main():
    csc_file = sys.argv[1]
    ions_annotation(csc_file)


if __name__ == "__main__":
    sys.exit(main())
