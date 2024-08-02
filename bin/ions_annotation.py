#!/usr/bin/env python

import argparse
import errno
import os
import sys
from pathlib import Path
import pandas as pd
import spectrum_utils as sus


def ions_annotation(csv_file):
    df = pd.read_csv(csv_file)

    def ion_annotation(row):
        usi_spectrum = sus.MsmsSpectrum(identifier=row["usi"], precursor_mz=row["exp_mass_to_charge"],
                                        precursor_charge=row["charge"],
                                        mz=row["mz_array"], intensity=row["intensity_array"],
                                        retention_time=row["retention_time"] / 60)
        peptide = row["peptidoform"].replace("(", "[").replace(")", "]").replace("UniMod", "UNIMOD")
        usi_spectrum.annotate_proforma(
            peptide,
            fragment_tol_mass=20,
            fragment_tol_mode="ppm",
            ion_types="by",
            neutral_losses={"NH3": -17.026549, "H2O": -18.010565},
        )
        annotation = usi_spectrum.annotation
        ions_matched = [str(n) for n in annotation if str(n) != "?"]
        return ions_matched

    df["ions_matched"] = df.apply(lambda row, col: ion_annotation(row), axis=1)
    df.to_csv(f"{Path(csv_file).stem}_ion_psm.csv", index=False)


def main():
    csc_file = sys.argv[1]
    ions_annotation(csc_file)


if __name__ == "__main__":
    sys.exit(main())
