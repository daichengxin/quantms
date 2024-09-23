#!/usr/bin/env python

import argparse
import errno
import os
import sys
from pathlib import Path
import pandas as pd
import spectrum_utils.spectrum as sus
from pandarallel import pandarallel

import os
os.environ['MEMORY_FS_ROOT'] = "/tmp"
os.environ['JOBLIB_TEMP_FOLDER']  = "/tmp"

def ions_annotation(csv_file, fragment_tol_mass, fragment_tol_mode):
    df = pd.read_csv(csv_file, sep="\t")

    def ion_annotation(row):
        mz_array = eval(row["mz_array"])
        intensity_array = eval(row["intensity_array"])
        usi_spectrum = sus.MsmsSpectrum(identifier=row["USI"], precursor_mz=row["exp_mass_to_charge"],
                                        precursor_charge=row["charge"],
                                        mz=mz_array, intensity=intensity_array,
                                        retention_time=row["retention_time"] / 60)
        usi_spectrum.annotate_proforma(
            row["mass_offset_proforma"],
            fragment_tol_mass=fragment_tol_mass,
            fragment_tol_mode=fragment_tol_mode,
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
                     "id_scores", "hit_rank", "reference_file_name", "scan_number",
                     "num_peaks", "search_engines", "collision_energy"], inplace=True)

    df["ions_matched"] = df.parallel_apply(lambda row: ion_annotation(row), axis=1)
    df.drop(columns=["mass_offset_proforma"], inplace=True)
    df.to_parquet(f"{Path(csv_file).stem}_ion_psm.parquet", index=False, compression="gzip")


def main():
    csc_file = sys.argv[1]
    fragment_tol_mass = float(sys.argv[2])
    fragment_tol_mode = sys.argv[3]
    pandarallel.initialize(progress_bar=True, nb_workers=10, use_memory_fs = False)
    ions_annotation(csc_file, fragment_tol_mass, fragment_tol_mode)


if __name__ == "__main__":
    sys.exit(main())
