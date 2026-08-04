"""
02_feature_matrix.py 

Purpose
-------
1. Remove non brain pixels
2. Filter lipid m/z range
3. TIC normalisation
4. Peak alignment
5. Final TIC normalisation
6. Feature matrix generation

Outputs
-------
feature_matrix_brain.csv
"""

import sys
import numpy as np
import pandas as pd

from pyimzml.ImzMLParser import ImzMLParser

from config import (
    PROCESSED_DATA,
    OUTPUTS,
    LIPID_MZ_MIN,
    LIPID_MZ_MAX
)

# =====================
# ARGUMENT
# =====================

if len(sys.argv) != 2:
    raise ValueError(
        "Usage: python 02_feature_matrix.py B1_T0"
    )

dataset_id = sys.argv[1]
 
# =======================
# GLOBAL M/Z BINS
# =======================

BIN_SIZE = 0.02

BINS = np.arange(
    LIPID_MZ_MIN,
    LIPID_MZ_MAX + BIN_SIZE,
    BIN_SIZE,
    dtype=np.float32
)

FEATURE_COLUMNS = [
    f"{b:.4f}"
    for b in BINS
]

# =======================
# TIC NORMALISATION
# =======================

def tic_normalize(intensities):
    """
    Pixel-wise TIC normalisation 

    Equivalent to:
        TIC = sum(hyperimage,3)
        hyperimageN = hyperimage ./ TIC
    in MATLAB.
    """

    intensities = np.asarray(intensities, dtype=np.float32)

    tic = intensities.sum(dtype=np.float64)

    if tic == 0:
        return np.zeros_like(intensities)

    return intensities / tic

# =======================
# M/Z FILTERING
# =======================

def filter_mz_range(mzs, intensities):

    mask = (
        (mzs >= LIPID_MZ_MIN) 
        &
        (mzs <= LIPID_MZ_MAX)
    )

    return (
        mzs[mask],
        intensities[mask]
    )


# =======================
# PEAK ALIGNMENT
# =======================

def align_peaks(mzs, intensities):

    aligned = np.zeros(
        len(BINS),
        dtype=np.float32
    )

    idx = np.digitize(
        mzs,
        BINS
    ) - 1

    valid = (
        (idx >= 0) 
        &
        (idx < len(BINS))
    )

    np.add.at(
        aligned,
        idx[valid],
        intensities[valid]
    )

    return aligned

# ================================
# BUILD BRAIN-ONLY FEATURE MATRIX
# ================================

dataset_folder = PROCESSED_DATA / dataset_id

for day in ["Day0", "Day5"]:

    rows = []

    folder = dataset_folder / day

    if not folder.exists():
        continue

    if not folder.is_dir():
        continue

    imzml_path = folder / "common_representation.imzML"
    mask_path = folder / "segmentation.npy"

    if (
        not imzml_path.exists()
        or
        not mask_path.exists()
    ):
        continue

    print(f"Processing {dataset_id} {day}")

    mask = np.load(mask_path)

    parser = ImzMLParser(str(imzml_path))

    for idx, (x, y, _) in enumerate(parser.coordinates):

        # ----------------------
        # Keep only brain pixels
        # ----------------------

        try:
            if mask[y - 1, x - 1] <= 0:
                continue
        except IndexError:
            continue

        # ----------------------
        # Read spectrum
        # ----------------------

        mzs, spectrum = parser.getspectrum(idx)

        mzs = np.asarray(mzs, dtype=np.float32)
        spectrum = np.asarray(spectrum, dtype=np.float32)

        # --------------------------
        # Keep only lipid m/z range
        # --------------------------

        mzs, spectrum = filter_mz_range(
            mzs,
            spectrum
        )

        # -----------------------------------
        # TIC normalise retained lipid peaks
        # -----------------------------------

        spectrum = tic_normalize(
            spectrum
        )

        # ----------------------
        # Peak alignment
        # ----------------------

        aligned = align_peaks(
            mzs,
            spectrum
        )

        # -------------------------
        # Final TIC normalisation
        # -------------------------

        aligned = tic_normalize(
            aligned
        )

        # ----------------------
        # Store feature vector
        # ----------------------

        row = {
            "sample_id": dataset_id,
            "day": day,
            "x": x,
            "y": y
        }

        row.update(
            dict(
                zip(
                    FEATURE_COLUMNS,
                    aligned
                )
            )
        )

        rows.append(row)

    # ===================================
    # BUILD FEATURE MATRIX
    # ===================================

    feature_matrix = pd.DataFrame(rows)

    OUTPUTS.mkdir(
        parents=True,
        exist_ok=True
    )

    output_dir = (
        OUTPUTS /
        "feature_matrices" /
        dataset_id /
        day
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    # ===================================
    # SAVE
    # ===================================

    feature_matrix.to_csv(
        output_dir /
        "feature_matrix_brain.csv",
        index=False
    )

    print(f"Saved {len(feature_matrix)} brain pixels")

print(f"\nFinished dataset {dataset_id}")