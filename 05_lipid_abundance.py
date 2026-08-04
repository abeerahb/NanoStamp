"""
05_lipid_annotation.py

Purpose
-------
1. Calculate mean abundance of representative lipids
2. Use feature_matrix_brain.csv and feature_matrix_rotated.csv
3. Add lipid identities from lipid_annotations.csv

Output
------
lipid_abundance.csv
"""

# Load modules
import sys
import numpy as np
import pandas as pd

from config import (
    OUTPUTS,
    REPRESENTATIVE_PEAKS
)

# =====================
# ARGUMENT
# =====================

if len(sys.argv) !=2:
    raise ValueError(
        "Usage: python 05_lipid_abundance.py B1_T0"
    )

dataset_id = sys.argv[1]

results = []

# =====================
# PROCESS DAY0 AND DAY5
# =====================

for day in ["Day0", "Day5"]:

    print(f"\nProcessing {dataset_id} {day}")


    # ---------------------
    # Locate feature matrix
    # ---------------------

    feature_dir = (
        OUTPUTS /
        "feature_matrices" /
        dataset_id /
        day
    )


    if day == "Day5":

        feature_file = (
            OUTPUTS /
            "matlab_2" /
            dataset_id /
            "Day5" /
            "feature_matrix_rotated.csv"
        )

    else:

        feature_file = (
            feature_dir /
            "feature_matrix_brain.csv"
        )


    if not feature_file.exists():

        print(
            f"Missing feature matrix: {feature_file}"
        )

        continue


    print(
        f"Using feature matrix: {feature_file}"
    )


    # =====================
    # LOAD FEATURE HEADERS
    # =====================

    feature_matrix = pd.read_csv(
        feature_file,
        nrows=1
    )

    # Clean MATLAB lipid column names
    feature_matrix.columns = [
        col.replace("x", "").replace("_", ".")
        if col.startswith("x") else col
        for col in feature_matrix.columns
    ]

    annotation_path = (
        OUTPUTS / 
        "lipid_annotation_brain" /
        dataset_id /
        day /
        "lipid_annotations.csv"
    )

    if (
        not feature_file.exists()
        or
        not annotation_path.exists()
    ):
        print(f"Missing files for {dataset_id} {day}")
        continue

    print(f"Processing {dataset_id} {day}")

    feature_matrix = pd.read_csv(feature_file)

    # Cleans MATLAB column names
    feature_matrix.columns = [
        str(col).replace("x", "").replace("_", ".")
        if str(col).startswith("x")
        else str(col)
        for col in feature_matrix.columns
    ]

    annotations = pd.read_csv(annotation_path)

    if annotations.empty:
        print(f"No annotations found for {dataset_id} {day}")
        continue

    metadata_cols = [
        "sample_id",
        "day",
        "x",
        "y"
    ]

    feature_cols = []

    feature_mzs = []

    for col in feature_matrix.columns:

        if col in metadata_cols:
            continue

        try:

            mz = float(col)

            feature_cols.append(col)

            feature_mzs.append(mz)

        except ValueError:

            continue

    feature_mzs = np.array(feature_mzs)

    print(f"Found {len(feature_cols)} lipid features")

    # =======================
    # REPRESENTATIVE PEAKS
    # =======================

    for target_mz in REPRESENTATIVE_PEAKS:

        closest_idx = np.argmin(
            np.abs(feature_mzs - target_mz)
        )

        feature_mz = feature_mzs[closest_idx]
        feature_col = feature_cols[closest_idx]

        annotation = annotations.iloc[
            (annotations["mz"] - feature_mz).abs().idxmin()
        ]

        mean_abundance = feature_matrix[
            feature_col
        ].mean()

        results.append({

            "sample_id": dataset_id,
            "day": day,
            
            "target_mz": target_mz,
            "feature_mz": feature_mz,

            "lipid_name":
                annotation["lipid_name"],

            "lipid_class":
                annotation["lipid_class"],
            
            "mean_abundance":
                mean_abundance
        })

# =======================
# SAVE
# =======================

results = pd.DataFrame(results)

output_dir = (
    OUTPUTS / 
    "lipid_abundance_brain" /
    dataset_id
)

output_dir.mkdir(
    parents=True,
    exist_ok=True
)

results.to_csv(
    output_dir /
    "lipid_abundance.csv",
    index=False
)

print(f"Saved {dataset_id} lipid abundance table")