"""
04_lipid_annotation.py

Purpose
-------
1. Annotate lipid features using feature matrices
2. Match m/z bins to lipid database
3. Save annotations for each dataset and day

Outputs
-------
lipid_annotations.csv
"""


# =====================
# IMPORT MODULES
# =====================

import sys
import numpy as np
import pandas as pd
from pathlib import Path

from config import OUTPUTS


# =====================
# ARGUMENT
# =====================

if len(sys.argv) != 2:
    raise ValueError(
        "Usage: python 04_lipid_annotation.py B1_T0"
    )

dataset_id = sys.argv[1]


# =====================
# LOAD LIPID DATABASE
# =====================

DB_PATH = Path(
    "/scratch/prj/crb_nanostamp/abeerah_lipidomics/NanoStamp/resources/lipid_database.csv"
)

if not DB_PATH.exists():
    raise FileNotFoundError(
        f"Lipid database not found:\n{DB_PATH}"
    )

lipid_db = pd.read_csv(DB_PATH)


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

    # =====================
    # EXTRACT M/Z FEATURES
    # =====================

    mz_columns = []

    for col in feature_matrix.columns:

        try:

            mz = float(col)

            mz_columns.append(mz)

        except ValueError:

            # ignore x,y columns etc.
            continue


    print(
        f"Found {len(mz_columns)} m/z features"
    )


    # =====================
    # ANNOTATE LIPIDS
    # =====================

    annotations = []


    for mz in mz_columns:


        matches = lipid_db[
            np.abs(
                lipid_db["mz"] - mz
            ) <= 0.01
        ]


        if matches.empty:

            continue


        best = matches.iloc[0]

        annotations.append({

            "mz": mz,

            "lipid_name":
                best["lipid_name"],

            "lipid_class":
                best["lipid_class"]

        })

    annotations = pd.DataFrame(
        annotations
    )

    # =====================
    # SAVE OUTPUT
    # =====================

    output_dir = (
        OUTPUTS /
        "lipid_annotation_brain" /
        dataset_id /
        day
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    output_file = (
        output_dir /
        "lipid_annotations.csv"
    )

    annotations.to_csv(
        output_file,
        index=False
    )

    print(
        f"Saved {len(annotations)} annotations:"
    )

    print(output_file)

print(
    f"\nFinished {dataset_id}"
)