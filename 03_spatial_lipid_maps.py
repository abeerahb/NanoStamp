"""
03_spatial_lipid_maps.py

Purpose of script
-----------------
1. Generate spatial lipid maps 
2. Use brain-only feature matrices
3. Compare Day0 vs Day5
4. Plot representative lipids only
"""

# Load modules
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from config import (
    OUTPUTS,
    REPRESENTATIVE_PEAKS
)

# =====================
# ARGUMENT
# =====================

if len(sys.argv) !=2:
    raise ValueError(
        "Usage: python 03_spatial_lipid_maps.py B1_T0"
    )

dataset_id = sys.argv[1]

# =====================
# PROCESS DAYS
# =====================

for day in ["Day0", "Day5"]:

    feature_dir = (
        OUTPUTS /
        "feature_matrices" /
        dataset_id /
        day
    )

    feature_path = (
        feature_dir /
        "feature_matrix_brain.csv"
    )

    if not feature_path.exists():
        print(f"Missing {feature_path}")
        continue

    print(f"Generating maps for {dataset_id} {day}")

    df = pd.read_csv(feature_path)

    feature_cols = [
        c
        for c in df.columns
        if c not in [
            "sample_id",
            "day",
            "x",
            "y"
        ]
    ]

    feature_mzs = np.array(
        [float(mz) for mz in feature_cols],
        dtype=np.float32
    )

    output_dir = (
        OUTPUTS /
        "spatial_maps_brain" /
        dataset_id /
        day
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    # ==================================
    # REPRESENTATIVE LIPIDS ONLY
    # ==================================

    for target_mz in REPRESENTATIVE_PEAKS:

        closest_index = np.argmin(
            np.abs(feature_mzs - target_mz)
        )

        actual_mz = feature_mzs[closest_index]
        lipid = feature_cols[closest_index]

        print(
            f"{target_mz:.2f} -> using {actual_mz:.4f}"
        )

        width = int(df["x"].max())
        height = int(df["y"].max())

        image = np.full(
            (height, width),
            np.nan
        )

        for _, row in df.iterrows():

            image[
                int(row["y"]) - 1,
                int(row["x"]) - 1
            ] = row[lipid]

        plt.figure(figsize=(6,6))

        plt.imshow(
            image,
            cmap="magma",
            origin="lower"
        )

        plt.colorbar(
            label="cbar.set_label('Absolute Intensity')"
        )

        plt.title(
            f"{dataset_id} {day}\n{target_mz:.2f} m/z"
        )

        plt.axis("off")

        plt.tight_layout()

        plt.savefig(
            output_dir /
            f"{target_mz:.2f}_closest_{actual_mz:.4f}.png",
            dpi=300
        )

        plt.close()

print("Spatial lipid maps complete")