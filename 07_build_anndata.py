"""
07_build_anndata.py

Purpose
-------
1. Build AnnData objects for Spateo
2. Use all lipid features
3. Store spatial coordinates
4. Store lipid annotations
5. Store per-pixel metadata

Output
------
B1_T0_Day0.h5ad
B1_T0_Day5.h5ad
"""

# Load modules
import sys
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt

from config import OUTPUTS

# =====================
# ARGUMENT
# =====================

if len(sys.argv) != 2:
    raise ValueError(
        "Usage: python 07_build_anndata.py B1_T0"
    )

dataset_id = sys.argv[1]

# =====================
# DATASET ORIENTATION
# =====================

orientation_map = {

    "B1_T0": -78.75,
    "B1_T1": -78.75,
    "B2_T0": -78.75,
    "B2_T1": -78.75,
    "B3_T0": -78.75,
    "B3_T1": -78.75

}


def correct_orientation(coords, angle_deg):

    coords = coords.copy()

    # centre coordinates
    centre = coords.mean(axis=0)

    centred = coords - centre

    theta = np.deg2rad(angle_deg)

    rotation_matrix = np.array(
        [
            [np.cos(theta), -np.sin(theta)],
            [np.sin(theta), np.cos(theta)]
        ]
    )

    rotated = centred @ rotation_matrix.T

    # shift back to positive coordinates
    rotated[:,0] -= rotated[:,0].min()
    rotated[:,1] -= rotated[:,1].min()

    return rotated

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

    feature_matrix = pd.read_csv(feature_file)

    # Clean MATLAB lipid column names
    new_columns = []

    for col in feature_matrix.columns:

        # Only MATLAB lipid columns look like x600_5100
        if col.startswith("x") and len(col) > 1 and col[1].isdigit():
            new_columns.append(col[1:].replace("_", "."))

        else:
            new_columns.append(col)

    feature_matrix.columns = new_columns

    annotation_path = (
        OUTPUTS /
        "lipid_annotation_brain" /
        dataset_id /
        day /
        "lipid_annotations.csv"
    )

    region_path = (
        OUTPUTS /
        "matlab_2" /
        dataset_id /
        day /
        "pixel_regions.csv"
    )

    if not feature_file.exists():
        print(f"Missing feature matrix: {feature_file}")
        continue

    if not annotation_path.exists():
        print(f"Missing annotations: {annotation_path}")
        continue

    if not region_path.exists():
        print(f"Missing region labels: {region_path}")
        continue

    # Load files
    annotations = pd.read_csv(annotation_path)
    regions = pd.read_csv(region_path)

    # =========================
    # FIX COORDINATE PRECISION
    # =========================
    # MATLAB rotation creates floating point coordinates.
    # Round before merging so feature matrix and segmentation match.

    feature_matrix["x"] = feature_matrix["x"].astype(float).round(6)
    feature_matrix["y"] = feature_matrix["y"].astype(float).round(6)

    regions["x"] = regions["x"].astype(float).round(6)
    regions["y"] = regions["y"].astype(float).round(6)


    print("\nCoordinate check:")
    print("Feature matrix:")
    print(feature_matrix[["x", "y"]].head())

    print("\nRegions:")
    print(regions[["x", "y"]].head())


    # =========================
    # Convert day column
    # =========================

    if feature_matrix["day"].dtype == object:

        feature_matrix["day"] = (
            feature_matrix["day"]
            .astype(str)
            .str.replace("Day", "", regex=False)
            .astype(int)
        )

    # =========================
    # Merge manual segmentation
    # =========================

    feature_matrix = feature_matrix.merge(
        regions,
        on=["x", "y"],
        how="left"
    )


    merged_check = feature_matrix.merge(
        regions,
        on=["x", "y"],
        how="inner"
    )

    print(
        "Matched pixels:",
        len(merged_check),
        "/",
        len(feature_matrix)
    )

    feature_matrix["region"] = (
        feature_matrix["region"]
        .fillna(0)
        .astype(int)
    )

    region_names = {
        0: "Background",
        1: "Tumour",
        2: "Grey",
        3: "White"
    }

    feature_matrix["region_name"] = (
        feature_matrix["region"]
        .map(region_names)
    )

    # Remove non-biological background pixels for Spateo
    feature_matrix = feature_matrix[
        feature_matrix["region_name"] != "Background"
    ].copy()

    # =========================
    # Identify lipid columns
    # =========================

    feature_cols = []

    for col in feature_matrix.columns:

        try:
            float(col)
            feature_cols.append(col)
        except ValueError:
            continue

    print(f"Found {len(feature_cols)} lipid features")
    
    # =========================
    # X matrix
    # =========================

    X = feature_matrix[
        feature_cols
    ].to_numpy(dtype=np.float32)

    # =========================
    # obs
    # =========================

    obs = feature_matrix[
        [
            "sample_id",
            "day",
            "x",
            "y",
            "region",
            "region_name"
        ]
    ].copy()

    print("\nRegion counts")

    print(obs["region_name"].value_counts())

    obs["patient"] = dataset_id.split("_")[0]
    obs["treatment"] = dataset_id.split("_")[1]
    
    obs["region_name"] = pd.Categorical(
    obs["region_name"],
    categories=[
        "Background",
        "Tumour",
        "Grey",
        "White"
        ]
    )   

    # Spateo fix
    
    obs.index = [f"pixel_{day}_{i}" for i in range(len(obs))]

    # =========================
    # Annotation prep
    # =========================

    annotations["mz"] = annotations["mz"].astype(float)

    annotation_mzs = annotations["mz"].to_numpy()

    # =========================
    # var (feature metadata)
    # =========================    

    var_rows = []

    for feature in feature_cols:

        mz = float(feature)

        lipid_name = "Unknown"
        lipid_class = "Unknown"

        if len(annotation_mzs) > 0:

            idx = np.argmin(
                np.abs(annotation_mzs - mz)
            )

            if np.abs(annotation_mzs[idx] - mz) <= 0.01:

                lipid_name = annotations.iloc[idx]["lipid_name"]
                lipid_class = annotations.iloc[idx]["lipid_class"]

        var_rows.append(
            {
                "mz": mz,
                "lipid_name": lipid_name,
                "lipid_class": lipid_class
            }
        )

    var = pd.DataFrame(var_rows)
    var.index = feature_cols

    # =========================
    # Check dimensions
    # =========================

    print(f"X shape   : {X.shape}")
    print(f"obs shape : {obs.shape}")
    print(f"var shape : {var.shape}")

    if X.shape[0] != len(obs):
        raise ValueError(
            "obs does not match number of pixels."
        )

    if X.shape[1] != len(var):
        raise ValueError(
            "var does not match number of lipid features"
        )

    # =========================
    # Build AnnData
    # =========================

    adata = ad.AnnData(
        X=X,
        obs=obs,
        var=var
    )

    # =========================
    # Spatial coordinates
    # =========================

    adata.obsm["spatial"] = (
        obs[
            ["x", "y"]
        ]
        .to_numpy(dtype=np.float32)
    )

    # =====================================
    # Apply spatial orientation correction
    # =====================================

    if day == "Day5":

        orientation = orientation_map.get(
            dataset_id,
            0
        )

    else:

        orientation = 0

    print(
        f"Applying rotation angle: {orientation} degrees"
    )

    adata.obsm["spatial"] = correct_orientation(
        adata.obsm["spatial"],
        orientation
    )

    # =====================
    # SAVE DIRECTORY
    # =====================

    output_dir = (
        OUTPUTS /
        "anndata_raw" /
        dataset_id
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    # =========================
    # SPATEO ROTATION VERIFY
    # =========================

    plt.scatter(
        adata.obsm["spatial"][:,0],
        adata.obsm["spatial"][:,1],
        s=2,
        c=adata.obs["region"]
    )

    plt.gca().invert_yaxis()
    plt.title(f"{dataset_id} {day} AnnData spatial check")

    plt.savefig(
        output_dir / f"{dataset_id}_{day}_spatial_check.png",
        dpi=300
    )

    plt.close()

    # =========================
    # Save
    # =========================

    adata.write_h5ad(
        output_dir /
        f"{dataset_id}_{day}.h5ad"
    )

    print(
        f"Saved {dataset_id}_{day}.h5ad"
    )

print("\nFinished building AnnData objects.")