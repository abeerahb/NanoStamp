"""
01_preprocessing.py

Purpose
-------
- Processes ONE biological dataset at a time. 
1. Read raw imzmML files
2. Apply lock-mass calibration
3. Save aligned imzML files
4. Crop tissue ROIs
5. Generate:
    - common representation
    - meaningful signal
    - segmentation
"""

# =====================
# LOAD MODULES 
# =====================

import sys
import pandas as pd

from processing import (
    aligned_representation, 
    process
)

from config import(
    RAW_DATA,
    ALIGNED_DATA,
    PROCESSED_DATA,
    METADATA_PATH,
    LOCK_MASS_PEAK,
    LOCK_MASS_TOL,
    MZ_START,
    MZ_END,
    MASS_RESOLUTION,
    REPRESENTATIVE_PEAKS
)

# =====================
# ARGUMENT
# =====================

if len(sys.argv) != 2:
    raise ValueError(
        "Usage: python 01_preprocessing.py B1_T0"
    )

dataset_id = sys.argv[1]

print(f"\nProcessing dataset: {dataset_id}")

# =====================
# LOAD METADATA
# =====================

print("Loading metadata...")

metadata_df = pd.read_csv(METADATA_PATH)

dataset_rows = metadata_df[
    metadata_df["sample_file_name"].str.contains(
        dataset_id,
        na=False
    )
]

if len(dataset_rows) == 0:
    raise ValueError(
        f"No metadata rows found for {dataset_id}"
    )

# =====================
# PROCESS EACH ROI
# =====================

for _, roi in dataset_rows.iterrows():

    sample_name = roi.sample_file_name

    if "Day0" in sample_name:
        day = "Day0"
    elif "Day5" in sample_name:
        day = "Day5"
    else:
        day = "Unknown"
    
    print(f"\n--- {sample_name} ---")

    # =====================
    # MASS LOCK CALIBRATION
    # =====================

    raw_path = RAW_DATA / f"{roi.file_name}.imzml"

    aligned_path = (
        ALIGNED_DATA /
        f"{roi.file_name}.imzML"
    )

    aligned_path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    print("Lock mass calibration...")

    aligned_representation(
        raw_path,
        aligned_path,
        LOCK_MASS_PEAK,
        LOCK_MASS_TOL
    )

    # =====================
    # ROI PROCESSING
    # =====================

    output_path = (
        PROCESSED_DATA /
        dataset_id /
        day
    )

    output_path.mkdir(
        parents=True,
        exist_ok=True
    )

    print("ROI processing...")

    process(
        aligned_path,
        output_path,
        roi.x_min,
        roi.x_max,
        roi.y_min,
        roi.y_max,
        MZ_START,
        MZ_END,
        MASS_RESOLUTION,
        REPRESENTATIVE_PEAKS
    )

print("\nFinished")
