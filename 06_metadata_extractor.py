"""
06_metadata_extractor.py

Purpose
-------
1. Combine metadata from all datasets
2. Add representative lipid annotations
3. Produce one metadata table for downstream analysis

Output
------
sample_metadata.csv
"""

# Load modules
import pandas as pd

from config import OUTPUTS

# ====================
# DATASETS
# ====================

DATASETS = [
    "B1_T0",
    "B1_T1",
    "B2_T0",
    "B2_T1",
    "B3_T0",
    "B3_T1",
]

metadata = []

# =====================
# PROCESS ALL DATASETS
# =====================

for dataset in DATASETS:

    for day in ["Day0", "Day5"]:

        abundance_path = (
            OUTPUTS /
            "lipid_abundance_brain" /
            dataset /
            "lipid_abundance.csv"
        )

        annotation_path = (
            OUTPUTS /
            "lipid_annotation_brain" /
            dataset /
            day /
            "lipid_annotations.csv"
        )

        if not abundance_path.exists():
            print(f"Missing{abundance_path}")
            continue

        abundance = pd.read_csv(abundance_path)

        if annotation_path.exists():
            annotations = pd.read_csv(annotation_path)
        else:
            annotations = None

        patient, treatment = dataset.split("_")

        day_number = int(day.replace("Day", ""))

        for _, row in abundance.iterrows():

            metadata.append({

                "sample_id": dataset,
                "patient": patient,
                "treatment": treatment,
                "day": day_number,

                "target_mz": row["target_mz"],
                "feature_mz": row["feature_mz"],

                "lipid_name": row["lipid_name"],
                "lipid_class": row["lipid_class"],

                "mean_abundance": row["mean_abundance"]

            })

# ====================
# SAVE
# ====================

metadata = pd.DataFrame(metadata)

metadata.to_csv(
    OUTPUTS / "sample_metadata_brain.csv",
    index=False
)

print("Metadata extraction complete")