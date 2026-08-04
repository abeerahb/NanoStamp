from pathlib import Path

# =======================
# PATHS
# =======================

DATASET_PATH = Path(
    "/scratch/prj/crb_nanostamp/abeerah_lipidomics/NanoStamp/data"
)

RAW_DATA = DATASET_PATH / "raw"

ALIGNED_DATA = DATASET_PATH / "aligned"

PROCESSED_DATA = DATASET_PATH / "processed"

OUTPUTS = DATASET_PATH / "outputs"

METADATA_PATH = DATASET_PATH / "metadata.csv"

# =======================
# PROCESSING SETTINGS
# =======================

MZ_START = 50
MZ_END = 1200

MASS_RESOLUTION = 0.025

LOCK_MASS_PEAK = 885.5498
LOCK_MASS_TOL = 0.3

LIPID_MZ_MIN = 600
LOCK_MASS_TOL = 0.3

LIPID_MZ_MIN = 600
LIPID_MZ_MAX = 900

REPRESENTATIVE_PEAKS = [
    600.51, # GM
    682.59, # Tumour
    888.62  # WM
]