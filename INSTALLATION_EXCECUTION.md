# NanoStamp Installation & Execution Guide

## 1. System Requirements
NanoStamp was developed using:
- Linux HPC environment
- Python v3.x
- Conda environment management
- MATLAB, required for manual segmentation step
- Spateo v1.1.1

## 2. Clone Repository

```text
git clone https://github.com/abeerahb/NanoStamp

cd NanoStamp
```

## 3. Create Environment

Using conda:

```text
conda env create -f nanostamp.yml
```

Activate:
```text
conda activate nanostamp
```

Check installation:
```text
python --version

pip list
```

## 4. Configure NanoStamp

Before running the pipeline, users must edit the script:
```text
config.py
```
This file controls:
### Project directory
```text
DATASET_PATH = Path(
    "/scratch/prj/crb_nanostamp/abeerah_lipidomics/NanoStamp/data"
)
```

### Input/output locations
```text
RAW_DATA = DATASET_PATH / "raw"
```

## 5. Dataset organisation
User data should include the following files in the below format:
```text
data/
  raw/
    DATASET_ID (e.g. B1_T0)/
      sample.imzml
      sample.ibd
```

## 6. Pipeline Execution
Run sequentially:
```text
01_preprocessing.py
02_feature_matrix.py
03_spatial_lipid_maps
MATLAB segmentation
04_lipid_annotation.py
05_lipid_abundance.py
06_metadata_extractor.py
07_build_anndata.py
08_spateo_analysis.py
```

## 7. Pipeline Output
Explain outputs:
```text
outputs/

processed/
  preprocessed MSI data

feature_matrices/
  spatial lipid matrices

lipid_annotations/
  annotated lipid features/

lipid_abundance/
  regional abundance tables

anndata/
  AnnData objects

spateo/
  aligned spatial models
```
## Troubleshooting

#### Missing Python packages e.g.
```text
ModuleNotFoundError
```

Please use the commands below:
```text
conda activate nanostamp

pip install -r requirements.txt
```

#### Path errors
Check:
```text
config.py
```
