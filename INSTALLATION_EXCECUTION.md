# NanoStamp Installation & Execution Guide

## 1. System Requirements
NanoStamp was developed using:
- Linux HPC environment
- Python v3.x
- Conda environment management
- MATLAB, required for manual segmentation step
- Spateo v1.1.1

## 2. Clone Repository

```bash
git clone https://github.com/abeerahb/NanoStamp

cd NanoStamp
```

## 3. Create Environment

Using conda:

```bash
conda env create -f nanostamp.yml
```

Activate:
```bash
conda activate nanostamp
```

Check installation:
```bash
python --version

pip list
```

## 4. Configure NanoStamp

Before running the pipeline, users must edit the script:
```python
config.py
```
This file controls:
### Project directory
```python
DATASET_PATH = Path(
    "/scratch/prj/crb_nanostamp/abeerah_lipidomics/NanoStamp/data"
)
```

### Input/output locations
```python
RAW_DATA = DATASET_PATH / "raw"
```

## 5. Dataset organisation
User data should include the following files in the below format:
```bash
data/
  raw/
    DATASET_ID (e.g. B1_T0)/
      sample.imzml
      sample.ibd
```

## 6. Pipeline Execution
Run sequentially:
```bash
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
```bash
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
```bash
ModuleNotFoundError
```

Please use the commands below:
```bash
conda activate nanostamp

pip install -r requirements.txt
```

#### Path errors
Check:
```text
config.py
```
