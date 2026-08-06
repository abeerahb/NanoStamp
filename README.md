# NanoStamp

Gliomas account for approximately 78-80% of malignant central nervous system tumours and remain associated with poor clinical outcomes despite advances in molecular profiling and targeted therapies. A major contributor to treatment resistance is metabolic reprogramming, whereby tumour cells alter lipid metabolism to support proliferation, membrane integrity, and survival. Although spatial lipidomics has improved understanding of lipid heterogeneity within gliomas, conventional approaches remain largely destructive and are unable to capture longitudinal molecular changes within living tissues. 

This study developed NanoStamp, an end-to-end computational pipeline for the processing and explanatory analysis of longitudinal nanoneedle-derived mass spectrometry imaging (MSI) datasets. The workflow was developed and evaluated using six longitudinal ex vivo murine glioblastoma DESI-MSI datasets (N = 6) and integrates spectral preprocessing, mass calibration, feature matrix generation, lipid annotation, spatial registration, AnnData construction, and spatiotemporal modelling using Spateo to reconstruct lipidomic changes across treatment timepoints. The pipeline was subsequently applied to investigate lipidomic responses to temozolomide (TMZ) treatment in longitudinal glioblastoma tissue samples. 

NanoStamp successfully processed and integrated longitudinal MSI datasets within a unified analytical framework while preserving spatial information and anatomical context. Exploratory longitudinal analyses identified region-specific changes in representative lipid-associated features following treatment, with distinct response patterns observed across grey matter, tumour, and white matter. Grey matter exhibited the greatest overall temporal lipidomic change, whereas tumour tissue demonstrated the most pronounced treatment-associated remodelling, supporting the concept that lipid responses to TMZ are spatially heterogenous. These findings suggest that lipid remodelling represents an important component of glioma adaptation to therapy. 

Overall, this study demonstrates the feasibility of integrating longitudinal nanoneedle-derived MSI datasets within a reproducible computational workflow. By combining spatial lipidomics with spatiotemporal modelling, NanoStamp provides a reproducible computational framework for investigating treatment-induced metabolic adaptation in living glioma tissues and establishes a foundation for future longitudinal spatial omics analyses and multimodal integration.


## Pipeline Steps:
1. Preprocessing
2. Feature matrix
3. Spatial lipid maps

-> MATLAB segmentation

4. Lipid annotation
5. Lipid abundance
6. Metadata extractor
7. Build AnnData
8. Spateo Analysis

NOTE: This repository also contains the environment yaml used, the link to the raw dataset for users to recreate results, as well a file containing extra scripts for reading outputs from the main pipeline, alongside scripts that can be run to keep as a reference guide for running Spateo.

#### This pipeline was built on a Linux High Performance Computer (HPC). Resultantly, all scripts are configured appropriately to run on Linux. Users will need to independently customise scripts to the interface they use accordingly. All customisations would be appreciated and if you wish to share this, please feel free to fork a branch!
