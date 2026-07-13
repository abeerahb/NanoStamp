# NanoStamp

Gliomas are a type of malignant brain tumour originating from glial cells in the brain and spinal cord representing 78-80% of all malignant CNS brain tumours. Despite advances in molecular and genetic research, gliomas pose ongoing challenges in treatment and management, due to the gliomas’ characteristic to reprogram a lipid’s metabolism to enhance cell membrane integrity and support abnormal cell proliferation. Spatial lipidomics is a particular type of spatial omics involving the analysis of a complete set of lipids within a biological system, looking at the structure, function and interaction of lipids. It has recently been used in research to analyse gliomas through the use of porous silicon nanoneedles allowing for, repeated sampling of live glioma specimens with minimal invasiveness. Our study builds on this research by developing a computational pipeline, known as Nanostamp to analyse nanoneedle-based living tissue datasets, with a focus on validating and integrating lipidomic and metabolic profiles. By reconstructing dynamic relationships across lipidomic activity through space and time, this project seeks to characterise the spatiotemporal evolution of molecular states in gliomas, ultimately improving our understanding of tumour progression and metabolic reprogramming in cancer (Mayo Clinic, 2019).

## Pipeline Steps:
1. Preprocessing
2. Feature matrix
3. Spatial lipid maps
4. Lipid annotation
5. Lipid abundance
6. Metadata extractor
7. Build AnnData
8. Spateo Analysis

NOTE: This repository also contains the environment yaml used, the link to the raw dataset for users to recreate results, as well a file containing extra scripts for reading outputs from the main pipeline, alongside scripts that can be run to keep as reference guide for running Spateo.

#### This pipeline was built on a Linux High Performance Computer (HPC). Resultantly, all scripts are configured appropriately to run on Linux. Users will need to independently customise scripts to the interface they use acoordingly. All customisations would be appreciated and if you wish to share this, please feel free to fork a branch!
