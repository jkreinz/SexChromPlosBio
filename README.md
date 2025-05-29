# Amaranthus tuberculatus Sex Chromosome Scripts - PLOS 2025

## Overview
Comprehensive genomic analysis of sex chromosome evolution in *Amaranthus tuberculatus*, focusing on the sex-linked region on Scaffold_1.

## Data Availability

### Analysis Scripts (GitHub)
All scripts used to process and plot all figures (main and sup) are available in this repository under `scripts/r`.
Data corresponding to these scripts are linked and detailed below in the Data section. 

Several scripts related to mapping read, calling SNPs, and assembling genomes are also available in `scripts/bash`.

### Data (Zenodo)
Complete dataset organized by analysis type:
**[![DOI](https://zenodo.org/badge/DOI/YOUR_DOI_HERE.svg)](https://doi.org/YOUR_DOI_HERE)**

The Zenodo archive contains a `data_by_script/` folder with:

data_by_script/
├── 01_GWAS/           # GWAS and population genetics data
├── 02_Comparative/    # GENESPACE comparative genomics data
├── 03_PCA/           # PCA results for inversions
├── 04_Depth/         # Sequencing depth matrices
├── 05_Phylogenetic/  # Phylogenetic trees and metadata
└── 06_Recombination/ # Recombination rate data

PLOS Bio also requires deposition of the numeric values underlying each plot, and those are available as well in the Zenodo archive in the folder `figure_data_export/`. The only exception to this is juicebox plots of the Hi-C maps and GENESPACE plots of synteny. We are not able to provide scripts for the former since these were manually curated, however, one can reproduce all GENESPACE plots with the scripts provided.  

## Quick Start

### 1. Download Data
```bash
# Download from Zenodo and extract
wget https://zenodo.org/record/YOUR_RECORD/files/data_by_script.zip
unzip data_by_script.zip
```

### 2. Set up environment
```
# Install required packages
source("config/requirements.R")

# Update file paths to point to your downloaded data
source("config/file_paths.R")
```

### 3. Run Analyses from Zenodo Deposited Datasets
```
# Run in order (each script uses its corresponding data folder)
source("scripts/01_GWAS_PopGen_Analyses.R")      # Uses data_by_script/01_GWAS/
source("scripts/02_Comparative_Genomic_Analyses.R") # Uses data_by_script/02_Comparative/
source("scripts/03_PCA_Structure_Analyses.R")    # Uses data_by_script/03_PCA/
source("scripts/04_Depth_Based_Analyses.R")      # Uses data_by_script/04_Depth/
source("scripts/05_Phylogenetic_Analysis.R")     # Uses data_by_script/05_Phylogenetic/
source("scripts/06_Recombination_Analysis.R")    # Uses data_by_script/06_Recombination/

```
