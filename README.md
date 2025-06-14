# Amaranthus tuberculatus Sex Chromosome Scripts - PLOS 2025

## Overview
Comprehensive genomic analysis of sex chromosome evolution in *Amaranthus tuberculatus*, focusing on the sex-linked region on Scaffold_1.

## Data Availability

### Analysis Scripts (GitHub)
All scripts used to process and plot all figures (main and sup) are available in this repository under `scripts/r`.
Data corresponding to these scripts are linked and detailed below in the Data section. 

Several scripts related to mapping read, calling SNPs, and assembling genomes are also available in `scripts/bash`.

### Data (Zenodo)
Datasets required to reproduce figures and related analyses are organized by analysis type:
**[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15556890.svg)](https://doi.org/10.5281/zenodo.15556890)**

The Zenodo archive contains a `data_by_script/` folder with:
```text
data_by_script/

├── 01_GWAS/           # GWAS and stats across the genome
├── 02_Comparative/    # GENESPACE comparative genomics 
├── 03_PCA/           # PCA results for inversions
├── 04_Depth/         # Depth based analyses
├── 05_Phylogenetic/  # Phylogenetic trees and metadata
└── 06_Recombination/ # Recombination rate data
```

PLOS Bio also requires deposition of the numeric values underlying each plot, and those are available as well in the Zenodo archive in the folder `figure_data_export/`. The only exception to this the plots resulting from GENESPACE and RepeatOBserver analyses, however, one can reproduce these analyses and related plots with the scripts provided.  

## Quick Start

### 1. Download Data
```bash
# Download from Zenodo and extract
wget https://zenodo.org/record/YOUR_RECORD/files/data_by_script.zip
unzip data_by_script.zip
```

### 2. Run Analyses from Zenodo Deposited Datasets
```
# Run in order (each script uses its corresponding data folder)
source("scripts/01_GWAS_PopGen_Analyses.R")      # Uses data_by_script/01_GWAS/
source("scripts/02_Comparative_Genomic_Analyses.R") # Uses data_by_script/02_Comparative/
source("scripts/03_PCA_Structure_Analyses.R")    # Uses data_by_script/03_PCA/
source("scripts/04_Depth_Based_Analyses.R")      # Uses data_by_script/04_Depth/
source("scripts/05_Phylogenetic_Analysis.R")     # Uses data_by_script/05_Phylogenetic/
source("scripts/06_Recombination_Analysis.R")    # Uses data_by_script/06_Recombination/
```
