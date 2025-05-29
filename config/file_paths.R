# =============================================================================
# File Paths Configuration - Amaranthus Sex Chromosome Analysis
# =============================================================================
# This file centralizes all data file paths for the analysis scripts.
# Users only need to update DATA_BASE_DIR to point to their downloaded data.

library(here)

# =============================================================================
# BASE DIRECTORIES - UPDATE THIS SECTION
# =============================================================================

# **USERS: Update this path to where you downloaded the Zenodo data**
# Default assumes data_by_script folder is in the same directory as this repo
DATA_BASE_DIR <- file.path(here::here(), "..", "data_by_script")

# Alternative paths (uncomment the one that matches your setup):
# DATA_BASE_DIR <- "~/Downloads/data_by_script"                    # Downloaded to Downloads
# DATA_BASE_DIR <- "/path/to/your/data_by_script"                  # Custom path
# DATA_BASE_DIR <- file.path(here::here(), "data_by_script")       # If you put data inside repo

# =============================================================================
# SCRIPT-SPECIFIC DIRECTORIES (Don't change these)
# =============================================================================

SCRIPT01_DIR <- file.path(DATA_BASE_DIR, "01_GWAS")
SCRIPT02_DIR <- file.path(DATA_BASE_DIR, "02_Comparative") 
SCRIPT03_DIR <- file.path(DATA_BASE_DIR, "03_PCA")
SCRIPT04_DIR <- file.path(DATA_BASE_DIR, "04_Depth")
SCRIPT05_DIR <- file.path(DATA_BASE_DIR, "05_Phylogenetic")
SCRIPT06_DIR <- file.path(DATA_BASE_DIR, "06_Recombination")

# Output directory for figures and results
OUTPUT_DIR <- file.path(here::here(), "outputs")
FIGURES_DIR <- file.path(OUTPUT_DIR, "figures")
TABLES_DIR <- file.path(OUTPUT_DIR, "tables")

# Create output directories if they don't exist
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# SCRIPT 01: GWAS AND POPULATION GENETICS FILES
# =============================================================================

SCRIPT01_FILES <- list(
  # Main GWAS results
  main_gwas = file.path(SCRIPT01_DIR, "commongarden_allfiltsnps_193_hap2_remvdrepeats_maf05.assoc.txt"),
  nomulti_gwas = file.path(SCRIPT01_DIR, "sexassocsnps_nomultimapping_recalled_metareorder.assoc.txt"),
  xy_gwas = file.path(SCRIPT01_DIR, "193_XYScaffolds_filteredsnps.assoc.txt"),
  
  # Genotype data
  tophit_genos = file.path(SCRIPT01_DIR, "tophit.012"),
  tophit_indv = file.path(SCRIPT01_DIR, "tophit.012.indv"),
  fam_file = file.path(SCRIPT01_DIR, "commongarden_allfiltsnps_193_hap2.fam"),
  raw_genos = file.path(SCRIPT01_DIR, "sexassocsnps_nomultimapping_recalled.raw"),
  bim_file = file.path(SCRIPT01_DIR, "sexassocsnps_nomultimapping_recalled.vcf.bim"),
  
  # Phenotype data
  sex_pheno = file.path(SCRIPT01_DIR, "sexpheno"),
  
  # Population genetics data
  fst_data = file.path(SCRIPT01_DIR, "male_female_1kb_fst.txt"),
  het_data = file.path(SCRIPT01_DIR, "male_female_hetprop.txt"),
  diversity_data = file.path(SCRIPT01_DIR, "0fold_allinds_allsites_100kb_pi.txt"),
  
  # Genomic features
  gene_density = file.path(SCRIPT01_DIR, "193_2_genes_100kbdensity.bed"),
  te_density = file.path(SCRIPT01_DIR, "193_2_TEs_100kbdensity.bed")
)

# =============================================================================
# SCRIPT 02: COMPARATIVE GENOMICS FILES
# =============================================================================

SCRIPT02_FILES <- list(
  # GENESPACE requires setup - these are the input directories
  bed_dir = file.path(SCRIPT02_DIR, "bed"),
  peptide_dir = file.path(SCRIPT02_DIR, "peptide")
)

# =============================================================================
# SCRIPT 03: PCA STRUCTURE FILES
# =============================================================================

SCRIPT03_FILES <- list(
  # Main PCA data
  plink_fam = file.path(SCRIPT03_DIR, "herb_commongarden_pnas_missing9_maf05.fam"),
  sample_info = file.path(SCRIPT03_DIR, "3waymerged_sampleinfo.txt"),
  ancestry_data = file.path(SCRIPT03_DIR, "commongarden_and_PNAS_metadata_ancestry.txt"),
  
  # Sex-linked region PCA
  scaf1_eigenval = file.path(SCRIPT03_DIR, "commongarden_allfiltsnps_193_hap2_remvdrepeats_outliersdropped_scaf1inversion.eigenval"),
  scaf1_eigenvec = file.path(SCRIPT03_DIR, "commongarden_allfiltsnps_193_hap2_remvdrepeats_outliersdropped_scaf1inversion.eigenvec")
)

# Function to get inversion file paths (inversion1-18)
get_inversion_files <- function(inversion_num) {
  list(
    eigenvec = file.path(SCRIPT03_DIR, paste0("inversion", inversion_num, ".eigenvec")),
    eigenval = file.path(SCRIPT03_DIR, paste0("inversion", inversion_num, ".eigenval"))
  )
}

# =============================================================================
# SCRIPT 04: DEPTH ANALYSIS FILES
# =============================================================================

SCRIPT04_FILES <- list(
  # Sample metadata
  pop_list = file.path(SCRIPT04_DIR, "herbcontemp_poplist.txt"),
  plink_fam = file.path(SCRIPT04_DIR, "herb_commongarden_pnas_missing9_maf05.fam"),
  sample_info = file.path(SCRIPT04_DIR, "3waymerged_sampleinfo.txt"),
  pca_data = file.path(SCRIPT04_DIR, "commongarden_allfiltsnps_193_hap2.eigenvec"),
  ancestry_data = file.path(SCRIPT04_DIR, "commongarden_and_PNAS_metadata_ancestry.txt"),
  mispheno_samples = file.path(SCRIPT04_DIR, "misphenoed_samples.txt"),
  
  # Depth matrices
  depth_100bp = file.path(SCRIPT04_DIR, "allsamples_numericorder_depthsby100bpwind_tab.txt"),
  depth_sex_hits = file.path(SCRIPT04_DIR, "allsamples_numericorder_sexhitsdepth_tab.txt"),
  
  # Phenotype data
  garden_phenotypes = file.path(SCRIPT04_DIR, "Common Garden Seq Stats - geno w pheno.csv")
)

# =============================================================================
# SCRIPT 05: PHYLOGENETIC FILES
# =============================================================================

SCRIPT05_FILES <- list(
  # Sample metadata (shared with other scripts)
  pop_list = file.path(SCRIPT05_DIR, "herbcontemp_poplist.txt"),
  plink_fam = file.path(SCRIPT05_DIR, "herb_commongarden_pnas_missing9_maf05.fam"), 
  sample_info = file.path(SCRIPT05_DIR, "3waymerged_sampleinfo.txt"),
  pca_data = file.path(SCRIPT05_DIR, "commongarden_allfiltsnps_193_hap2.eigenvec"),
  ancestry_data = file.path(SCRIPT05_DIR, "commongarden_and_PNAS_metadata_ancestry.txt"),
  
  # Phylogenetic tree
  tree_file = file.path(SCRIPT05_DIR, "commongarden_allsites_193_hap2_hemi.min4.phy.treefile")
)

# =============================================================================
# SCRIPT 06: RECOMBINATION FILES
# =============================================================================

SCRIPT06_FILES <- list(
  # Recombination rate data
  y_recomb_rates = file.path(SCRIPT06_DIR, "chromo_malesonly.Scaffold_1.results.txt")
)

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

# Check if all required files exist for a given script
check_script_files <- function(script_num) {
  script_files <- switch(script_num,
    "01" = SCRIPT01_FILES,
    "02" = SCRIPT02_FILES, 
    "03" = SCRIPT03_FILES,
    "04" = SCRIPT04_FILES,
    "05" = SCRIPT05_FILES,
    "06" = SCRIPT06_FILES,
  )
  
  if (is.null(script_files)) {
    cat("Invalid script number. Use 01-06.\n")
    return(FALSE)
  }
  
  missing_files <- c()
  for (name in names(script_files)) {
    file_path <- script_files[[name]]
    if (!file.exists(file_path)) {
      missing_files <- c(missing_files, paste0(name, ": ", basename(file_path)))
    }
  }
  
  if (length(missing_files) > 0) {
    cat("Script", script_num, "- Missing files:\n")
    cat(paste("  ✗", missing_files, collapse = "\n"))
    cat("\n")
    return(FALSE)
  } else {
    cat("Script", script_num, "- All files found! ✓\n")
    return(TRUE)
  }
}

# Check all files for all scripts
check_all_files <- function() {
  cat("Checking data files for all scripts...\n")
  cat("=====================================\n")
  
  all_good <- TRUE
  for (script_num in c("01", "02", "03", "04", "05", "06")) {
    script_good <- check_script_files(script_num)
    all_good <- all_good && script_good
  }
  
  if (all_good) {
    cat("\n🎉 All files found! Ready to run analyses.\n")
  } else {
    cat("\n❌ Some files are missing. Please check the paths above.\n")
    cat("💡 Tip: Update DATA_BASE_DIR at the top of this file.\n")
  }
  
  return(all_good)
}

# List the directory structure that should exist
show_expected_structure <- function() {
  cat("Expected data directory structure:\n")
  cat("=================================\n")
  cat("data_by_script/\n")
  cat("├── 01_GWAS/           # GWAS results, genotypes, population genetics\n")
  cat("├── 02_Comparative/    # GENESPACE bed/ and peptide/ folders\n")
  cat("├── 03_PCA/           # PCA eigenvectors/values, inversion1-18 files\n")
  cat("├── 04_Depth/         # Depth matrices, sample metadata, phenotypes\n")
  cat("├── 05_Phylogenetic/  # Tree files and sample metadata\n")
  cat("└── 06_Recombination/ # Recombination rate data\n")
}

# =============================================================================
# INITIALIZATION MESSAGE
# =============================================================================

cat("📁 File paths configuration loaded!\n")
cat("Base data directory:", DATA_BASE_DIR, "\n")

if (!dir.exists(DATA_BASE_DIR)) {
  cat("⚠️  Data directory not found!\n")
  cat("💡 Please:\n")
  cat("   1. Download data from Zenodo\n") 
  cat("   2. Update DATA_BASE_DIR at the top of this file\n")
  cat("   3. Run check_all_files() to verify\n")
} else {
  cat("✓ Data directory exists\n")
  cat("💡 Run check_all_files() to verify all required files are present\n")
}

cat("\nUseful functions:\n")
cat("- check_all_files()           # Check if all data files exist\n")
cat("- check_script_files('01')    # Check files for specific script\n") 
cat("- show_expected_structure()   # Show expected folder structure\n")
