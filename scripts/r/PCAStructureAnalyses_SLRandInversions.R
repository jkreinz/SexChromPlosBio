# PCA Analyses for SLR and for Inversions segregating within A. Tuberculatus, based on coordinates genespace analyses
# Author: Julia Kreiner
# Date: May 27
# Description: Automated analysis of PCA patterns across multiple inversions
#              with statistical testing for PC1, PC2, and PC3

# Load required libraries
library(data.table)
library(ggplot2)
library(dplyr)
library(cowplot)
library(car)
library(viridis)

# =============================================================================
# DATA PREPARATION (Run once before loop)
# =============================================================================

# Read PLINK order file
plink <- read.table("~/herb_commongarden_pnas_missing9_maf05.fam")
names(plink)[1] <- "sample"
plink$order <- seq(1, nrow(plink))

# Read comprehensive sample information
allmeta <- read.table("~/3waymerged_sampleinfo.txt", sep = "\t", header = TRUE)

# Join with PLINK order
both2 <- inner_join(plink, allmeta, by = "sample")

# Read metadata with ancestry information
meta_ordered <- read.table("~/Documents/UToronto_PhD/Papers/Evolution Preadapt Submission/Evolution2021-Preadapt/commongarden_and_PNAS_metadata_ancestry.txt", 
                           header = TRUE)
names(meta_ordered)[1] <- "sample"


############################

e_val<-read.table("commongarden_allfiltsnps_193_hap2_remvdrepeats_outliersdropped_scaf1inversion.eigenval")
contemp_herb<-fread("~/commongarden_allfiltsnps_193_hap2_remvdrepeats_outliersdropped_scaf1inversion.eigenvec")

names(contemp_herb)[2]<-"sample"
contemp_herb$sample<-as.character(contemp_herb$sample)
contemp_herb_pca<-inner_join(contemp_herb, allmeta, by="sample")
contemp_herb_pca_anc<-inner_join(contemp_herb_pca,meta_ordered,by="sample")

eval_props<-e_val$V1/sum(e_val$V1)


p1<-contemp_herb_pca_anc %>% #filter(year==2019) %>%
  ggplot(aes(PC1, PC2, color=sex,shape=sex)) +
  geom_point(alpha=.9, size=2) +
  theme_bw() +
  #lims(y=c(-0.05,0.15)) +
  labs(x=paste("PC1 (", round(eval_props[1],digits = 2)*100, "%)",sep=""),
       y=paste("PC2 (", round(eval_props[2],digits = 2)*100, "%)",sep=""),
       color="var. rudis\nancestry\nproportion") +
  scale_color_viridis_d(option="magma",end = .92) + guides(color=F, shape=F) #+
theme(plot.margin=margin(0.2,0.2,0.2,0.2,"cm")) 
p1


p4<-contemp_herb_pca_anc %>% #filter(year==2019) %>%
  ggplot(aes(PC1, PC2, color=K1,shape=sex)) +
  geom_point(alpha=.9, size=2) +
  theme_bw() +
  #lims(y=c(-0.05,0.15)) +
  labs(x=paste("PC1 (", round(eval_props[1],digits = 2)*100, "%)",sep=""),
       y=paste("PC2 (", round(eval_props[2],digits = 2)*100, "%)",sep=""),
       color="var. rudis\nancestry\nproportion") +
  scale_color_viridis_c() + guides(color=F, shape=F) #+

#Figure 3B
plot_grid(p1,p4, nrow=2)



# =============================================================================
# INVERSION ANALYSIS LOOP
# =============================================================================

# Initialize storage for results
pca_plots_list <- list()
anova_results <- list()
significant_terms <- data.frame()
marginal_terms <- data.frame()

# Loop through inversions 1-18
for (inversion_num in 1:18) {
  
  cat("\n=== PROCESSING INVERSION", inversion_num, "===\n")
  
  # Construct file paths
  eigenvec_file <- paste0("inversion", inversion_num, ".eigenvec")
  eigenval_file <- paste0("inversion", inversion_num, ".eigenval")
  
  # Check if files exist before processing
  if (!file.exists(eigenvec_file) || !file.exists(eigenval_file)) {
    cat("Warning: Files for inversion", inversion_num, "not found. Skipping...\n")
    next
  }
  
  # Read PCA results
  contemp_herb <- fread(eigenvec_file)
  names(contemp_herb)[2] <- "sample"
  contemp_herb$sample <- as.character(contemp_herb$sample)
  
  # Join with metadata
  contemp_herb_pca <- inner_join(contemp_herb, allmeta, by = "sample")
  contemp_herb_pca_anc <- inner_join(contemp_herb_pca, meta_ordered, by = "sample")
  
  # Read eigenvalues and calculate variance explained
  e_val <- read.table(eigenval_file)
  eval_props <- e_val$V1 / sum(e_val$V1)
  
  # ==========================================================================
  # CREATE PCA PLOTS
  # ==========================================================================
  
  # PC1 vs PC2 (by sex)
  p1 <- contemp_herb_pca_anc %>%
    ggplot(aes(PC1, PC2, color = sex, shape = sex)) +
    geom_point(alpha = .9, size = 2) +
    theme_bw() +
    labs(x = paste("PC1 (", round(eval_props[1], digits = 2) * 100, "%)", sep = ""),
         y = paste("PC2 (", round(eval_props[2], digits = 2) * 100, "%)", sep = "")) +
    scale_color_viridis_d(option = "magma", end = .92) +
    guides(color = FALSE, shape = FALSE)
  
  # PC1 vs PC3 (by sex)
  p2 <- contemp_herb_pca_anc %>%
    ggplot(aes(PC1, PC3, color = sex, shape = sex)) +
    geom_point(alpha = .9, size = 2) +
    theme_bw() +
    labs(x = paste("PC1 (", round(eval_props[1], digits = 2) * 100, "%)", sep = ""),
         y = paste("PC3 (", round(eval_props[3], digits = 2) * 100, "%)", sep = "")) +
    scale_color_viridis_d(option = "magma", end = .92) +
    guides(color = FALSE, shape = FALSE)
  
  # PC2 vs PC3 (by sex)
  p3 <- contemp_herb_pca_anc %>%
    ggplot(aes(PC2, PC3, color = sex, shape = sex)) +
    geom_point(alpha = .9, size = 2) +
    theme_bw() +
    labs(x = paste("PC2 (", round(eval_props[2], digits = 2) * 100, "%)", sep = ""),
         y = paste("PC3 (", round(eval_props[3], digits = 2) * 100, "%)", sep = "")) +
    scale_color_viridis_d(option = "magma", end = .92) +
    guides(color = FALSE, shape = FALSE)
  
  # PC1 vs PC2 (by ancestry - K1)
  p4 <- contemp_herb_pca_anc %>%
    ggplot(aes(PC1, PC2, color = K1, shape = sex)) +
    geom_point(alpha = .9, size = 2) +
    theme_bw() +
    labs(x = paste("PC1 (", round(eval_props[1], digits = 2) * 100, "%)", sep = ""),
         y = paste("PC2 (", round(eval_props[2], digits = 2) * 100, "%)", sep = ""),
         color = "var. rudis\nancestry\nproportion") +
    scale_color_viridis_c() +
    guides(color = FALSE, shape = FALSE)
  
  # Combine plots
  combined_plot <- plot_grid(p1, NULL, p2, p3, nrow = 2)
  pca_plots_list[[paste0("inversion_", inversion_num)]] <- combined_plot
  
  # ==========================================================================
  # STATISTICAL ANALYSIS (ANOVA) - PC1, PC2, PC3
  # ==========================================================================
  
  # Check if required columns exist
  required_cols <- c("sex", "lat", "long", "env", "K1")
  missing_cols <- required_cols[!required_cols %in% names(contemp_herb_pca_anc)]
  
  if (length(missing_cols) > 0) {
    cat("Warning: Missing columns for inversion", inversion_num, ":", 
        paste(missing_cols, collapse = ", "), "\n")
    next
  }
  
  # Perform ANOVA for each PC
  tryCatch({
    anova_pc1 <- Anova(lm(data = contemp_herb_pca_anc, 
                          PC1 ~ sex + lat + long + env + K1), type = 2)
    anova_pc2 <- Anova(lm(data = contemp_herb_pca_anc, 
                          PC2 ~ sex + lat + long + env + K1), type = 2)
    anova_pc3 <- Anova(lm(data = contemp_herb_pca_anc, 
                          PC3 ~ sex + lat + long + env + K1), type = 2)
    
    # Store ANOVA results
    anova_results[[paste0("inversion_", inversion_num)]] <- list(
      PC1 = anova_pc1, PC2 = anova_pc2, PC3 = anova_pc3)
    
    # Process results for PC1, PC2, PC3
    pc_names <- c("PC1", "PC2", "PC3")
    anova_list <- list(anova_pc1, anova_pc2, anova_pc3)
    
    for (i in 1:3) {
      pc_name <- pc_names[i]
      anova_result <- anova_list[[i]]
      
      # Find significant (p < 0.05) and marginally significant (0.05 <= p < 0.1) terms
      sig_rows <- which(anova_result$`Pr(>F)` < 0.05)
      marginal_rows <- which(anova_result$`Pr(>F)` >= 0.05 & anova_result$`Pr(>F)` < 0.1)
      
      # Print results for this PC
      cat(pc_name, "results:\n")
      
      if (length(sig_rows) > 0) {
        sig_terms <- data.frame(
          Inversion = inversion_num, PC = pc_name,
          Term = rownames(anova_result)[sig_rows],
          F_value = anova_result$`F value`[sig_rows],
          P_value = anova_result$`Pr(>F)`[sig_rows],
          Significance = ifelse(anova_result$`Pr(>F)`[sig_rows] < 0.001, "***",
                                ifelse(anova_result$`Pr(>F)`[sig_rows] < 0.01, "**", "*"))
        )
        significant_terms <- rbind(significant_terms, sig_terms)
        cat("  Significant (p < 0.05):\n")
        print(sig_terms[, c("Term", "F_value", "P_value", "Significance")])
      }
      
      if (length(marginal_rows) > 0) {
        marg_terms <- data.frame(
          Inversion = inversion_num, PC = pc_name,
          Term = rownames(anova_result)[marginal_rows],
          F_value = anova_result$`F value`[marginal_rows],
          P_value = anova_result$`Pr(>F)`[marginal_rows],
          Significance = "."
        )
        marginal_terms <- rbind(marginal_terms, marg_terms)
        cat("  Marginally significant (0.05 ≤ p < 0.1):\n")
        print(marg_terms[, c("Term", "F_value", "P_value", "Significance")])
      }
      
      if (length(sig_rows) == 0 && length(marginal_rows) == 0) {
        cat("  No significant or marginally significant factors\n")
      }
      cat("\n")
    }
    
  }, error = function(e) {
    cat("Error in ANOVA for inversion", inversion_num, ":", e$message, "\n")
  })
  
  # Save plots
  ggsave(filename = paste0("inversion_", inversion_num, "_pca_plots.png"),
         plot = combined_plot, width = 12, height = 8, dpi = 300)
  ggsave(filename = paste0("inversion_", inversion_num, "_ancestry_plot.png"),
         plot = p4, width = 8, height = 6, dpi = 300)
  
  cat("Completed inversion", inversion_num, "\n")
}

# =============================================================================
# SUMMARY RESULTS
# =============================================================================

cat("\n=== FINAL SUMMARY ===\n")

# Save results
if (nrow(significant_terms) > 0) {
  write.csv(significant_terms, "significant_terms_all_inversions.csv", row.names = FALSE)
  cat("Significant terms saved to: significant_terms_all_inversions.csv\n")
}

if (nrow(marginal_terms) > 0) {
  write.csv(marginal_terms, "marginal_terms_all_inversions.csv", row.names = FALSE)
  cat("Marginal terms saved to: marginal_terms_all_inversions.csv\n")
}

# Create summary by factor
if (nrow(significant_terms) > 0 || nrow(marginal_terms) > 0) {
  all_terms <- rbind(
    if (nrow(significant_terms) > 0) significant_terms %>% mutate(Type = "Significant") else data.frame(),
    if (nrow(marginal_terms) > 0) marginal_terms %>% mutate(Type = "Marginal") else data.frame()
  )
  
  summary_by_inversion <- all_terms %>%
    group_by(Inversion, PC) %>%
    summarise(
      Significant_count = sum(Type == "Significant"),
      Marginal_count = sum(Type == "Marginal"),
      Terms = paste(Term, collapse = ", "),
      .groups = 'drop'
    ) %>%
    arrange(Inversion, PC)
  
  write.csv(summary_by_inversion, "summary_by_inversion_and_PC.csv", row.names = FALSE)
  cat("Summary by inversion and PC saved to: summary_by_inversion_and_PC.csv\n")
}

save(pca_plots_list, anova_results, significant_terms, marginal_terms,
     file = "inversion_analysis_results.RData")

# note that Inversion 5 is the one directly upstream of the SLR
# plot "inversion_5_pca_plots.png --> Sup Fig 17
