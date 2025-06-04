# GENESPACE Orthology Analysis for Amaranthus Genus
# Author: Julia Kreiner
# Date: May 27th
# Description: Comparative genomics analysis using GENESPACE to identify synteny,
#              orthologous groups, and chromosomal inversions across Amaranthus species
#              and haplotypes. Includes visualization of synteny patterns and statistical
#              analysis of orthologous group distributions.

# =============================================================================
# SETUP AND INITIALIZATION
# =============================================================================

# NOTE: Before running this script:
# 1. Open terminal and activate orthofinder conda environment
# 2. Start RStudio from that environment
# Commands:
# conda activate orthofinder
# open -na rstudio

library(GENESPACE)

# Initialize GENESPACE parameters
# Set working directory and path to MCScanX software
gpar <- init_genespace(
  wd = "data_byscript/02_Comparative/genespace/",     # Working directory with input files
  path2mcscanx = "~/MCScanX/"                       # Path to MCScanX installation
)

# Run the main GENESPACE analysis
# This performs orthology detection, synteny analysis, and prepares data structures
out <- run_genespace(gpar, overwrite = FALSE)       # Set overwrite=TRUE to rerun analysis

# Set synteny parameters for downstream analyses
gpar <- set_syntenyParams(out)

# Define chromosomes that need to be flipped due to complete inversions
# This corrects for assemblies where some chromosomes are in reverse orientation
invchr <- data.frame(
  genome = c("Acrutensus", "Acrutensus", "Acrutensus", "Acrutensus", "Acrutensus", "Acrutensus",
             "Atricolor", "Atricolor", "Atricolor", "Atricolor", "Atricolor",
             "Apalm", "Apalm", "Apalm"),
  chr = c("chr02B", "chr01", "chr02A", "chr09", "chr08", "chr04",          # A. cruentus inversions
          "Chr16", "Chr4", "Chr5", "Chr12", "Chr14",                       # A. tricolor inversions  
          "Scaffold_6", "Scaffold_8", "Scaffold_14"                        # A. palmeri inversions
  )
)

# Genome order for visualization:
# "Walpole_X","Walpole_Y","Nat_X","Nat_Y","Rudis_X","Rudis_Y"

# =============================================================================
# GENOME-WIDE SYNTENY VISUALIZATION
# =============================================================================

# FIGURE 2A
# Create full genome synteny plot (riparian plot)
# Shows synteny relationships across all chromosomes for all species/haplotypes
ripd <- plot_riparian(
  gsParam = out,
  useRegions = FALSE,                               # Use full chromosomes, not specific regions
  invertTheseChrs = invchr,                        # Apply chromosome inversions
  backgroundColor = NULL,                           # Default background
  genomeIDs = c("Atub_Nune5_hap2", "Atub_Nune5_hap1",    # Nune population haplotypes
                "Atub_Nat_hap2", "Atub_Nat_hap1",         # Natural population haplotypes  
                "Atub_193_hap2", "Atub_193_hap1",         # Walpole (193) population haplotypes
                "Atricolor", "Apalm", "Acrutensus"),       # Outgroup species
  refGenome = "Atub_193_hap2",                     # Reference genome for alignment
  inversionColor = "white"                          # Color for inverted regions
)

# =============================================================================
# SEX-LINKED REGION (SLR) FOCUSED ANALYSIS
# =============================================================================

# FIGURE 2C
# Define region of interest: Scaffold_1 containing the sex-linked region
roi <- data.frame(
  genome = c("Atub_193_hap2"),
  chr = c("Scaffold_1"), 
  start = 0, 
  end = 60000000                                    # 60 Mb region (full scaffold)
)

# Create focused synteny plot for sex-linked region
ripd <- plot_riparian(
  gsParam = out,
  useRegions = FALSE,
  reorderBySynteny = TRUE,                          # Reorder chromosomes by synteny strength
  invertTheseChrs = invchr, 
  highlightBed = roi,                               # Highlight the sex-linked region
  backgroundColor = NULL,
  
  genomeIDs = c("Atub_Nune5_hap2","Atub_Nune5_hap1",
                "Atub_Nat_hap2", "Atub_Nat_hap1",
                "Atub_193_hap2", "Atub_193_hap1",
                "Atricolor", "Apalm", "Acrutensus"),
  refGenome = "Atub_193_hap2",
  inversionColor = "white"
)

# =============================================================================
# HAPLOTYPE-SPECIFIC SYNTENY ANALYSIS
# =============================================================================

# ---- NUNE REFERENCE PHASED HAPLOTYPES ----
# Compare synteny between Nune haplotype 1 and 2 on Scaffold_1
roi <- data.frame(
  genome = c("Atub_Nune5_hap2"), 
  chr = c("Scaffold_1")
)

# Query syntenic hits for the region of interest
qreturn <- query_hits(gsParam = out, bed = roi, synOnly = TRUE)

# Extract synteny between Nune haplotypes
hap1 <- subset(qreturn[["Atub_Nune5_hap2, Scaffold_1: 0-Inf"]], genome2 == "Atub_Nune5_hap1")

# Visualize synteny hits colored by syntenic blocks
q1 <- gghits(hap1, useOrder = FALSE, colorByBlks = TRUE)
#y on the X axis, X on the Y axis

# Check unique syntenic block IDs
unique(hap1$blkID)

# ---- WALPOLE (193) PHASED HAPLOTYPES ----
# Compare synteny between Walpole haplotype 1 and 2 on Scaffold_1
roi <- data.frame(
  genome = c("Atub_193_hap2"),
  chr = c("Scaffold_1")
)

# Extract synteny between Walpole haplotypes
hap1 <- subset(qreturn[["Atub_193_hap2, Scaffold_1: 0-Inf"]], genome2 == "Atub_193_hap1")

# Visualize and analyze syntenic blocks
q1 <- gghits(hap1, useOrder = FALSE, colorByBlks = TRUE)
#y on the X axis, X on the Y axis

print(q1)
unique(hap1$blkID)

# ---- NATURAL (NAT) PHASED HAPLOTYPES ----
# Compare synteny between Natural population haplotypes on Scaffold_1
roi <- data.frame(
  genome = c("Atub_Nat_hap2"),
  chr = c("Scaffold_1")
  # Commented regions show specific coordinate ranges that could be analyzed:
  # start = 3.8e07, end = 4.7e07,                  # 38-47 Mb region
  # start = c(40.5e06, 44.5e06, 46e06),           # Multiple specific regions
  # end = c(43.5e06, 46e06, 47e06),
  # color = c("green", "yellow", "orange")          # Color coding for regions
)

# Query and extract synteny between Natural haplotypes
qreturn <- query_hits(gsParam = out, bed = roi, synOnly = TRUE)
hap1 <- subset(qreturn[["Atub_Nat_hap2, Scaffold_1: 0-Inf"]], genome2 == "Atub_Nat_hap1")

# Visualize syntenic blocks
q1 <- gghits(hap1, useOrder = FALSE, colorByBlks = TRUE)
unique(hap1$blkID)
#y on the X axis, X on the Y axis

# =================================
### FIGURE 1C - Focal Y vs X dotplot
# =================================

roi <- data.frame(
  genome = c("Atub_193_hap2"), 
  chr = c("Scaffold_1"), 
  start= 3.25e07, end=4.7e07)
qreturn <- query_hits(gsParam = out, bed = roi, synOnly = TRUE)
hap1 <- subset(qreturn[["Atub_193_hap2, Scaffold_1: 32500000-4.7e+07"]], genome2 == "Atub_193_hap1")

gghits2 <- function(hits, colorByBlks = TRUE, alpha = ifelse(colorByBlks, 
                                                             1, 0.25), useOrder = TRUE, minScore = 0, minGenes2plot = 0, pointSize = 2) 
{
  ofID1 <- ofID2 <- sameOg <- ngene1 <- ngene2 <- ord1 <- ord2 <- blkID <- inBuffer <- rnd2 <- rnd1 <- n <- isArrayRep2 <- isArrayRep1 <- chr1 <- noAnchor <- bitScore <- quantile <- chr2 <- sameOG <- isAnchor <- start1 <- start2 <- x <- y <- NULL
  tp <- data.table(hits)
  
  if (colorByBlks) {
    hc <- subset(tp, !is.na(blkID) & isAnchor)
  } else {
    hc <- subset(tp, bitScore > minScore)
  }
  
  ng1 <- as.integer(uniqueN(hc$ofID1))
  ng2 <- as.integer(uniqueN(hc$ofID2))
  
  hc[, `:=`(ngene1, uniqueN(ofID1[!noAnchor & isArrayRep1], na.rm = TRUE)), by = "chr1"]
  hc[, `:=`(ngene2, uniqueN(ofID2[!noAnchor & isArrayRep2], na.rm = TRUE)), by = "chr2"]
  hc <- subset(hc, ngene1 > minGenes2plot & ngene2 > minGenes2plot)
  
  if (useOrder) {
    hc[, `:=`(x = ord1, y = ord2)]
    xlab <- "query gene rank order position"
    ylab <- "target gene rank order position"
  } else {
    hc[, `:=`(x = start1 / 1e+06, y = start2 / 1e+06)]
    xlab <- "Y Chromosome (Mb) gene position"
    ylab <- "X Chromosome (Mb)\n gene position"
  }
  
  if (!colorByBlks) {
    # Define a custom color palette
    colorPalette <- c("blue", "red")  # Alternating colors
    
    p <- ggplot(hc, aes(x = x, y = y)) + 
      geom_point(alpha = alpha, size = pointSize, aes(color = factor(row_number() %% 2))) +  # Alternate colors based on row number
      scale_color_manual(values = blkCols, 
                         guide = "none") +  # Apply custom color palette
      scale_x_continuous(expand = c(0, 0), breaks = pretty(hc$x, n = 10)) + 
      scale_y_continuous(expand = c(0, 0), breaks = seq(floor(min(hc$y)), ceiling(max(hc$y)), by = 1)) +  # Increase by 1
      facet_grid(genome2 + chr2 ~ genome1 + chr1, scales = "free", space = "free", as.table = FALSE) + 
      labs(x = xlab, y = ylab) + 
      theme_bw() +
      theme(
        panel.grid.minor = element_blank(), 
        # panel.grid.major = element_line(color = "gray80", size = 0.5, linetype = 2),  # Adjust color for dashed lines
        panel.spacing = unit(0.1, "mm"), 
        axis.ticks = element_line(color = "black"),  # Make axis ticks visible
        axis.line = element_line(color = "black", size = 0.5),  # Make axis lines visible
        strip.background = element_blank(), 
        axis.text.y = element_text(family = "Helvetica", size = 8, color = "black"), 
        axis.text.x = element_text(family = "Helvetica", size = 8, color = "black"),  # Make x-axis text color black
        axis.title.x = element_text(family = "Helvetica", size = 12, color = "black"),  # Set x-axis title color to black
        axis.title.y = element_text(family = "Helvetica", size = 12, color = "black"),  # Set y-axis title color to black
        plot.title = element_blank(),
        strip.text = element_blank(),
        plot.margin = unit(c(0, 0.0, 0, 0), "cm")
      ) +       theme_bw() 
  } else {
    blkCols <- sample(gs_colors(uniqueN(hc$blkID)))
    p <- ggplot(hc, aes(x = x, y = y, col = blkID)) + 
      geom_point(alpha = alpha, size = pointSize) + 
      scale_color_manual(values = c("black","black","black","#db8d40","#db8d40"), guide = "none") + 
      scale_x_continuous(expand = c(0, 0), breaks = pretty(hc$x, n = 10)) + 
      scale_y_continuous(expand = c(0, 0), breaks = seq(floor(min(hc$y)), ceiling(max(hc$y)), by = 1)) +  # Increase by 1
      facet_grid(genome2 + chr2 ~ genome1 + chr1, scales = "free", space = "free", as.table = FALSE) + 
      labs(x = xlab, y = ylab) + 
      theme_bw() +
      theme(
        strip.background = element_blank(), 
        strip.text.x = element_blank(),
        strip.text = element_blank(),
        plot.title = element_blank()
        
      ) 
  }
  
  print(p)
}

q1<-gghits2(hap1, useOrder = F,colorByBlks = T) #FIG 1C

# =============================================================================
# ORTHOLOGOUS GROUP ANALYSIS ACROSS ALL SCAFFOLDS
# Figure 2B
# =============================================================================

library(tidyverse)

# Define all scaffolds for comprehensive orthology analysis
roi <- data.frame(
  genome = c("Atub_193_hap2"),
  chr = c("Scaffold_1", "Scaffold_2", "Scaffold_3", "Scaffold_4", "Scaffold_5", "Scaffold_6",
          "Scaffold_7", "Scaffold_8", "Scaffold_9", "Scaffold_10", "Scaffold_11", "Scaffold_12",
          "Scaffold_13", "Scaffold_14", "Scaffold_15", "Scaffold_16")
)

# Query syntenic hits across all scaffolds
qreturn <- query_hits(gsParam = out, bed = roi, synOnly = TRUE)

# Combine all results into single data frame
concatenated_df <- do.call(rbind, qreturn)
head(concatenated_df)

# Parse to examine within-scaffold synteny (Scaffold_1 to Scaffold_1)
check <- concatenated_df %>% 
  filter(chr1 == "Scaffold_1" & chr2 == "Scaffold_1")

# Get unique orthologous groups (remove duplicates by block ID)
check_uniqOG <- check[!duplicated(check$blkID), ]

# =============================================================================
# ORTHOLOGOUS GROUP COUNTING AND ANALYSIS
# =============================================================================

# Count orthologous groups by genome comparison and chromosome
OGcounts <- concatenated_df %>%
  group_by(genome1, genome2, chr1) %>%
  summarise(
    OGcount = length(unique(blkID)),               # Count unique orthologous groups
    length = max(end1),                            # Chromosome/scaffold length
    .groups = 'drop'
  )

# Preliminary visualization (excluding outgroup species)
# Shows distribution of orthologous groups across chromosomes
OGcounts %>% 
  filter(genome2 != "Apalm" & genome2 != "Acrutensus" & 
           genome2 != "Atricolor" & genome2 != "Atub_193_hap2") %>%
  ggplot(aes(OGcount, reorder(chr1, length), color = genome2)) +
  geom_point() +
  theme_bw()

# Calculate summary statistics by chromosome
summary_data <- OGcounts %>%
  group_by(chr1) %>%
  summarize(
    OGcount_median = median(OGcount),
    OGcount_upper = quantile(OGcount, 0.75),      # 75th percentile
    OGcount_lower = quantile(OGcount, 0.25),      # 25th percentile
    .groups = 'drop'
  )
print(summary_data)

# Overall median orthologous group count (excluding outgroups and self-comparisons)
OGcounts %>% 
  filter(genome2 != "Apalm" & genome2 != "Acrutensus" & 
           genome2 != "Atricolor" & genome2 != "Atub_193_hap2") %>%
  summarize(OGcount_median = median(OGcount))

# Extract numeric chromosome identifiers for proper ordering
OGcounts2 <- OGcounts %>% 
  mutate(chr_numeric = as.numeric(str_extract(chr1, "\\d+")))

# =============================================================================
# FIGURE 2B: ORTHOLOGOUS GROUP DISTRIBUTION VISUALIZATION
# =============================================================================

# Create publication-quality boxplot showing orthologous group counts by chromosome
OGcounts2 %>%
  filter(genome2 != "Apalm" & genome2 != "Acrutensus" & 
           genome2 != "Atricolor" & genome2 != "Atub_193_hap2") %>%
  ggplot() +
  geom_boxplot(aes(y = as.factor(chr_numeric), x = OGcount)) +           # Boxplot by chromosome
  geom_jitter(aes(y = as.factor(chr_numeric), x = OGcount, color = genome2), 
              width = 0.1) +                                              # Individual points colored by genome
  theme_bw() +
  labs(y = "Chromosome", x = "Count of\nOrthologous Groups") +
  scale_color_viridis_d(option = "E") +                                  # Viridis color palette
  coord_flip() +                                                         # Flip coordinates
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# =============================================================================
# STATISTICAL ANALYSIS OF ORTHOLOGOUS GROUP VARIATION
# =============================================================================

library(car)

# Filter data for statistical analysis (exclude outgroups and self-comparisons)
OGless <- OGcounts %>% 
  filter(genome2 != "Apalm" & genome2 != "Acrutensus" & 
           genome2 != "Atricolor" & genome2 != "Atub_193_hap2")

# Linear model testing effects of chromosome, length, and genome on orthologous group counts
# Tests whether certain chromosomes have systematically different numbers of orthologous groups
lm1_reduced <- lm(data = OGless, OGcount ~ chr1 + genome2) #genome2 is paired assembly, length is chr length, chr1 is scaffold
lm1 <- lm(data = OGless, OGcount ~ chr1 + length + genome2) #genome2 is paired assembly, length is chr length, chr1 is scaffold

# Model summary and ANOVA
Anova(lm1, type = 2)                               # Type II ANOVA for balanced design
summary(lm1)
Anova(lm1_reduced, type = 2)                               # Type II ANOVA for balanced design
summary(lm1_reduced)


# =============================================================================
# POST-HOC ANALYSIS: CHROMOSOME-SPECIFIC DIFFERENCES
# =============================================================================

library(lsmeans)

# Calculate least squares means for each chromosome (corrected for scaffold length)
# This shows whether Scaffold_1 (sex chromosome) differs from autosomes
# after accounting for chromosome length differences
test <- data.frame(lsmeans(lm1, "chr1"))

# Visualize chromosome-specific means with confidence intervals
ggplot(data = test, aes(lsmean, chr1)) +
  geom_point() +                                   # Point estimates
  theme_bw() +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL)) +    # Confidence intervals
  labs(x = "Least Squares Mean\n(Orthologous Group Count)", 
       y = "Chromosome/Scaffold",
       title = "Orthologous Group Counts by Chromosome\n(Adjusted for Length)")

# ==============================================================
# Test for gene convservation across SLR regions on Scaffold 1
# ==============================================================

#query pan genes
roi <- data.frame(
  genome = c("Atub_193_hap2"), 
  chr = c("Scaffold_1"))

test <- query_pangenes(
  gsParam = out, bed = roi,transform = T)

pangen_df<-test$`Atub_193_hap2, Scaffold_1: 0-Inf`
pangen_df_atub<-pangen_df[grep("Atub_193_hap2",pangen_df$genome) ,]
pangen_df_atub_scaf1<-pangen_df_atub[grep("Scaffold_1",pangen_df_atub$chr) ,]

#########
#library(data.table)
#
## Convert both datasets to data.table if not already
#setDT(gwa_bon)
#setDT(pangen_df_atub)
#
## Set keys for efficient joining
#setkey(gwa_bon, chr, ps)
#setkey(pangen_df_atub, chr, start, end)
#
#
## Perform overlap join - find SNPs that fall within gene boundaries
## This finds all SNPs where: gene_start <= SNP_position <= gene_end on same chromosome
#overlaps <- foverlaps(
#  gwa_bon[, .(chr, start = ps, end = ps, rs, ps, af, beta, se, p_bon)],  # SNP intervals (point intervals)
#  pangen_df_atub,                                                        # Gene intervals (all columns kept)
#  by.x = c("chr", "start", "end"),
#  by.y = c("chr", "start", "end"),
#  type = "within"  # SNP position must be within gene boundaries
#)
#
## Remove rows where no overlap was found (NA values) and keep only gene columns
#overlaps_clean <- overlaps[!is.na(pgID)]
#
## Get unique gene IDs that have overlapping SNPs
#genes_with_snps <- unique(overlaps_clean$pgID)
#
## Filter original gene dataframe to only genes with overlapping SNPs
#genes_filtered <- pangen_df_atub[pgID %in% genes_with_snps]
#
## Result: unique gene records that have at least one overlapping significant SNP
#result <- unique(genes_filtered)


#coordinates of SLR  including upstream and downstream of main tract on scaffold 1) end and start
pangen_df_atub_sdr<-pangen_df_atub[pangen_df_atub$start < 43548508 & pangen_df_atub$end > 40.5e06 & pangen_df_atub$chr == "Scaffold_1" |
                                     pangen_df_atub$start < 33647320 & pangen_df_atub$end > 32285336 & pangen_df_atub$chr == "Scaffold_1"|
                                   pangen_df_atub$start < 44592427 & pangen_df_atub$end > 44444604 & pangen_df_atub$chr == "Scaffold_1",]
names(pangen_df_atub_sdr)[13:18]<-c("Walpole_X","Walpole_Y","Nat_X","Nat_Y","Rudis_X","Rudis_Y")
#count number of missing genes
# Assuming your data frame is called 'df'
count_non_blanks <- function(row) {
  sum(row != "character(0)")
}


# consider only other Ys
names(pangen_df_atub_sdr)[13:18]
columns_to_check <- c(16,18)
#columns_to_check <- c(13)

Ys<-pangen_df_atub_sdr[,..columns_to_check]
# Apply the function row-wise
count_of_non_blanks <- apply(pangen_df_atub_sdr[,..columns_to_check], 1, count_non_blanks)
sum(count_of_non_blanks<2)
sum(count_of_non_blanks==2)
count_df<-as.data.frame(count_of_non_blanks)

SL_pan<-ggplot(data=count_df, aes(count_of_non_blanks)) +
  geom_bar(color="black",fill="lightgrey") +
  theme_bw() +
  labs(x="Number of alternate Y assemblies\nwhere focal gene is present",
       title="Chromosome 1 Sex-Linked Genes")


# Define counts
(SLR_core <- sum(count_of_non_blanks==2))
(SLR_total <-  length(count_of_non_blanks))
SLR_core

propcore = (SLR_core <- sum(count_of_non_blanks==2))/(SLR_total <-  length(count_of_non_blanks))
propcore

# do the same but for genome wide
count_of_non_blanks <- apply(pangen_df_atub[,c(16,18)], 1, count_non_blanks)
sum(count_of_non_blanks<2)
sum(count_of_non_blanks==2)

count_df<-as.data.frame(count_of_non_blanks)
head(count_df)

gw_pan<-ggplot(data=count_df, aes(count_of_non_blanks)) +
  geom_bar(color="black",fill="lightgrey") +
  theme_bw() +
  labs(x="Number of alternate Y assemblies\nwhere focal gene is present", title="Genome-Wide")

(Genome_core <- sum(count_of_non_blanks==2))
(Genome_total <-  length(count_of_non_blanks))

propcore_Genome = (Genome_core <- sum(count_of_non_blanks==2))/(Genome_total <-  length(count_of_non_blanks))
propcore_Genome

# Proportion test (two-sample test for equality of proportions)
prop.test(c(SLR_core, Genome_core), c(SLR_total, Genome_total), correct=FALSE)


#SUP Figure 18
plot_grid(SL_pan, gw_pan, nrow=1,ncol=2)

# ================================================================
#Gain of genes on the Y? Lets compare our focal Y to all other Xs
# ================================================================

#query pan genes
roi <- data.frame(
  genome = c("Atub_193_hap2"), 
  chr = c("Scaffold_1"), 
  # start=30e06, end=50e06,
  #start=c(30e06),end=c(60e06), 
  color=c("yellow"))

test <- query_pangenes(
  gsParam = out, bed = roi,transform = T)

pangen_df<-test$`Atub_193_hap2, Scaffold_1: 0-Inf`
pangen_df_atub<-pangen_df[grep("Atub_193_hap2",pangen_df$genome) ,]
pangen_df_atub_scaf1<-pangen_df_atub[grep("Scaffold_1",pangen_df_atub$chr) ,]
pangen_df_atub_sdr<-pangen_df_atub[pangen_df_atub$start < 43548508 & pangen_df_atub$end > 40.5e06 & pangen_df_atub$chr == "Scaffold_1" |
                                     pangen_df_atub$start < 33647320 & pangen_df_atub$end > 32285336 & pangen_df_atub$chr == "Scaffold_1"|
                                     pangen_df_atub$start < 44592427 & pangen_df_atub$end > 44444604 & pangen_df_atub$chr == "Scaffold_1",]

#pangen_df_atub_sdr<-pangen_df_atub[pangen_df_atub$start < 44585719 & pangen_df_atub$start > 40.5e06 & pangen_df_atub$chr == "Scaffold_1",]
names(pangen_df_atub_scaf1)[13:18]<-c("Walpole_X","Walpole_Y","Nat_X","Nat_Y","Rudis_X","Rudis_Y")
#count number of missing genes
# Assuming your data frame is called 'df'
count_non_blanks <- function(row) {
  sum(row != "character(0)")
}

# Columns to consider
names(pangen_df_atub_scaf1)[13:18]
columns_to_check <- c(13,15,17)

# Apply the function row-wise
count_of_non_blanks <- apply(pangen_df_atub_sdr[,..columns_to_check], 1, count_non_blanks)
hist(count_of_non_blanks,nclass=3, xlab="Number of X assemblies\nwhere Y gene is present",main="")
count_df<-as.data.frame(count_of_non_blanks)

(YgenesabsentfromX=sum(count_of_non_blanks==0))
(allY=length(count_of_non_blanks))

YgenesabsentfromX/allY
