# Recombination Rate Analysis for X and Y Chromosomes/Scaffolds
# Author: Julia Kreiner
# Date: May 27, 2025
# Description: Analysis of recombination landscapes across sex chromosomes

# Load required libraries
library(data.table)   # Fast data reading and manipulation
library(tidyverse)    # Data manipulation and visualization toolkit

# =============================================================================
# Y CHROMOSOME/SCAFFOLD ANALYSIS
# =============================================================================

# Read Y chromosome recombination data (males only)
recomb_Y <- fread("data_by_script/06_Recombination/chromo_malesonly.Scaffold_1.results.txt", sep = "\t")
head(recomb_Y)

# Standardize column names
names(recomb_Y) <- c("Scaf", "Loci", "Mean_rho", "Median_rho", "L95", "U95")

# Explore distribution of Y chromosome recombination rates
hist(recomb_Y$Mean_rho)

# Create windowed analysis for Y chromosome
winrecomb_Y <- recomb_Y %>%
  mutate(win = floor(Loci/100)) %>%
  group_by(win) %>%
  dplyr::summarise(mean_winrho = mean(Mean_rho))

# Plot Y chromosome recombination landscape with region boundaries
male_landscape <- ggplot(data = winrecomb_Y, aes(win/10, mean_winrho)) +
  geom_point() +
  geom_smooth(span = .1) +
  labs(y = "Mean 4NeR\n in 100kb Windows",
       x = "Y Genomic Position (Mb)") +
  theme_bw() +
  ylim(0, 125) +                          # Consistent y-axis with X chromosome
  geom_vline(xintercept = 40.5, lty = "dashed") +   # Mark region of interest
  geom_vline(xintercept = 43.5, lty = "dashed")     # boundaries
male_landscape

# Calculate cumulative recombination for Y chromosome genetic map
recomb_Y$cumulativerho <- cumsum(recomb_Y$Mean_rho)

# Create Y chromosome genetic map with multiple region markers
malemap <- recomb_Y %>%
  mutate(win = floor(Loci/100)) %>%
  group_by(win) %>%
  dplyr::summarise(mean_cumrho = mean(cumulativerho)) %>%
  ggplot(aes(win/10, mean_cumrho)) +
  geom_point() +
  geom_smooth(span = .1) +
  labs(y = "Rho",
       x = "Y Genomic Position (Mb)") +
  theme_bw() +
  # Primary region of interest (black dashed lines)
  geom_vline(xintercept = 40.5, lty = "dashed") +
  geom_vline(xintercept = 43.5, lty = "dashed") +
  # Secondary regions of interest (grey dashed lines)
  geom_vline(xintercept = 27, lty = "dashed", color = "darkgrey") +
  geom_vline(xintercept = 34, lty = "dashed", color = "darkgrey")
malemap

# =============================================================================
# DETAILED ANALYSIS OF SPECIFIC REGIONS
# =============================================================================

# Calculate overall mean recombination rate for Y chromosome
mean(recomb_Y$Mean_rho)

# Create high-resolution zoom plot of region of interest (38-47 Mb)
zoommap <- recomb_Y %>%
  filter(Loci > 3800 & Loci < 4700) %>%    # Filter to region of interest
  ggplot(aes(Loci/100, Mean_rho)) +
  geom_point(alpha = .4) +                 # Semi-transparent points
  geom_smooth(span = .1, method = "loess") +
  labs(y = "Mean 4NeR\n(Effective Recombination Rate)",
       x = "Y Genomic Position (Mb)") +
  theme_bw() +
  geom_vline(xintercept = 40.5, lty = "dashed") +
  geom_vline(xintercept = 43.5, lty = "dashed")

# Note: Commented line shows potential for multi-panel plot
#plot_grid(depthplot, zoommap, ncol = 1, align = "v", axis = "l")

##############
#comparison of recombination rates around the SLR
library(dplyr)

# Calculate mean and median recombination rates inside the primary region of interest
# Region: 40.5-43.5 Mb (40,500-43,500 bp)
mean_rho_inside <- recomb_Y %>%
  filter(Loci > 40500 & Loci < 43500) %>%
  summarise(mean_rho = mean(Mean_rho, na.rm = TRUE),
            median_rho = median(Mean_rho, na.rm = TRUE))

# Calculate mean and median recombination rates outside regions of interest
# Excludes both primary region (40.5-43.5 Mb) and secondary regions (27-34 Mb)
# Uses flanking regions: <37 Mb or >47 Mb, but excludes 27-34 Mb region
mean_rho_outside <- recomb_Y %>%
  filter(Loci <= 37000 | Loci >= 47000) %>%    # Outside primary region
  filter(Loci < 27000 | Loci > 34000) %>%      # Exclude secondary region
  summarise(mean_rho = mean(Mean_rho, na.rm = TRUE),
            median_rho = median(Mean_rho, na.rm = TRUE))

# Display results
print("Recombination rates inside region of interest:")
print(mean_rho_inside)

print("Recombination rates outside regions of interest:")
print(mean_rho_outside)

# Calculate fold-change in recombination rate
print("Fold-change (inside/outside):")
print(mean_rho_inside / mean_rho_outside)
