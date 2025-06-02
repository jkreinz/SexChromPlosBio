# ==============================================================================
# Depth Based-Analyses - Figure 1B, Figure 3C, Sup Fig XXX, 
# ==============================================================================
# Load required libraries
library(tidyverse)
library(data.table)
library(dplyr)
library(cowplot)
# ==============================================================================
# FIGURE 1B: Data Preparation and Sample Metadata
# ==============================================================================

# Read sample list and PLINK files
all <- read.table("~/herbcontemp_poplist.txt", 
                  na.strings = c("", "NA", "?"), sep = "\t", header = F)
names(all)[1] <- "sample"

plink <- read.table("~/herb_commongarden_pnas_missing9_maf05.fam")
names(plink)[1] <- c("sample")
plink$order <- seq(1, nrow(plink))

# Join sample metadata
both <- inner_join(plink, all, by = "sample")

# Read additional sample information
all2 <- read.table("~/3waymerged_sampleinfo.txt", sep = "\t", header = T)
both2 <- inner_join(plink, all2, by = "sample")

# Read PCA results from common garden analysis
contemp_herb <- fread("commongarden_allfiltsnps_193_hap2.eigenvec")
names(contemp_herb)[2] <- "sample"
contemp_herb$sample <- as.character(contemp_herb$sample)

inner_join(contemp_herb, all2, by = "sample")

# Assign dataset categories based on collection year
all2$dataset <- "PNAS"
all2$dataset[all2$year == 2019] <- "Paired"
all2$dataset[all2$year < 2015] <- "Herbarium"

# Combine PCA and metadata
contemp_herb_pca <- inner_join(contemp_herb, all2, by = "sample")

# Read ordered metadata with ancestry information
meta_ordered <- read.table(
  "~/Documents/UToronto_PhD/Papers/Evolution Preadapt Submission/Evolution2021-Preadapt/commongarden_and_PNAS_metadata_ancestry.txt",
  header = T
)
names(meta_ordered)[1] <- "sample"

# Add ancestry information to PCA data
contemp_herb_pca_anc <- inner_join(contemp_herb_pca, meta_ordered, by = "sample")

# Read misphenotyped samples
mispheno <- read.table("misphenoed_samples.txt")
names(mispheno) <- c("order", "population", "sample")
mispheno$sample <- as.character(mispheno$sample)

# ==============================================================================
# DEPTH ANALYSIS: Sex chromosome hemizygosity detection
# ==============================================================================

# Read depth data for 100bp windows
depthmat <- fread("~/allsamples_numericorder_depthsby100bpwind_tab.txt")
names(depthmat) <- c("scaf", "start", "end", str_sort(contemp_herb_pca_anc$sample, numeric = F))

# Calculate individual depth statistics for hemizygous region
# Hemizygous region: 40.5-45 Mb (putative sex chromosome region)
# Background region: outside 40.5-45 Mb
ind_means_hemiz <- data.frame(
  sample = names(depthmat[depthmat$start > 40500000 & depthmat$start < 45000000, -c(1:3)]),
  # Median depth in hemizygous region (40.5-43.6 Mb)
  ind_meandepth = miscTools::colMedians(depthmat[depthmat$start > 40500000 & depthmat$start < 43600000, -c(1:3)]),
  # Median depth in background scaffold regions
  scaf1depth_ = miscTools::colMedians(depthmat[depthmat$start < 40500000 | depthmat$start > 45000000, -c(1:3)])
)

# Reshape depth data for analysis region (32.5-47 Mb)
long_alldepth <- melt(depthmat[start > 32500000 & start < 47000000, ], 
                      id.vars = c("scaf", "start", "end"))
names(long_alldepth)[4] <- "sample"

# ==============================================================================
# SEX DETERMINATION AND GROUPING
# ==============================================================================


# Combine depth statistics with PCA and ancestry data
contemp_herb_pca_anc_hemidepth <- inner_join(ind_means_hemiz, contemp_herb_pca_anc, by = "sample")
contemp_herb_pca_anc_hemidepth_misgroup <- left_join(contemp_herb_pca_anc_hemidepth, mispheno, by = "sample")

# Classify samples based on sex and depth patterns
contemp_herb_pca_anc_hemidepth_misgroup$grouping <- NA

# Good matches: expected depth patterns for sex
contemp_herb_pca_anc_hemidepth_misgroup$grouping[
  is.na(contemp_herb_pca_anc_hemidepth_misgroup$population) & 
    contemp_herb_pca_anc_hemidepth_misgroup$sex == "M"
] <- "goodmale"

contemp_herb_pca_anc_hemidepth_misgroup$grouping[
  is.na(contemp_herb_pca_anc_hemidepth_misgroup$population) & 
    contemp_herb_pca_anc_hemidepth_misgroup$sex == "F"
] <- "goodfemale"

# Mismatches: unexpected depth patterns for assigned sex
# Females with high depth in hemizygous region (>= 7.5x)
contemp_herb_pca_anc_hemidepth_misgroup$grouping[
  contemp_herb_pca_anc_hemidepth_misgroup$sex == "F" & 
    contemp_herb_pca_anc_hemidepth_misgroup$ind_meandepth >= 7.5
] <- "wierdfemale"

# Males with low depth in hemizygous region (<= 7.5x)
contemp_herb_pca_anc_hemidepth_misgroup$grouping[
  contemp_herb_pca_anc_hemidepth_misgroup$sex == "M" & 
    contemp_herb_pca_anc_hemidepth_misgroup$ind_meandepth <= 7.5
] <- "wierdmale"

# Convert to factor with specific level order
contemp_herb_pca_anc_hemidepth_misgroup$grouping <- as.factor(contemp_herb_pca_anc_hemidepth_misgroup$grouping)
contemp_herb_pca_anc_hemidepth_misgroup$grouping <- factor(
  contemp_herb_pca_anc_hemidepth_misgroup$grouping, 
  levels = c("goodfemale", "wierdfemale", "wierdmale", "goodmale")
)


# ==============================================================================
# SEX MISMATCH STATISTICAL TEST
# ==============================================================================

# Test for excess of mismatched males vs females
table(contemp_herb_pca_anc_hemidepth_misgroup$grouping)

# Define counts for two-proportion test
female_mismatch <- 5   # Number of mismatched females
female_total <- 94     # Total females
male_mismatch <- 17    # Number of mismatched males  
male_total <- 89       # Total males

# Two-proportion test comparing mismatch rates between sexes
prop.test(c(female_mismatch, male_mismatch), c(female_total, male_total), correct = FALSE)


# ==============================================================================
# ANCESTRY GROUPING
# ==============================================================================

# Classify samples by ancestry based on admixture proportions (K1)
contemp_herb_pca_anc_hemidepth_misgroup$ancestry_grouping <- NA

# Predominantly var. tuberculatus (K1 <= 0.25)
contemp_herb_pca_anc_hemidepth_misgroup[
  contemp_herb_pca_anc_hemidepth_misgroup$K1 <= 0.25,
]$ancestry_grouping <- "predom. var. tuberculatus"

# Predominantly var. rudis (K1 >= 0.75)
contemp_herb_pca_anc_hemidepth_misgroup[
  contemp_herb_pca_anc_hemidepth_misgroup$K1 >= 0.75,
]$ancestry_grouping <- "predom. var. rudis"

# Admixed (0.25 < K1 < 0.75)
contemp_herb_pca_anc_hemidepth_misgroup[
  contemp_herb_pca_anc_hemidepth_misgroup$K1 < 0.75 & 
    contemp_herb_pca_anc_hemidepth_misgroup$K1 > 0.25,
]$ancestry_grouping <- "admixed"

# ==============================================================================
# PREPARE DATA FOR DEPTH VISUALIZATION
# ==============================================================================

# Combine depth data with metadata
long_alldepth_meta <- inner_join(long_alldepth, contemp_herb_pca_anc_hemidepth_misgroup, by = "sample")
long_alldepth_meta_newhemi <- inner_join(
  long_alldepth_meta[, c(1:5, 35)], 
  ind_means_hemiz, 
  by = "sample"
)

# Summarize depth data by 200kb windows
longall_summary <- long_alldepth_meta_newhemi %>%
  mutate(win = floor(end / 200000)) %>%  # Create 200kb windows
  group_by(scaf, win, sample) %>%
  summarise(
    mean_value = median(value),           # Median depth per window
    scafdepth = first(scaf1depth_),      # Background scaffold depth
    grouping = first(grouping),          # Sex by depth grouping
    .groups = 'drop'
  )

# ==============================================================================
# FIGURE 1C: Scaled Depth Plot Across Sex Chromosome Region
# ==============================================================================

p11 <- longall_summary %>%
  ggplot(aes(win / 10, mean_value / scafdepth, group = sample, color = grouping)) +
  geom_line(alpha = .5) +
  coord_cartesian(ylim = c(0, 2.2)) +
  labs(y = "Scaled Depth", x = "Genomic Position (Mb)") +
  theme_bw() +
  geom_hline(yintercept = 1) +  # Expected depth ratio for diploid regions
  guides(colour = guide_legend(override.aes = list(lwd = 2, alpha = .8))) +
  scale_color_viridis_d(
    option = "magma", 
    end = .92,
    labels = c('Female', 'Mismatched Female', 'Mismatched Male', 'Male', "")
  ) +
  scale_x_continuous(
    expand = c(0, 0), 
    breaks = seq(floor(min(longall_summary$win / 10)), ceiling(max(longall_summary$win / 10)), by = 1)
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  guides(color = "none")

p11

plot_grid(p1,p11,nrow=,ncol=1,align = "hv",axis = "lr",
          rel_heights = c(2,2),labels = "AUTO")


# ==============================================================================
# FIGURE 3C: ANCESTRY vs SCALED DEPTH BY SEX GROUPING
# ==============================================================================

contemp_herb_pca_anc_hemidepth_misgroup$scaleddepth<-contemp_herb_pca_anc_hemidepth_misgroup$ind_meandepth/contemp_herb_pca_anc_hemidepth_misgroup$scaf1depth_

q2 <- contemp_herb_pca_anc_hemidepth_misgroup %>% 
  filter(!is.na(grouping)) %>%
  ggplot(aes(K1, scaleddepth, color = grouping)) +
  geom_point() +
  geom_smooth(method = "lm", alpha = .3) +
  theme_bw() +
  labs(
    x = "Proportion var. rudis Ancestry",
    y = "Scaled Depth in Male Hemizygous Region"
  ) +
  scale_color_viridis_d(option = "magma", end = .92)

q2

# ==============================================================================
# MULTIVARIATE ANALYSIS: FACTORS AFFECTING SCALED DEPTH
# ==============================================================================

# Comprehensive model including ancestry, environment, and geography
contemp_depth <- lm(
  data = contemp_herb_pca_anc_hemidepth_misgroup, 
  scaleddepth ~ grouping * K1 + env + lat + long
) 

# Type II ANOVA for main effects and interactions
Anova(contemp_depth, type = 2)  # Main text result
Anova(contemp_depth, type = 3)  # Need tp update!!! Main text result


# ==============================================================================
# SEX CHROMOSOME COPY NUMBER DEVIATION ANALYSIS
# ==============================================================================


# Calculate mean scaled depth by sex grouping (for reference values)
means_scaled <- contemp_herb_pca_anc_hemidepth_misgroup %>% 
  group_by(grouping) %>% 
  summarise(meangroup = mean(scaleddepth), .groups = 'drop')

# Extract reference means for good matches
mean_goodfemale <- means_scaled$meangroup[means_scaled$grouping == 'goodfemale'][1]
mean_goodmale <- means_scaled$meangroup[means_scaled$grouping == 'goodmale'][1]

# Error checking for required values
if(length(mean_goodfemale) == 0 || length(mean_goodmale) == 0) {
  stop("Required mean values not found in means_scaled.")
}

# Assign expected depth values based on sex grouping
contemp_herb_pca_anc_hemidepth_misgroup$sex2 <- NA

# Assign expected female depth to all female groups
contemp_herb_pca_anc_hemidepth_misgroup$sex2[
  grepl("female", contemp_herb_pca_anc_hemidepth_misgroup$grouping, ignore.case = TRUE)
] <- mean_goodfemale

# Assign expected male depth to all male groups
contemp_herb_pca_anc_hemidepth_misgroup$sex2[
  !grepl("female", contemp_herb_pca_anc_hemidepth_misgroup$grouping, ignore.case = TRUE)
] <- mean_goodmale

# Calculate deviation from expected copy number
contemp_herb_pca_anc_hemidepth_misgroup$distancefromopt <- 
  contemp_herb_pca_anc_hemidepth_misgroup$scaleddepth - 
  contemp_herb_pca_anc_hemidepth_misgroup$sex2

# Visualize distribution of deviations
hist(abs(contemp_herb_pca_anc_hemidepth_misgroup$distancefromopt), 
     breaks = 50, main = "Distribution of Copy Number Deviations")

# Linear model testing geographic and environmental effects on copy number deviation
lm2 <- lm(data = contemp_herb_pca_anc_hemidepth_misgroup, 
          distancefromopt ~ lat + long + state + env + K1)

# Analysis of variance for model terms
Anova(lm2)
summary(lm2)


# ==============================================================================
# GEOGRAPHIC MAPPING of SEX by DEPTH in hemizygous region
# ==============================================================================

# Load required libraries for mapping and visualization
library(dplyr)
library(ggmap)
library(scatterpie)
library(wesanderson)
library(PNWColors)
library(MetBrewer)
library(maps)
library(mapdata)
library(maptools)
library(scales)
library(raster)
library(sf)
library(tidyverse)
library(ggrepel)
library(rnaturalearth)
library(rnaturalearthhires)
library(ggthemes)
library(ggnewscale)
library(viridis)
# Load ocean color palette
#devtools::install_github("kwstat/pals")   
library(pals)

# ==============================================================================
# DATA PREPARATION AND SAMPLE METADATA
# ==============================================================================

# Read sample information with geographic coordinates
cg <- read.table("~/3waymerged_sampleinfo.txt", 
                 sep = "\t", na.strings = c("", "NA"), header = T)

# Assign dataset categories based on collection year
cg$dataset <- "PNAS"
cg$dataset[cg$year == 2019] <- "Paired"
cg$dataset[cg$year < 2015] <- "Herbarium"

# Get population-level coordinates and sample counts by environment
coord_byenv <- cg %>% 
  group_by(env, long, lat, dataset) %>% 
  summarise(n = n(), state = first(state), .groups = 'drop')

# ==============================================================================
# GEOGRAPHIC MAP BASE PREPARATION
# ==============================================================================

# Define states and provinces to include in map
states <- c("New York", "Pennsylvania", "Maryland", "West Virginia", "Virginia", 
            "Kentucky", "Ohio", "Michigan", "Indiana", "Illinois", "Wisconsin", 
            "Minnesota", "Iowa", "Missouri", "Kansas", "Nebraska", "South Dakota", 
            "North Carolina", "Tennessee", "Mississippi", "Oklahoma",
            "Lake Michigan", "Lake Ontario", "Lake Superior")
provinces <- c("Ontario")

# Get administrative boundaries
us <- getData("GADM", country = "USA", level = 1)
canada <- getData("GADM", country = "CAN", level = 1)
us.states <- us[us$NAME_1 %in% states, ]
ca.provinces <- canada[canada$NAME_1 %in% provinces, ]

# Get natural features (lakes, ocean, rivers)
lakes <- rnaturalearth::ne_download(scale = 10, type = 'lakes', category = 'physical') %>%
  sf::st_as_sf(crs = 4269)

ocean <- rnaturalearth::ne_download(scale = 10, type = 'ocean', category = 'physical') %>%
  sf::st_as_sf(crs = 4269)

rivers <- rnaturalearth::ne_download(scale = 10, type = 'rivers_lake_centerlines', 
                                     category = 'physical') %>%
  sf::st_as_sf(crs = 4269)

# Prepare state labels with adjusted positions
in_sf <- ne_states(geounit = "United States of America", returnclass = "sf")
uss <- st_as_sf(us.states) %>%
  mutate(
    lon = map_dbl(geometry, ~st_centroid(.x)[[1]]),
    lat = map_dbl(geometry, ~st_centroid(.x)[[2]])
  )

# Adjust Michigan label position and add Ontario
uss[uss$NAME_1 == "Michigan", ]$lon <- -85.0554
uss[uss$NAME_1 == "Michigan", ]$lat <- 44.00902
uss <- uss %>% add_row(NAME_1 = "Ontario", lon = -78.5554, lat = 45)

# ==============================================================================
# ENVIRONMENT VARIABLE STANDARDIZATION
# ==============================================================================

# Rename environmental variables for better legends
coord_byenv$env[coord_byenv$env == "Ag"] <- "Agricultural"
coord_byenv$env[coord_byenv$env == "Nat"] <- "Natural"

cg$env[cg$env == "Ag"] <- "Agricultural"
cg$env[cg$env == "Nat"] <- "Natural"
cg$env[cg$env == "Dist"] <- "Disturbed"

# Convert to factor with specific level order
coord_byenv$env <- as.factor(coord_byenv$env)
cg$env <- as.factor(cg$env)
cg$env <- factor(cg$env, levels = c("Agricultural", "Natural", "Disturbed"))

# ==============================================================================
# CREATE BASE MAP
# ==============================================================================

plain1 <- ggplot() +
  # Add state and province boundaries
  geom_path(data = us.states, aes(x = long, y = lat, group = group)) +
  geom_path(data = ca.provinces, aes(x = long, y = lat, group = group), fill = "white") +
  coord_map() +
  # Add natural features
  geom_sf(data = lakes, mapping = aes(geometry = geometry), color = "black") +
  geom_sf(data = ocean, mapping = aes(geometry = geometry), color = "black") +
  geom_sf(data = rivers, mapping = aes(geometry = geometry), 
          color = "grey50", alpha = .75) +
  # Set map extent to focus on study region
  scale_x_continuous(limits = c(-97, -82)) +
  scale_y_continuous(limits = c(38, 43)) +
  coord_sf(xlim = c(-96, -82)) +
  # Add state/province labels
  geom_text(data = uss, aes(x = lon, y = lat, label = NAME_1), 
            cex = 3, alpha = .8) +
  theme_map()


# ==============================================================================
# GEOGRAPHIC VISUALIZATION: COPY NUMBER DEVIATIONS
# ==============================================================================

# Create map showing copy number deviations across sampling locations
map_copynum <- plain1 +
  geom_jitter(
    data = contemp_herb_pca_anc_hemidepth_misgroup,
    aes(y = lat, x = long, color = distancefromopt, shape = env),
    fill = "black", alpha = 1, height = .25, width = .2, cex = 2.5
  ) +
  scale_shape_manual(values = c("circle", "triangle")) +
  scale_color_gradientn(colours = ocean.balance(100), guide = "colourbar") +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(
    colour = "Difference in Copy Number\nfrom Matched-Sex",
    shape = "Collection\nHabitat"
  ) +
  theme(legend.position = "top")

map_copynum


# ========================================================
#Investigation of Life History Phenotypes by sex groupings
# ========================================================

phenos<-fread("~/Downloads/Common Garden Seq Stats - geno w pheno.csv")

head(contemp_herb_pca_anc_hemidepth)
head(phenos)

names(phenos)[1]<-"sample"
phenos$sample<-as.character(phenos$sample)
groups_wphenos<-inner_join(contemp_herb_pca_anc_hemidepth_misgroup,phenos,by="sample")

groups_wphenos$grouping<-factor(groups_wphenos$grouping,levels=c("goodmale","wierdmale","goodfemale","wierdfemale"))
groups_wphenos$matchstatus <- ifelse(grepl("wierd", groups_wphenos$grouping), "mismatched", "matched")

lm1<-lm(data=groups_wphenos, Dry.Biomass ~ sex + matchstatus + lat + long + K1 + env + Treatment )
lm1_sex<-lm(data=groups_wphenos, Dry.Biomass ~ sex  + lat + long + K1 + env + Treatment )
AIC(lm1,lm1_sex)

lm1<-lm(data=groups_wphenos, DaystoFlower ~ sex + matchstatus + lat + long + K1 + env + Treatment )
lm1_sex<-lm(data=groups_wphenos, DaystoFlower ~ sex  + lat + long + K1 + env + Treatment )
AIC(lm1,lm1_sex)


#for plotting only
lm1<-lm(data=groups_wphenos, Dry.Biomass ~ grouping + lat + long + K1 + env + Treatment )
Anova(lm2,type=2)
library(lsmeans)
out1<-as.data.frame(lsmeans(lm2, ~grouping))

p1<-ggplot() + 
  geom_point(data=out2, aes( grouping,lsmean)) + 
  geom_errorbar(data=out2,aes(ymax=upper.CL, ymin=lower.CL, x=grouping)) +
  labs(y="Days to Flowering") +
  geom_jitter(data=groups_wphenos, aes( grouping,DaystoFlower,color=Treatment)) +
  theme_bw()


lm2<-lm(data=groups_wphenos, DaystoFlower ~ grouping + lat + long + K1 + env + Treatment )
Anova(lm2,type=2)
library(lsmeans)
out2<-as.data.frame(lsmeans(lm2, ~grouping))

p2<-ggplot() + 
  geom_point(data=out2, aes( grouping,lsmean)) + 
  geom_errorbar(data=out2,aes(ymax=upper.CL, ymin=lower.CL, x=grouping)) +
  labs(y="Days to Flowering") +
  geom_jitter(data=groups_wphenos, aes( grouping,DaystoFlower,color=Treatment)) +
  theme_bw()

#SUP FIGURE 21
plot_grid(p1,p2,ncol=1,nrow=2)


# ========================================================
# Depth GWAS - Figure 4 analyses and plotting
# ========================================================

###below files read and in processed earlier in this script
#depthmat<-fread("allsamples_numericorder_depthsby100bpwind_tab.txt",sep="\t",sep2 = " ")
#depthmat[1:10,1:10]
#library(stringr)
#names(depthmat)<-c("scaf","start","end",str_sort(contemp_herb_pca_anc$sample, numeric = F))

#GWAS
samp<-data.frame(sample=names(depthmat)[-c(1:3)])
samp_meta<-inner_join(samp,contemp_herb_pca_anc,by="sample")
means_depth<-colMeans(depthmat[depthmat$start < 40500000 | depthmat$start > 47000000,-c(1:3)],na.rm = T)

#run sex GWAS based on gene copy number
pval<-list()
tval<-list()
slope<-list()
slope_error<-list()
temp<-data.frame("geno"=NA, "pheno"="NA")
b0<-list()

#depthmat_MADS represents a matrix of depth at a given locus across samples
#by scaling by mean depth outside of our focal region, we get an estimate of Copy number (Scaled/Normalized Depth)
depthmat_MADS<-depthmat[depthmat$start > 40000000 & depthmat$start < 45000000,]

for (i in 1:(nrow(depthmat_MADS))) {
  temp<-data.frame(geno=as.numeric(depthmat_MADS[i,-c(1:3)]/means_depth), 
                   pheno=samp_meta$sex)
  
  temp[is.na(temp) | temp=="Inf"] = NA
  temp_comp<-temp[complete.cases(temp),]
  
  lm1<-lm(data=temp_comp, geno ~ pheno)
  temp2<-as.data.frame(summary(lm1)$coefficients)
  
  pval[[i]]<-temp2$`Pr(>|t|)`[2]
  tval[[i]]<-temp2$`t value`[2]
  slope[[i]]<-temp2$Estimate[2]
  slope_error[[i]]<-temp2$`Std. Error`[2]
  b0[[i]]<-temp2$Estimate[1]
  
}

#merge results into a data frame
results<-data.frame(scaf=depthmat_MADS$scaf, start=depthmat_MADS$start, end=depthmat_MADS$end, pval=unlist(pval),tval=unlist(tval),slope=unlist(slope),slope_error=unlist(slope_error))
hist(-log10(results$pval))
results$scaf <- gsub("Scaffold_","",results$scaf)

#perform multiple test correction
results$fdr<-p.adjust(results$pval, method = "fdr") 
results_fdrsig<-results[results$fdr < 0.05,]

results$bon<-p.adjust(results$pval, method = "bonferroni") 
results_bonsig<-results[results$bon < 0.05,]

#plot
#FIGURE 4
results$scaf<-as.numeric(results$scaf)
z1<-ggplot(data=results, aes(start/1000000, -log10(pval))) +
  geom_point(alpha=.4) +
  #acet_grid(~scaf, scales = "free_x",space = "free_x") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  xlab("Genomic Position (Mb)") +
  geom_hline(yintercept = -log10(max(na.omit(results_fdrsig$pval))),lty="dashed") +
  geom_hline(yintercept = -log10(max(na.omit(results_bonsig$pval))),lty="dashed") +
  geom_vline(xintercept = 42.3208) +
  geom_vline(xintercept = 42.3246) +
  theme_bw()

z1 #figure 4A

contemp_herb_pca_anc_hemidepth_misgroup

wierdmale_depth <- subset(depthmat, select = na.omit(contemp_herb_pca_anc_hemidepth_misgroup[contemp_herb_pca_anc_hemidepth_misgroup$grouping == "wierdmale",]$sample))
wierdfemale_depth <- subset(depthmat, select = na.omit(contemp_herb_pca_anc_hemidepth_misgroup[contemp_herb_pca_anc_hemidepth_misgroup$grouping == "wierdfemale",]$sample))
goodfemales_depth <- subset(depthmat, select = na.omit(contemp_herb_pca_anc_hemidepth_misgroup[contemp_herb_pca_anc_hemidepth_misgroup$grouping == "goodfemale",]$sample))
goodmales_depth <- subset(depthmat, select = na.omit(contemp_herb_pca_anc_hemidepth_misgroup[contemp_herb_pca_anc_hemidepth_misgroup$grouping == "goodmale",]$sample))

library(matrixStats)
depth_summary<-data.frame(scaf=depthmat$scaf,
                          start=depthmat$start,
                          end=depthmat$end,
                          wierdmale=rowMeans(as.matrix(wierdmale_depth)),
                          wierdfemale=rowMeans(as.matrix(wierdfemale_depth)),
                          goodmale=rowMeans(as.matrix(goodmales_depth)),
                          goodfemale=rowMeans(as.matrix(goodfemales_depth)),
                          wierdmale_var=rowVars(as.matrix(wierdmale_depth)),
                          wierdfemale_var=rowVars(as.matrix(wierdfemale_depth)),
                          goodmale_var=rowVars(as.matrix(goodmales_depth)),
                          goodfemale_var=rowVars(as.matrix(goodfemales_depth)))


library(reshape2) 
long_depthsum<-melt(depth_summary, id.vars=c("scaf","start","end"))

medians_depth<-long_depthsum[long_depthsum$start < 40500000 | long_depthsum$start > 45000000,] %>%
  filter(!grepl("_var", variable)) %>%
  group_by(variable) %>%
  dplyr:::summarise(mean_value=mean(value))

#figure 4B upper
long_depthsum_scaled<-inner_join(long_depthsum,medians_depth,by="variable")
long_depthsum_scaled$scaled<-long_depthsum_scaled$value/long_depthsum_scaled$mean_value

fig4bupper<-long_depthsum_scaled %>% filter(start>42290000 & start < 42340000) %>%
  filter(variable == "goodmale" | variable == "goodfemale") %>%
  # filter(win>4229 & win < 4234) %>%
  ggplot( aes(start/1000000,scaled,color=variable)) +
  geom_point(alpha=.3) +
  geom_smooth(method="loess",se=T,span=.2) +
  ylim(0,10) +
  #coord_cartesian(ylim=c(0,1.7)) +
  labs(y="Scaled Depth (Copy Number)", x="Genomic Position (Mb)") +
  theme_bw() + #guides(color="none") +
  geom_hline(yintercept = 1, lty="dashed") +
  geom_vline(xintercept = 42.3208) +
  geom_vline(xintercept = 42.3246) +
  scale_color_viridis_d(option="magma",end = .92,direction=-1) +
  guides(color = "none")


fig4blower<-results %>% filter(start>42290000 & start < 42340000) %>%
  ggplot(aes(start/1000000, -log10(pval))) +
  geom_point(alpha=.4) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  xlab("Genomic Position (Mb)") +
  geom_hline(yintercept = -log10(max(na.omit(results_fdrsig$pval))),lty="dashed") +
  geom_hline(yintercept = -log10(max(na.omit(results_bonsig$pval))),lty="dashed") +
  geom_vline(xintercept = 42.3208) +
  geom_vline(xintercept = 42.3246) +
  theme_bw()

plot_grid(fig4bupper, fig4blower, ncol = 1, nrow = 2, align = "hv" )

#########
tophit<-depthmat[depthmat$start ==42327700,]
tophip_df<-data.frame(depth=as.numeric(t(tophit)[-c(1:3)]), sex=samp_meta$sex, sample=samp_meta$sample)
tophip_df2<-inner_join(tophip_df,contemp_herb_pca_anc,by="sample")
means_depthdf<-data.frame(sample=names(means_depth), mean_depthscaf1=means_depth)
tophip_df3<-inner_join(tophip_df2,means_depthdf,by="sample")
tophip_df3$scaled_depth<-tophip_df3$depth/tophip_df3$mean_depthscaf1

summary(lm(tophip_df3$scaled_depth ~ tophip_df3$sex.x))

#FIGURE 4C
tophip_df3 %>%
  filter(!is.na(sex.x)) %>%
  ggplot(aes( sex.x, scaled_depth, color=sex.x)) +
  geom_jitter() +
  theme_classic() +
  scale_color_viridis_d(option="magma",end = .92) +
  geom_boxplot(outlier.shape = NA)


#########
#Sup Figure 12
library(reshape2)
library(stringr)


depthmat<-fread("allsamples_numericorder_sexhitsdepth_tab.txt",sep="\t",sep2 = " ")
depthmat[1:10,1:10]
names(depthmat)<-c("scaf","start","end",str_sort(contemp_herb_pca_anc$sample, numeric = F))

longdepthsex<-melt(depthmat,id.vars = c("scaf","start","end"))
head(longdepthsex)
names(longdepthsex)[4]<-"sample"
longdepthsex_meta<-inner_join(longdepthsex, contemp_herb_pca_anc[,c(2,14)],by="sample")
longdepthsex_meta_scaf1depth<-inner_join(longdepthsex_meta, ind_means_hemiz[,c(1,3)],by="sample")
longdepthsex_meta_scaf1depth$scaleddepth<-longdepthsex_meta_scaf1depth$value/longdepthsex_meta_scaf1depth$scaf1depth_

depth_diff<-longdepthsex_meta_scaf1depth %>% filter(sex != "NA") %>%
  group_by(scaf,start, sex) %>%
  summarise(mean_value=mean(scaleddepth)) %>% 
  spread(sex, mean_value) %>%
  mutate(difference = M / F)

depth_bysex<-longdepthsex_meta_scaf1depth %>% filter(sex != "NA") %>%
  group_by(scaf,end, sex) %>%
  summarise(mean_value=mean(scaleddepth)) %>%
  spread(sex, mean_value)

head(depth_bysex)
depth_bysex$targetscaf<-NA
depth_bysex[depth_bysex$scaf == "Scaffold_1",]$targetscaf<-"TRUE"

z1<-ggplot(data=depth_bysex, aes(M,F,color=targetscaf)) +
  geom_point(alpha=.3) + 
  geom_abline(slope = 1,intercept = 0) +
  # lims(y=c(0,120), x=c(0,170)) +
  theme_bw() +
  guides(color="none") +
  labs(y="Female Coverage",x="Male Coverage")

z2<-ggplot(data=depth_bysex, aes(M/F,fill=targetscaf)) +
  geom_histogram() +
  # lims(y=c(0,120), x=c(0,170)) +
  theme_bw() +
  labs(x="Male:Female Coverage")

library(cowplot)
plot_grid(z1,z2, rel_widths =c(1,1.3)) #SUP FIGURE 12


# ==============================================================================
# Plot median depth in hemizygous region in herbarium dataset: Sup Figure 20
# =============================================================================


# Read depth data for 100bp windows
depthmat <- fread("~/allsamples_numericorder_depthsby100bpwind_herbarium_tab.txt")
samp<-read.table("herborder.txt")
names(depthmat) <- c("scaf", "start", "end", samp$V1)

# Calculate individual depth statistics for hemizygous region
# Hemizygous region: 40.5-45 Mb (putative sex chromosome region)
# Background region: outside 40.5-45 Mb
ind_means_hemiz <- data.frame(
  sample = names(depthmat[depthmat$start > 40500000 & depthmat$start < 45000000, -c(1:3)]),
  # Median depth in hemizygous region (40.5-43.6 Mb)
  ind_meandepth = miscTools::colMedians(depthmat[depthmat$start > 40500000 & depthmat$start < 43600000, -c(1:3)]),
  # Median depth in background scaffold regions
  scaf1depth_ = miscTools::colMedians(depthmat[depthmat$start < 40500000 | depthmat$start > 45000000, -c(1:3)])
)

ind_means_hemiz$sample<-gsub("HBO","HB0",ind_means_hemiz$sample)

herb_depth<-inner_join(all2[all2$year < 2015,],ind_means_hemiz,by="sample" )

#SUP FIGURE 20
write.csv(herb_depth, file="figure_data_export/SupFig20_dataset.csv")

herb_depth %>% filter(sex == "F" | sex == "M") %>%
  ggplot(aes(long,ind_meandepth/scaf1depth_, color=sex)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(y="Median Scaled Depth in Male Hemizygous Region", x="Longitude") +
  theme_bw()


