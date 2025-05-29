# ==============================================================================
# GWAS and Pop Gen Analyses of sex-linkage (Figure 1A, Sup Figure 1-3)
# ==============================================================================

# Load required libraries
library(tidyverse)
library(data.table)
library(dplyr)
library(cowplot)


# ==============================================================================
# FIGURE 1A: GWAS Manhattan Plot
# ==============================================================================

#read in data
gwa<-fread("commongarden_allfiltsnps_193_hap2_remvdrepeats_maf05.assoc.txt")

#calculate p value thresholds
gwa$FDR_corr<-p.adjust(gwa$p_wald,method = "fdr")
gwa$p_bon<-p.adjust(gwa$p_wald,method = "bonferroni")
gwa_FDR<-gwa[gwa$FDR_corr < 0.05,]
max(gwa_FDR$p_wald)
gwa_bon<-gwa[gwa$p_bon < 0.05,]
max(gwa_bon$p_wald)


# Create Manhattan plot for GWAS results on Scaffold_1
p1 <- gwa %>% 
  filter(chr == "Scaffold_1") %>%    
  filter(ps > 32500000 & ps < 47000000) %>%  # Focus on 30-47 Mb region
  ggplot(aes(ps / 1000000, -log10(p_wald), color = abs(beta))) +
  geom_point(alpha = .5) +
  theme_bw() +
  xlab("Position in Mb") +
  scale_x_continuous(
    expand = c(0, 0), 
    breaks = seq(floor(min(gwa$ps / 1000000)), ceiling(max(gwa$ps / 1000000)), by = 1)
  ) +
  # Add Bonferroni correction threshold line
  geom_hline(yintercept = -log10(max(gwa_bon$p_wald)), 
             lty = "dashed", lwd = .25, color = "black") +
  # Black to light grey color scale for effect sizes
  scale_color_gradient(high = "black", low = "grey80") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(0.25, 0.5, 0.25, 1.25), "cm")
  )

p1 #FIGURE 1A

# GWAS plot across the genome
#SUP FIGURE 1
p1a <- gwa %>% 
  mutate(chr_numeric = as.numeric(str_extract(chr, "\\d+"))) %>%
  filter(chr_numeric >= 1 & chr_numeric <= 16) %>%
  ggplot(aes(ps / 1000000, -log10(p_wald), color = abs(beta))) +
  geom_point(alpha = .5) +
  theme_bw() +
  xlab("Position in Mb") +
  scale_x_continuous(
    expand = c(0, 0), 
    breaks = seq(floor(min(gwa$ps / 1000000)), ceiling(max(gwa$ps / 1000000)), by = 1)
  ) +
  # Add Bonferroni correction threshold line
  geom_hline(yintercept = -log10(max(gwa_bon$p_wald)), 
             lty = "dashed", lwd = .25, color = "black") +
  # Black to light grey color scale for effect sizes
  scale_color_gradient(high = "black", low = "grey80") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(0.25, 0.5, 0.25, 1.25), "cm")
  )

p1a

###look at GWAS results from re-calling SNPs of sex-linked regions after removing multimapping reads
####SUP FIGURE 2
#read in data
gwa<-fread("~/sexassocsnps_nomultimapping_recalled_metareorder.assoc.txt")

#calculate p value thresholds
gwa$FDR_corr<-p.adjust(gwa$p_wald,method = "fdr")
gwa$p_bon<-p.adjust(gwa$p_wald,method = "bonferroni")
gwa_FDR<-gwa[gwa$FDR_corr < 0.05,]
max(gwa_FDR$p_wald)
gwa_bon<-gwa[gwa$p_bon < 0.05,]
max(gwa_bon$p_wald)

p1.1 <- gwa %>% 
  mutate(chr_numeric = as.numeric(str_extract(chr, "\\d+"))) %>%
  filter(chr_numeric >= 1 & chr_numeric <= 16) %>%
  ggplot(aes(ps / 1000000, -log10(p_wald), color = abs(beta))) +
  geom_point(alpha = .5) +
  theme_bw() +
  xlab("Position in Mb") +
  scale_x_continuous(
    expand = c(0, 0), 
    breaks = seq(floor(min(gwa$ps / 1000000)), ceiling(max(gwa$ps / 1000000)), by = 1)
  ) +
  # Add Bonferroni correction threshold line
  geom_hline(yintercept = -log10(max(gwa_bon$p_wald)), 
             lty = "dashed", lwd = .25, color = "black") +
  # Black to light grey color scale for effect sizes
  scale_color_gradient(high = "black", low = "grey80") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(0.25, 0.5, 0.25, 1.25), "cm")
  ) +
  facet_grid(~chr)

p1.1

###what about GWAS results from competitive mapping of the X and Y?
###SUP FIGURE 3
#read in data
gwa<-fread("~/193_XYScaffolds_filteredsnps.assoc.txt")

#calculate p value thresholds
gwa$FDR_corr<-p.adjust(gwa$p_wald,method = "fdr")
gwa$p_bon<-p.adjust(gwa$p_wald,method = "bonferroni")
gwa_FDR<-gwa[gwa$FDR_corr < 0.05,]
max(gwa_FDR$p_wald)
gwa_bon<-gwa[gwa$p_bon < 0.05,]
max(gwa_bon$p_wald)

p1.2 <- gwa %>% 
  mutate(chr_numeric = as.numeric(str_extract(chr, "\\d+"))) %>%
  filter(chr_numeric >= 1 & chr_numeric <= 16) %>%
  ggplot(aes(ps / 1000000, -log10(p_wald), color = abs(beta))) +
  geom_point(alpha = .5) +
  theme_bw() +
  xlab("Position in Mb") +
  scale_x_continuous(
    expand = c(0, 0), 
    breaks = seq(floor(min(gwa$ps / 1000000)), ceiling(max(gwa$ps / 1000000)), by = 1)
  ) +
  # Add Bonferroni correction threshold line
  geom_hline(yintercept = -log10(max(gwa_bon$p_wald)), 
             lty = "dashed", lwd = .25, color = "black") +
  # Black to light grey color scale for effect sizes
  scale_color_gradient(high = "black", low = "grey80") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(0.25, 0.5, 0.25, 1.25), "cm")
  ) +
  facet_grid(~chr)

p1.2

###SUP FIGURE 8
#what does the genotype-to-phenotype map look like for the top SNP?

tophit<- read.table("~/tophit.012",header = F)
th_ind<-read.table("~/tophit.012.indv",header = F)
tophit$ind<-th_ind$V1
th_meta<- fread("~/commongarden_allfiltsnps_193_hap2.fam")
th_meta$ind<-paste(th_meta$V1,th_meta$V2,sep="_")
sex<-read.table("sexpheno")

tophit_all<-inner_join(tophit,th_meta,by="ind")
tophit_all$sex<-sex$V1

tophit_all$sex2<-NA
tophit_all$sex2[tophit_all$sex=="2"]<-"F"
tophit_all$sex2[tophit_all$sex=="1"]<-"M"

lm1<-lm(data=tophit_all, V2.x ~ sex)
Anova(lm1)
summary(lm1)

tophit_all %>% filter(!is.na(sex2)) %>%
  ggplot(aes( as.factor(V2.x), sex2)) +
  geom_jitter(height = .1, width = .1, alpha = .5) +
  ylab("Sex") +
  xlab("Genotype") +
  theme_bw()

###Sup Figure 11
#Lasso regresssion of top sex-linked SNPs, using the GWAS results where multimapping reads have been removed

genos<-read.table("sexassocsnps_nomultimapping_recalled.raw",header=T) #012 data of sex linked snps
snppos<-read.table("sexassocsnps_nomultimapping_recalled.vcf.bim")
names(genos)<-c(names(genos)[1:6],paste(snppos$V1,snppos$V4,sep = "_"))

library(car)
summary(lm(genos$Scaffold_1_43312648 ~ genos$PHENOTYPE))
summary(lm(genos$Scaffold_1_41377895 ~ genos$PHENOTYPE))
summary(lm(genos$Scaffold_1_33443179 ~ genos$PHENOTYPE))

#install.packages("glmnet")
library(glmnet)

y<-genos$PHENOTYPE
x<-genos[,7:ncol(genos)]
xt<-t(x)
xt<-xt[complete.cases(xt),]
x_comp<-t(xt)
#perform k-fold cross-validation to find optimal lambda value
cv_model <- cv.glmnet(x_comp, y, alpha = 1)

#find optimal lambda value that minimizes test MSE
best_lambda <- cv_model$lambda.min
best_lambda
plot(cv_model)

#find coefficients of best model
best_model <- glmnet(x_comp, y, alpha = 1, lambda = best_lambda)
coef(best_model)

lasso_gwas<-coef(best_model)

# extract non-zero coefficients properly
ind_snps <- lasso_gwas[lasso_gwas[,1] != 0, 1]
ind_snps <- ind_snps[ind_snps != 0]  # Remove any actual zeros

# Get the names of non-zero coefficients
nonzero_names <- names(ind_snps)

# Create the effect_check dataframe
effect_check <- data.frame(
  coeff = as.numeric(ind_snps), 
  val = nonzero_names
)

# First create the chr column
effect_check$chr <- gsub("^(Scaffold_[0-9]+).*", "\\1", effect_check$val)

# Then create the region column
effect_check$region <- case_when(
  effect_check$chr == "Scaffold_1" & 
    as.numeric(gsub(".*_", "", effect_check$val)) >= 40500000 & 
    as.numeric(gsub(".*_", "", effect_check$val)) <= 43500000 ~ "Scaffold_1_target_region",
  effect_check$chr == "Scaffold_1" ~ "Scaffold_1_outside_region",
  TRUE ~ effect_check$chr
)

table(effect_check$region)

# Set up a 1x2 grid (1 row, 2 columns)
par(mfrow = c(1, 2))
plot(cv_model)  # First plot
hist(effect_check$coeff[-1], breaks=20, xlab="Effect Size", main="")  # SUP FIGURE 11
# Reset to single plot layout
par(mfrow = c(1, 1))

rsq <- best_model$dev.ratio
cat("R-squared:", rsq, "\n")

# ======================
# Het and Fst along the Genome
# ======================

fst_all<-fread("male_female_1kb_fst.txt")
head(fst_all)
hist(fst_all[fst_all$chromosome != "Scaffold_1",]$avg_wc_fst)

p3<-fst_all %>% filter(chromosome == "Scaffold_1") %>%  filter(window_pos_1 > 32000000 & window_pos_1 < 47000000) %>%
  ggplot(aes(window_pos_1/1000000,avg_wc_fst)) +
  geom_point(alpha=.2) +
  #facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
  #facet_wrap( ~ chr, ncol = 1, nrow=2) +
  theme_classic() +
  #xlim(33.425,33.46) +
  xlab("Position in Mb") +
  ylab("Weir Cockerham FST") +
  scale_x_continuous(expand=c(0,0)) +
  geom_hline(yintercept = quantile(fst_all[fst_all$chromosome != "Scaffold_1",]$avg_wc_fst, probs=c(0.025)),color="lightgrey",lty="dashed") +
  geom_hline(yintercept = quantile(fst_all[fst_all$chromosome != "Scaffold_1",]$avg_wc_fst, probs=c(0.975)),color="lightgrey",lty="dashed") +
  # geom_ribbon(aes(ymax=quantile(fst_all[fst_all$chromosome != "Scaffold_1",]$avg_wc_fst, probs=c(0.975)),
  #                ymin=quantile(fst_all[fst_all$chromosome != "Scaffold_1",]$avg_wc_fst, probs=c(0.025))),fill="white",alpha=.4) +
  scale_color_viridis_c(option="magma") +
  theme_classic() 

p3

#stats for fst within the SLR
fst_all %>% filter(chromosome == "Scaffold_1") %>%  filter(window_pos_1 > 41500000 & window_pos_1 < 43500000) %>% 
  summarise(mean_Fst=mean(avg_wc_fst), median=median(avg_wc_fst))

######

malefemale_hetprop<-fread("male_female_hetprop.txt")
head(malefemale_hetprop)
names(malefemale_hetprop) <- c("chr","ps","male","male_n","female","female_n")

library(reshape2)
library(tidyverse)
wind_het<-malefemale_hetprop %>% mutate(win=floor(ps/1000)) %>% group_by(chr,win) %>%
  summarise(male_hetprop=mean(male), female_hetprop=mean(female))
head(wind_het)
#hetprop_long<-melt(malefemale_hetprop,id.vars = c("chr","ps"))

p2<-wind_het %>% filter(chr == "Scaffold_1") %>% 
  filter(win > 32000 & win < 47000) %>%
  #mutate(chr_numeric = as.numeric(str_extract(chr, "\\d+"))) %>%
  #filter(chr_numeric >= 1 & chr_numeric <= 16) %>%
  ggplot(aes(win/100,male_hetprop-female_hetprop)) +
  geom_point(alpha=.2) +
  #facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
  # facet_grid( ~ chr_numeric, scales = "free_x") +
  theme_classic() +
  #ylim(0,8) +
  xlab("Position in Mb") +
  ylab("Male:Female\nHeterozygosity") +
  scale_x_continuous(expand=c(0,0)) +
  geom_hline(yintercept=quantile(wind_het$male_hetprop-wind_het$female_hetprop,probs = c(0.025)),color="grey70",lty="dashed") +
  geom_hline(yintercept=quantile(wind_het$male_hetprop-wind_het$female_hetprop,probs = c(0.975)),color="grey70",lty="dashed") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

p2
plot_grid(p2, p3, ncol = 1, nrow=2, align = "hv") #SUP FIGURE 7

#stats for fst within the SLR
wind_het %>% filter(chromosome == "Scaffold_1") %>%  filter(window_pos_1 > 41500000 & window_pos_1 < 43500000) %>% 
  summarise(mean_Fst=mean(avg_wc_fst), median=median(avg_wc_fst))


###SUP FIGURE 5
#lets just take a look at densities across the genome

genes<-read.table("193_2_genes_100kbdensity.bed")
#TEs<-fread("allrepeats_densities_10kb.bed")
TEs<-read.table("193_2_TEs_100kbdensity.bed")

names(genes)<-c("chr","start","end","count","density")
names(TEs)<-c("chr","start","end","count","density")
div<-read.table("0fold_allinds_allsites_100kb_pi.txt",header = T)

head(div)
head(genes)


pTE<-TEs %>% filter(chr == "Scaffold_1") %>%
  mutate(chr_numeric = as.numeric(str_extract(chr, "\\d+"))) %>%
  #filter(chr_numeric >= 1 & chr_numeric <= 16) %>%
  ggplot(aes(end/1000000,density)) +
  geom_point(alpha=.15) +
  #facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
  facet_grid( ~ chr_numeric, scales = "free_x") +
  theme_classic() +
  #coord_cartesian(ylim=c(0,0.3)) +
  xlab("100Kb Window") +
  ylab("TE Density") +
  geom_vline(xintercept = 28.5,lwd=4, alpha=.2) +
  geom_vline(xintercept = 42, lwd=5, alpha=.2) +
  geom_smooth(span=.1, method="loess",alpha=.3) +
  # xlim(37,47)
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 

TEs %>% filter(chr == "Scaffold_1") %>%
  mutate(chr_numeric = as.numeric(str_extract(chr, "\\d+"))) %>%
  filter(chr_numeric >= 1 & chr_numeric <= 16) %>%
  ggplot(aes(end/1000000,density)) +
  geom_col() +
  #facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
  facet_grid( ~ chr_numeric, scales = "free_x") +
  theme_classic() +
  #coord_cartesian(ylim=c(0,0.3)) +
  xlab("100Kb Window") +
  ylab("TE Density") +
  geom_vline(xintercept = 42.325) 

pTE

pGene <- genes %>% filter(chr == "Scaffold_1") %>% #filter(start>40000000 & end < 45000000) %>%
  mutate(chr_numeric = as.numeric(str_extract(chr, "\\d+"))) %>%
  filter(chr_numeric >= 1 & chr_numeric <= 16) %>%
  ggplot(aes(end/100000,density)) +
  geom_point(alpha=.15) +
  #facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
  facet_grid( ~ chr_numeric, scales = "free_x") +
  theme_classic() +
  #ylim(0,8) +
  xlab("100Kb Window") +
  ylab("Gene Density") +
  geom_vline(xintercept = 285,lwd=4, alpha=.2) +
  geom_vline(xintercept = 420, lwd=5, alpha=.2) +
  geom_smooth(span=.1) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())



pDIV<-div %>% filter(chromosome == "Scaffold_1") %>% #filter(start>40000000 & end < 45000000) %>%
  mutate(chr_numeric = as.numeric(str_extract(chromosome, "\\d+"))) %>%
  filter(chr_numeric >= 1 & chr_numeric <= 16) %>%
  ggplot(aes((((window_pos_2-50000)/100000)/10),avg_pi)) +
  geom_point(alpha=.15) +
  #facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
  facet_grid( ~ chr_numeric, scales = "free_x") +
  theme_classic() +
  coord_cartesian(ylim=c(0,0.06)) +
  xlab("Genomic Position (Mb)") +
  ylab("0 Fold Diversity") +
  geom_vline(xintercept = 28.5,lwd=4, alpha=.2) +
  geom_vline(xintercept = 42.0, lwd=5, alpha=.2) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_smooth(span=.1) 

#pDIV

library(cowplot)
plot_grid(pGene,pTE,pDIV,nrow = 3,align = "hv") #Sup FIGURE 5

########


