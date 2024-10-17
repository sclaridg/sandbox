## ---------------------------
##
## Script name:
##
## Purpose of script:
##
## Author: Sally Claridge
##
## Date Created: 26 October 2022
##
## Copyright (c) Sally Claridge, 2022
## Email: sally.claridge@icahn.mssm.edu
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

options(scipen = 6, digits = 4)

## ---------------------------

## load up the packages we will need:

pkgs <- c("tidyverse", "reshape2", "magrittr", "ggpubr")

check <- sapply(pkgs, require, warn.conflicts = TRUE, character.only = TRUE)
if(any(!check)){
  pkgs.missing <- pkgs[!check]
  install.packages(pkgs.missing)
  check <- sapply(pkgs.missing, require, warn.conflicts = TRUE, character.only = TRUE)
}

## ---------------------------

drug <- read.delim(file = "./Drug_sensitivity_AUC_(PRISM_Repurposing_Secondary_Screen)_19Q4.csv", header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)
exp <- read.delim(file = "./Expression_22Q2_Public_subsetted (EGFR).csv", header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)

# Remove BRD IDs in parentheses from column names
colnames(drug) <- gsub(pattern = "\\ \\(.*", replacement = "", x = colnames(drug))

# Convert wide-format drug dataframe to long format where all drug responses are collapsed into "drug" and "AUC" columns
drug_long <- melt(drug, id.vars = colnames(drug)[1:6], measure.vars = colnames(drug)[7:ncol(drug)], variable.name = "drug", value.name = "AUC")

# Merge gene expression and drug response data
df_all <- merge(drug_long, exp, by = colnames(drug)[1:6])

# Filter for drugs of interest in the drugs2filt vector
drugs2filt <- c("dovitinib", "trametinib")
df_filt <- df_all[df_all$drug %in% drugs2filt,]

# Make scatter plots
plot <- ggplot(data = df_filt, mapping = aes(x = EGFR, y = AUC, color = drug)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~ drug) +
  labs(x = "EGFR expression", y = "PRISM AUC") +
  coord_cartesian(ylim = c(0, 1.1)) +
  scale_y_continuous(breaks = seq(0, 1.1, by = 0.1), labels = seq(0, 1.1, by = 0.1)) +
  scale_x_continuous(breaks = seq(0, 10, by = 1), labels = seq(0, 10, by = 1)) +
  geom_smooth(method = "lm", color = "black", lwd = 0.5) +
  stat_cor(method = "spearman", mapping = aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.y = 1.075, label.x = 0.1) +
  theme_linedraw() +
  theme(legend.position = "none")

plot
