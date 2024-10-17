## ---------------------------
##
## Script name: AUCHeatmap.R
##
## Purpose of script:
##
## Author: Sally Claridge
##
## Date Created: 19 September 2023
##
## Copyright (c) Sally Claridge, 2023
## Email: sally.claridge@icahn.mssm.edu
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

## Global variables

auc_df <- read.delim(file = "./BenjaminHopkins/20230317_AUC_merged_51lines_lessdrugs_noNAs for SHalini with aaronson name changed.csv", header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)

auc_df$Target <- NULL # Do not need the Target column; if auc file doesn't have this column, comment out this line

# Required for lineage annotation (add/remove cell lines as necessary)
# Only cell lines in this file will be included in the heatmap
# Column names must be "cell_line_name", and "lineage"
lineages <- read.delim(file = "./BenjaminHopkins/screened_lineages.csv", header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)

# Adjust width and height of heatmaps according to cell lines and drugs included
# For 35 cell lines, use 10 (accounts for height of drug names)
# For 117 drugs, use 22 (accounts for width of cell lines names and width of legend)
width_inches <- 22
height_inches <- 10

output_path <- "./BenjaminHopkins/" # Must end in `/`

# Optional: Only include given drugs in heatmap
# Input CSV file must have one column named "Drug"
# Un-comment the following two lines if you want to filter for specific drugs and run whole script in order
# drug_filter <- read.delim(file = "./drug_list.csv", header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)
# auc_df <- filter(auc_df, Drug %in% drug_filter$Drug)

## DO NOT EDIT BELOW THIS LINE ---------------------------

## load up the packages we will need:

pkgs <- c("tidyverse", "reshape2", "magrittr", "ComplexHeatmap", "circlize", "viridis")

check <- sapply(pkgs, require, warn.conflicts = TRUE, character.only = TRUE)
if(any(!check)){
  pkgs.missing <- pkgs[!check]
  install.packages(pkgs.missing)
  check <- sapply(pkgs.missing, require, warn.conflicts = TRUE, character.only = TRUE)
}

## ---------------------------

# Assign "Drug" names to rownames
rownames(auc_df) <- auc_df$Drug
# Delete Drug column
auc_df$Drug <- NULL
auc <- auc_df %>% dplyr::select(lineages$cell_line_name)
auc <- auc %>% select(order(colnames(auc)))


# Raw AUC matrix
auc_raw <- t(auc)

# z by drug matrix
auc_ZxDrug <- as.data.frame(apply(auc, 1, scale))
rownames(auc_ZxDrug) <- colnames(auc)
auc_ZxDrug <- as.matrix(auc_ZxDrug)

# z by CCL matrix
auc_ZxCCL <- as.data.frame(apply(auc, 2, scale))
auc_ZxCCL <- t(auc_ZxCCL)
colnames(auc_ZxCCL) <- rownames(auc)

# Set annotation
lin_col_values <- data.frame(color = viridis::viridis(length(unique(lineages$lineage))),
                             lineage = unique(lineages$lineage))
lineages <- merge(lineages, lin_col_values, by = "lineage")
lineages <- with(lineages,  lineages[order(cell_line_name), ])

ann <- data.frame(Lineage = lineages$lineage)
color_list <- c(lineages$color)
names(color_list) <- lineages$lineage

colours <- list("Lineage" = color_list)
rowAnn <- HeatmapAnnotation(df = ann,
                            which = 'row',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'))


# Make heatmap objects
# Resource for ComplexHeatmap package options: https://jokergoo.github.io/ComplexHeatmap-reference/book/
# Look into the clustering options and see what customizations you can make
ht_raw <- Heatmap(matrix = auc_raw, # The data
                  name = "AUC", # The title that will appear over the legend
                  col = circlize::colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c("navy", "mediumblue", "white", "red", "darkred")), # Change heatmap colors
                  row_title = NA,
                  column_title = NA,
                  show_row_dend = FALSE,
                  show_column_dend = FALSE,
                  cluster_rows = TRUE, cluster_columns = TRUE)

drg_max <- pmax(abs(floor(range(auc_ZxDrug)[1])), round(range(auc_ZxDrug)[2]))
ht_drg <- Heatmap(matrix = auc_ZxDrug, # The data
                  name = "AUC z-score", # The title that will appear over the legend
                  col = circlize::colorRamp2(c(-drg_max, -drg_max/2, 0, drg_max/2, drg_max), c("navy", "mediumblue", "white", "red", "darkred")), # Change heatmap colors
                  row_title = NA,
                  column_title = NA,
                  right_annotation = rowAnn,
                  show_row_dend = FALSE,
                  show_column_dend = FALSE,
                  cluster_rows = TRUE, cluster_columns = TRUE)

CCL_max <- pmax(abs(floor(range(auc_ZxCCL)[1])), round(range(auc_ZxCCL)[2]))
ht_CCL <- Heatmap(matrix = auc_ZxCCL, # The data
                  name = "AUC z-score", # The title that will appear over the legend
                  col = circlize::colorRamp2(c(-CCL_max, -CCL_max/2, 0, CCL_max/2, CCL_max), c("navy", "mediumblue", "white", "red", "darkred")), # Change heatmap colors
                  row_title = NA,
                  column_title = NA,
                  right_annotation = rowAnn,
                  show_row_dend = FALSE,
                  show_column_dend = FALSE,
                  cluster_rows = TRUE, cluster_columns = TRUE)


png(paste0(output_path, "heatmap_AUC_raw.png"), height = height_inches, width = width_inches, units = "in", res = 500)
draw(ht_raw)
dev.off()

png(paste0(output_path, "heatmap_AUC_ZxDrug.png"), height = height_inches, width = width_inches, units = "in", res = 500)
draw(ht_drg)
dev.off()

png(paste0(output_path, "heatmap_AUC_ZxCCL.png"), height = height_inches, width = width_inches, units = "in", res = 500)
draw(ht_CCL)
dev.off()