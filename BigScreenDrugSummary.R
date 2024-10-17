## ---------------------------
##
## Script name: BigScreenDrugSummary.R
##
## Purpose of script: For the single drug in the provided CSV file of
## processed_data, make a heatmap of the AUC values across all cell lines and
## overlay all dose response curves on a single plot.
##
## Author: Sally Claridge
##
## Date Created: 08 November 2021
##
## Copyright (c) Sally Claridge, 2021
## Email: sally.claridge@icahn.mssm.edu
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

## Global variables you edit

drug_name <- "KH-4-43"
# Must match drug name in the processed_data and AUC files

processed_data_path <- "./data_output/archive/SallyClaridge/20210413_KH-4-43_ovarian_processed_data_truncated.csv"
# Must have these case-sensitive columns: concentreation, drug, norm, cell_line

AUC_path <- "./data_output/growth_metrics_merged/20210407_AUC_merged.csv"
# Can be z-scores or raw AUC (heatmap will say AUC on legend for both cases)

output_directory <- "./data_output/archive/SallyClaridge/"
# Must end with "/"

## Do not edit below this line

## ---------------------------

options(scipen = 6, digits = 4)

## ---------------------------

## load up the packages we will need:

pkgs <- c("tidyverse", "reshape2", "magrittr", "drc", "ComplexHeatmap")

check <- sapply(pkgs, require, warn.conflicts = TRUE, character.only = TRUE)
if(any(!check)){
  pkgs.missing <- pkgs[!check]
  install.packages(pkgs.missing)
  check <- sapply(pkgs.missing, require, warn.conflicts = TRUE, character.only = TRUE)
}

## ---------------------------

## load up our functions into memory

source("./scripts/MakeDrugDataframe.R")

## ---------------------------

## Overlay dose response curves

## Read in data
dataset_init <- read.delim(file = processed_data_path, header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)
dataset_processed <- MakeDrugDataframe(dataset = dataset_init, drugname = drug_name)

## Make plot

plot <- ggplot(data = dataset_processed, mapping = aes(x = concentration, y = norm, color = cell_line)) +
  geom_point() +
  stat_smooth(method = "drm",
              se = FALSE,
              method.args = list(fct = L.4(names = c("Slope",
                                                     "Lower Limit",
                                                     "Upper Limit",
                                                     "ED50")))) +
  scale_y_continuous(breaks = seq(0, 1.4, by = 0.2)) +
  coord_cartesian(ylim = c(0, 1.2)) +
  scale_x_log10() +
  theme_bw() +
  labs(x = "Concentration (uM)",
       y = "Relative intensity (CellTiter-Glo)",
       title = paste0("Drug: ", drug_name),
       color = "Cell line")

## Save plot

ggsave(plot = plot, filename = paste0(output_directory, format(Sys.Date(), "%Y%m%d"), "_", drug_name, "_all_DRCs.png"),
       device = "png", dpi = 500, width = 5, height = 4, units = "in")

## ---------------------------

## Make heatmap

## Read in data

z <- read.delim(file = AUC_path, header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)

## Format data

z_drug <- z[z$Drug == drug_name, ]
z_drug$Target <- NULL
z_drug <- data.frame(z_drug[, -1], row.names = z_drug[, 1])
z_drug <- data.frame(t(z_drug))
z_drug$cell_line_name <- rownames(z_drug)

## Read in lineages for annotation

lin <- read.delim(file = "./data_raw/masters/master_aaronson_line_lineages_withA2780.csv", header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)

## Format lineate annotation

z_drug <- merge(z_drug, lin, by = "cell_line_name")
rownames(z_drug) <- z_drug$cell_line_name

## Convert to matrix

annot <- matrix(z_drug[, 3])
colnames(annot) <- "Lineage"
rownames(annot) <- rownames(z_drug)
body <- matrix(z_drug[, 2])
colnames(body) <- "KH-4-43"
rownames(body) <- rownames(z_drug)

## Make heatmap

lin_row_annot <- rowAnnotation(Lineage = annot[,1])
ht <- Heatmap(matrix = body, name = "AUC", row_title = "Cell line", column_title = "", right_annotation = lin_row_annot)

## Save heatmap

png(paste0(output_directory, format(Sys.Date(), "%Y%m%d"), "_", drug_name, "_heatmap.png"), height = 10, width = 4, units = "in", res = 500)
draw(ht)
dev.off()

## ---------------------------

## Print R session information

sessionInfo()
