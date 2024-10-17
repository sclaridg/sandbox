## ---------------------------
##
## Script name: MakeDRC.R
##
## Purpose of script: Make single dose response curve plot object.
##
## Author: Sally Claridge
##
## Date Created: 27 July 2021
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

options(scipen = 6, digits = 4)

## ---------------------------

## load up the packages we will need:

pkgs <- c("tidyverse", "reshape2", "magrittr", "drc")

check <- sapply(pkgs, require, warn.conflicts = TRUE, character.only = TRUE)
if(any(!check)){
  pkgs.missing <- pkgs[!check]
  install.packages(pkgs.missing)
  check <- sapply(pkgs.missing, require, warn.conflicts = TRUE, character.only = TRUE)
}

## ---------------------------

MakeDRC <- function(dataset, drugname, cell_line, is_combo, DMSO_or_TX, combo_drug, combo_drug_conc, whichIC, endpoint_type) {
  # Make plot object of dose response curves from
  # normalized CellTiter-Glo readout in long format

  require(drc)
  require(ggplot2)
  require(tidyverse)

  # Check file formatting and filter for drug of interest
  dataset_drug <- MakeDrugDataframe(dataset, drugname)

  # Make dose response curve object and calculate ED30 and ED50 from dose response curve object
  growth_metrics <- ComputeGrowthMetrics(dataset_drug)

  # Make plot title
  if(!is_combo) {
    plot_title <- paste0("Cell line: ", cell_line, "\nDrug: ", drugname)
  }
  if(is_combo) {
    if(DMSO_or_TX == "TX") {
      if(whichIC == "IC50") {
        plot_title <- paste0("Cell line: ", cell_line, "\nDrug: ", drugname, " and ", combo_drug, " (IC50 = ", combo_drug_conc, ")")
      }
      if(whichIC == "IC30") {
        plot_title <- paste0("Cell line: ", cell_line, "\nDrug: ", drugname, " and ", combo_drug, " (IC30 = ", combo_drug_conc, ")")
      }
    }
    if(DMSO_or_TX == "DMSO") {
      plot_title <- paste0("Cell line: ", cell_line, "\nDrug: ", drugname, " and DMSO")
    }
  }
  # Make plot object of dose response curve
  plot <- ggplot(data = dataset_drug) +
    geom_point(mapping = aes(x = concentration, y = norm)) +
    stat_smooth(mapping = aes(x = concentration, y = norm),
                method = "drm",
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
         y = ifelse(endpoint_type == "CTG", "Relative intensity (CellTiter-Glo)", ifelse(endpoint_type == "confluence", "Relative confluence", "Arbitrary units")),
         title = plot_title) +
    theme(legend.position = "none")

  # Check if ED30/50 and relative intensity values make sense
  growth_metrics2 <- CheckIC3050(dataset_drug, growth_metrics)
  ic30 <- growth_metrics2[[1]]
  ic30_label <- ifelse(is.numeric(ic30), paste0("IC30: ", round(ic30, digits = 2), " uM"), paste0("IC30: ", ic30))
  ic50 <- growth_metrics2[[2]]
  ic50_label <- ifelse(is.numeric(ic50), paste0("IC50: ", round(ic50, digits = 2), " uM"), paste0("IC50: ", ic50))
  auc <- growth_metrics2[[3]]
  auc_label <- paste0("AUC: ", round(auc, digits = 2))

  # Print appropriate IC30/50 values on plot object
  plot <- plot +
    geom_text(y = 1.20, x = log10(7), mapping = aes(label = auc_label)) +
    geom_text(y = 1.15, x = log10(7), mapping = aes(label = ic30_label)) +
    geom_text(y = 1.10, x = log10(7), mapping = aes(label = ic50_label))

  # Return plot object
  return(plot)
}