## ---------------------------
##
## Script name: process_DRC.R
##
## Purpose of script:
##
## Author: Sally Claridge
##
## Date Created: 06 April 2021
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

ctg <- read.delim(file = "./data_output/archive/SallyClaridge/20220131_CTG_IGROV1_OVCAR8_MEKi_UPSi_processed_data.csv", header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)
drug_name <- "33-11"
cell_line_name <- "OVCAR8"

outpath <- "./data_output/archive/SallyClaridge/20220131_CTG_IGROV1_OVCAR8_MEKi_UPSi_DRCs"

## ---------------------------

set.seed(1234)
options(scipen = 6, digits = 4)

## ---------------------------

## load up the packages we will need:

pkgs <- c("tidyverse", "reshape2", "magrittr", "ComplexHeatmap", "gghighlight", "plotly", "ggrepel", "drc")

check <- sapply(pkgs, require, warn.conflicts = TRUE, character.only = TRUE)
if(any(!check)){
  pkgs.missing <- pkgs[!check]
  install.packages(pkgs.missing)
  check <- sapply(pkgs.missing, require, warn.conflicts = TRUE, character.only = TRUE)
}

## ---------------------------

source("scripts/unfactorize.R")

## ---------------------------

ComputeGrowthMetrics <- function(dataset) {
  require(PharmacoGx)
  require(drc)

  # Make dose response curve object
  drc <- drm(formula = norm ~ log10(concentration),
             curveid = cell_line,
             data = dataset,
             fct = L.4(names = c("hill", "lower", "upper", "ED50")))

  # Calculate ED30 and ED50 from dose response curve object
  ic <- data.frame(ED(drc, c(50, 70),
                      interval = "delta",
                      type = "absolute",
                      display = FALSE))
  ic50 <- 10 ^ ic$Estimate[1]
  ic30 <- 10 ^ ic$Estimate[2]

  aac <- computeAUC(concentration = dataset$concentration,
                    viability = dataset$norm,
                    viability_as_pct = FALSE,
                    trunc = FALSE,
                    area.type = "Fitted")
  auc <- 1 - aac

  growth_metrics <- list(ic30, ic50, auc)

  return(growth_metrics)
}

## ---------------------------

CheckIC3050 <- function(dataset, growth_metrics) {
  # Check if ED30/50 and relative intensity values make sense
  test10 <- dplyr::filter(dataset, concentration == max(dataset$concentration))
  test0 <- dplyr::filter(dataset, concentration == min(dataset$concentration))
  avg_10um <- mean(test10$norm)
  avg_0um <- mean(test0$norm)

  ic30 <- growth_metrics[[1]]
  ic50 <- growth_metrics[[2]]

  # Assign +, -, or numerical value as the IC30 and IC50
  # if(avg_0um > 0.7 & avg_10um < 0.7 & ic30 <= 10 & ic30 >= 0.04) { ic30 <- ic30 }
  # if(avg_0um <= 0.7 & avg_10um <= 0.7 | ic30 < 0.04) { ic30 <- "-" }
  # if(avg_0um > 0.7 & avg_10um >= 0.7 | ic30 > 10) { ic30 <- "+" }
  # if(avg_0um > 0.5 & avg_10um < 0.5 & ic50 <= 10 & ic50 >= 0.04) { ic50 <- ic50 }
  # if(avg_0um <= 0.5 & avg_10um <= 0.5 | ic50 < 0.04) { ic50 <- "-" }
  # if(avg_0um > 0.5 & avg_10um >= 0.5 | ic50 > 10) { ic50 <- "+" }

  growth_metrics2 <- list("IC30" = ic30, "IC50" = ic50, "AUC" = growth_metrics[[3]])

  return(growth_metrics2)
}

## ---------------------------

MakeDRC <- function(dataset, drugname, cell_line, growth_metrics) {
  # Make plot object of dose response curves from
  # normalized CellTiter-Glo readout in long format

  require(drc)
  require(ggplot2)
  require(tidyverse)

  # Make plot title
  plot_title <- paste0("Cell line: ", cell_line, "\nDrug: ", drugname)

  # Make plot object of dose response curve
  plot <- ggplot(data = dataset) +
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
         y = "Relative intensity (CellTiter-Glo)",
         title = plot_title) +
    theme(legend.position = "none")

  ic30 <- growth_metrics[[1]]
  ic30_label <- ifelse(is.numeric(ic30), paste0("IC30: ", round(ic30, digits = 3), " uM"), paste0("IC30: ", ic30))
  ic50 <- growth_metrics[[2]]
  ic50_label <- ifelse(is.numeric(ic50), paste0("IC50: ", round(ic50, digits = 3), " uM"), paste0("IC50: ", ic50))
  auc <- growth_metrics[[3]]
  auc_label <- paste0("AUC: ", round(auc, digits = 3))

  # Print appropriate IC30/50 values on plot object
  plot <- plot +
    geom_text(y = 1.20, x = log10(7), mapping = aes(label = auc_label)) +
    geom_text(y = 1.15, x = log10(7), mapping = aes(label = ic30_label)) +
    geom_text(y = 1.10, x = log10(7), mapping = aes(label = ic50_label))
  # plot <- plot +
  #   geom_text(y = 1.20, x = 1, mapping = aes(label = auc_label)) +
  #   geom_text(y = 1.15, x = 1, mapping = aes(label = ic30_label)) +
  #   geom_text(y = 1.10, x = 1, mapping = aes(label = ic50_label))

  # Return plot object
  return(plot)
}

## ---------------------------

ctg2 <- filter(ctg, drug == drug_name)
ctg3 <- filter(ctg2, cell_line == cell_line_name)

metrics <- ComputeGrowthMetrics(ctg3)
metrics2 <- CheckIC3050(dataset = ctg3, growth_metrics = metrics)

plot_ctg <- MakeDRC(dataset = ctg3, drugname = drug_name, cell_line = cell_line_name, growth_metrics = metrics2)

plot_name <- paste0(format(Sys.Date(), "%Y%m%d"), "_", cell_line_name, "_", drug_name, "_DRC.png")

ggsave(plot = plot_ctg, filename = plot_name, path = outpath, device = "png", dpi = 300, width = 7, height = 7.5, units = "in", )

