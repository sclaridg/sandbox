## ---------------------------
##
## Script name: ComputeGrowthMetrics.R
##
## Purpose of script: Compute IC30, IC50, and AUC from a dose response.
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

pkgs <- c("tidyverse", "reshape2", "magrittr", "PharmacoGx", "drc")

check <- sapply(pkgs, require, warn.conflicts = TRUE, character.only = TRUE)
if(any(!check)){
  pkgs.missing <- pkgs[!check]
  install.packages(pkgs.missing)
  check <- sapply(pkgs.missing, require, warn.conflicts = TRUE, character.only = TRUE)
}

## ---------------------------

ComputeGrowthMetrics <- function(dataset_drug) {
  require(PharmacoGx)
  require(drc)

  print(paste0("Computing ", dataset_drug$drug[1], " growth metrics."))

  # Make dose response curve object
  drc <- drm(formula = norm ~ log10(concentration),
             curveid = drug,
             data = dataset_drug,
             fct = L.4(names = c("hill", "lower", "upper", "ED50")))

  # Calculate ED30 and ED50 from dose response curve object
  ic <- data.frame(ED(drc, c(50, 70),
                      interval = "delta",
                      type = "absolute",
                      display = FALSE))
  ic50 <- 10 ^ ic$Estimate[1]
  ic30 <- 10 ^ ic$Estimate[2]

  aac <- computeAUC(concentration = dataset_drug$concentration,
                    viability = dataset_drug$norm,
                    viability_as_pct = FALSE,
                    trunc = FALSE,
                    area.type = "Fitted")
  auc <- 1 - aac

  growth_metrics <- list(ic30, ic50, auc)

  return(growth_metrics)
}