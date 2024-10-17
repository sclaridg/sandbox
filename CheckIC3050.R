## ---------------------------
##
## Script name: CheckIC3050
##
## Purpose of script: Check validity of IC30/50 values and replace with +/- accordingly.
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

pkgs <- c("tidyverse", "reshape2", "magrittr")

check <- sapply(pkgs, require, warn.conflicts = TRUE, character.only = TRUE)
if(any(!check)){
  pkgs.missing <- pkgs[!check]
  install.packages(pkgs.missing)
  check <- sapply(pkgs.missing, require, warn.conflicts = TRUE, character.only = TRUE)
}


## ---------------------------

CheckIC3050 <- function(dataset_drug, growth_metrics) {
  # Check if ED30/50 and relative intensity values make sense
  test10 <- dplyr::filter(dataset_drug, concentration == max(dataset_drug$concentration))
  test0 <- dplyr::filter(dataset_drug, concentration == min(dataset_drug$concentration))
  avg_10um <- mean(test10$norm)
  avg_0um <- mean(test0$norm)

  ic30 <- growth_metrics[[1]]
  ic50 <- growth_metrics[[2]]

  # Assign +, -, or numerical value as the IC30 and IC50
  if(avg_0um > 0.7 & avg_10um < 0.7 & ic30 <= 10 & ic30 >= 0.04) { ic30 <- ic30 }
  if(avg_0um <= 0.7 & avg_10um <= 0.7 | ic30 < 0.04) { ic30 <- "-" }
  if(avg_0um > 0.7 & avg_10um >= 0.7 | ic30 > 10) { ic30 <- "+" }
  if(avg_0um > 0.5 & avg_10um < 0.5 & ic50 <= 10 & ic50 >= 0.04) { ic50 <- ic50 }
  if(avg_0um <= 0.5 & avg_10um <= 0.5 | ic50 < 0.04) { ic50 <- "-" }
  if(avg_0um > 0.5 & avg_10um >= 0.5 | ic50 > 10) { ic50 <- "+" }

  growth_metrics2 <- list(ic30, ic50, growth_metrics[[3]])

  return(growth_metrics2)
}