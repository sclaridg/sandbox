## ---------------------------
##
## Script name: DilutionSeries
##
## Purpose of script: Create a vector of concentration values for a dilution seroes starting from a max value, decreasing at a specified dilution
##
## Author: Sally Claridge
##
## Date Created: 02 July 2020
##
## Copyright (c) Sally Claridge, 2020
## Email: sally.claridge@icahn.mssm.edu
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

## set working directory for Mac and PC

DilutionSeries <- function(max, dilution_factor, dilution_n, low_to_high) {
  conc <- c()
  for(i in 0:(dilution_n - 1)) { conc <- c(conc, max/dilution_factor^i) }
  conc <- formatC(conc, digits = 4, format = "f")

  if(low_to_high) { conc <- rev(conc)}

  return(conc)
}

## ---------------------------