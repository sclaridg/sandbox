## ---------------------------
##
## Script name: FilterDEGsByPathway.R
##
## Purpose of script: Filter differential expression datasets by genes in a given pathway.
##
## Author: Sally Claridge
##
## Date Created: 26 July 2023
##
## Copyright (c) Sally Claridge, 2023
## Email: sally.claridge@icahn.mssm.edu
##
## ---------------------------
##
## Notes:
##
##
## EDIT ---------------------------

input_filepath_pathway <- "./GO_term_summary_20230726_143700.csv"
input_filepath_DEGs <- "./20230323_JZP_KPA_Differential_expression_analysis_table.csv"
output_filepath <- "./20230726_test.csv"

## DO NOT EDIT BELOW THIS LINE ---------------------------

options(scipen = 6, digits = 4)

## Install and load the packages we will need:
pkgs <- c("tidyverse", "reshape2", "magrittr")

check <- sapply(pkgs, require, warn.conflicts = TRUE, character.only = TRUE)
if(any(!check)){
  pkgs.missing <- pkgs[!check]
  install.packages(pkgs.missing)
  check <- sapply(pkgs.missing, require, warn.conflicts = TRUE, character.only = TRUE)
}

# Read in pathway file from respective online database and differential expression tale from Genewiz
pathway <- read.delim(file = input_filepath_pathway, header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)
degs <- read.delim(file = input_filepath_DEGs, header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)

# Filter degs for only genes in the pathway of interest
degs_filt <- filter(degs, Gene.name %in% pathway$Symbol)

# Write dataframe to file
write.csv(file = output_filepath, x = degs_filt)