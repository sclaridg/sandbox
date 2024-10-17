## ---------------------------
##
## Script name: hallmark_pathway_analysis.R
##
## Purpose of script:
##
## Author: Sally Claridge
##
## Date Created: 17 October 2022
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

pkgs <- c("tidyverse", "reshape2", "magrittr", "msigdbr", "org.Hs.eg.db")

check <- sapply(pkgs, require, warn.conflicts = TRUE, character.only = TRUE)
if(any(!check)){
  pkgs.missing <- pkgs[!check]
  install.packages(pkgs.missing)
  check <- sapply(pkgs.missing, require, warn.conflicts = TRUE, character.only = TRUE)
}

## ---------------------------

hallmark_gene_sets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_pw <- split(x = hallmark_gene_sets$gene_symbol, f = hallmark_gene_sets$gs_name)

pw <- as.character(hallmark_pw$HALLMARK_P53_PATHWAY)

overlap <- read.delim(file = "./RNAseq_processed/20221012_KH-4-43_3CCL/DEGS_3CCL_overlapping.csv", header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)

overlap_up <- subset(x = overlap, direction == "Up")

intersect(pw, overlap_up$Gene)
