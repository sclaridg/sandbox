## ---------------------------
##
## Script name: MapGeneIDs.R
##
## Purpose of script: For a dataframe with a column of ENSEMBL IDs,
##                    add a Gene.name column with gene symbols and
##                    write CSV to file.
##
## Author: Sally Claridge
##
## Date Created: 09 May 2023
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

input_filepath <- "./BenjaminHopkins/normalized_counts.csv"
output_filepath <- "./BenjaminHopkins/normalized_counts_genenames.csv"

## DO NOT EDIT BELOW THIS LINE ---------------------------

options(scipen = 6, digits = 4)

## install packages and load libraries
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)

## read in dataframe
df <- read.delim(file = input_filepath, header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)

## rename first column
colnames(df)[1] <-  "ENSEMBL_id"

# add Gene.name column, mapping ENSEMBL IDs to gene symbols
df$Gene.name <- AnnotationDbi::mapIds(org.Mm.eg.db,
                                      keys = df$ENSEMBL_id,
                                      column = "SYMBOL",
                                      keytype = "ENSEMBL",
                                      multiVals = "first")
## write dataframe to file
write.csv(file = output_filepath, x = df)
