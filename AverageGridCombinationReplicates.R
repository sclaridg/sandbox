## ---------------------------
##
## Script name: AverageGridCombo
##
## Purpose of script: Average the three replicates from combination grid screens.
##
## Author: Sally Claridge
##
## Date Created: 17 December 2020
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

source("./scripts/unfactorize.R")

pkgs <- c("tidyverse", "magrittr")
check <- sapply(pkgs, require, warn.conflicts = TRUE, character.only = TRUE)
if(any(!check)){
  pkgs.missing <- pkgs[!check]
  install.packages(pkgs.missing)
  check <- sapply(pkgs.missing, require, warn.conflicts = TRUE, character.only = TRUE)
}

ReadData <- function(file) {
  file <-  read.delim(file,
                      header = FALSE,
                      sep = ",",
                      check.names = FALSE,
                      stringsAsFactors = FALSE)

  start <- which(file[, 1] == "A") - 1
  cols <- file[start, ]
  colnames(file) <- cols
  file <- file[(start + 1):(which(file[, 1] == "P")), ]

  rows <- file[, 1]
  file[, 1] <- NULL
  file <- file %>% mutate_all(type.convert)

  rownames(file) <- rows

  return(file)
}

path <- "./data_raw/archive/SallyClaridge/20220325_HA336_grids/"

files <- paste0(path, list.files(path = path))
data_raw <- lapply(files, ReadData)
data_avg <- Reduce(`+`, data_raw) / length(data_raw)
write.csv(data_avg, file = paste0(path, "average.csv"))

















