## ---------------------------
##
## Script name: CTG_DRC_processing.R
##
## Purpose of script: Compute growth metrics for and make plot of the dose response of normalized CellTiter-Glo data.
##
## Author: Sally Claridge
##
## Date Created: 07 August 2024
##
## Copyright (c) Sally Claridge, 2024
## Email: sally.e.claridge@gmail.com
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

set.seed(1234)
options(scipen = 6, digits = 4)

## ---------------------------

# Uncomment code in this section to install the following packages, if needed

# install.packages("tidyverse")
# install.packages("reshape2")
# install.packages("magrittr")
# install.packages("drc")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("PharmacoGx")

## ---------------------------

# When/if installed, load the packages

pkgs <- c("tidyverse", "reshape2", "magrittr", "PharmacoGx", "drc")

check <- sapply(pkgs, require, warn.conflicts = TRUE, character.only = TRUE)
if(any(!check)){
  pkgs.missing <- pkgs[!check]
  install.packages(pkgs.missing)
  check <- sapply(pkgs.missing, require, warn.conflicts = TRUE, character.only = TRUE)
}

## ---------------------------

# Input CSV must have the following case-sensitive column names: drug, cell_line, concentration, norm
df <- read.delim(file = "path/to/input_data.csv", header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)

# Define output directory where the dose resopnse plot will be saved
output_dir <- "path/to/save/output/plot/"

# Subset df by cell line and drug
cell_line_name <- "CELL_LINE"
drugname <- "DRUG"
df <- subset(df, drug == drugname & cell_line == cell_line_name)

# If you want to add a blood plasma concentration line, run the following line uncomment the geom_vline and annotate calls in the plotting function
# Else bpc object is ignored
bpc <- 2.68

# Calculate AUC

# Function returns the AUC (Area Under the drug response Curve) given concentration and viability as input, normalized by the concentration range of the experiment.
# The area returned is the response (1-Viablility) area, i.e. area under the curve when the response curve is plotted on a log10 concentration scale, with high AUC implying high sensitivity to the drug
# We reverse this
aac <- PharmacoGx::computeAUC(concentration = df$concentration,
                  viability = df$norm,
                  viability_as_pct = FALSE,
                  trunc = FALSE,
                  area.type = "Fitted")
auc <- 1 - aac

# Calculate IC30 and IC50
drc <- drc::drm(formula = norm ~ log10(concentration),
                curveid = drug,
                data = df,
                fct = L.4(names = c("hill", "lower", "upper", "ED50")))

ic <- data.frame(drc::ED(drc, c(50, 70),
                         interval = "delta",
                         type = "absolute",
                         display = FALSE))
ic50 <- 10 ^ ic$Estimate[1]
ic30 <- 10 ^ ic$Estimate[2]

# Make text objects for plotting
title_text <- paste0("Cell line: ", cell_line_name, "\nDrug: ", drugname)
subtitle_text <- paste0("AUC: ", round(auc, digits = 3), ", IC30: ", round(ic30, digits = 3), " uM, IC50: ", round(ic50, digits = 3), " uM")

# Make plot object
plot <- ggplot(data = df, mapping = aes(x = concentration, y = norm)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "drm",
              se = FALSE,
              method.args = list(fct = L.4(names = c("Slope",
                                                     "Lower Limit",
                                                     "Upper Limit",
                                                     "ED50")))) +
  scale_y_continuous(breaks = seq(0, 20.0, by = 0.2)) +
  coord_cartesian(ylim = c(0, round(max(df$norm), digits = 1))) +
  # coord_cartesian(ylim = c(0, 1.2)) + # uncomment if you want to control the y-axis dimensions, this call will overwrite the previous coord_cartesian call
  scale_x_log10() +
  theme_bw() +
  labs(x = "Concentration (nM)",
       y = "Relative Intensity (CellTiter-Glo)",
       title = title_text,
       subtitle = subtitle_text) +
  # geom_vline(xintercept = bpc, lty = 2, alpha = 2) +
  # annotate(x = bpc, #change this variable to change position along x-axis
  #          y = 0.65, geom = "label",
  #          label = paste0("BPC = ", bpc, " uM"),
  #          size = 5, fill = "white", alpha = 0.7) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        plot.title = element_text(size = 22),
        plot.subtitle = element_text(size = 15))

plot

# Save plot to output directory
# ggsave(plot = plot,
#        filename = paste0(output_dir, "/", Sys.Date(), "_", cell_line_name, "_", drugname,"_DRC.png"),
#        device = "png",
#        dpi = 500,
#        width = 6,
#        height = 6,
#        units = "in")


## ---------------------------

sessionInfo()

## ---------------------------