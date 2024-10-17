## ---------------------------
##
## Script name: VennDiagramGSEA.R
##
## Purpose of script:
##
## Author: Sally Claridge
##
## Date Created: 31 July 2023
##
## Copyright (c) Sally Claridge, 2023
## Email: sally.claridge@icahn.mssm.edu
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

## load up the packages we will need:

pkgs <- c("tidyverse", "reshape2", "magrittr", "ggvenn", "ggVennDiagram", "RColorBrewer")

check <- sapply(pkgs, require, warn.conflicts = TRUE, character.only = TRUE)
if(any(!check)){
  pkgs.missing <- pkgs[!check]
  install.packages(pkgs.missing)
  check <- sapply(pkgs.missing, require, warn.conflicts = TRUE, character.only = TRUE)
}

## Compare 2 Cell Lines ---------------------------

# START EDITING

gsea1 <- read.delim(file = "./1.csv", header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)
cell_line1 <- "1" # cell line name

gsea2 <- read.delim(file = "./2.csv", header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)
cell_line2 <- "2"  # cell line name

output_path <- "./" # Must end in `/`

# STOP EDITING

gsea1_signif <- filter(gsea1, padj < 0.05)$terms
gsea2_signif <- filter(gsea2, padj < 0.05)$terms

gsea_list <- list(gsea1_signif, gsea2_signif)
names(gsea_list) <- c(cell_line1, cell_line2)

intersection <- as.data.frame(intersect(gsea1_signif, gsea2_signif))
colnames(intersection) <- "intersection"

venn_plot <- ggvenn(gsea_list,
       show_elements = FALSE,
       label_sep = "\n",
       fill_color = brewer.pal(name = "Set1", n = 3),
       show_percentage = FALSE,
       auto_scale = TRUE,
       stroke_size = 1,
       text_size = 6,
       set_name_size = 6)

ggsave(plot = venn_plot,
       filename = paste0(cell_line1, "_", cell_line2, "_venn.png"),
       path = output_path,
       dpi = 300, device = "png", units = "in", height = 5, width = 5)

write.csv(x = intersection,
          file = paste0(output_path, cell_line1, "_", cell_line2, "_intersection.csv"))

## Compare 6 Cell Lines ---------------------------

# START EDITING

gsea1 <- read.delim(file = "./1.csv", header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)
cell_line1 <- "1"  # cell line name

gsea2 <- read.delim(file = "./2.csv", header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)
cell_line2 <- "2"  # cell line name

gsea3 <- read.delim(file = "./3.csv", header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)
cell_line3 <- "3"  # cell line name

gsea4 <- read.delim(file = "./4.csv", header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)
cell_line4 <- "4"  # cell line name

gsea5 <- read.delim(file = "./5.csv", header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)
cell_line5 <- "5"  # cell line name

gsea6 <- read.delim(file = "./6.csv", header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)
cell_line6 <- "6"  # cell line name

output_path <- "./" # Must end in `/`

# STOP EDITING

gsea1_signif <- filter(gsea1, padj < 0.05)$terms
gsea2_signif <- filter(gsea2, padj < 0.05)$terms
gsea3_signif <- filter(gsea3, padj < 0.05)$terms
gsea4_signif <- filter(gsea4, padj < 0.05)$terms
gsea5_signif <- filter(gsea5, padj < 0.05)$terms
gsea6_signif <- filter(gsea2, padj < 0.05)$terms

gsea_list <- list(gsea1_signif, gsea2_signif, gsea3_signif, gsea4_signif, gsea5_signif, gsea6_signif)
names(gsea_list) <- c(cell_line1, cell_line2, cell_line3, cell_line4, cell_line5, cell_line6)

intersection <- as.data.frame(Reduce(intersect, gsea_list))
colnames(intersection) <- "intersection"

venn_plot <- ggVennDiagram(x = gsea_list,
                           category.names = names(gsea_list),
                           edge_size = 4,
                           label = "count",
                           label_alpha = 0) +
  scale_fill_distiller(palette = "Blues") +
  scale_color_manual(values = c("red", "orange", "yellow", "limegreen", "purple", "magenta"))

ggsave(plot = venn_plot,
       filename = paste0(cell_line1, "_", cell_line2, "_", cell_line3, "_", cell_line4, "_", cell_line5, "_", cell_line6, "_venn.png"),
       path = output_path,
       dpi = 300, device = "png", units = "in", height = 6, width = 6)

write.csv(x = intersection,
          file = paste0(output_path, cell_line1, "_", cell_line2, "_", cell_line3, "_", cell_line4, "_", cell_line5, "_", cell_line6, "_intersection.csv"))