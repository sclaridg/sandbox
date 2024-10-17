## ---------------------------
##
## Script name: differential_expression_analysis.R
##
## Purpose of script: Analyze differential gene expression data provided by Genewiz.
##
## Author: Sally Claridge
##
## Date Created: 26 September 2022
##
## Copyright (c) Sally Claridge, 2022
## Email: sally.claridge@icahn.mssm.edu
##
## ---------------------------
##
## Notes:
##
##
## Install Packages ------------------------------------------------------------
# Install these packages the first time you ever run this script
# Comment out these install lines for subsequent runs
# install.packages("tidyverse")
# BiocManager::install("ComplexHeatmap")
# install.packages("ggplot2")
# BiocManager::install("org.Mm.eg.db")
# install.packages("msigdbr")
# BiocManager::install("fgsea")

## Read in data ----------------------------------------------------------------
# Make global variables that will be used to add appropriate labels to plots
cell_line_name <- "A2780"
drug_name <- "KH-4-43"

# Read in DEGs
degs <- read.delim(file = "./Genewiz/20200831_Genewiz_A2780_KH-4-43_Differential_expression_analysis_table.csv", header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)
# Next line is for A2780 +/- KH-4-43 or cisplatin only
degs$log2FoldChange <- -degs$log2FoldChange

# Provide path to output directory where you want all files saved to, e.g. "./RNAseq/"
output_path <- "./RNAseq_processed/"
# These next two lines of code will make a dated directory within the output_path folder
outdir <- paste0(output_path, format(Sys.Date(), "%Y%m%d"), "_", cell_line_name,"_", drug_name)
if(!dir.exists(file.path(outdir))) { dir.create(file.path(outdir), recursive = TRUE) }

## Process DEGs ----------------------------------------------------------------
# Assign Up, Down, Not to differential expression based on log2FC and adjusted p-value
degs$diffexpressed <- "Not"
# if log2Foldchange > 1 and pvalue < 0.05, set as "UP"
degs$diffexpressed[degs$log2FoldChange > 1 & degs$padj < 0.05] <- "Up"
# if log2Foldchange < -1 and pvalue < 0.05, set as "DOWN"
degs$diffexpressed[degs$log2FoldChange < -1 & degs$padj < 0.05] <- "Down"

## Make Heatmap ----------------------------------------------------------------
library(ComplexHeatmap)

ht_mat <- subset(degs, padj < 0.05)
ht_mat <- ht_mat[, c(1, 5:10)]
rownames(ht_mat) <- ht_mat$ID
ht_mat$ID <- NULL
ht_mat <- as.matrix(ht_mat)

ht <- Heatmap(matrix = ht_mat,
              name = "Differential Gene Expression",
              row_title = "Genes",
              column_title = paste0(cell_line_name, ": ", drug_name, " treated vs. control"),
              show_row_names = FALSE)

# To save the heatmap as a png, run the next three lines all together
# png(filename = paste0(outdir, "/heatmap.png"), height = 15, width = 5, units = "in", res = 500)
# draw(ht)
# dev.off()

## Make Volcano Plot -----------------------------------------------------------
library(ggplot2)

volcano_plot <- ggplot(data = degs) +
  geom_point(mapping = aes(x = log2FoldChange, y = -log(padj), col = diffexpressed), alpha = 0.5) +
  scale_color_manual(values = c("blue", "black", "red")) +
  labs(y = expression(-log[10](p[adj])),
       x = expression(log[2](Fold~Change)),
       title = paste0(cell_line_name, ": ", drug_name, " treated vs. control"),
       color = "Differential expression status")  +
  theme_light() +
  theme(legend.position = "bottom")
# volcano_plot

# ggsave(plot = volcano_plot, filename = "volcano.png", path = outdir, dpi = 300, device = "png", units = "in", height = 5, width = 6)

## GSEA: Get Pathways ----------------------------------------------------------
# MOUSE
# Convert gene names to ENTREZ IDs
# library(org.Mm.eg.db)
#
# degs$ENTREZID_id <- mapIds(x = org.Mm.eg.db,
#                            keys = degs$ID,
#                            column = "ENTREZID",
#                            keytype = "ENSEMBL",
#                            multiVals = "first")
# ## Remove DEG rows with the same ENTREZ_id (the data loss isn't very significant)
# library(tidyverse)
# degs_filt <- degs %>% group_by(ENTREZID_id) %>% count() %>% filter(n != 1)
# degs <- subset(degs, !(ENTREZID_id %in% degs_filt$ENTREZID_id))
# rownames(degs) <- degs$ENTREZID_id
#
# # Retrieve Hallmark gene sets for Mus musculus
# library(msigdbr)
# # Note: You can set the category flag to be different thigns to grab various gene set types
# # Brief intro here https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html
# hallmark_gene_sets <- msigdbr(species = "Mus musculus", category = "H")
# hallmark_pw <- split(x = hallmark_gene_sets$entrez_gene, f = hallmark_gene_sets$gs_name)

#HUMAN
# Convert gene names to ENTREZ IDs
library(org.Hs.eg.db)
degs$ENTREZID_id <- mapIds(org.Hs.eg.db,
                           keys = degs$ID,
                           column = "ENTREZID",
                           keytype = "ENSEMBL",
                           multiVals = "first")
library(tidyverse)
degs_filt <- degs %>% group_by(ENTREZID_id) %>% count() %>% filter(n != 1)
degs <- subset(degs, !(ENTREZID_id %in% degs_filt$ENTREZID_id))
rownames(degs) <- degs$ENTREZID_id
library(msigdbr)
# Note: You can set the category flag to be different thigns to grab various gene set types
# Brief intro here https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html
hallmark_gene_sets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_pw <- split(x = hallmark_gene_sets$entrez_gene, f = hallmark_gene_sets$gs_name)

# Pathway databases
# data(carta.hs) # BioCarta
# data(kegg.sets.hs) # Kyoto Encyclopedia of Genes and Genomes
# data(go.sets.hs) # Gene Ontology
# reactome_pw <- reactomePathways(degs$ENTREZID_id) # Reactome
hallmark_gene_sets <- msigdbr(species = "human", category = "H") # Hallmark
hallmark_pw <- split(x = hallmark_gene_sets$human_entrez_gene, f = hallmark_gene_sets$gs_name)




## GSEA: Rank Genes ------------------------------------------------------------
# Rank order DEGs to use with the fgsea package
degs_ordered <- degs[order(degs$log2FoldChange, decreasing = TRUE),]
degs_ranked_FC <- degs_ordered$log2FoldChange
names(degs_ranked_FC) <- degs_ordered$ENTREZID_id

# Optional: Can look at enrichment plots for specific pathways by name
library(fgsea)

# values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#e31a1c", "#fb9a99", "#fdbf6f")
# breaks = c("COv362", "COv434", "SkOv3", "HA336", "IGROV1", "OvCar8", "A2780")
p53 <- plotEnrichment(hallmark_pw[["HALLMARK_P53_PATHWAY"]], degs_ranked_FC, gseaParam = 1) +
  geom_line(color = "#fdbf6f", lwd = 1.5) +
  labs(y = NULL, x = NULL) +
  theme(text = element_text(size = 20))
ggsave(plot = p53, filename = "p53_enrichment.png", path = outdir, dpi = 300, device = "png", units = "in", height = 6, width = 10)

kras <- plotEnrichment(hallmark_pw[["HALLMARK_KRAS_SIGNALING_UP"]], degs_ranked_FC, gseaParam = 1) +
  geom_line(color = "#fdbf6f", lwd = 1.5) +
  labs(y = NULL, x = NULL) +
  theme(text = element_text(size = 20))
ggsave(plot = kras, filename = "kras_enrichment.png", path = outdir, dpi = 300, device = "png", units = "in", height = 6, width = 10)

e2f <- plotEnrichment(hallmark_pw[["HALLMARK_E2F_TARGETS"]], degs_ranked_FC, gseaParam = 1) +
  geom_line(color = "#fdbf6f", lwd = 1.5) +
  labs(y = NULL, x = NULL) +
  theme(text = element_text(size = 20))
ggsave(plot = e2f, filename = "e2f_enrichment.png", path = outdir, dpi = 300, device = "png", units = "in", height = 6, width = 10)

g2m <- plotEnrichment(hallmark_pw[["HALLMARK_G2M_CHECKPOINT"]], degs_ranked_FC, gseaParam = 1) +
  geom_line(color = "#fdbf6f", lwd = 1.5) +
  labs(y = NULL, x = NULL) +
  theme(text = element_text(size = 20))
ggsave(plot = g2m, filename = "g2m_enrichment.png", path = outdir, dpi = 300, device = "png", units = "in", height = 6, width = 10)



 ## GSEA: Run and Process Output ------------------------------------------------
# Run GSEA with the fgsea package
hallmark_res <- fgsea(pathways = hallmark_pw, degs_ranked_FC, maxSize = 500)
hallmark_res$signif_label <- hallmark_res$signif_label <- ifelse(hallmark_res$padj < 0.05, "Adjusted p-value < 0.05", "Not significant")
hallmark_res$leadingEdge <- NULL
write.csv(x = hallmark_res, file = paste0(outdir, "/fgsea_hallmark_all_results.csv"), row.names = FALSE)

# Can collapse pathways to isolate the independent pathways
## This may take a while to run if you change the category flag above in the msigdbr function call
collapsed_pws <- collapsePathways(fgseaRes = hallmark_res[order(padj)][padj < 0.05],
                                  pathways = hallmark_pw,
                                  stats = degs_ranked_FC)
collapsed <- hallmark_res[pathway %in% collapsed_pws$mainPathways]
collapsed$signif_label <- ifelse(collapsed$padj < 0.05, "*", "")
write.csv(x = collapsed, file = paste0(outdir, "/fgsea_hallmark_collapsed_pathways_results.csv"), row.names = FALSE)

# Barplot of of enriched gene sets
hallmark_res$pathway <- gsub("HALLMARK_", "", hallmark_res$pathway)
gsea_plot <- ggplot(data = hallmark_res, aes(x = reorder(pathway, NES), y = NES)) +
  geom_bar(aes(fill = signif_label), stat = "identity") +
  scale_fill_manual(breaks = c("Adjusted p-value < 0.05", "Not significant"), values = c("darkorchid1", "turquoise3"), na.value = NA) +
  scale_x_discrete(position = "top") +
  coord_flip() +
  labs(x = "Hallmark Gene Sets", y = "Magnitude of gene-set level changes (Normalized Enrichment Score)", fill = "Significance",
       title = paste0(cell_line_name, ": ", drug_name, " treated vs. control")) +
  theme_light() +
  theme(legend.position = "top")
gsea_plot

ggsave(plot = gsea_plot, filename = "gsea.png", path = outdir, dpi = 300, device = "png", units = "in", height = 13, width = 8)