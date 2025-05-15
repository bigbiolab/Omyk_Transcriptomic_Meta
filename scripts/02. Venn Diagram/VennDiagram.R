# Common differentially expressed genes across stress related studies
# Author: Muntasim Fuad

# 00. Load required packages
library(tidyverse)
library(VennDiagram)
library(ggVennDiagram)
library(openxlsx)

# 01. Load datasets 
PRJNA559610 <- read.csv("outputs/DESeq2/PRJNA559610.csv")
PRJNA1000995 <- read.csv("outputs/DESeq2/PRJNA1000995.csv")
PRJNA1064938 <- read.csv("outputs/DESeq2/PRJNA1064938.csv")
PRJNA1092638 <- read.csv("outputs/DESeq2/PRJNA1092638.csv")

# 02. Filter significantly expressed genes
PRJNA559610 <- PRJNA559610 |> filter(padj < 0.05 & abs(log2FoldChange) > 1)
PRJNA1000995 <- PRJNA1000995 |> filter(padj < 0.05 & abs(log2FoldChange) > 1)
PRJNA1064938 <- PRJNA1064938 |> filter(padj < 0.05 & abs(log2FoldChange) > 1)
PRJNA1092638 <- PRJNA1092638 |> filter(padj < 0.05 & abs(log2FoldChange) > 1)

# 03. Create a list of gene sets
gene_sets <- list(
  PRJNA559610 = PRJNA559610$Gene_ID,
  PRJNA1000995 = PRJNA1000995$Gene_ID,
  PRJNA1064938 = PRJNA1064938$Gene_ID,
  PRJNA1092638 = PRJNA1092638$Gene_ID
)


# Create Venn diagram
png("figures/Venn Diagram/Venn_diagram.png", units = "in", width = 6, 
    height = 6, res = 2000)

venn.plot <- venn.diagram(
  x = gene_sets,
  filename = NULL,
  fill = c("#99d8c9", "#addd8e", "#bcbddc", "#fec44f"),
  alpha = 0.6,
  cex = 0.8,
  cat.cex = 1,
  cat.col = c("#2ca25f", "#31a354", "#756bb1", "#e6550d"),
  cat.pos = c(-20, 20, -20, 20),
  cat.dist = c(0.22, 0.22, 0.1, 0.09),
  margin = 0.1,
  lwd = 2,
  disable.logging = TRUE
)

grid.draw(venn.plot)
dev.off()

# 06. Retrieve intersection of three gene sets
intersect <- process_region_data(Venn(gene_sets))

# 07. Extract genes that are common across all studies
intersect_genes <- intersect |> 
  filter(id == "1/2/3/4") |> 
  pull(item) |>  as.data.frame()

colnames(intersect_genes) <- "Gene_ID"

write.csv(intersect_genes, "outputs/Venn Diagram/intersect_genes.csv", row.names = FALSE)
