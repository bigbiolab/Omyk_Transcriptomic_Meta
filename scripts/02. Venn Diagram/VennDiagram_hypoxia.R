# Common differentially expressed genes across hypoxia studies
# Author: Muntasim Fuad

# 00. Load required packages
library(tidyverse)
library(ggVennDiagram)
library(openxlsx)

# 01. Load datasets 
PRJNA1064938 <- read.csv("outputs/DESeq2/PRJNA1064938.csv")
PRJNA1000995 <- read.csv("outputs/DESeq2/PRJNA1000995.csv")

# 02. Filter significantly expressed genes
PRJNA1064938 <- PRJNA1064938 |> filter(padj < 0.05 & abs(log2FoldChange) > 1)
PRJNA1000995 <- PRJNA1000995 |> filter(padj < 0.05 & abs(log2FoldChange) > 1)

# 02. Create a list of gene sets
gene_sets <- list(PRJNA1064938 = PRJNA1064938$Gene_ID,
                  PRJNA1000995 = PRJNA1000995$Gene_ID)


# Create Venn diagram
png("figures/Venn Diagram/Venn_diagram_hypoxia.png", units = "in", width = 6, 
    height = 6, res = 2000)

venn.plot <- venn.diagram(
  x = gene_sets,
  filename = NULL,
  fill = c("#bcbddc", "#addd8e"),
  alpha = 0.6,
  cex = 0.8,
  cat.cex = 1,
  cat.col = c("#756bb1", "#2ca25f"),
  cat.pos = c(0, 0),
  cat.dist = c(0.02 , 0.02),
  margin = 0.1,
  lwd = 2,
  scaled = FALSE,
  disable.logging = TRUE
)

grid.draw(venn.plot)
dev.off()

# 06. Retrieve intersection of three gene sets
intersect <- process_region_data(Venn(gene_sets))

# 07. Extract genes that are common across all studies
intersect_genes <- intersect |> 
  filter(id == "1/2") |> 
  pull(item) |>  as.data.frame()

colnames(intersect_genes) <- "Gene_ID"

write.csv(intersect_genes, "outputs/Venn Diagram/intersect_genes_hypoxia.csv", row.names = FALSE)
