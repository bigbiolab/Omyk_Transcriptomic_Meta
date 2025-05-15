# KEGG Enrichment Analysis
# Author: Muntasim Fuad

# Load required packages
library(tidyverse)
library(clusterProfiler)
library(biomaRt)
library(hrbrthemes)
library(openxlsx)

conflicted:: conflict_prefer("select", "dplyr")
conflicted::conflicts_prefer(clusterProfiler::filter)

set.seed(123)

# Load results
meta_result <- read.csv("outputs/Meta-Analysis/Meta_combined_heat_studies.csv")


# Set up biomaRt for annotation
mart <- useMart("ensembl", dataset = "omykiss_gene_ensembl")

# Fetch gene annotations
entrez_geneList<- getBM(
  attributes = c("entrezgene_id", "ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = meta_result$Gene_ID,
  mart = mart
)

# KEGG enrichment analysis
KEGG <- enrichKEGG(
  entrez_geneList$entrezgene_id,
  organism = "omy",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = NULL,
  minGSSize = 10,
  maxGSSize = 1000,
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
)

# Extract KEGG results
kegg_results <- KEGG@result

# Extract top 30 KEGG results 
kegg_data <- kegg_results |> 
  filter(p.adjust < 0.05) |> 
  arrange(desc(RichFactor))


# Save the KEGG results to an Excel file
write.xlsx(kegg_data,
           file = "outputs/Enrichment Analysis/KEGG_heat_stress.xlsx")
