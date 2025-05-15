# 00. Load required pacakges
library(tidyverse)
library(biomaRt)
library(openxlsx)
library(conflicted)

meta_result <- read.xlsx("outputs/Meta-Analysis/Filtered_Meta_DEGs.xlsx") |> 
  select(-c(Gene_Symbol, Gene_Description))

# Set up biomaRt for annotation
mart <- useMart("ensembl", dataset = "omykiss_gene_ensembl")

# Fetch gene annotations
annotations <- getBM(
  attributes = c("ensembl_gene_id",
                 "external_gene_name",
                 "description",
                 "omykiss_paralog_ensembl_gene",
                 "omykiss_paralog_associated_gene_name",
                 "omykiss_paralog_subtype"),
  filters = "ensembl_gene_id",
  values = meta_result$Gene_ID,
  mart = mart
)


# 14. Merge annotations directly with the combined results
annotated_results <- merge(meta_result, 
                           annotations, by.x = "Gene_ID", 
                           by.y = "ensembl_gene_id", 
                           all.x = TRUE)

write.xlsx(annotated_results, "outputs/Meta-Analysis/Filtered_Meta_DEGs_homologs.xlsx")
