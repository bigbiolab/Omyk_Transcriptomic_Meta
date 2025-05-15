# Meta-Analysis of RNA-seq data using MetaVolcanoR
# Author: Muntasim Fuad

# 00. Load required pacakges
library(MetaVolcanoR) 
library(tidyverse)
library(biomaRt)
library(openxlsx)

# Data pre-processing

# 01. Define datasets
bioprojects <- c("PRJNA559610",
                 "PRJNA1000995", 
                 "PRJNA1064938",
                 "PRJNA1092638")

# 02. Generate file paths
results <- paste0("outputs/DESeq2/", bioprojects, ".csv")

# 03. Read all files into a list and name them properly
studies <- lapply(results, read.csv)
names(studies) <- bioprojects

# ------------------------------------------------------------------------------
# 04. Meta-Analysis

# Random Effect Model
meta_degs_rem <- rem_mv(diffexp= studies,
                        pcriteria='padj',
                        foldchangecol= "log2FoldChange",
                        genenamecol= "Gene_ID",
                        geneidcol= NULL,
                        collaps= TRUE,
                        vcol= "lfcse", 
                        cvar=FALSE,
                        metathr= 0.01,
                        jobname= "MetaVolcano",
                        outputfolder= "figures/Meta-Analysis/", 
                        draw= 'PDF',
                        ncores=1)

meta_results <- meta_degs_rem@metaresult

# 12. Set up biomaRt for annotation
mart <- useMart("ensembl", dataset = "omykiss_gene_ensembl")

# 13. Fetch gene annotations
annotations <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "description"),
  filters = "ensembl_gene_id",
  values = meta_results$Gene_ID,
  mart = mart
)

# 14. Merge annotations directly with the combined results
annotated_results <- merge(meta_results, 
                           annotations, by.x = "Gene_ID", 
                           by.y = "ensembl_gene_id", 
                           all.x = TRUE)

annotated_results <- annotated_results |> 
  rename(Gene_Symbol = external_gene_name, 
         Gene_Description = description)

write.xlsx(annotated_results, "outputs/Meta-Analysis/Random Effect Model.xlsx")
# ------------------------------------------------------------------------------
# Load meta analysis results
annotated_results <- read.xlsx("outputs/Meta-Analysis/Random Effect Model.xlsx")

# 05. Filter out DEGs
# Filtering based on statistical significance
significant_genes <- annotated_results |> 
  filter(randomP < 0.05)

# Filtering based reliable effect size estimates
key_genes <- significant_genes |> 
  filter(abs(randomSummary) >= 1)

# Save filtered results
write.xlsx(key_genes,
           "outputs/Meta-Analysis/Filtered_Meta_DEGs.xlsx")

# ------------------------------------------------------------------------------
# Data pre-processing

# 01. Define datasets
bioprojects <- c("PRJNA559610",
                 "PRJNA1000995", 
                 "PRJNA1064938",
                 "PRJNA1092638")

# 02. Generate file paths
results <- paste0("outputs/DESeq2/", bioprojects, ".csv")

# 03. Read all files into a list and name them properly
studies <- lapply(results, read.csv)
names(studies) <- bioprojects

# Combining-approach: Mean
meta_degs_comb <- combining_mv(diffexp= studies,
                               pcriteria='padj',
                               foldchangecol= "log2FoldChange",
                               genenamecol= "Gene_ID",
                               metafc= 'Mean',
                               metathr= 0.01, 
                               collaps= TRUE,
                               jobname= "MetaVolcano",
                               outputfolder= "figures/Meta-Analysis/",
                               draw= 'PDF')

write.xlsx(meta_degs_comb@metaresult, "outputs/Meta-Analysis/Meta_combining-approach_mean.xlsx")

# Load Venn Diagram results
intersect <- read.csv("outputs/Venn Diagram/intersect_genes.csv")

combined_meta <- meta_degs_comb@metaresult

intersects <- combined_meta |> 
  filter(Gene_ID %in% intersect$Gene_ID)

# 12. Set up biomaRt for annotation
mart <- useMart("ensembl", dataset = "omykiss_gene_ensembl")

# 13. Fetch gene annotations
annotations <- getBM(
  attributes = c("ensembl_gene_id",
                 "external_gene_name",
                 "description"),
  filters = "ensembl_gene_id",
  values = intersects$Gene_ID,
  mart = mart
)

# 14. Merge annotations directly with the combined results
intersec_results <- merge(intersects, 
                           annotations, 
                          by.x = "Gene_ID", 
                          by.y = "ensembl_gene_id", 
                           all.x = TRUE)

intersec_results <- intersec_results |>
  rename(`Gene Symbol` = external_gene_name, 
         `Gene Description` = description)

write.csv(intersec_results,
          "outputs/Meta-Analysis/Meta_combining-approach_mean_intersect_genes.csv",
          row.names = FALSE)
# ------------------------------------------------------------------------------
# Data pre-processing

# 01. Define datasets
bioprojects <- c("PRJNA559610",
                 "PRJNA1092638")

# 02. Generate file paths
heat_results <- paste0("outputs/DESeq2/", bioprojects, ".csv")

# 03. Read all files into a list and name them properly
heat_studies <- lapply(heat_results , read.csv)
names(heat_studies) <- bioprojects

# Combining-approach: Mean
meta_degs_comb <- combining_mv(heat_studies,
                               pcriteria='padj',
                               foldchangecol= "log2FoldChange",
                               genenamecol= "Gene_ID",
                               metafc= 'Mean',
                               metathr= 0.01, 
                               collaps= TRUE,
                               jobname= "MetaVolcano",
                               outputfolder= "figures/Meta-Analysis/",
                               draw= 'PDF')

heat_meta <- meta_degs_comb@metaresult

# 12. Set up biomaRt for annotation
mart <- useMart("ensembl", dataset = "omykiss_gene_ensembl")

# 13. Fetch gene annotations
annotations <- getBM(
  attributes = c("ensembl_gene_id",
                 "external_gene_name",
                 "description"),
  filters = "ensembl_gene_id",
  values = heat_meta$Gene_ID,
  mart = mart
)

# 14. Merge annotations directly with the combined results
heat_results <- merge(heat_meta, 
                          annotations, 
                          by.x = "Gene_ID", 
                          by.y = "ensembl_gene_id", 
                          all.x = TRUE)

write.csv(heat_results,
          "outputs/Meta-Analysis/Meta_combined_heat_studies.csv",
          row.names = FALSE)
# ------------------------------------------------------------------------------
# Data pre-processing

# 01. Define datasets
bioprojects <- c("PRJNA1000995",
                 "PRJNA1064938")

# 02. Generate file paths
hypoxia <- paste0("outputs/DESeq2/", bioprojects, ".csv")

# 03. Read all files into a list and name them properly
hypoxia_studies <- lapply(hypoxia , read.csv)
names(hypoxia_studies) <- bioprojects

# Combining-approach: Mean
meta_degs_comb <- combining_mv(hypoxia_studies,
                               pcriteria='padj',
                               foldchangecol= "log2FoldChange",
                               genenamecol= "Gene_ID",
                               metafc= 'Mean',
                               metathr= 0.01, 
                               collaps= TRUE,
                               jobname= "MetaVolcano",
                               outputfolder= "figures/Meta-Analysis/",
                               draw= 'PDF')

hypoxia_meta <- meta_degs_comb@metaresult

# 12. Set up biomaRt for annotation
mart <- useMart("ensembl", dataset = "omykiss_gene_ensembl")

# 13. Fetch gene annotations
annotations <- getBM(
  attributes = c("ensembl_gene_id",
                 "external_gene_name",
                 "description"),
  filters = "ensembl_gene_id",
  values = hypoxia_meta$Gene_ID,
  mart = mart
)

# 14. Merge annotations directly with the combined results
hypoxia_results <- merge(hypoxia_meta, 
                      annotations, 
                      by.x = "Gene_ID", 
                      by.y = "ensembl_gene_id", 
                      all.x = TRUE)

write.csv(hypoxia_results,
          "outputs/Meta-Analysis/Meta_combined_hypoxia_studies.csv",
          row.names = FALSE)
  
