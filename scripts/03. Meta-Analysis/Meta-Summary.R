# Retrieve summary of meta-analysis result
# Author: Muntasim Fuad

# Load necessary libraries
library(tidyverse)
library(openxlsx)

# Load data
meta_result <- read.xlsx("outputs/Meta-Analysis/Random Effect Model.xlsx")

# Rename columns 
meta_result <- meta_result |> 
  select(Gene_Symbol, randomSummary, randomP, Gene_Description) |> 
  rename(log2FC = randomSummary, P.Value = randomP)

# Define significance thresholds
meta_result <- meta_result %>%
  mutate(Significance = case_when(
    P.Value < 0.05 & log2FC > 1 ~ "Up",
    P.Value < 0.05 & log2FC < -1 ~ "Down",
    TRUE ~ "NoSignificant"
  ))

# Count the number of up-regulated and down-regulated genes
gene_stats <- meta_result %>% 
  filter(Significance != "NoSignificant") %>% 
  group_by(Significance) %>% 
  summarise(
    Count = n(),
    Percentage = n() / nrow(meta_result) * 100
  )

gene_stats

# Extract annotated differentially expressed genes 
annotated_genes <- meta_result |> mutate(Gene_Symbol = na_if(Gene_Symbol, "")) |> 
  filter(!is.na(Gene_Symbol)) |> filter(!Significance == "NoSignificant")

# Export results
write.xlsx(annotated_genes, 
           "outputs/Meta-Analysis/Filtered_Meta_DEGs_annotated_only.xlsx")

# ------------------------------------------------------------------------------
# Load filtered results
filtered_results <- read.xlsx("outputs/Meta-Analysis/Filtered_Meta_DEGs.xlsx") |> 
                                select(c(1,13,14))

stress_genes <- filtered_results[grepl("heat", filtered_results$Gene_Description), ]

immune_genes <- filtered_results[grepl("immune", filtered_results$Gene_Description), ]

hypoxia_genes <- filtered_results[grepl("hypoxia", filtered_results$Gene_Description), ]

damage_genes <- filtered_results[grepl("damage", filtered_results$Gene_Description), ] 

# ------------------------------------------------------------------------------
# Load filtered results
meta_results <- read.xlsx("outputs/Meta-Analysis/Filtered_Meta_DEGs.xlsx") |> 
  select(c(1,13,14))

# Import GeneCards search results
heat_stress <- read.csv("data/GeneCards/GeneCards-SearchResults_heat_stress.csv")
hypoxia <- read.csv("data/GeneCards/GeneCards-SearchResults_hypoxia.csv")
house_keeping <- read.csv("data/GeneCards/GeneCards-SearchResults_house_keeping_genes.csv")
# Merge all datasets
genecards <-
  rbind(heat_stress[, 1:2]) %>%
  rbind(hypoxia[, 1:2]) |> distinct()

key_genes <- meta_results |> filter(toupper(Gene_Symbol) %in% genecards$Gene.Symbol)

# Exclude house keeping genes
key_genes <- key_genes |> 
  filter(!toupper(Gene_Symbol) %in% house_keeping$Gene.Symbol)

key_genes <- rbind(key_genes,hypoxia_genes,damage_genes,
                   immune_genes, stress_genes) |> distinct()

# Export results
write.xlsx(key_genes, 
           "outputs/Meta-Analysis/Meta_DEGs_key_genes.xlsx")
# ------------------------------------------------------------------------------
# Add regulation in the meta_analysis result
meta_key_results <-  read.xlsx("outputs/Meta-Analysis/Filtered_Meta_DEGs.xlsx")

# Define genes based on regulation
meta_key_results <- meta_key_results |>
  select(Gene_ID,Gene_Symbol,Gene_Description, log2FoldChange = randomSummary,
         adj.P = randomP) |> 
  mutate(Regulation = ifelse(log2FoldChange > 0, "UP", "DOWN"))

# Export results
write.xlsx(meta_key_results, 
           "outputs/Meta-Analysis/Meta_DEGs_regulation.xlsx")

# ------------------------------------------------------------------------------
# Load filtered results
meta_results <- read.xlsx("outputs/Meta-Analysis/Filtered_Meta_DEGs.xlsx",
                          na.strings = "") |> 
  select(c(1,3,6,13,14))

meta_results <- meta_results |> 
  rename(log2FoldChange = randomSummary, 
         adjusted.P.Value = randomP) |> 
  relocate(Gene_ID, Gene_Symbol, Gene_Description, log2FoldChange, adjusted.P.Value) |> 
  drop_na(Gene_Symbol) |> 
  mutate(Regulation = ifelse(log2FoldChange > 0, "UP", "DOWN"))


# Import GeneCards search results
heat_stress <- read.csv("data/GeneCards/GeneCards-SearchResults_heat_stress.csv")
hypoxia <- read.csv("data/GeneCards/GeneCards-SearchResults_hypoxia.csv")

# Merge all datasets
heat_stress_meta_analysis <- meta_results |> 
  filter(toupper(Gene_Symbol) %in% heat_stress$Gene.Symbol)

hypoxia_meta_analysis <- meta_results |> 
  filter(toupper(Gene_Symbol) %in% hypoxia$Gene.Symbol)

# Export results
write.xlsx(heat_stress_meta_analysis, 
           "outputs/Meta-Analysis/Meta_DEGs_heat_stress_genes.xlsx")

write.xlsx(hypoxia_meta_analysis,
           "outputs/Meta-Analysis/Meta_DEGs_hypoxia_genes.xlsx")

# ------------------------------------------------------------------------------
# Load Venn Diagram results
intersect <- read.csv("outputs/Venn Diagram/intersect_genes.csv")
heat_stress <- read.csv("outputs/Venn Diagram/intersect_heat_genes.csv")
hypoxia <- read.csv("outputs/Venn Diagram/intersect_genes_hypoxia.csv")

# Load all DGE results
PRJNA559610 <- read.csv("outputs/DESeq2/PRJNA559610.csv")
PRJNA1000995 <- read.csv("outputs/DESeq2/PRJNA1000995.csv")
PRJNA1064938 <- read.csv("outputs/DESeq2/PRJNA1064938.csv")
PRJNA1092638 <- read.csv("outputs/DESeq2/PRJNA1092638.csv")


combined_results <- rbind(PRJNA559610, PRJNA1000995, 
                          PRJNA1064938, PRJNA1092638) |>  
  mutate(Regulation = ifelse(log2FoldChange > 0, "UP", "DOWN"))

# Merge all datasets
heat_stress_intersect <- combined_results |> 
  filter(Gene_ID %in% heat_stress$Gene_ID)

hypoxia_intersect <- combined_results |> 
  filter(Gene_ID %in% hypoxia$Gene_ID)

# Load meta results
combined_meta <- read.xlsx("outputs/Meta-Analysis/Meta_combining-approach_mean.xlsx")

intersects <- combined_meta |> 
  filter(Gene_ID %in% intersect$Gene_ID)

# Export results
write.csv(intersects,
          "outputs/Venn Diagram/intersect_genes_meta_analysis.csv",
          row.names = FALSE)
write.csv(heat_stress_intersect,
            "outputs/Venn Diagram/intersect_heat_genes_meta_analysis.csv",
          row.names = FALSE)
write.csv(hypoxia_intersect,
          "outputs/Venn Diagram/intersect_genes_hypoxia_meta_analysis.csv",
          row.names = FALSE)

            