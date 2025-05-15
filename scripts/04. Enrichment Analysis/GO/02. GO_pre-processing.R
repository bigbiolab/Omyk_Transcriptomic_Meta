# Data pre-processing GO enrichment analysis using g:Profiler
# Author: Muntasim Fuad

# Load required packages
library(tidyverse)
library(ggh4x)
library(reshape2)
library(stringr)
library(officer)
library(openxlsx)

# Load results from g:Profiler
gprofiler <- read.csv("outputs/Enrichment Analysis/GO/gProfiler_omykiss_GO.csv")

# Load meta-analysis results
meta_results <- read.xlsx("outputs/Meta-Analysis/Filtered_Meta_DEGs.xlsx")

# Filter results by category (BP,MF, CC)
bp <- gprofiler |> filter(source == "GO:BP")
mf <- gprofiler |> filter(source == "GO:MF")
cc <- gprofiler |> filter(source == "GO:CC")

# Filter out necessary columns
bp <- bp |> select(term_name, term_id, intersection_size, intersection)
mf <- mf |> select(term_name, term_id, intersection_size, intersection)
cc <- cc |> select(term_name, term_id, intersection_size, intersection) 

# Biological Process (BP)
# Split intersections column and unnest into long format
formatted_bp <- bp |> 
  mutate(Gene_ID = strsplit(as.character(intersection), ",")) |> 
  unnest(Gene_ID ) |> 
  select(Gene_ID , GO_ID = term_id, GO_Term = term_name)

# Merge formatted data with meta-analysis results
merged_bp <- formatted_bp |> 
  inner_join(meta_results, by = "Gene_ID",
             relationship = "many-to-many") |> 
  select(Gene_ID , GO_ID, GO_Term, log2FoldChange = randomSummary)

# Define genes based on regulation
merged_bp <- merged_bp |>
  mutate(Regulation = ifelse(log2FoldChange > 0, "UP", "DOWN"))

# Group by biological process and regulation, then count the occurrences
bp_counts <- merged_bp |> 
  group_by(GO_Term, Regulation) |> 
  summarize(Count = n(), .groups = 'drop')

# Export BP count data
write.xlsx(bp_counts,"outputs/Enrichment Analysis/GO/BP_counts.xlsx")

# Molecular Function (MF)
# Split intersections column and unnest into long format
formatted_mf <- mf |> 
  mutate(Gene_ID = strsplit(as.character(intersection), ",")) |> 
  unnest(Gene_ID ) |> 
  select(Gene_ID , GO_ID = term_id, GO_Term = term_name)

# Merge formatted data with meta-analysis results
merged_mf <- formatted_mf |> 
  inner_join(meta_results, by = "Gene_ID",
             relationship = "many-to-many") |> 
  select(Gene_ID, GO_ID, GO_Term, log2FoldChange = randomSummary)

# Define genes based on regulation
merged_mf <- merged_mf |>
  mutate(Regulation = ifelse(log2FoldChange > 0, "UP", "DOWN"))

# Group by molecular function and regulation, then count the occurrences
mf_counts <- merged_mf |> 
  group_by(GO_Term, Regulation) |> 
  summarize(Count = n(), .groups = 'drop')

# Export MF count data
write.xlsx(mf_counts,"outputs/Enrichment Analysis/GO/MF_counts.xlsx")

# Cellular Component (CC)
# Split intersections column and unnest into long format
formatted_cc <- cc |>
  mutate(Gene_ID  = strsplit(as.character(intersection), ",")) |> 
  unnest(Gene_ID) |> 
  select(Gene_ID, GO_ID = term_id, GO_Term = term_name)

# Merge formatted data with meta-analysis results
merged_cc <- formatted_cc |> 
  inner_join(meta_results, by = "Gene_ID",
             relationship = "many-to-many") |> 
  select(Gene_ID, GO_ID, GO_Term, log2FoldChange = randomSummary)

# Define genes based on regulation
merged_cc <- merged_cc |>
  mutate(Regulation = ifelse(log2FoldChange > 0, "UP", "DOWN"))

# Group by cellular component and regulation, then count the occurrences
cc_counts <- merged_cc |> 
  group_by(GO_Term, Regulation) |> 
  summarize(Count = n(), .groups = 'drop')

# Export CC count data
write.xlsx(cc_counts,"outputs/Enrichment Analysis/GO/CC_counts.xlsx")
