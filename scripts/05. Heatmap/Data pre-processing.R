# Data pre-processing for heatmap
# Author: Muntasim Fuad

# Load required packages
library(tidyverse)
library(biomaRt)

conflicted::conflicts_prefer(dplyr::select)

# Fuction to clean the data
clean_fpkm_data <- function(fpkm_df) {
  # Function to clean a single column
  clean_column <- function(col) {
    # Extract the numeric value (everything after last space)
    as.numeric(sub("^.*\\s", "", col))
  }
  
  # Clean gene IDs - remove "gene:" prefix and anything after the ID
  gene_ids <- sub("^gene:(ENSOMYG\\d+).*", "\\1", fpkm_df[[1]])
  
  # Clean all columns except the first
  cleaned_values <- as.data.frame(lapply(fpkm_df[-1], clean_column))
  
  # Simplify column names by removing "_FPKM"
  colnames(cleaned_values) <- sub("_FPKM$", "", colnames(cleaned_values))
  
  # Create new data frame with cleaned Gene_ID column
  cleaned_df <- data.frame(
    Gene_ID = gene_ids,
    cleaned_values,
    stringsAsFactors = FALSE,
    check.names = TRUE  # Ensures valid column names
  )
  
  # Remove duplicate rows (keeping first occurrence)
  cleaned_df <- cleaned_df[!duplicated(cleaned_df$Gene_ID), ]
  
  # Remove rows with NA gene IDs
  cleaned_df <- cleaned_df[complete.cases(cleaned_df$Gene_ID), ]
  
  return(cleaned_df)
}

# Load and clean all datasets
PRJNA559610 <- read.delim("data/FPKM_StringTie/PRJNA559610_FPKM.txt")
PRJNA1000995 <- read.delim("data/FPKM_StringTie/PRJNA1000995_FPKM.txt")
PRJNA1064938 <- read.delim("data/FPKM_StringTie/PRJNA1064938_FPKM.txt")
PRJNA1092638 <- read.delim("data/FPKM_StringTie/PRJNA1092638_FPKM.txt")

PRJNA559610_clean <- clean_fpkm_data(PRJNA559610)
PRJNA1000995_clean <- clean_fpkm_data(PRJNA1000995)
PRJNA1064938_clean <- clean_fpkm_data(PRJNA1064938)
PRJNA1092638_clean <- clean_fpkm_data(PRJNA1092638)


# Combine all heat stress datasets
heat_stress <- PRJNA559610_clean %>%
  full_join(PRJNA1092638_clean, by = "Gene_ID") %>%
  mutate(across(where(is.numeric), ~ coalesce(., 0)))

# Set up biomaRt for annotation
mart <- useMart("ensembl", dataset = "omykiss_gene_ensembl")

# Fetch gene annotations
annotations <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "description"),
  filters = "ensembl_gene_id",
  values = heat_stress$Gene_ID,
  mart = mart
)

# Merge annotations with the combined results
heat_stress <- merge(heat_stress, 
                   annotations, by.x = "Gene_ID", 
                   by.y = "ensembl_gene_id", 
                   all.x = TRUE)
  
# ───────────────────────────────────────────────────────────
# Filter to genes with ≥ 5 non-zero samples across ALL numeric columns
# ───────────────────────────────────────────────────────────
heat_stress <- heat_stress %>%
  # count non-zero entries per row among numeric columns
  filter(
    rowSums(
      select(., where(is.numeric)) != 0,
      na.rm = TRUE
    ) >= 5
  )

# Relocate and rename
heat_stress <- heat_stress |>
  relocate(Gene_ID, external_gene_name, description) |>
  rename(Gene_Symbol = external_gene_name) |>
  select(-description)

# ------------------------------------------------------------------------------
# Combine all hypoxia datasets
hypoxia <- PRJNA1064938_clean %>%
  full_join(PRJNA1000995_clean, by = "Gene_ID") %>%
  mutate(across(where(is.numeric), ~ coalesce(., 0)))

# Fetch gene annotations
annotations <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "description"),
  filters = "ensembl_gene_id",
  values = hypoxia$Gene_ID,
  mart = mart
)

# Merge annotations with the combined results
hypoxia <- merge(hypoxia, 
                  annotations, by.x = "Gene_ID", 
                  by.y = "ensembl_gene_id", 
                  all.x = TRUE)

# ───────────────────────────────────────────────────────────
# Filter to genes with ≥ 5 non-zero samples across ALL numeric columns
# ───────────────────────────────────────────────────────────
hypoxia <- hypoxia %>%
  # count non-zero entries per row among numeric columns
  filter(
    rowSums(
      select(., where(is.numeric)) != 0,
      na.rm = TRUE
    ) >= 10
  )

# Relocate and rename
hypoxia <- hypoxia |>
  relocate(Gene_ID, external_gene_name, description) |>
  rename(Gene_Symbol = external_gene_name) |>
  select(-description)


# Save the cleaned data
dir.create("outputs/Heatmap", showWarnings = FALSE)
write.csv(heat_stress, "outputs/Heatmap/heatstress_cleaned_fpkm.csv",
          row.names = FALSE)
write.csv(hypoxia, "outputs/Heatmap/hypoxia_cleaned_fpkm.csv", 
          row.names = FALSE)
