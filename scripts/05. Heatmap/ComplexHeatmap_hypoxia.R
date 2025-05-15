# Data Pre-processing for heatmap
# Author: Muntasim Fuad

# Load required packages
library(tidyverse)
library(fuzzyjoin)
library(ComplexHeatmap)
library(circlize)
library(openxlsx)

conflicted::conflicts_prefer(dplyr::filter)

# Load data for heatmap
hypoxia_data <- read.csv("outputs/Heatmap/hypoxia_cleaned_fpkm.csv",
                          na.strings ="")
glimpse(hypoxia_data)

# Import key genes
hypoxia_genes <- read.csv("outputs/Meta-Analysis/Meta_combined_hypoxia_studies.csv",
                      na.strings = "")

key_genes <- hypoxia_genes |> drop_na(external_gene_name) |> 
  filter(abs(metafc) >= 1.4) |> select(Gene_ID,external_gene_name) |> 
  rename(Gene.ID = Gene_ID,
         Gene.Symbol = external_gene_name)

intersect_genes <- read.csv("outputs/Meta-Analysis/Meta_combining-approach_mean_intersect_genes.csv")

intersect_genes <-  intersect_genes |> select(Gene.ID,
                                              Gene.Symbol)

key_genes <- rbind(key_genes, intersect_genes) |> distinct()


# Filter the combined data to include only key genes
combined_data <- hypoxia_data |>
  filter(Gene_ID %in% key_genes$Gene.ID)

# Export the combined data
write.csv(combined_data, "outputs/Heatmap/hypoxia_key_genes.csv",
          row.names = FALSE)

# ------------------------------------------------------------------------------
# Load data
combined_data <- read.csv("outputs/Heatmap/hypoxia_key_genes.csv")

combined_data <- combined_data |>
  group_by(Gene_Symbol) |>
  summarise(across(where(is.numeric), mean)) |>
  column_to_rownames("Gene_Symbol")

# Load metadata
PRJNA1000995_meta <- read.csv("data/Raw Counts/metadata/PRJNA1000995.csv")
PRJNA1064938_meta <- read.csv("data/Raw Counts/metadata/PRJNA1064938.csv")

# Combine all the metadata
metadata_combined <- rbind(
  PRJNA1000995_meta,
  PRJNA1064938_meta
)

# Apply log2 transformation
combined_data <- log2(combined_data + 1)

# Convert to matrix
combined_matrix <- as.matrix(combined_data)

# Prepare column annotations with specific order
sample_annot <- metadata_combined %>%
  filter(Samples %in% colnames(combined_matrix)) %>%
  column_to_rownames("Samples") |> 
  mutate(Condition = factor(Conditions, 
                            levels = c("control", "hypoxia"),
                            ordered = TRUE)) %>% 
  arrange(Conditions)

# Reorder the matrix columns to match the annotation order
combined_matrix <- combined_matrix[, rownames(sample_annot)]

# Define colors for each condition
condition_colors <- c("control" = "#a1dab4", 
                      "hypoxia" = "#9e9ac8")

# Create the annotation object
column_ha <- HeatmapAnnotation(
  Condition = sample_annot$Condition,
  col = list(Condition = condition_colors),
  annotation_name_gp = gpar(fontsize = 14),
  show_annotation_name = TRUE,
  annotation_legend_param = list(
    title = "Condition",
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 10)
  )
)

# Define the color mapping function
col_fun <- colorRamp2(
  c(c(0, 2, 4, 6, 8)),
  c("#0868ac", "#9ecae1",
    "#ffffb2", "#fee391",
    "#fe9929"))


png("figures/Heatmaps/hypoxia.png", width = 7000, height = 4800, res = 600)

ht <- Heatmap(
  combined_matrix,
  name = "Log2 (FPKM +1)",
  col = col_fun,
  top_annotation = column_ha,
  
  # Heatmap parameters
  row_names_gp = gpar(fontsize = 14),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  show_row_names = TRUE,
  
  # Clustering
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "complete",
  
  # Add black grid lines for each cell
  rect_gp = gpar(col = "black", lwd = 0.5),  # Adjust `lwd` for line thickness
  
  # Heatmap legend
  heatmap_legend_param = list(
    at = c(0, 2, 4, 6, 8),  # Explicitly set breakpoints
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 10)
  )
)

# Draw the heatmap
draw(ht, heatmap_legend_side = "right")

dev.off()

