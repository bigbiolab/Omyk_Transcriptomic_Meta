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
heat_data <- read.csv("outputs/Heatmap/heatstress_cleaned_fpkm.csv",
                         na.strings ="")
glimpse(heat_data)

# Import key genes
heat_genes <- read.csv("outputs/Meta-Analysis/Meta_combined_heat_studies.csv",
                          na.strings = "")

key_genes <- heat_genes |> drop_na(external_gene_name) |> 
  filter(abs(metafc) >= 4.5) |> select(Gene_ID,external_gene_name) |> 
  rename(Gene.ID = Gene_ID,
         Gene.Symbol = external_gene_name)

intersect_genes <- read.csv("outputs/Meta-Analysis/Meta_combining-approach_mean_intersect_genes.csv")

intersect_genes <-  intersect_genes |> select(Gene.ID,
                                              Gene.Symbol)

key_genes <- rbind(key_genes, intersect_genes) |> distinct()


# Filter the combined data to include only key genes
combined_data <- heat_data |>
  filter(Gene_ID %in% key_genes$Gene.ID)

# Export the combined data
write.csv(combined_data, "outputs/Heatmap/heat_key_genes.csv",
          row.names = FALSE)

# ------------------------------------------------------------------------------
# Load data
combined_data <- read.csv("outputs/Heatmap/heat_key_genes.csv")

combined_data <- combined_data |>
  group_by(Gene_Symbol) |>
  summarise(across(where(is.numeric), mean)) |>
  column_to_rownames("Gene_Symbol")

# Load metadata
PRJNA559610_meta <- read.csv("data/Raw Counts/metadata/PRJNA559610.csv")
PRJNA1092638_meta <- read.csv("data/Raw Counts/metadata/PRJNA1092638.csv")

# Combine all the metadata
metadata_combined <- rbind(
  PRJNA559610_meta,
  PRJNA1092638_meta
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
                            levels = c("control", "heat stress"),
                            ordered = TRUE)) %>% 
  arrange(Conditions)

# Reorder the matrix columns to match the annotation order
combined_matrix <- combined_matrix[, rownames(sample_annot)]

# Define colors for each condition
condition_colors <- c("control" = "#a1dab4", 
                      "heat stress" = "#fb6a4a")

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
  c(c(0, 2, 4, 6, 8, 10, 12)),
  c("#0868ac", "#9ecae1",
    "#ffffb2", "#fee391",
    "#fe9929", "#f03b20",
    "#bd0026"
    )
  )


png("figures/Heatmaps/Heat Stress.png", width = 7000, height = 4800, res = 600)

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
    at = c(0, 2, 4, 6, 8, 10, 12),  # Explicitly set breakpoints
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 10)
  )
)

# Draw the heatmap
draw(ht, heatmap_legend_side = "right")

dev.off()

