# Volcano plot visualizing meta-analysis results
# Author: Muntasim Fuad

# Load necessary libraries
library(tidyverse)
library(ggrepel)
library(openxlsx)
library(conflicted)

# Load data
meta_result <- read.csv("outputs/Meta-Analysis/Meta_combined_heat_studies.csv",
                        na.strings = "")

# Rename columns 
meta_result <- meta_result |> 
  rename(Gene_Symbol = external_gene_name,
    log2FC = metafc, P.Value = metap) |> drop_na()

# Define significance thresholds
meta_result <- meta_result %>%
  mutate(Significance = case_when(
    P.Value < 0.05 & log2FC > 1 ~ "Up",
    P.Value < 0.05 & log2FC < -1 ~ "Down",
    TRUE ~ "NoSignificant"
  ))

# Extract top differentially expressed genes for annotation
top_genes <- meta_result |> mutate(Gene_Symbol = na_if(Gene_Symbol, "")) |> 
  filter(!is.na(Gene_Symbol)) |> 
  filter(P.Value < 0.05 & abs(log2FC) > 1)|> 
  slice_max(order_by = abs(log2FC), n = 100)


# Plot the volcano plot with border, threshold lines, updated legend, and gene names
volcano <- ggplot(meta_result, aes(x = log2FC, y = -log10(P.Value), color = Significance)) +
  geom_point(alpha = 0.8,
             size = 2.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  scale_color_manual(values = c("#2c7fb8", "#636363", "#e34a33")) +
  theme_minimal() +
  labs(title = "",
       x = "log2 (Fold Change)",
       y = "-log10 (adjusted P-value)") +
  theme(legend.title = element_text(size = 15, face = "bold" ),
        legend.text = element_text(size = 14),
        legend.position = "right",
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(colour = "black", 
                                   hjust = 1, size = 12),
        axis.text.y = element_text(colour = "black", 
                                   hjust = 1, size = 12),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  geom_text_repel(data = top_genes, aes(label = Gene_Symbol),
                  colour = "black",
                  size = 3,
                  max.overlaps = 5,
                  direction = "both",
                  max.time = 5,
                  force = 5,
                  force_pull = 5,
                  point.padding = 0.5,
                  seed = 40
  ) +
  scale_x_continuous(limits = c(-8, 8), 
                     breaks = seq(-8, 8, by = 1))

# Save the volcano plot
ggsave(filename = ("figures/Meta-Analysis/volcano plot_heat.png"), 
       plot = volcano, 
       width = 12, 
       height = 10, 
       dpi = 600,
       bg = "white")

