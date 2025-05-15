# biological function enrichment analysis using g:Profiler
# Author: Muntasim Fuad

# Load required packages
library(tidyverse)
library(openxlsx)

# Load Gene Ontology data
bp <- read.xlsx("outputs/Enrichment Analysis/GO/BP_counts.xlsx")

# Add category
bp$Category <- "Biological Process"

# Filter out top results
# Number of pathways to show 
n <- 20
# - Biological Process
bp <- bp |> 
  group_by(GO_Term) |> 
  mutate(Total = sum(Count)) |> 
  ungroup()

top_terms_bp <- bp |> 
  distinct(GO_Term, Total) |> 
  arrange(desc(Total)) |> 
  slice_head(n = n) |> 
  pull(GO_Term)

bp <- bp |> 
  filter(GO_Term %in% top_terms_bp)

# Reorder GO_Term based on the sum of Count for each level
bp <- bp |> 
  group_by(GO_Term) |> 
  ungroup() |> 
  mutate(GO_Term = reorder(GO_Term, Total))

# Create the plot
go_bp <- ggplot(bp, aes(x = Count, y = GO_Term, fill = Regulation)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +
  facet_wrap(~Category, strip.position = "right") +
  theme(
    axis.text.x = element_text(colour = "black",
                               hjust = 1, vjust = 1,
                               size = 45),
    axis.title.x = element_text(size = 45),
    axis.text.y = element_text(colour = "black", 
                               hjust = 1,
                               size = 50),
    strip.background.y = element_rect(fill = "#bae4b3",colour = "black",
                                      linewidth = 2.5), 
    strip.text = element_text(size = 45, color = "black"),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 125),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 2.5),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.height = unit(2, "cm"),
    legend.key.width = unit(2, "cm"),
    legend.text = element_text(size = 45),
    legend.title = element_text(size = 45)
  ) +
  labs(x = "Number of genes", y = "") +
  scale_fill_manual(values = c("UP" = "#fe9929", "DOWN" = "#6a51a3"))


# Export the plot
ggsave(filename = "figures/Enrichment Analysis/GO/Biological Process.png",
       plot = go_bp,
       width = 40,
       height = 25,
       units = "in",
       dpi = 600,
       bg = "white")

