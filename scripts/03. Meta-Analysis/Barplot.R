# 
# Author: Muntasim Fuad

# load required packages 

library(tidyverse)
library(ggpubr)
library(openxlsx)

# load data
data <- read.xlsx("outputs/Meta-Analysis/top 40 genes.xlsx")

# Create the plot
plot <- ggplot(data, aes(x = reorder(Gene_Symbol, -log2FoldChange), 
                 y = log2FoldChange, 
                 fill = Regulation)) + 
  geom_bar(stat = "identity", 
           width = 0.8) +
  scale_fill_manual(values = c("UP" = "#e34a33", "DOWN" = "#2c7fb8"))  +
labs(x = "", y = "log2(Fold Change)")+
  theme(axis.title.y = element_text(size = 18,
                                    face = "plain",
                                    colour = "black"),
        axis.text.x = element_text(angle = 45,
                                   size = 14,
                                   hjust = 1,
                                   vjust = 1,
                                   face = "bold",
                                   colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black",
                                   face = "bold"),
        axis.ticks.y = element_blank(),
        legend.position = "top",
        panel.border = element_rect(linewidth = 1,
                                    fill = NA),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.height = unit(0.8, "cm"),
        legend.key.width = unit(0.8, "cm"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)) +
  geom_hline(yintercept = 0, color = "black") +
  scale_y_continuous(limits = c(-10, 22), 
                     breaks = seq(-10, 20, by = 5))

# Save the plot
ggsave(filename = "figures//Meta-Analysis/Barplot.png",
       plot = plot,
       width = 16, 
       height = 8,
       dpi = 600,
       bg = "white")
