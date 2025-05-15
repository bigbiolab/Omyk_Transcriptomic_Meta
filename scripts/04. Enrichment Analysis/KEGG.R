# KEGG Enrichment Analysis
# Author: Muntasim Fuad

# Load required packages
library(tidyverse)
library(clusterProfiler)
library(biomaRt)
library(hrbrthemes)
library(openxlsx)

conflicted:: conflict_prefer("select", "dplyr")
conflicted::conflicts_prefer(clusterProfiler::filter)

set.seed(123)

# Load results
meta_result <- read.xlsx("outputs/Meta-Analysis/Random Effect Model.xlsx")


# Set up biomaRt for annotation
mart <- useMart("ensembl", dataset = "omykiss_gene_ensembl")

# Fetch gene annotations
entrez_geneList<- getBM(
  attributes = c("entrezgene_id", "ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = meta_result$Gene_ID,
  mart = mart
)

# KEGG enrichment analysis
KEGG <- enrichKEGG(
  entrez_geneList$entrezgene_id,
  organism = "omy",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = NULL,
  minGSSize = 10,
  maxGSSize = 1000,
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
)

# Extract KEGG results
kegg_results <- KEGG@result

# Extract top 30 KEGG results 
kegg_data <- kegg_results |> 
  filter(p.adjust < 0.05) |> 
  arrange(desc(RichFactor)) |> 
  head(22) 

#. Generate color gradients for CC
pal <- c("#bdd7e7","#6baed6", "#3182bd", "#08519c")

# Generate dot plot
kegg_dot <- ggplot(kegg_data, aes(x = RichFactor, 
                                  y = reorder(Description, RichFactor),
                                  size = Count, 
                                  color = p.adjust)) +
  geom_point() +
  labs(title = "KEGG Pathway Enrichment",
       x = "Rich Factor",
       y = "",
       color = "Adjusted P-value",
       size = "Gene Count") +
  theme_ipsum_pub(base_size = 30, 
                  grid_col = "#bdbdbd",
                  grid = "X", 
                  plot_title_margin = 15,
                  axis_title_size = 20,
                  axis_title_just = "c") +
  theme(legend.position = "right",
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 35),
        legend.key.height = unit(1.2, "cm"),
        legend.key.width = unit(0.8, "cm"), 
        panel.border = element_rect(colour = "black", fill = NA,
                                    linewidth = 2),
        axis.text.y = element_text(colour = "black", 
                                   hjust = 1, size = 35),
        axis.text.x = element_text(size = 35),
        axis.title.x = element_text(size = 35),
        plot.title = element_text(size = 40)) +
  scale_color_gradientn(
    colors = pal,
    limits = c(0, 0.05),
    breaks = c(0.04, 0.03, 0.02, 0.01),
    labels = c("0.04", "0.03", "0.02", "0.01"),
    values = scales::rescale(c(0.04, 0.03, 0.02, 0.01))
  ) +
  scale_size_continuous(range = c(5, 15))

# Save the dot plot
ggsave("figures/Enrichment Analysis/KEGG_dotplot.png", 
       plot = kegg_dot, 
       width = 28, 
       height = 20, 
       dpi = 600,
       units = "in",
       bg = "white")

# Save the KEGG results to an Excel file
write.xlsx(kegg_data,
           file = "outputs/Enrichment Analysis/KEGG.xlsx")

# ------------------------------------------------------------------------------
browseKEGG(KEGG, "omy04218")
browseKEGG(KEGG, "omy04141")
browseKEGG(KEGG, "omy04136")
