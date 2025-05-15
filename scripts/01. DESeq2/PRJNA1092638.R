# Differential gene expression analysis
# Author: Muntasim Fuad

# BioProject: PRJNA1092638

# 00. Load required packages
library(DESeq2)
library(tidyverse)
library(ggrepel)
library(openxlsx)


# 01. Load raw counts
PRJNA704549 <- read.xlsx("data/Raw Counts/PRJNA1092638.xlsx")
glimpse(PRJNA704549)

# 02.Load metadata
metadata <- read.csv("data/Raw Counts/metadata/PRJNA1092638.csv")
glimpse(metadata)

# 03.Create a matrix &  add gene ids as row names
count_data <- PRJNA704549 |> 
  column_to_rownames("Geneid") |> 
  as.data.frame()

count_data <- count_data |>
  mutate_all(as.numeric) |> 
  as.matrix()

# 04. Match metadata with count data
metadata <- metadata |> 
  filter(Samples %in% colnames(count_data)) |> 
  arrange(match(Samples, colnames(count_data)))

# 05. Prepare Sample information
colData <- data.frame( condition = as.factor(metadata$Conditions), 
                       row.names = colnames(count_data))

# 06. Create DESeq2 data set object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = colData,
                              design = ~ condition)

# 07. filter any counts less than 10
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# 08. Run DESeq2
dds <- DESeq(dds)

# 09. Analyze the data set and compile result
res <- results(dds)
res$Gene_ID <- sub("gene:", "", rownames(res))

# 10. Export results to a CSV file
write.csv(res, "outputs/DESeq2/PRJNA1092638.csv", 
          row.names = FALSE)
# ------------------------------------------------------------------------------
# 18. Get Normalized counts
normalized_counts <- counts(dds, normalized=TRUE) |> as.data.frame()

normalized_counts$Gene_ID <- sub("gene:", "", rownames(normalized_counts))
rownames(normalized_counts) <- NULL

# Save normalized counts to a file
write.csv(normalized_counts, "outputs/DESeq2/Normalized Counts/PRJNA1092638_normalized_counts.csv",
          row.names = FALSE)
# ------------------------------------------------------------------------------
# Perform variance-stabilizing transformation (VST) for PCA
vsd <- vst(dds, blind = FALSE)

# Perform PCA on the transformed data
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)

# Calculate the percentage of variance explained by each PC
percent_var <- round(100 * attr(pca_data, "percentVar"))

# Define custom colors for each conditions
custom_colors <- c("control" = "#1c9099", 
                   "heat stress" = "#de2d26")

# Plot PCA with custom styling
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_text_repel(aes(label = name),
                  box.padding = 0.5, 
                  max.overlaps = Inf,
                  size = 3) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("PRJNA1092638") +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),  # Add border
    legend.position = "right",
    aspect.ratio = 1
  ) +
  guides(color = guide_legend(title = "condition"))

# Export the plot
ggsave(filename = "figures/PCA/PRJNA1092638.png",
       plot = pca_plot,
       width = 5,
       height = 5,
       units = "in",
       bg = "white",
       dpi = 600)
