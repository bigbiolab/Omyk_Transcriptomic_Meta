# Gene Ontology Enrichment Analysis using g:Profiler
# Author: Muntasim Fuad

# Load required packages
library(gprofiler2)
library(tidyverse)
library(openxlsx)

# Load meta-analysis results 
meta_results <- read.xlsx("outputs/Meta-Analysis/Filtered_Meta_DEGs.xlsx")

# Run gprofiler2
gostres <- gost(query = meta_results$Gene_ID,
                organism = "omykiss",
                ordered_query = FALSE,
                multi_query = FALSE,
                significant = TRUE,
                exclude_iea = FALSE,
                measure_underrepresentation = FALSE,
                evcodes = TRUE,
                user_threshold = 0.05,
                correction_method = "fdr",
                domain_scope = c("annotated", "known",
                                 "custom", "custom_annotated"),
                custom_bg = NULL,
                numeric_ns = "",
                sources = NULL,
                as_short_link = FALSE,
                highlight = TRUE)

# ------------------------------------------------------------------------------
# Run, if as_short_link = FALSE
# --------------------------------

result <- gostres$result

result[] <- lapply(result, function(column) {
  if (is.list(column)) {
    return(sapply(column, toString))  # Convert list elements to strings
  } else {
    return(column)  # Keep other columns as is
  }
})

write.csv(result , "outputs/Enrichment Analysis/GO/gProfiler_omykiss_GO.csv",
          row.names = FALSE)
