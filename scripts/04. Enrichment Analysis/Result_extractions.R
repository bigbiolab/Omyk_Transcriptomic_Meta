# Table for KEGG and GO enrichment results
# Author: Muntasim Fuad

# Load required packages
library(tidyverse)
library(openxlsx)
library(officer)

# ------------------------------------------------------------------------------
# Export KEGG & GO  enrichment results
# Load enrichment results
gprofiler <- read.csv("outputs/Enrichment Analysis/GO/gProfiler_omykiss_GO.csv")
kegg <- read.xlsx("outputs/Enrichment Analysis/KEGG.xlsx")

kegg$geneID <- gsub("/", ", ", kegg$geneID)


# Filter results by category (BP,MF, CC)
bp_data <- gprofiler |> filter(source == "GO:BP")
mf_data <- gprofiler |> filter(source == "GO:MF")
cc_data <- gprofiler |> filter(source == "GO:CC")

# Filter out top 20 results
top_kegg <- kegg |> arrange(desc(Count)) |> select(Description, 
                                                   `KEGG ID` = ID,GeneRatio,
                                                   RichFactor,
                                                   `Adjust P-value` = p.adjust,
                                                   Count)

top_bp <- bp_data |> arrange(desc(intersection_size)) |> 
  slice_head(n = 10) |> 
  select(`GO Source` = source,
         `GO Term` = term_name,
         `GO ID` = term_id, 
         `Term Size` = `term_size`,
         `Query Size` = `query_size`,
         `Adjust P-value` = p_value,
         Count = intersection_size)

top_mf <- mf_data |> arrange(desc(intersection_size)) |> 
  slice_head(n = 10)|>
  select(`GO Source` = source,
         `GO Term` = term_name,
         `GO ID` = term_id, 
         `Term Size` = `term_size`,
         `Query Size` = `query_size`,
         `Adjust P-value` = p_value,
         Count = intersection_size)

top_cc <- cc_data |> arrange(desc(intersection_size)) |> 
  slice_head(n = 10)|>
  select(`GO Source` = source,
         `GO Term` = term_name,
         `GO ID` = term_id, 
         `Term Size` = `term_size`,
         `Query Size` = `query_size`,
         `Adjust P-value` = p_value,
         Count = intersection_size)

# Combined all GO terms
combined_go <- rbind(top_bp, top_mf, top_cc) |> ungroup()


# Export the results

# Create a Word document and export results
doc <- read_docx()

# Add KEGG results
doc <- doc %>%
  body_add_par("KEGG Enrichment Analysis Results", style = "heading 1") %>%
  body_add_table(value = top_kegg, style = "table_template")

# Add GO Biological Process (BP) results
doc <- doc %>%
  body_add_par("Gene Ontology Enrichment Analysis Results", style = "heading 1") %>%
  body_add_table(value = combined_go, style = "table_template")

# Save the Word document
print(doc, target = "tables/Enrichment_Results.docx")
