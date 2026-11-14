# 07. KEGG Pathway Enrichment Analysis

setwd("C:/Users/Admin/Desktop/Bioinformatics")

library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# Load mapped gene IDs
gene_mapping <- read.csv("results/gene_id_mapping.csv")

# Extract unique Entrez IDs
gene_ids <- unique(as.character(gene_mapping$ENTREZID))

# Perform KEGG pathway enrichment analysis
kegg_results <- enrichKEGG(
  gene = gene_ids,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20
)

# Convert results into a data frame
kegg_table <- as.data.frame(kegg_results)

# Display enrichment results
print(head(kegg_table))

# Save KEGG enrichment results
write.csv(
  kegg_table,
  "results/KEGG_enrichment_results.csv",
  row.names = FALSE
)

# Create KEGG pathway visualization
if (nrow(kegg_table) > 0) {
  
  kegg_plot <- dotplot(
    kegg_results,
    showCategory = 15,
    title = "KEGG Pathway Enrichment"
  )
  
  print(kegg_plot)
  
  ggsave(
    "results/KEGG_enrichment_dotplot.png",
    plot = kegg_plot,
    width = 10,
    height = 8,
    dpi = 300
  )
  
} else {
  print("No significantly enriched KEGG pathways found.")
}