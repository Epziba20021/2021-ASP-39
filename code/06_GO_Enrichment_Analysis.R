if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))

BiocManager::install("DOSE", ask = FALSE, update = FALSE)

BiocManager::install("enrichplot", ask = FALSE, update = FALSE)

packages <- c("clusterProfiler", "org.Hs.eg.db", "DOSE", "enrichplot", "GO.db")

sapply(packages, requireNamespace, quietly = TRUE)


setwd("C:/Users/Admin/Desktop/Bioinformatics")

datasets <- readRDS("results/imported_datasets.rds")

head(read.csv("results/significant_gene_frequency.csv"))




# 06. GO Enrichment Analysis

setwd("C:/Users/Admin/Desktop/Bioinformatics")

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

# Load significant gene frequency results
gene_data <- read.csv("results/significant_gene_frequency.csv")

# Extract unique gene symbols
gene_symbols <- unique(gene_data$Gene)

# Convert gene symbols into Entrez IDs
gene_mapping <- bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Remove duplicate Entrez IDs
gene_mapping <- gene_mapping %>%
  distinct(ENTREZID, .keep_all = TRUE)

# Display mapping results
print(gene_mapping)

# Save gene mapping results
write.csv(
  gene_mapping,
  "results/gene_id_mapping.csv",
  row.names = FALSE
)

names(datasets)




# Perform Gene Ontology enrichment analysis

go_results <- enrichGO(
  gene = gene_mapping$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

# Display enrichment results
go_table <- as.data.frame(go_results)

print(head(go_table))

# Save GO enrichment results
write.csv(
  go_table,
  "results/GO_enrichment_results.csv",
  row.names = FALSE
)




# Visualize top enriched GO terms

library(enrichplot)
library(ggplot2)

go_plot <- dotplot(
  go_results,
  showCategory = 15,
  title = "GO Biological Process Enrichment"
)

print(go_plot)

ggsave(
  "results/GO_enrichment_dotplot.png",
  plot = go_plot,
  width = 10,
  height = 8,
  dpi = 300
)