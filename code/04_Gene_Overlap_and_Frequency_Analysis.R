# 04. GENE OVERLAP AND FREQUENCY ANALYSIS

# Load required packages
library(dplyr)
library(ggplot2)

# 2. Set working directory
setwd("C:/Users/Admin/Desktop/Bioinformatics")

# 3. Load imported datasets
datasets <- readRDS("results/imported_datasets.rds")

# 4. Fix duplicate column names
datasets <- lapply(datasets, function(data) {
  names(data) <- make.unique(names(data))
  data
})


# ==========================================================
# PART A: COMMON SIGNIFICANT GENE ANALYSIS
# ==========================================================

# 5. Extract significant genes from each dataset

significant_gene_lists <- lapply(datasets, function(data) {
  
  data %>%
    filter(
      !is.na(GeneName),
      !is.na(`q-value`),
      `q-value` <= 0.05
    ) %>%
    pull(GeneName) %>%
    unique()
  
})


# Display significant gene counts
cat("\n===== PART A: COMMON SIGNIFICANT GENE ANALYSIS =====\n")

cat("\nSignificant genes in each dataset:\n")

for (i in seq_along(significant_gene_lists)) {
  cat(names(datasets)[i], ":",
      length(significant_gene_lists[[i]]), "\n")
}


# Identify genes common to all five datasets
common_genes <- Reduce(
  intersect,
  significant_gene_lists
)


# Display common genes
cat("\nCommon significant genes:\n")

print(common_genes)

cat("\nTotal common significant genes:\n")

print(length(common_genes))


# Create result table
common_gene_table <- data.frame(
  GeneName = common_genes
)


# Save common genes
write.csv(
  common_gene_table,
  "results/common_genes.csv",
  row.names = FALSE
)


# OVERALL GENE FREQUENCY ANALYSIS
cat("\n===== PART B: OVERALL GENE FREQUENCY ANALYSIS =====\n")

# Extract all unique genes from each dataset
gene_lists <- lapply(datasets, function(data) {
  unique(na.omit(data$GeneName))
})

# Count gene frequency across datasets
all_genes <- unlist(gene_lists)

gene_frequency <- as.data.frame(table(all_genes))

colnames(gene_frequency) <- c("Gene", "Dataset_Count")


# Sort genes by frequency
gene_frequency <- gene_frequency %>%
  arrange(desc(Dataset_Count), Gene)


# Display top 20 frequently occurring genes
cat("\nTop 20 frequently occurring genes:\n")

print(head(gene_frequency, 20))


# Save overall gene frequency results
write.csv(
  gene_frequency,
  "results/overall_gene_frequency.csv",
  row.names = FALSE
)

# Create top 15 gene frequency visualization

top_genes <- head(gene_frequency, 15)

p <- ggplot(
  top_genes,
  aes(
    x = reorder(Gene, Dataset_Count),
    y = Dataset_Count
  )
) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top 15 Frequently Occurring Genes",
    x = "Gene",
    y = "Number of Datasets"
  ) +
  theme_minimal()


# Display graph
print(p)


# Save graph
ggsave(
  "results/gene_frequency.png",
  plot = p,
  width = 9,
  height = 7,
  dpi = 300
)

# 19. Completion message
cat("\n============================================\n")
cat("GENE OVERLAP AND FREQUENCY ANALYSIS COMPLETED!\n")
cat("============================================\n")

cat("\nCompleted analyses:\n")
cat("1. Common Significant Gene Analysis\n")
cat("2. Overall Gene Frequency Analysis\n")