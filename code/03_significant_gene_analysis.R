# 03 - SIGNIFICANT GENE ANALYSIS

# 1. Load required packages

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

# SIGNIFICANT GENE COUNT ANALYSIS
# Count significant genes in each dataset

sig_counts <- data.frame(
  Strain = names(datasets),
  
  Significant_Genes = sapply(datasets, function(data) {
    sum(
      !is.na(data[["q-value"]]) &
        data[["q-value"]] <= 0.05
    )
  }),
  
  Total_Genes = sapply(datasets, nrow)
)

# Calculate percentage of significant genes
sig_counts$Percentage <- round(
  (sig_counts$Significant_Genes /
     sig_counts$Total_Genes) * 100,
  2
)

# Display results
cat("\nSIGNIFICANT GENE COUNTS:\n")
print(sig_counts)

# Create comparison bar chart
sig_plot <- ggplot(
  sig_counts,
  aes(x = Strain, y = Significant_Genes)
) +
  geom_col() +
  geom_text(
    aes(label = Significant_Genes),
    vjust = -0.5
  ) +
  labs(
    title = "Significant Gene Counts Across VZV Datasets",
    x = "VZV Dataset",
    y = "Number of Significant Genes"
  ) +
  theme_minimal()

# Display graph
print(sig_plot)

# Save significant gene count table
write.csv(
  sig_counts,
  "results/significant_gene_counts.csv",
  row.names = FALSE
)

# Save graph
ggsave(
  "results/significant_gene_counts.png",
  plot = sig_plot,
  width = 9,
  height = 6,
  dpi = 300
)

# SIGNIFICANT GENE FREQUENCY ANALYSIS
#Extract significant genes from each dataset

significant_gene_lists <- lapply(datasets, function(data) {
  
  data %>%
    filter(
      !is.na(GeneName),
      GeneName != "",
      !is.na(`q-value`),
      `q-value` <= 0.05
    ) %>%
    pull(GeneName) %>%
    unique()
  
})

# Combine significant genes

unique_genes <- unique(
  unlist(significant_gene_lists)
)

# Calculate gene frequency across datasets
gene_frequency <- data.frame(
  Gene = unique_genes
)

gene_frequency$Strain_Count <- sapply(
  gene_frequency$Gene,
  function(g) {
    
    sum(sapply(
      significant_gene_lists,
      function(x) g %in% x
    ))
    
  }
)


# Sort genes by frequency
gene_frequency <- gene_frequency %>%
  arrange(desc(Strain_Count), Gene)

# Display gene frequency table
cat("\nSIGNIFICANT GENE FREQUENCY:\n")
print(gene_frequency)


# Display total unique significant genes
cat("\nTotal unique significant genes:\n")
print(nrow(gene_frequency))


# Display genes found in all five datasets
cat("\nGenes found in all five datasets:\n")

print(
  gene_frequency %>%
    filter(Strain_Count == length(datasets))
)

# Save gene frequency table
write.csv(
  gene_frequency,
  "results/significant_gene_frequency.csv",
  row.names = FALSE
)



# IMPORTANT SIGNIFICANT GENE ANALYSIS

# 20. Display top 20 significant genes
cat("\nTop 20 Significant Genes:\n")
print(head(gene_frequency, 20))

# Select genes found in at least 4 datasets
important_genes <- gene_frequency %>%
  filter(Strain_Count >= 4) %>%
  arrange(desc(Strain_Count), Gene)

# Display important significant genes

cat("\nIMPORTANT SIGNIFICANT GENES:\n")
print(important_genes)

# Count important significant genes
cat(
  "\nTotal Important Significant Genes:",
  nrow(important_genes),
  "\n"
)

# Save important significant genes
write.csv(
  important_genes,
  "results/important_significant_genes.csv",
  row.names = FALSE
)

# COMPLETION MESSAGE
cat("\n====================================")
cat("\nSIGNIFICANT GENE ANALYSIS COMPLETED!")
cat("\n====================================")

cat("\n\nCompleted analyses:")
cat("\n1. Significant Gene Count Analysis")
cat("\n2. Significant Gene Frequency Analysis")
cat("\n3. Important Significant Gene Analysis\n")
