
# 03 - SIGNIFICANT GENE ANALYSIS

library(dplyr)
library(ggplot2)

setwd("C:/Users/Admin/Desktop/Bioinformatics")

datasets <- readRDS("results/imported_datasets.rds")

# Fix duplicate column names
datasets <- lapply(datasets, function(data) {
  names(data) <- make.unique(names(data))
  data
})

dir.create("results", showWarnings = FALSE)

# Extract unique significant genes from each dataset
significant_gene_lists <- lapply(datasets, function(data) {
  
  gene_names <- trimws(as.character(data[["GeneName"]]))
  q_values <- as.numeric(data[["q-value"]])
  
  valid <- !is.na(gene_names) &
    gene_names != "" &
    !is.na(q_values) &
    q_values <= 0.05
  
  unique(gene_names[valid])
})

# Count unique significant genes and total unique genes
sig_counts <- data.frame(
  Strain = names(datasets),
  
  Significant_Genes = sapply(
    significant_gene_lists,
    length
  ),
  
  Total_Genes = sapply(datasets, function(data) {
    genes <- trimws(as.character(data[["GeneName"]]))
    
    length(unique(
      genes[!is.na(genes) & genes != ""]
    ))
  })
)

sig_counts$Percentage <- round(
  sig_counts$Significant_Genes /
    sig_counts$Total_Genes * 100,
  2
)

cat("\nSIGNIFICANT GENE COUNTS:\n")
print(sig_counts)

# Plot significant gene counts
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

print(sig_plot)

# Save significant gene counts
write.csv(
  sig_counts,
  "results/significant_gene_counts.csv",
  row.names = FALSE
)

ggsave(
  "results/significant_gene_counts.png",
  plot = sig_plot,
  width = 9,
  height = 6,
  dpi = 300
)

# Plot significant gene percentages
percentage_plot <- ggplot(
  sig_counts,
  aes(x = Strain, y = Percentage)
) +
  geom_col() +
  geom_text(
    aes(label = paste0(Percentage, "%")),
    vjust = -0.5
  ) +
  labs(
    title = "Percentage of Significant Genes Across VZV Datasets",
    x = "Dataset",
    y = "Significant Genes (%)"
  ) +
  theme_minimal() +
  ylim(0, 60)

print(percentage_plot)

ggsave(
  "results/significant_gene_percentage.png",
  plot = percentage_plot,
  width = 9,
  height = 6,
  dpi = 300
)

# Calculate gene frequency across datasets
unique_genes <- unique(unlist(significant_gene_lists))

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

gene_frequency <- gene_frequency %>%
  arrange(desc(Strain_Count), Gene)

cat("\nSIGNIFICANT GENE FREQUENCY:\n")
print(gene_frequency)

cat("\nTotal unique significant genes:\n")
print(nrow(gene_frequency))

cat("\nGenes found in all five datasets:\n")
print(
  gene_frequency %>%
    filter(Strain_Count == length(datasets))
)

write.csv(
  gene_frequency,
  "results/significant_gene_frequency.csv",
  row.names = FALSE
)

# Identify important genes found in at least four datasets
cat("\nTop 20 Significant Genes:\n")
print(head(gene_frequency, 20))

important_genes <- gene_frequency %>%
  filter(Strain_Count >= 4) %>%
  arrange(desc(Strain_Count), Gene)

cat("\nIMPORTANT SIGNIFICANT GENES:\n")
print(important_genes)

cat(
  "\nTotal Important Significant Genes:",
  nrow(important_genes),
  "\n"
)

write.csv(
  important_genes,
  "results/important_significant_genes.csv",
  row.names = FALSE
)

cat("\nSIGNIFICANT GENE ANALYSIS COMPLETED!\n")
