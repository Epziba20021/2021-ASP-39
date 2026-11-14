# 05. Gene Visualization

library(ggplot2)
library(dplyr)

# Set working directory
setwd("C:/Users/Admin/Desktop/Bioinformatics")

# Load significant gene frequency results
gene_frequency <- read.csv("results/significant_gene_frequency.csv")

# Sort genes by frequency
gene_frequency <- gene_frequency %>%
  arrange(desc(Strain_Count), Gene)

# Divide genes into two groups
total_genes <- nrow(gene_frequency)
split_point <- ceiling(total_genes / 2)

genes_1 <- gene_frequency %>%
  slice(1:split_point)

genes_2 <- gene_frequency %>%
  slice((split_point + 1):total_genes)

# Create the first graph
plot_1 <- ggplot(
  genes_1,
  aes(
    x = reorder(Gene, Strain_Count),
    y = Strain_Count
  )
) +
  geom_col() +
  geom_text(
    aes(label = Strain_Count),
    hjust = -0.3,
    size = 3.5
  ) +
  coord_flip() +
  labs(
    title = "Significant Gene Frequency Across VZV Datasets (1 of 2)",
    x = "Gene Names",
    y = "Number of Datasets"
  ) +
  scale_y_continuous(
    breaks = 0:5,
    limits = c(0, 5.5)
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(size = 12, face = "bold")
  )

# Display and save the first graph
print(plot_1)

ggsave(
  "results/gene_frequency_1.png",
  plot = plot_1,
  width = 10,
  height = 9,
  dpi = 300
)

# Create the second graph
plot_2 <- ggplot(
  genes_2,
  aes(
    x = reorder(Gene, Strain_Count),
    y = Strain_Count
  )
) +
  geom_col() +
  geom_text(
    aes(label = Strain_Count),
    hjust = -0.3,
    size = 3.5
  ) +
  coord_flip() +
  labs(
    title = "Significant Gene Frequency Across VZV Datasets (2 of 2)",
    x = "Gene Names",
    y = "Number of Datasets"
  ) +
  scale_y_continuous(
    breaks = 0:5,
    limits = c(0, 5.5)
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(size = 12, face = "bold")
  )

# Display and save the second graph
print(plot_2)

ggsave(
  "results/gene_frequency_2.png",
  plot = plot_2,
  width = 10,
  height = 9,
  dpi = 300
)

# Display completion summary
cat("\nGENE VISUALIZATION COMPLETED SUCCESSFULLY!\n")
cat("Total significant genes:", total_genes, "\n")
cat("Genes in Graph 1:", nrow(genes_1), "\n")
cat("Genes in Graph 2:", nrow(genes_2), "\n")