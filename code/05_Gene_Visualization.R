
# 05 - GENE VISUALIZATION

library(ggplot2)
library(dplyr)

setwd("C:/Users/Admin/Desktop/Bioinformatics")

# Load significant gene frequency results
gene_frequency <- read.csv(
  "results/significant_gene_frequency.csv"
)

gene_frequency <- gene_frequency %>%
  arrange(desc(Strain_Count), Gene)

# Load corrected significant gene counts
sig_counts <- read.csv(
  "results/significant_gene_counts.csv"
)

# Visualize significant gene percentages
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

# Divide genes into two groups
total_genes <- nrow(gene_frequency)
split_point <- ceiling(total_genes / 2)

genes_1 <- gene_frequency %>%
  slice_head(n = split_point)

genes_2 <- gene_frequency %>%
  slice_tail(n = total_genes - split_point)

# Create the first gene frequency graph
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

print(plot_1)

ggsave(
  "results/gene_frequency_1.png",
  plot = plot_1,
  width = 10,
  height = 9,
  dpi = 300
)

# Create the second gene frequency graph
if (nrow(genes_2) > 0) {
  
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
  
  print(plot_2)
  
  ggsave(
    "results/gene_frequency_2.png",
    plot = plot_2,
    width = 10,
    height = 9,
    dpi = 300
  )
}

# Display completion summary
cat("\nGENE VISUALIZATION COMPLETED SUCCESSFULLY!\n")
cat("Total significant genes:", total_genes, "\n")
cat("Genes in Graph 1:", nrow(genes_1), "\n")
cat("Genes in Graph 2:", nrow(genes_2), "\n")
