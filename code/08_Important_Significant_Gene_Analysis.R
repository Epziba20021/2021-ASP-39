# 08. Important Significant Gene Analysis

setwd("C:/Users/Admin/Desktop/Bioinformatics")

library(dplyr)
library(ggplot2)

# Load important significant genes
important_genes <- read.csv(
  "results/important_significant_genes.csv"
)

# Display important genes
print(important_genes)

# Sort genes by frequency
important_genes <- important_genes %>%
  arrange(desc(Strain_Count))

# Create gene frequency visualization
gene_plot <- ggplot(
  important_genes,
  aes(
    x = reorder(Gene, Strain_Count),
    y = Strain_Count
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Important Significant Genes Across Five Datasets",
    x = "Gene",
    y = "Number of Datasets"
  ) +
  theme_minimal()

print(gene_plot)

# Save visualization
ggsave(
  "results/important_gene_frequency.png",
  plot = gene_plot,
  width = 9,
  height = 7,
  dpi = 300
)

# Save sorted gene data
write.csv(
  important_genes,
  "results/important_significant_genes_sorted.csv",
  row.names = FALSE
)