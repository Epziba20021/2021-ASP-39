
# 09 - FINAL PROJECT SUMMARY

setwd("C:/Users/Admin/Desktop/Bioinformatics")

# Load analysis results
gene_counts <- read.csv("results/significant_gene_counts.csv")
important_genes <- read.csv("results/important_significant_genes_sorted.csv")
common_genes <- read.csv("results/common_genes.csv")
go_results <- read.csv("results/GO_enrichment_results.csv")
kegg_results <- read.csv("results/KEGG_enrichment_results.csv")

# Calculate summary statistics
total_important_genes <- nrow(important_genes)
total_common_genes <- nrow(common_genes)
total_go_terms <- nrow(go_results)
total_kegg_pathways <- nrow(kegg_results)

# Identify top enriched GO terms
top_go <- head(
  go_results[order(go_results$p.adjust), "Description"],
  5
)

# Identify top enriched KEGG pathways
top_kegg <- head(
  kegg_results[order(kegg_results$p.adjust), "Description"],
  5
)

# Create final project summary
summary_text <- c(
  
  "CHICKENPOX (VARICELLA) AND VARICELLA-ZOSTER VIRUS (VZV)",
  "",
  "INTRODUCTION",
  "",
  "Chickenpox, also known as varicella, is a contagious infectious disease caused by the Varicella-Zoster Virus (VZV). VZV is a member of the Herpesviridae family. After the initial infection, the virus can remain latent in the human body and may reactivate later in life, causing herpes zoster (shingles).",
  "",
  "During VZV infection, human host cells respond by activating different genes involved in antiviral defence and immune responses. Studying these host gene expression changes helps researchers understand how human cells respond to VZV infection.",
  "",
  "VZV HOST GENE EXPRESSION ANALYSIS",
  "",
  "PROJECT SUMMARY",
  "",
  "This project analyses five VZV-related host gene expression datasets obtained from the NCBI Gene Expression Omnibus (GEO), study GSE121385.",
  "",
  paste("Total datasets analysed:", nrow(gene_counts)),
  paste("Total unique significant genes:", length(unique(unlist(
    lapply(
      list(
        read.csv("results/significant_gene_frequency.csv")$Gene
      ),
      identity
    )
  )))),
  paste("Important significant genes:", total_important_genes),
  paste("Common significant genes:", total_common_genes),
  paste("Enriched GO terms:", total_go_terms),
  paste("Enriched KEGG pathways:", total_kegg_pathways),
  "",
  "SIGNIFICANT GENE ANALYSIS",
  "",
  capture.output(print(gene_counts)),
  "",
  "COMMON SIGNIFICANT GENES",
  paste(common_genes$GeneName, collapse = ", "),
  "",
  "IMPORTANT SIGNIFICANT GENES",
  paste(important_genes$Gene, collapse = ", "),
  "",
  "TOP GO BIOLOGICAL PROCESSES",
  top_go,
  "",
  "TOP KEGG PATHWAYS",
  top_kegg,
  "",
  "CONCLUSION",
  "",
  "This project identifies significant host gene expression changes associated with VZV-related datasets. The analysis highlights common and frequently occurring genes, including MX1, ISG15 and MX2. Gene Ontology and KEGG enrichment analyses provide further information about biological processes and pathways associated with the identified genes.",
  "",
  "The findings contribute to understanding host antiviral responses during Varicella-Zoster Virus infection."
)

# Display summary
cat(summary_text, sep = "\n")

# Save final summary report
writeLines(
  summary_text,
  "report/final_project_summary.txt"
)

cat("\nFINAL PROJECT SUMMARY CREATED SUCCESSFULLY!\n")
