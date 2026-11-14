# 02. DATA EXPLORATION AND QUALITY CHECKING

# Load required package
library(dplyr)


# Set working directory
setwd("C:/Users/Admin/Desktop/Bioinformatics")


# Load imported datasets
datasets <- readRDS("results/imported_datasets.rds")

#Explore dataset dimensions
cat("\n========== DATASET DIMENSIONS ==========\n")

for (name in names(datasets)) {
  
  cat("\nDataset:", name, "\n")
  
  print(dim(datasets[[name]]))
  
}


# Display column names
cat("\n========== COLUMN NAMES ==========\n")

for (name in names(datasets)) {
  
  cat("\nDataset:", name, "\n")
  
  print(names(datasets[[name]]))
  
}


# Check dataset structure
cat("\n========== DATASET STRUCTURE ==========\n")

for (name in names(datasets)) {
  
  cat("\nDataset:", name, "\n")
  
  str(datasets[[name]])
  
}


# Summary of numerical values
cat("\n========== NUMERICAL SUMMARY ==========\n")

for (name in names(datasets)) {
  
  cat("\nDataset:", name, "\n")
  
  print(summary(datasets[[name]]))
  
}


# Check missing values
cat("\n========== MISSING VALUES ==========\n")

for (name in names(datasets)) {
  
  cat("\nDataset:", name, "\n")
  
  print(colSums(is.na(datasets[[name]])))
  
}


# Check duplicate genes
cat("\n========== DUPLICATE GENES ==========\n")

for (name in names(datasets)) {
  
  cat("\nDataset:", name, "\n")
  
  print(sum(duplicated(datasets[[name]]$GeneName)))
  
}


# Count significant genes using q-value
cat("\n========== SIGNIFICANT GENES (q-value <= 0.05) ==========\n")

for (name in names(datasets)) {
  
  count <- sum(
    datasets[[name]][["q-value"]] <= 0.05,
    na.rm = TRUE
  )
  
  cat(name, ":", count, "\n")
  
}


# Count genes using p-value
cat("\n========== SIGNIFICANT GENES (p-value <= 0.05) ==========\n")

for (name in names(datasets)) {
  
  count <- sum(
    datasets[[name]][["p-value"]] <= 0.05,
    na.rm = TRUE
  )
  
  cat(name, ":", count, "\n")
  
}


# Count UP and DOWN regulated genes
cat("\n========== UP AND DOWN REGULATED GENES ==========\n")

for (name in names(datasets)) {
  
  cat("\nDataset:", name, "\n")
  
  print(table(datasets[[name]][["DEG:UP_DOWN"]], useNA = "ifany"))
  
}


cat("\nData exploration and quality checking completed!\n")