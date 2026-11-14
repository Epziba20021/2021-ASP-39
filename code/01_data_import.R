# DATA IMPORT

library(readxl)
library(dplyr)
library(ggplot2)

# Set project directory
setwd("C:/Users/Admin/Desktop/Bioinformatics")

# Find Excel files
files <- list.files(
  "data",
  pattern = "\\.xlsx$",
  full.names = TRUE
)

print(files)

# Check five files
if (length(files) != 5) {
  stop("Expected five Excel files in the data folder.")
}

# Import datasets
yc01_low <- read_excel(files[1], .name_repair = "minimal")
yc01_high <- read_excel(files[2], .name_repair = "minimal")
ca <- read_excel(files[3], .name_repair = "minimal")
cf <- read_excel(files[4], .name_repair = "minimal")
varilrix <- read_excel(files[5], .name_repair = "minimal")

# Combine datasets
datasets <- list(
  YC01_Low = yc01_low,
  YC01_High = yc01_high,
  Suduvax_CA = ca,
  Suduvax_CF = cf,
  VarilRix = varilrix
)

# Display dataset dimensions
for (name in names(datasets)) {
  cat("\nDataset:", name, "\n")
  print(dim(datasets[[name]]))
}

# Save imported datasets
saveRDS(
  datasets,
  "results/imported_datasets.rds"
)

cat("\nAll five datasets imported successfully!\n")

list.files("results")