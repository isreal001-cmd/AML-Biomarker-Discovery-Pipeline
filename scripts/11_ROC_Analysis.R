###############################################################
# AML Biomarker Discovery Pipeline
# Script: 11_ROC_Analysis.R
# Purpose:
# ROC analysis of top biomarker candidates
###############################################################

# rm(list = ls())

source("config.R")
source("functions/plotting_functions.R")

###############################################################
# Load packages
###############################################################

packages <- c(
  "readr",
  "pROC",
  "ggplot2",
  "dplyr"
)

for(pkg in packages){
  
  if(!require(pkg, character.only = TRUE)){
    
    install.packages(pkg)
    
    library(pkg, character.only = TRUE)
    
  }
  
}

###############################################################
# Load data
###############################################################

vsd <- readRDS(
  "results/vst/vsd.rds"
)

metadata <- read_csv(
  "data/metadata/sample_metadata.csv",
  show_col_types = FALSE
)

biomarkers <- read_csv(
  "results/biomarkers/Top20_Biomarkers.csv",
  show_col_types = FALSE
)

cat("---------------------------------------\n")
cat("Input Files Loaded Successfully\n")
cat("---------------------------------------\n\n")

###############################################################
# Expression matrix
###############################################################

expr <- assay(vsd)

###############################################################
# Prepare phenotype
###############################################################

group <- ifelse(
  metadata$Condition == "AML",
  1,
  0
)

###############################################################
# Output folder
###############################################################

dir.create(
  "results/roc",
  recursive = TRUE,
  showWarnings = FALSE
)

###############################################################
# ROC analysis
###############################################################

roc_results <- data.frame()

for(gene in biomarkers$GeneID){
  
  values <- expr[gene, ]
  
  roc_obj <- roc(
    response = group,
    predictor = values,
    quiet = TRUE
  )
  
  auc_value <- as.numeric(
    auc(roc_obj)
  )
  
  roc_results <- rbind(
    roc_results,
    data.frame(
      GeneID = gene,
      AUC = auc_value
    )
  )
  
}

roc_results <- roc_results %>%
  
  arrange(desc(AUC))

###############################################################
# Save AUC table
###############################################################

write.csv(
  roc_results,
  "results/roc/ROC_AUC_Table.csv",
  row.names = FALSE
)

###############################################################
# Plot best ROC curve
###############################################################

best_gene <- roc_results$GeneID[1]

best_values <- expr[best_gene, ]

best_roc <- roc(
  response = group,
  predictor = best_values,
  quiet = TRUE
)

png(
  "results/roc/Best_Biomarker_ROC.png",
  width = 1800,
  height = 1600,
  res = 300
)

plot(
  best_roc,
  col = AML_RED,
  lwd = 3,
  main = paste(
    "ROC Curve:",
    best_gene
  )
)

abline(
  a = 0,
  b = 1,
  lty = 2
)

dev.off()

pdf(
  "results/roc/Best_Biomarker_ROC.pdf",
  width = 6,
  height = 6
)

plot(
  best_roc,
  col = AML_RED,
  lwd = 3,
  main = paste(
    "ROC Curve:",
    best_gene
  )
)

abline(
  a = 0,
  b = 1,
  lty = 2
)

dev.off()

###############################################################
# Summary
###############################################################

cat("---------------------------------------\n")
cat("ROC Analysis Completed Successfully\n")
cat("---------------------------------------\n\n")

cat("Best Biomarker:\n")
cat(best_gene, "\n")

cat("\nHighest AUC:\n")
print(max(roc_results$AUC))

cat("\nFiles created:\n")

cat("results/roc/ROC_AUC_Table.csv\n")
cat("results/roc/Best_Biomarker_ROC.png\n")
cat("results/roc/Best_Biomarker_ROC.pdf\n")