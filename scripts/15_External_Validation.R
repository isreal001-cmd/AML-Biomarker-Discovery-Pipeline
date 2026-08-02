###############################################################
# AML Biomarker Discovery Pipeline
# Script: 15_External_Validation.R
#
# Purpose:
# External validation of candidate biomarkers using
# TCGA, GEO, GTEx or user-supplied datasets.
###############################################################

# rm(list = ls())

###############################################################
# Load Configuration
###############################################################

source("config.R")
source("functions/validation_functions.R")

###############################################################
# Required Packages
###############################################################

packages <- c(
  "readr",
  "dplyr"
)

for(pkg in packages){
  
  if(!require(pkg, character.only = TRUE)){
    
    install.packages(pkg)
    
    library(pkg, character.only = TRUE)
    
  }
  
}

###############################################################
# Choose Dataset
###############################################################

validation_dataset <- "TCGA"

###############################################################
# Load Dataset
###############################################################

cat("---------------------------------------\n")
cat("Loading Validation Dataset\n")
cat("---------------------------------------\n\n")

dataset <- load_validation_dataset(
  validation_dataset
)

expr <- dataset$expression
metadata <- dataset$metadata

###############################################################
# Load Biomarkers
###############################################################

biomarkers <- load_biomarkers()

###############################################################
# Match Biomarkers
###############################################################

matched <- match_biomarkers(
  expr,
  biomarkers
)

###############################################################
# Detection Rate
###############################################################

detection <- calculate_detection_rate(
  matched
)

###############################################################
# Expression Statistics
###############################################################

statistics <- calculate_expression_statistics(
  matched
)

###############################################################
# Dataset Type
###############################################################

dataset_type <- detect_dataset_type(
  metadata
)

###############################################################
# Output Folder
###############################################################

dir.create(
  "results/validation",
  recursive = TRUE,
  showWarnings = FALSE
)

###############################################################
# Save Tables
###############################################################

write_csv(
  
  matched,
  
  "results/validation/Matched_Biomarkers.csv"
  
)

write_csv(
  
  detection,
  
  "results/validation/Biomarker_Detection.csv"
  
)

write_csv(
  
  statistics,
  
  "results/validation/Expression_Statistics.csv"
  
)

###############################################################
# Validation Summary
###############################################################

summary_table <- data.frame(
  
  Dataset = validation_dataset,
  
  TotalGenes = nrow(expr),
  
  BiomarkersMatched = nrow(matched),
  
  Samples = ncol(expr)-1,
  
  DatasetType = dataset_type
  
)

write_csv(
  
  summary_table,
  
  "results/validation/Validation_Summary.csv"
  
)

###############################################################
# Validation Report
###############################################################

report <- c(
  
  "AML Biomarker Discovery Pipeline",
  
  "",
  
  paste("Date:", Sys.time()),
  
  "",
  
  paste("Validation Dataset:", validation_dataset),
  
  paste("Dataset Type:", dataset_type),
  
  "",
  
  paste("Genes:", nrow(expr)),
  
  paste("Samples:", ncol(expr)-1),
  
  paste("Matched Biomarkers:", nrow(matched)),
  
  "",
  
  if(dataset_type=="single_group"){
    
    "Only one biological group detected."
    
  }else{
    
    "Two biological groups detected."
    
  },
  
  "",
  
  "Validation completed successfully."
  
)

writeLines(
  
  report,
  
  "results/validation/Validation_Report.txt"
  
)

###############################################################
# Console Output
###############################################################

cat("---------------------------------------\n")
cat("External Validation Completed\n")
cat("---------------------------------------\n\n")

cat("Dataset :", validation_dataset,"\n")
cat("Samples :", ncol(expr)-1,"\n")
cat("Genes :", nrow(expr),"\n")
cat("Matched Biomarkers :", nrow(matched),"\n")
cat("Dataset Type :", dataset_type,"\n\n")

cat("Files created:\n")

cat("results/validation/Matched_Biomarkers.csv\n")

cat("results/validation/Biomarker_Detection.csv\n")

cat("results/validation/Expression_Statistics.csv\n")

cat("results/validation/Validation_Summary.csv\n")

cat("results/validation/Validation_Report.txt\n")