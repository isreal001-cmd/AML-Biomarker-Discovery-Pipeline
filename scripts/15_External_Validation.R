###############################################################
# AML Biomarker Discovery Pipeline
# Script: 15_External_Validation.R
#
# Purpose:
# Validate candidate AML biomarkers using independent
# external datasets (TCGA, GEO, GTEx or User dataset)
#
# Author:
# Isreal Oluwafemi Abiodun
###############################################################

###############################################################
# Initialise Environment
###############################################################

# rm(list = ls())

source("config.R")
source("functions/validation_functions.R")

cat("=========================================\n")
cat("Stage 15 : External Biomarker Validation\n")
cat("=========================================\n\n")

###############################################################
# Load Required Packages
###############################################################

required_packages <- c(
  
  "readr",
  
  "dplyr"
  
)

install_and_load_packages(required_packages)

###############################################################
# Select Validation Dataset
###############################################################

validation_dataset <- "TCGA"

###############################################################
# Load Validation Dataset
###############################################################

cat("Loading validation dataset...\n")

dataset <- load_validation_dataset(
  validation_dataset
)

expression <- dataset$expression

metadata <- dataset$metadata

###############################################################
# Load Biomarkers
###############################################################

cat("Loading biomarker candidates...\n")

biomarkers <- load_biomarkers()

###############################################################
# Match Biomarkers
###############################################################

matched <- match_biomarkers(
  
  expression,
  
  biomarkers
  
)

###############################################################
# Calculate Detection Rate
###############################################################

cat("Calculating detection rate...\n")

detection <- calculate_detection_rate(
  matched
)

###############################################################
# Expression Statistics
###############################################################

cat("Calculating expression statistics...\n")

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
# Validation Summary
###############################################################

summary_table <- validation_summary(
  
  expression,
  
  metadata,
  
  matched
  
)

summary_table$Dataset <- validation_dataset

###############################################################
# Create Output Directory
###############################################################

output_dir <- file.path(
  
  RESULTS_DIR,
  
  "validation"
  
)

dir.create(
  
  output_dir,
  
  recursive = TRUE,
  
  showWarnings = FALSE
  
)

###############################################################
# Save Tables
###############################################################

write.csv(
  
  matched,
  
  file.path(
    
    output_dir,
    
    "Matched_Biomarkers.csv"
    
  ),
  
  row.names = FALSE
  
)

write.csv(
  
  detection,
  
  file.path(
    
    output_dir,
    
    "Biomarker_Detection.csv"
    
  ),
  
  row.names = FALSE
  
)

write.csv(
  
  statistics,
  
  file.path(
    
    output_dir,
    
    "Expression_Statistics.csv"
    
  ),
  
  row.names = FALSE
  
)

write.csv(
  
  summary_table,
  
  file.path(
    
    output_dir,
    
    "Validation_Summary.csv"
    
  ),
  
  row.names = FALSE
  
)

###############################################################
# Generate Validation Report
###############################################################

report <- c(
  
  "AML Biomarker Discovery Pipeline",
  
  "",
  
  paste("Pipeline Version :", PIPELINE_VERSION),
  
  paste("Author :", AUTHOR),
  
  paste("Date :", Sys.time()),
  
  "",
  
  "External Validation Summary",
  
  "",
  
  paste("Validation Dataset :", validation_dataset),
  
  paste("Dataset Type :", dataset_type),
  
  "",
  
  paste("Total Genes :", nrow(expression)),
  
  paste("Total Samples :", ncol(expression) - 1),
  
  paste("Candidate Biomarkers :", nrow(biomarkers)),
  
  paste("Matched Biomarkers :", nrow(matched)),
  
  "",
  
  "Detection Summary",
  
  "",
  
  paste(
    
    "Average Detection Rate :",
    
    round(
      
      mean(detection$DetectionRate),
      
      3
      
    )
    
  ),
  
  paste(
    
    "Minimum Detection Rate :",
    
    round(
      
      min(detection$DetectionRate),
      
      3
      
    )
    
  ),
  
  paste(
    
    "Maximum Detection Rate :",
    
    round(
      
      max(detection$DetectionRate),
      
      3
      
    )
    
  ),
  
  "",
  
  if(dataset_type == "single_group"){
    
    c(
      
      "Interpretation",
      
      "",
      
      "Only one biological group is available.",
      
      "Differential validation cannot be performed.",
      
      "Only biomarker presence and expression",
      
      "statistics are reported."
      
    )
    
  }else{
    
    c(
      
      "Interpretation",
      
      "",
      
      "Multiple biological groups detected.",
      
      "Dataset is suitable for comparative",
      
      "validation analysis."
      
    )
    
  },
  
  "",
  
  "Validation completed successfully."
  
)

writeLines(
  
  report,
  
  file.path(
    
    output_dir,
    
    "Validation_Report.txt"
    
  )
  
)

###############################################################
# Console Summary
###############################################################

cat("\n")
cat("=========================================\n")
cat("External Validation Completed Successfully\n")
cat("=========================================\n\n")

cat("Validation Dataset :", validation_dataset, "\n")
cat("Dataset Type       :", dataset_type, "\n")
cat("Genes             :", nrow(expression), "\n")
cat("Samples           :", ncol(expression) - 1, "\n")
cat("Matched Biomarkers:", nrow(matched), "\n\n")

cat("Output Files\n\n")

cat("Matched_Biomarkers.csv\n")
cat("Biomarker_Detection.csv\n")
cat("Expression_Statistics.csv\n")
cat("Validation_Summary.csv\n")
cat("Validation_Report.txt\n")

list.files(
  "data",
  recursive = TRUE,
  full.names = TRUE
)

list.files(
  "results/validation",
  recursive = TRUE,
  full.names = TRUE
)

