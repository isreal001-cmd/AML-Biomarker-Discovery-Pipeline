##############################################################
# AML Biomarker Discovery Pipeline
#
# Script:
# 02_prepare_metadata.R
#
# Purpose:
# Generate sample metadata required for DESeq2 analysis
# from the merged featureCounts count matrix.
#
# Author:
# Isreal Oluwafemi Abiodun
#
# Version:
# 1.0.0
##############################################################

##############################################################
# Load Global Configuration
##############################################################

source("config.R")

cat("\n")
cat("=========================================\n")
cat("Stage 2 : Prepare Sample Metadata\n")
cat("=========================================\n\n")

start_time <- Sys.time()

##############################################################
# Define Input and Output Files
##############################################################

INPUT_FILE <- file.path(
  DATA_DIR,
  "processed",
  "merged_counts.csv"
)

METADATA_DIR <- file.path(
  DATA_DIR,
  "metadata"
)

OUTPUT_FILE <- file.path(
  METADATA_DIR,
  "sample_metadata.csv"
)

dir.create(
  METADATA_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

##############################################################
# Validate Input
##############################################################

if (!file.exists(INPUT_FILE)) {
  
  stop(
    "Merged count matrix not found.\nRun 01_merge_featureCounts.R first."
  )
  
}

##############################################################
# Read Count Matrix
##############################################################

counts <- readr::read_csv(
  
  INPUT_FILE,
  
  show_col_types = FALSE
  
)

##############################################################
# Extract Sample Names
##############################################################

sample_names <- colnames(counts)[-1]

metadata <- data.frame(
  
  Sample = sample_names,
  
  stringsAsFactors = FALSE
  
)

##############################################################
# Assign Biological Condition
##############################################################

metadata$Condition <- ifelse(
  
  stringr::str_detect(
    metadata$Sample,
    "^AML"
  ),
  
  "AML",
  
  "Healthy"
  
)

##############################################################
# Assign Sequencing Platform
##############################################################

metadata$Platform <- ifelse(
  
  metadata$Condition == "AML",
  
  "NextSeq550",
  
  "HiSeq1000"
  
)

##############################################################
# Convert Variables to Factors
##############################################################

metadata$Condition <- factor(
  
  metadata$Condition,
  
  levels = c(
    "Healthy",
    "AML"
  )
  
)

metadata$Platform <- factor(
  
  metadata$Platform
  
)

##############################################################
# Quality Control
##############################################################

duplicate_samples <- sum(
  
  duplicated(metadata$Sample)
  
)

cat("Metadata Summary\n")
cat("-----------------------------------------\n\n")

cat("Total Samples      :", nrow(metadata), "\n")

cat("Duplicate Samples  :", duplicate_samples, "\n\n")

cat("Condition Counts\n")
print(table(metadata$Condition))

cat("\nPlatform Counts\n")
print(table(metadata$Platform))

##############################################################
# Save Metadata
##############################################################

readr::write_csv(
  
  metadata,
  
  OUTPUT_FILE
  
)

##############################################################
# Summary
##############################################################

end_time <- Sys.time()

cat("\n")
cat("=========================================\n")
cat("Stage 2 Completed Successfully\n")
cat("=========================================\n\n")

cat("Metadata File     :", OUTPUT_FILE, "\n")

cat("Samples           :", nrow(metadata), "\n")

cat("Healthy Samples   :",
    
    sum(metadata$Condition == "Healthy"),
    
    "\n")

cat("AML Samples       :",
    
    sum(metadata$Condition == "AML"),
    
    "\n")

cat("Time Elapsed      :",
    
    round(end_time - start_time, 2),
    
    "\n\n")