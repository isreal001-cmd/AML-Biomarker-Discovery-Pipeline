##############################################################
# AML Biomarker Discovery Pipeline
#
# Script:
# 04_VST_normalization.R
#
# Purpose:
# Perform Variance Stabilizing Transformation (VST)
# on the DESeq2 dataset for downstream visualization
# and machine learning analyses.
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
cat("Stage 4 : Variance Stabilizing Transformation\n")
cat("=========================================\n\n")

start_time <- Sys.time()

##############################################################
# Define Directories
##############################################################

VST_DIR <- file.path(
  RESULTS_DIR,
  "vst"
)

dir.create(
  VST_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

##############################################################
# Input File
##############################################################

DDS_FILE <- file.path(
  RESULTS_DIR,
  "deseq2",
  "dds.rds"
)

##############################################################
# Output Files
##############################################################

VSD_FILE <- file.path(
  VST_DIR,
  "vsd.rds"
)

VST_MATRIX_FILE <- file.path(
  VST_DIR,
  "vst_expression_matrix.csv"
)

##############################################################
# Validate Input
##############################################################

if (!file.exists(DDS_FILE)) {
  
  stop(
    "DESeq2 object not found.\nRun 03_DESeq2_analysis.R first."
  )
  
}

##############################################################
# Load DESeq2 Object
##############################################################

cat("Loading DESeq2 object...\n")

dds <- readRDS(
  DDS_FILE
)

cat("DESeq2 Object Successfully Loaded\n\n")

##############################################################
# Perform Variance Stabilizing Transformation
##############################################################

cat("Running Variance Stabilizing Transformation...\n\n")

vsd <- DESeq2::vst(
  
  dds,
  
  blind = FALSE
  
)

cat("VST Completed Successfully\n")

##############################################################
# Extract Expression Matrix
##############################################################

vst_matrix <- SummarizedExperiment::assay(vsd)

cat("\n")
cat("-----------------------------------------\n")
cat("VST Matrix Summary\n")
cat("-----------------------------------------\n")

cat(
  "Genes   :",
  nrow(vst_matrix),
  "\n"
)

cat(
  "Samples :",
  ncol(vst_matrix),
  "\n\n"
)

##############################################################
# Save Outputs
##############################################################

saveRDS(
  
  vsd,
  
  VSD_FILE
  
)

readr::write_csv(
  
  as.data.frame(vst_matrix),
  
  VST_MATRIX_FILE
  
)

##############################################################
# Summary
##############################################################

end_time <- Sys.time()

cat("\n")
cat("=========================================\n")
cat("Stage 4 Completed Successfully\n")
cat("=========================================\n\n")

cat(
  "VST Object        :",
  VSD_FILE,
  "\n"
)

cat(
  "Expression Matrix :",
  VST_MATRIX_FILE,
  "\n\n"
)

cat(
  "Genes             :",
  nrow(vst_matrix),
  "\n"
)

cat(
  "Samples           :",
  ncol(vst_matrix),
  "\n\n"
)

cat(
  "Time Elapsed      :",
  round(end_time - start_time, 2),
  "\n\n"
)