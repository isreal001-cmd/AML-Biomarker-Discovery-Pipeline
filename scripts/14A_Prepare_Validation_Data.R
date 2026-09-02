##############################################################
# AML Biomarker Discovery Pipeline
#
# Script:
# 14A_Prepare_Validation_Data.R
#
# Purpose:
# Verify that external validation datasets are available
# and correctly organised before downstream validation.
#
# Author:
# Isreal Oluwafemi Abiodun
#
# Version:
# 1.0.0
##############################################################

##############################################################
# Load Configuration
##############################################################

source("config.R")

cat("\n")
cat("=========================================\n")
cat("Stage 14A : Validation Dataset Assessment\n")
cat("=========================================\n\n")

start_time <- Sys.time()

##############################################################
# Validation Directories
##############################################################

VALIDATION_DIR <- file.path(
  "data",
  "validation"
)

VALIDATION_RESULTS_DIR <- file.path(
  RESULTS_DIR,
  "validation"
)

dir.create(
  VALIDATION_RESULTS_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

validation_dirs <- c(
  
  file.path(VALIDATION_DIR),
  
  file.path(VALIDATION_DIR, "TCGA"),
  
  file.path(VALIDATION_DIR, "GEO"),
  
  file.path(VALIDATION_DIR, "GTEx"),
  
  file.path(VALIDATION_DIR, "User")
  
)

for (directory in validation_dirs) {
  
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
}

##############################################################
# Expected Validation Files
##############################################################

expected <- data.frame(
  
  Dataset = c(
    
    "TCGA",
    
    "GEO",
    
    "GTEx",
    
    "User"
    
  ),
  
  Expression = c(
    
    file.path(VALIDATION_DIR, "TCGA", "expression.csv"),
    
    file.path(VALIDATION_DIR, "GEO", "expression.csv"),
    
    file.path(VALIDATION_DIR, "GTEx", "expression.csv"),
    
    file.path(VALIDATION_DIR, "User", "expression.csv")
    
  ),
  
  Metadata = c(
    
    file.path(VALIDATION_DIR, "TCGA", "metadata.csv"),
    
    file.path(VALIDATION_DIR, "GEO", "metadata.csv"),
    
    file.path(VALIDATION_DIR, "GTEx", "metadata.csv"),
    
    file.path(VALIDATION_DIR, "User", "metadata.csv")
    
  ),
  
  stringsAsFactors = FALSE
  
)

##############################################################
# Assess Dataset Availability
##############################################################

status <- data.frame(
  
  Dataset = expected$Dataset,
  
  Expression = file.exists(expected$Expression),
  
  Metadata = file.exists(expected$Metadata)
  
)

status$Ready <- status$Expression & status$Metadata

##############################################################
# Save Status Report
##############################################################

STATUS_FILE <- file.path(
  
  VALIDATION_RESULTS_DIR,
  
  "Validation_Dataset_Status.csv"
  
)

readr::write_csv(
  
  status,
  
  STATUS_FILE
  
)

##############################################################
# Console Summary
##############################################################

cat("Validation Dataset Status\n")
cat("-----------------------------------------\n\n")

print(status)

cat("\n")

cat(
  "Datasets Ready :",
  sum(status$Ready),
  "of",
  nrow(status),
  "\n\n"
)

cat(
  "Status Report :",
  STATUS_FILE,
  "\n\n"
)

##############################################################
# Completion
##############################################################

end_time <- Sys.time()

cat("=========================================\n")
cat("Stage 14A Completed Successfully\n")
cat("=========================================\n\n")

cat(
  "Time Elapsed :",
  round(end_time - start_time, 2),
  "\n\n"
)