###############################################################
# AML Biomarker Discovery Pipeline
# Script: 14A_Prepare_Validation_Data.R
# Purpose:
# Verify that external validation datasets are available
# and correctly organized before downstream validation.
###############################################################

# rm(list = ls())

###############################################################
# Load Configuration
###############################################################

source("config.R")

###############################################################
# Create Validation Directories
###############################################################

validation_dirs <- c(
  
  "data/validation",
  
  "data/validation/TCGA",
  
  "data/validation/GEO",
  
  "data/validation/GTEx",
  
  "data/validation/User"
  
)

for(dir in validation_dirs){
  
  dir.create(
    dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
}

###############################################################
# Expected Files
###############################################################

expected <- data.frame(
  
  Dataset = c(
    "TCGA",
    "GEO",
    "GTEx",
    "User"
  ),
  
  Expression = c(
    "data/validation/TCGA/expression.csv",
    "data/validation/GEO/expression.csv",
    "data/validation/GTEx/expression.csv",
    "data/validation/User/expression.csv"
  ),
  
  Metadata = c(
    "data/validation/TCGA/metadata.csv",
    "data/validation/GEO/metadata.csv",
    "data/validation/GTEx/metadata.csv",
    "data/validation/User/metadata.csv"
  ),
  
  stringsAsFactors = FALSE
  
)

###############################################################
# Check Availability
###############################################################

status <- data.frame()

for(i in seq_len(nrow(expected))){
  
  expr_exists <- file.exists(expected$Expression[i])
  
  meta_exists <- file.exists(expected$Metadata[i])
  
  status <- rbind(
    
    status,
    
    data.frame(
      
      Dataset = expected$Dataset[i],
      
      Expression = expr_exists,
      
      Metadata = meta_exists,
      
      Ready = expr_exists & meta_exists
      
    )
    
  )
  
}

###############################################################
# Save Report
###############################################################

dir.create(
  "results/validation",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  
  status,
  
  "results/validation/Validation_Dataset_Status.csv",
  
  row.names = FALSE
  
)

###############################################################
# Console Output
###############################################################

cat("---------------------------------------\n")
cat("Validation Dataset Assessment\n")
cat("---------------------------------------\n\n")

print(status)

cat("\n---------------------------------------\n")

ready <- sum(status$Ready)

cat("Datasets Ready:", ready, "of", nrow(status), "\n")

cat("---------------------------------------\n")

cat("\nFile created:\n")
cat("results/validation/Validation_Dataset_Status.csv\n")