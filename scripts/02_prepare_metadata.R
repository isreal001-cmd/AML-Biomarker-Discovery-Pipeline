###############################################################
# AML Biomarker Discovery Pipeline
# Stage 2: Prepare Sample Metadata
#
# Author: Isreal Oluwafemi Abiodun
#
# Description:
# Creates the metadata table required for DESeq2 from the
# merged featureCounts matrix.
###############################################################

#--------------------------------------------------------------
# 1. Clear Workspace
#--------------------------------------------------------------

# rm(list = ls())

#--------------------------------------------------------------
# 2. Load Packages
#--------------------------------------------------------------

packages <- c("readr", "dplyr", "stringr")

for(pkg in packages){
  
  if(!require(pkg, character.only = TRUE)){
    
    install.packages(pkg)
    
    library(pkg, character.only = TRUE)
    
  }
  
}

#--------------------------------------------------------------
# 3. Define Directories
#--------------------------------------------------------------

input_file <- "data/processed/merged_counts.csv"

metadata_dir <- "data/metadata"

dir.create(metadata_dir,
           recursive = TRUE,
           showWarnings = FALSE)

#--------------------------------------------------------------
# 4. Read Count Matrix
#--------------------------------------------------------------

counts <- read_csv(
  input_file,
  show_col_types = FALSE
)

#--------------------------------------------------------------
# 5. Extract Sample Names
#--------------------------------------------------------------

sample_names <- colnames(counts)[-1]

metadata <- data.frame(
  
  Sample = sample_names,
  
  stringsAsFactors = FALSE
  
)

#--------------------------------------------------------------
# 6. Assign Biological Condition
#--------------------------------------------------------------

metadata$Condition <- ifelse(
  
  str_detect(metadata$Sample, "^AML"),
  
  "AML",
  
  "Healthy"
  
)

#--------------------------------------------------------------
# 7. Assign Sequencing Platform
#--------------------------------------------------------------

metadata$Platform <- ifelse(
  
  metadata$Condition == "AML",
  
  "NextSeq550",
  
  "HiSeq1000"
  
)

#--------------------------------------------------------------
# 8. Convert to Factors
#--------------------------------------------------------------

metadata$Condition <- factor(
  
  metadata$Condition,
  
  levels = c("Healthy","AML")
  
)

metadata$Platform <- factor(
  
  metadata$Platform
  
)

#--------------------------------------------------------------
# 9. Quality Checks
#--------------------------------------------------------------

cat("---------------------------------------\n")

cat("Metadata Summary\n")

cat("---------------------------------------\n\n")

cat("Total Samples :", nrow(metadata), "\n\n")

cat("Condition Counts\n")

print(table(metadata$Condition))

cat("\nPlatform Counts\n")

print(table(metadata$Platform))

cat("\nDuplicate Samples : ")

print(sum(duplicated(metadata$Sample)))

#--------------------------------------------------------------
# 10. Save Metadata
#--------------------------------------------------------------

write_csv(
  
  metadata,
  
  file.path(
    
    metadata_dir,
    
    "sample_metadata.csv"
    
  )
  
)

cat("\n")

cat("---------------------------------------\n")

cat("Metadata Successfully Generated\n")

cat("---------------------------------------\n")

cat("\nSaved to:\n")

cat("data/metadata/sample_metadata.csv\n")
