##############################################################
# AML Biomarker Discovery Pipeline
#
# Script:
# 01_merge_featureCounts.R
#
# Purpose:
# Merge Galaxy featureCounts output files into a single
# count matrix for downstream differential expression analysis.
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
cat("Stage 1 : Merge featureCounts Files\n")
cat("=========================================\n\n")

start_time <- Sys.time()

##############################################################
# Input Directories
##############################################################

COUNT_DIR <- file.path(DATA_DIR, "Counts")

PROCESSED_DIR <- file.path(DATA_DIR, "processed")

dir.create(
  PROCESSED_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

##############################################################
# Locate Count Files
##############################################################

count_files <- list.files(
  
  path = COUNT_DIR,
  
  pattern = "\\.tabular$",
  
  full.names = TRUE
  
)

cat("Searching featureCounts files...\n")

cat("Files Detected :", length(count_files), "\n\n")

if (length(count_files) == 0) {
  
  stop(
    "No featureCounts (.tabular) files found in data/Counts."
  )
  
}

##############################################################
# Read First Count File
##############################################################

merged_counts <- readr::read_tsv(
  
  count_files[1],
  
  show_col_types = FALSE
  
)

##############################################################
# Validate First File
##############################################################

if (ncol(merged_counts) != 2) {
  
  stop(
    "First featureCounts file must contain exactly two columns."
  )
  
}

if (colnames(merged_counts)[1] != "Geneid") {
  
  stop(
    "First column must be named 'Geneid'."
  )
  
}

##############################################################
# Merge Remaining Files
##############################################################

if (length(count_files) > 1) {
  
  for (i in 2:length(count_files)) {
    
    cat(
      
      "[",
      
      i,
      
      "/",
      
      length(count_files),
      
      "] ",
      
      basename(count_files[i]),
      
      "\n",
      
      sep = ""
      
    )
    
    temp <- readr::read_tsv(
      
      count_files[i],
      
      show_col_types = FALSE
      
    )
    
    if (ncol(temp) != 2) {
      
      stop(
        
        basename(count_files[i]),
        
        " does not contain exactly two columns."
        
      )
      
    }
    
    merged_counts <- dplyr::full_join(
      
      merged_counts,
      
      temp,
      
      by = "Geneid"
      
    )
    
  }
  
}

##############################################################
# Quality Control
##############################################################

duplicate_genes <- duplicated(merged_counts$Geneid)

missing_values <- sum(is.na(merged_counts))

cat("\n")
cat("Duplicate Gene IDs :", sum(duplicate_genes), "\n")

cat("Missing Values     :", missing_values, "\n")

##############################################################
# Save Merged Count Matrix
##############################################################

output_file <- file.path(
  
  PROCESSED_DIR,
  
  "merged_counts.csv"
  
)

readr::write_csv(
  
  merged_counts,
  
  output_file
  
)

##############################################################
# Save Log File
##############################################################

log_file <- file.path(
  
  LOG_DIR,
  
  "merge_log.txt"
  
)

sink(log_file)

cat("AML Biomarker Discovery Pipeline\n\n")

cat("Pipeline Version : ", PIPELINE_VERSION, "\n")

cat("Date             : ", Sys.time(), "\n\n")

cat("Files Merged     : ", length(count_files), "\n")

cat("Genes            : ", nrow(merged_counts), "\n")

cat("Samples          : ", ncol(merged_counts) - 1, "\n")

cat("Duplicate Genes  : ", sum(duplicate_genes), "\n")

cat("Missing Values   : ", missing_values, "\n")

sink()

##############################################################
# Summary
##############################################################

end_time <- Sys.time()

cat("\n")
cat("=========================================\n")
cat("Stage 1 Completed Successfully\n")
cat("=========================================\n\n")

cat("Genes             :", nrow(merged_counts), "\n")

cat("Samples           :", ncol(merged_counts) - 1, "\n")

cat("Output            :", output_file, "\n")

cat("Log               :", log_file, "\n")

cat("Time Elapsed      :", round(end_time - start_time, 2), "\n\n")