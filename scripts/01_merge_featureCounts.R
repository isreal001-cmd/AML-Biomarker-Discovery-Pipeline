###############################################################
# AML Biomarker Discovery Pipeline
# Stage 1: Merge Galaxy featureCounts Outputs
#
# Author: Isreal Oluwafemi Abiodun
# Description:
#   Merges all Galaxy featureCounts (.tabular) files into a
#   single count matrix suitable for DESeq2.
###############################################################

# rm(list = ls())

##############################
# Load Required Packages
##############################

packages <- c("readr", "dplyr", "purrr")

for(pkg in packages){
  
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
  
}

##############################
# Define Project Directories
##############################

count_dir <- "data/counts"

output_dir <- "data/processed"

log_dir <- "results/logs"

dir.create(output_dir,
           recursive = TRUE,
           showWarnings = FALSE)

dir.create(log_dir,
           recursive = TRUE,
           showWarnings = FALSE)

##############################
# Locate Count Files
##############################

count_files <- list.files(
  path = count_dir,
  pattern = "\\.tabular$",
  full.names = TRUE
)

cat("---------------------------------------\n")
cat("Galaxy featureCounts Merge Pipeline\n")
cat("---------------------------------------\n")

cat("Files detected:", length(count_files), "\n")

if(length(count_files) == 0){
  
  stop("No .tabular files found in data/counts")
  
}

##############################
# Read First File
##############################

merged_counts <- read_tsv(
  count_files[1],
  show_col_types = FALSE
)

##############################
# Validate First File
##############################

if(ncol(merged_counts) != 2){
  
  stop("First file does not contain exactly 2 columns.")
  
}

if(colnames(merged_counts)[1] != "Geneid"){
  
  stop("First column must be named 'Geneid'.")
  
}

##############################
# Merge Remaining Files
##############################

for(i in 2:length(count_files)){
  
  temp <- read_tsv(
    count_files[i],
    show_col_types = FALSE
  )
  
  if(ncol(temp) != 2){
    
    stop(
      paste(
        basename(count_files[i]),
        "does not contain exactly 2 columns."
      )
    )
    
  }
  
  merged_counts <- full_join(
    merged_counts,
    temp,
    by = "Geneid"
  )
  
}

##############################
# Check Duplicate Genes
##############################

duplicates <- duplicated(merged_counts$Geneid)

cat("Duplicate Gene IDs:", sum(duplicates), "\n")

##############################
# Missing Values
##############################

na_values <- sum(is.na(merged_counts))

cat("Missing Values:", na_values, "\n")

##############################
# Save Count Matrix
##############################

write_csv(
  merged_counts,
  file.path(output_dir,
            "merged_counts.csv")
)

##############################
# Save Merge Log
##############################

log_file <- file.path(
  log_dir,
  "merge_log.txt"
)

sink(log_file)

cat("AML Biomarker Discovery Pipeline\n\n")

cat("Merge Date:\n")

print(Sys.time())

cat("\n")

cat("Files merged:\n")

print(length(count_files))

cat("\n")

cat("Genes:\n")

print(nrow(merged_counts))

cat("\n")

cat("Samples:\n")

print(ncol(merged_counts)-1)

cat("\n")

cat("Duplicate Genes:\n")

print(sum(duplicates))

cat("\n")

cat("Missing Values:\n")

print(na_values)

sink()

##############################
# Summary
##############################

cat("\n")

cat("---------------------------------------\n")

cat("Merge Completed Successfully\n")

cat("---------------------------------------\n")

cat("Genes :", nrow(merged_counts), "\n")

cat("Samples :", ncol(merged_counts)-1, "\n")

cat("\n")

cat("Merged count matrix saved to:\n")

cat("data/processed/merged_counts.csv\n")

cat("\n")

cat("Log saved to:\n")

cat("results/logs/merge_log.txt\n")
