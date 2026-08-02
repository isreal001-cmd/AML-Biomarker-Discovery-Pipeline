###############################################################
# AML Biomarker Discovery Pipeline
# Script: 14B_Download_TCGA_Data.R
#
# Purpose:
# Convert downloaded TCGA-LAML STAR count files into a
# standardized validation dataset.
#
# Input:
# GDCdata/TCGA-LAML/Transcriptome_Profiling/
# Gene_Expression_Quantification/
#
# Output:
# data/validation/TCGA/
# ├── expression.csv
# └── metadata.csv
###############################################################

# rm(list = ls())

###############################################################
# Load Configuration
###############################################################

source("config.R")

###############################################################
# Required Packages
###############################################################

packages <- c(
  "readr",
  "dplyr",
  "stringr",
  "purrr",
  "tibble"
)

for(pkg in packages){
  
  if(!require(pkg, character.only = TRUE)){
    
    install.packages(pkg)
    
    library(pkg, character.only = TRUE)
    
  }
  
}

###############################################################
# Input Folder
###############################################################

input_dir <-
  "GDCdata/TCGA-LAML/Transcriptome_Profiling/Gene_Expression_Quantification"

###############################################################
# Locate STAR Count Files
###############################################################

count_files <- list.files(
  
  input_dir,
  
  pattern = "augmented_star_gene_counts.tsv$",
  
  recursive = TRUE,
  
  full.names = TRUE
  
)

cat("---------------------------------------\n")
cat("TCGA STAR Files Located\n")
cat("---------------------------------------\n")
cat("Files:", length(count_files), "\n\n")

if(length(count_files) == 0){
  
  stop("No STAR count files found.")
  
}

###############################################################
# Read First File
###############################################################

first <- read_tsv(
  
  count_files[1],
  
  comment = "#",
  
  show_col_types = FALSE
  
)

###############################################################
# Keep Genes Only
###############################################################

first <- first %>%
  
  filter(str_detect(gene_id, "^ENSG"))

###############################################################
# Initialize Expression Matrix
###############################################################

expression <- tibble(
  
  GeneID = first$gene_id
  
)

###############################################################
# Initialize Metadata
###############################################################

metadata <- tibble(
  
  SampleID = character(),
  
  File = character()
  
)

###############################################################
# Read All Samples
###############################################################

cat("Reading expression files...\n\n")

for(i in seq_along(count_files)){
  
  file <- count_files[i]
  
  sample <- basename(dirname(file))
  
  cat(i, "/", length(count_files), ":", sample, "\n")
  
  dat <- read_tsv(
    
    file,
    
    comment = "#",
    
    show_col_types = FALSE
    
  )
  
  dat <- dat %>%
    
    filter(str_detect(gene_id, "^ENSG"))
  
  expression[[sample]] <- dat$unstranded
  
  metadata <- bind_rows(
    
    metadata,
    
    tibble(
      
      SampleID = sample,
      
      File = basename(file)
      
    )
    
  )
  
}

###############################################################
# Save Files
###############################################################

dir.create(
  
  "data/validation/TCGA",
  
  recursive = TRUE,
  
  showWarnings = FALSE
  
)

write_csv(
  
  expression,
  
  "data/validation/TCGA/expression.csv"
  
)

write_csv(
  
  metadata,
  
  "data/validation/TCGA/metadata.csv"
  
)

###############################################################
# Summary
###############################################################

cat("\n---------------------------------------\n")
cat("TCGA Dataset Prepared Successfully\n")
cat("---------------------------------------\n\n")

cat("Genes   :", nrow(expression), "\n")
cat("Samples :", ncol(expression) - 1, "\n\n")

cat("Files created:\n")
cat("data/validation/TCGA/expression.csv\n")
cat("data/validation/TCGA/metadata.csv\n")