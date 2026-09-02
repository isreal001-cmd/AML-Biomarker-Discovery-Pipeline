##############################################################
# AML Biomarker Discovery Pipeline
#
# Script:
# 14B_Download_TCGA_Data.R
#
# Purpose:
# Convert downloaded TCGA-LAML STAR count files into a
# standardized validation dataset.
##############################################################

##############################################################
# Load Configuration
##############################################################

source("config.R")

cat("\n")
cat("=========================================\n")
cat("Stage 14B : TCGA Validation Dataset Preparation\n")
cat("=========================================\n\n")

start_time <- Sys.time()

##############################################################
# Required Packages
##############################################################

required_packages <- c(
  "readr",
  "dplyr",
  "stringr",
  "purrr",
  "tibble"
)

install_and_load_packages(required_packages)

##############################################################
# Define Directories
##############################################################

INPUT_DIR <- file.path(
  "GDCdata",
  "TCGA-LAML",
  "Transcriptome_Profiling",
  "Gene_Expression_Quantification"
)

OUTPUT_DIR <- file.path(
  "data",
  "validation",
  "TCGA"
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

##############################################################
# Locate STAR Count Files
##############################################################

count_files <- list.files(
  path = INPUT_DIR,
  pattern = "augmented_star_gene_counts.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

cat("Searching for TCGA STAR Count Files...\n\n")
cat("Files Found :", length(count_files), "\n\n")

if(length(count_files) == 0){
  stop("No TCGA STAR count files were found.")
}

##############################################################
# Read First File
##############################################################

first <- readr::read_tsv(
  count_files[1],
  comment = "#",
  show_col_types = FALSE
)

cat("Detected columns:\n")
print(names(first))
cat("\n")

required_columns <- c(
  "gene_id",
  "unstranded"
)

missing_columns <- setdiff(
  required_columns,
  names(first)
)

if(length(missing_columns) > 0){
  
  stop(
    paste(
      "Missing required columns:",
      paste(missing_columns, collapse = ", ")
    )
  )
  
}

##############################################################
# Keep Gene Rows Only
##############################################################

first <- first %>%
  dplyr::filter(
    stringr::str_detect(
      gene_id,
      "^ENSG"
    )
  )

##############################################################
# Initialise Objects
##############################################################

expression <- tibble::tibble(
  GeneID = first$gene_id
)

metadata <- tibble::tibble(
  SampleID = character(),
  File = character()
)

##############################################################
# Read All Samples
##############################################################

cat("Reading Expression Files...\n\n")

for(i in seq_along(count_files)){
  
  file <- count_files[i]
  
  sample_id <- basename(dirname(file))
  
  cat(
    "[",
    i,
    "/",
    length(count_files),
    "] ",
    sample_id,
    "\n",
    sep = ""
  )
  
  dat <- readr::read_tsv(
    file,
    comment = "#",
    show_col_types = FALSE
  )
  
  ##########################################################
  # Validate Columns
  ##########################################################
  
  if(!all(required_columns %in% names(dat))){
    
    stop(
      paste(
        "Required columns missing in:",
        sample_id
      )
    )
    
  }
  
  ##########################################################
  # Keep Genes Only
  ##########################################################
  
  dat <- dat %>%
    dplyr::filter(
      stringr::str_detect(
        gene_id,
        "^ENSG"
      )
    )
  
  ##########################################################
  # Validate Gene Count
  ##########################################################
  
  if(nrow(dat) != nrow(expression)){
    
    stop(
      paste(
        "Gene count mismatch detected in:",
        sample_id
      )
    )
    
  }
  
  ##########################################################
  # Validate Gene Order
  ##########################################################
  
  if(!identical(
    expression$GeneID,
    dat$gene_id
  )){
    
    stop(
      paste(
        "Gene order mismatch detected in:",
        sample_id
      )
    )
    
  }
  
  ##########################################################
  # Append Counts
  ##########################################################
  
  expression <- dplyr::bind_cols(
    
    expression,
    
    tibble::tibble(
      !!sample_id := dat$unstranded
    )
    
  )
  
  ##########################################################
  # Metadata
  ##########################################################
  
  metadata <- dplyr::bind_rows(
    
    metadata,
    
    tibble::tibble(
      
      SampleID = sample_id,
      
      File = basename(file)
      
    )
    
  )
  
}

##############################################################
# Remove Duplicate GeneID Column
##############################################################

expression <- expression %>%
  dplyr::select(
    !matches("^GeneID\\.\\.\\.")
  )

##############################################################
# Save Files
##############################################################

expression_file <- file.path(
  OUTPUT_DIR,
  "expression.csv"
)

metadata_file <- file.path(
  OUTPUT_DIR,
  "metadata.csv"
)

readr::write_csv(
  expression,
  expression_file
)

readr::write_csv(
  metadata,
  metadata_file
)

##############################################################
# Summary
##############################################################

cat("\n")
cat("---------------------------------------\n")
cat("TCGA Dataset Prepared Successfully\n")
cat("---------------------------------------\n\n")

cat("Genes   :", nrow(expression), "\n")
cat("Samples :", ncol(expression) - 1, "\n\n")

cat("Expression Matrix :", expression_file, "\n")
cat("Metadata          :", metadata_file, "\n\n")

##############################################################
# Completion
##############################################################

end_time <- Sys.time()

cat("=========================================\n")
cat("Stage 14B Completed Successfully\n")
cat("=========================================\n\n")

cat(
  "Time Elapsed :",
  round(as.numeric(difftime(end_time, start_time, units = "secs")), 2),
  "seconds\n\n"
)