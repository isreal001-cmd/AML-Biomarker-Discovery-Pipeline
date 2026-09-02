##############################################################
# AML Biomarker Discovery Pipeline
#
# Script:
# 11_ROC_Analysis.R
#
# Purpose:
# Evaluate diagnostic performance of top biomarker
# candidates using Receiver Operating Characteristic (ROC)
# analysis.
#
# Author:
# Isreal Oluwafemi Abiodun
#
# Version:
# 1.1.0
##############################################################


##############################################################
# Load Configuration
##############################################################

source("config.R")
source("functions/plotting_functions.R")


##############################################################
# Header
##############################################################

cat("\n")
cat("=========================================\n")
cat("Stage 11 : ROC Analysis\n")
cat("=========================================\n\n")

start_time <- Sys.time()


##############################################################
# Directories
##############################################################

ROC_DIR <- file.path(
  RESULTS_DIR,
  "roc"
)

dir.create(
  ROC_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


##############################################################
# Input Files
##############################################################

VSD_FILE <- file.path(
  RESULTS_DIR,
  "vst",
  "vsd.rds"
)

METADATA_FILE <- file.path(
  "data",
  "metadata",
  "sample_metadata.csv"
)

BIOMARKER_FILE <- file.path(
  RESULTS_DIR,
  "biomarkers",
  "Top20_Biomarkers.csv"
)


##############################################################
# Validate Inputs
##############################################################

required_files <- c(
  VSD_FILE,
  METADATA_FILE,
  BIOMARKER_FILE
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  
  stop(
    paste(
      "Missing input file(s):",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
  
}


##############################################################
# Load Data
##############################################################

cat("Loading input files...\n")

vsd <- readRDS(
  VSD_FILE
)

metadata <- readr::read_csv(
  METADATA_FILE,
  show_col_types = FALSE
)

biomarkers <- readr::read_csv(
  BIOMARKER_FILE,
  show_col_types = FALSE
)

cat(
  "Candidate Biomarkers :",
  nrow(biomarkers),
  "\n\n"
)


##############################################################
# Expression Matrix
##############################################################

expr <- SummarizedExperiment::assay(
  vsd
)

cat(
  "Expression Matrix Dimensions :",
  nrow(expr),
  "genes x",
  ncol(expr),
  "samples\n\n"
)


##############################################################
# Validate Metadata
##############################################################

if (!"Condition" %in% names(metadata)) {
  
  stop(
    paste(
      "Condition column not found in metadata.",
      "\nAvailable columns:",
      paste(
        names(metadata),
        collapse = ", "
      )
    )
  )
  
}


##############################################################
# Validate Sample Alignment
##############################################################

if (
  ncol(expr) != nrow(metadata)
) {
  
  stop(
    paste(
      "Expression matrix and metadata have different",
      "numbers of samples.",
      "\nExpression samples:",
      ncol(expr),
      "\nMetadata samples:",
      nrow(metadata)
    )
  )
  
}


##############################################################
# Binary Phenotype
##############################################################

group <- ifelse(
  metadata$Condition == "AML",
  1,
  ifelse(
    metadata$Condition %in% c(
      "Healthy",
      "Non_AML",
      "Non AML"
    ),
    0,
    NA
  )
)


##############################################################
# Validate Phenotype
##############################################################

if (
  any(
    is.na(group)
  )
) {
  
  warning(
    paste(
      "Some samples have unrecognized phenotype values.",
      "\nValues detected:",
      paste(
        unique(metadata$Condition),
        collapse = ", "
      )
    )
  )
  
}


if (
  length(
    unique(
      group[
        !is.na(group)
      ]
    )
  ) < 2
) {
  
  stop(
    "ROC analysis requires both AML and non AML samples."
  )
  
}


##############################################################
# Check Biomarker Column
##############################################################

if (
  !"GeneID" %in% names(biomarkers)
) {
  
  stop(
    paste(
      "GeneID column not found in biomarker file.",
      "\nAvailable columns:",
      paste(
        names(biomarkers),
        collapse = ", "
      )
    )
  )
  
}


##############################################################
# Standardize Ensembl Identifiers
##############################################################

cat(
  "Standardizing Ensembl gene identifiers...\n"
)

expr_gene_ids <- rownames(expr)

expr_gene_ids_clean <- sub(
  "\\..*$",
  "",
  expr_gene_ids
)

biomarker_gene_ids_clean <- sub(
  "\\..*$",
  "",
  biomarkers$GeneID
)


##############################################################
# Match Biomarkers to Expression Matrix
##############################################################

biomarker_rows <- match(
  biomarker_gene_ids_clean,
  expr_gene_ids_clean
)


##############################################################
# Matching Summary
##############################################################

cat(
  "Candidate biomarkers :",
  length(biomarker_gene_ids_clean),
  "\n"
)

cat(
  "Matched biomarkers   :",
  sum(
    !is.na(
      biomarker_rows
    )
  ),
  "\n"
)

cat(
  "Unmatched biomarkers :",
  sum(
    is.na(
      biomarker_rows
    )
  ),
  "\n\n"
)


##############################################################
# Stop if No Biomarkers Match
##############################################################

if (
  all(
    is.na(
      biomarker_rows
    )
  )
) {
  
  stop(
    paste(
      "No biomarkers matched the VST expression matrix.",
      "\nThis may indicate an identifier mismatch."
    )
  )
  
}


##############################################################
# ROC Analysis
##############################################################

roc_results <- data.frame(
  GeneID = character(),
  GeneSymbol = character(),
  AUC = numeric(),
  stringsAsFactors = FALSE
)

cat(
  "Running ROC analysis...\n\n"
)


for (
  i in seq_along(
    biomarkers$GeneID
  )
) {
  
  ############################################################
  # Skip unmatched biomarkers
  ############################################################
  
  if (
    is.na(
      biomarker_rows[i]
    )
  ) {
    
    warning(
      paste(
        biomarkers$GeneID[i],
        "not found in expression matrix."
      )
    )
    
    next
    
  }
  
  
  ############################################################
  # Get actual expression matrix row
  ############################################################
  
  expression_row <- biomarker_rows[i]
  
  
  ############################################################
  # Extract expression values
  ############################################################
  
  values <- as.numeric(
    expr[
      expression_row,
    ]
  )
  
  
  ############################################################
  # Validate Expression Values
  ############################################################
  
  if (
    length(values) != length(group)
  ) {
    
    warning(
      paste(
        biomarkers$GeneID[i],
        "has an incorrect number of expression values."
      )
    )
    
    next
    
  }
  
  
  if (
    any(
      !is.finite(
        values
      )
    )
  ) {
    
    warning(
      paste(
        biomarkers$GeneID[i],
        "contains invalid expression values."
      )
    )
    
    next
    
  }
  
  
  ############################################################
  # Calculate ROC
  ############################################################
  
  roc_obj <- pROC::roc(
    response = group,
    predictor = values,
    quiet = TRUE
  )
  
  
  ############################################################
  # Calculate AUC
  ############################################################
  
  auc_value <- as.numeric(
    pROC::auc(
      roc_obj
    )
  )
  
  
  ############################################################
  # Get Gene Symbol if Available
  ############################################################
  
  gene_symbol <- NA_character_
  
  if (
    "GeneSymbol" %in% names(biomarkers)
  ) {
    
    gene_symbol <- as.character(
      biomarkers$GeneSymbol[i]
    )
    
  }
  
  
  ############################################################
  # Add Result
  ############################################################
  
  roc_results <- rbind(
    
    roc_results,
    
    data.frame(
      GeneID = biomarkers$GeneID[i],
      GeneSymbol = gene_symbol,
      AUC = auc_value,
      stringsAsFactors = FALSE
    )
    
  )
  
}


##############################################################
# Validate ROC Results
##############################################################

if (
  nrow(roc_results) == 0
) {
  
  stop(
    "ROC analysis produced no valid biomarker results."
  )
  
}


##############################################################
# Rank Biomarkers
##############################################################

roc_results <- roc_results |>
  
  dplyr::arrange(
    dplyr::desc(
      AUC
    )
  )


##############################################################
# Output Files
##############################################################

AUC_TABLE <- file.path(
  ROC_DIR,
  "ROC_AUC_Table.csv"
)

PNG_FILE <- file.path(
  ROC_DIR,
  "Best_Biomarker_ROC.png"
)

PDF_FILE <- file.path(
  ROC_DIR,
  "Best_Biomarker_ROC.pdf"
)


##############################################################
# Save ROC Table
##############################################################

readr::write_csv(
  roc_results,
  AUC_TABLE
)


##############################################################
# Best Biomarker
##############################################################

best_gene <- roc_results$GeneID[1]

best_gene_clean <- sub(
  "\\..*$",
  "",
  best_gene
)


##############################################################
# Locate Best Biomarker in Expression Matrix
##############################################################

best_index <- match(
  best_gene_clean,
  expr_gene_ids_clean
)

if (
  is.na(best_index)
) {
  
  stop(
    paste(
      "Best biomarker could not be located",
      "in the expression matrix."
    )
  )
  
}


##############################################################
# Extract Best Biomarker Expression
##############################################################

best_values <- as.numeric(
  expr[
    best_index,
  ]
)


##############################################################
# Calculate Best Biomarker ROC
##############################################################

best_roc <- pROC::roc(
  response = group,
  predictor = best_values,
  quiet = TRUE
)


##############################################################
# Best Biomarker AUC
##############################################################

best_auc <- as.numeric(
  pROC::auc(
    best_roc
  )
)


##############################################################
# Plot PNG
##############################################################

png(
  PNG_FILE,
  width = 1800,
  height = 1600,
  res = FIG_DPI
)

plot(
  best_roc,
  col = AML_RED,
  lwd = 3,
  main = paste(
    "ROC Curve:",
    best_gene
  )
)

abline(
  a = 0,
  b = 1,
  lty = 2
)

legend(
  "bottomright",
  legend = paste0(
    "AUC = ",
    round(
      best_auc,
      4
    )
  ),
  bty = "n"
)

dev.off()


##############################################################
# Plot PDF
##############################################################

pdf(
  PDF_FILE,
  width = 6,
  height = 6
)

plot(
  best_roc,
  col = AML_RED,
  lwd = 3,
  main = paste(
    "ROC Curve:",
    best_gene
  )
)

abline(
  a = 0,
  b = 1,
  lty = 2
)

legend(
  "bottomright",
  legend = paste0(
    "AUC = ",
    round(
      best_auc,
      4
    )
  ),
  bty = "n"
)

dev.off()


##############################################################
# Summary
##############################################################

end_time <- Sys.time()

cat("\n")
cat("=========================================\n")
cat("Stage 11 Completed Successfully\n")
cat("=========================================\n\n")

cat(
  "Candidate Biomarkers :",
  nrow(biomarkers),
  "\n"
)

cat(
  "Matched Biomarkers   :",
  sum(
    !is.na(
      biomarker_rows
    )
  ),
  "\n"
)

cat(
  "ROC Results          :",
  nrow(roc_results),
  "\n\n"
)

cat(
  "Best Biomarker       :",
  best_gene,
  "\n"
)

cat(
  "Best Gene Symbol     :",
  roc_results$GeneSymbol[1],
  "\n"
)

cat(
  "Highest AUC          :",
  round(
    best_auc,
    4
  ),
  "\n\n"
)

cat(
  "AUC Table            :",
  AUC_TABLE,
  "\n"
)

cat(
  "PNG Figure            :",
  PNG_FILE,
  "\n"
)

cat(
  "PDF Figure            :",
  PDF_FILE,
  "\n\n"
)

cat(
  "Time Elapsed         :",
  round(
    end_time - start_time,
    2
  ),
  "\n\n"
)