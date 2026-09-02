##############################################################
# AML Biomarker Discovery Pipeline
#
# Script:
# 14C_External_Validation_ROC.R
#
# Purpose:
# External validation of discovered AML biomarkers using
# the independent GEO dataset GSE13159.
#
# Analysis:
# 1. Load discovered biomarkers
# 2. Read GSE13159 series matrix
# 3. Reconstruct GEO sample metadata
# 4. Identify AML and non AML samples
# 5. Extract GPL570 expression matrix
# 6. Annotate Affymetrix probes
# 7. Map discovered Ensembl biomarkers to probes
# 8. Calculate biomarker signature scores
# 9. Perform ROC analysis
# 10. Calculate AUC and 95% CI
# 11. Identify optimal threshold
# 12. Calculate sensitivity and specificity
# 13. Generate publication ready ROC plot
#
# Author:
# Isreal Oluwafemi Abiodun
#
# Version:
# 2.0.0
##############################################################


##############################################################
# Load Configuration
##############################################################

source("config.R")


##############################################################
# Required Packages
##############################################################

required_packages <- c(
  "readr",
  "dplyr",
  "stringr",
  "pROC",
  "AnnotationDbi",
  "hgu133plus2.db"
)


for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    stop(
      paste0(
        "Required package '",
        pkg,
        "' is not installed.\n",
        "Please install it before running this script."
      )
    )
    
  }
  
}


library(readr)
library(dplyr)
library(stringr)
library(pROC)
library(AnnotationDbi)
library(hgu133plus2.db)


##############################################################
# Header
##############################################################

cat("\n")
cat("=========================================\n")
cat("Stage 14C : External Validation ROC\n")
cat("=========================================\n\n")

start_time <- Sys.time()


##############################################################
# Directories
##############################################################

GEO_ROOT <- file.path(
  DATA_DIR,
  "validation",
  "GEO"
)


GEO_RAW_DIR <- file.path(
  GEO_ROOT,
  "raw"
)


GEO_PROCESSED_DIR <- file.path(
  GEO_ROOT,
  "processed"
)


GEO_METADATA_DIR <- file.path(
  GEO_ROOT,
  "metadata"
)


VALIDATION_DIR <- file.path(
  RESULTS_DIR,
  "validation"
)


ROC_DIR <- file.path(
  VALIDATION_DIR,
  "ROC"
)


dir.create(
  GEO_PROCESSED_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  GEO_METADATA_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  ROC_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


##############################################################
# Input Files
##############################################################

BIOMARKER_FILE <- file.path(
  RESULTS_DIR,
  "biomarkers",
  "Ranked_Biomarkers.csv"
)


GEO_FILE <- file.path(
  GEO_RAW_DIR,
  "GSE13159_series_matrix.txt.gz"
)


if (!file.exists(BIOMARKER_FILE)) {
  
  stop(
    paste0(
      "Ranked biomarker file not found:\n",
      BIOMARKER_FILE
    )
  )
  
}


if (!file.exists(GEO_FILE)) {
  
  stop(
    paste0(
      "GSE13159 GEO file not found:\n",
      GEO_FILE
    )
  )
  
}


##############################################################
# Output Files
##############################################################

METADATA_RAW_FILE <- file.path(
  GEO_METADATA_DIR,
  "GSE13159_sample_metadata_raw.csv"
)


METADATA_PROCESSED_FILE <- file.path(
  GEO_PROCESSED_DIR,
  "GSE13159_sample_metadata_processed.csv"
)


EXPRESSION_FILE <- file.path(
  GEO_PROCESSED_DIR,
  "GSE13159_expression_matrix.csv"
)


MAPPING_FILE <- file.path(
  GEO_PROCESSED_DIR,
  "GSE13159_biomarker_probe_mapping.csv"
)


SIGNATURE_FILE <- file.path(
  VALIDATION_DIR,
  "GSE13159_biomarker_signature_scores.csv"
)


ROC_RESULTS_FILE <- file.path(
  ROC_DIR,
  "GSE13159_ROC_results.csv"
)


ROC_PLOT_FILE <- file.path(
  ROC_DIR,
  "GSE13159_ROC_curve.png"
)


SUMMARY_FILE <- file.path(
  VALIDATION_DIR,
  "Stage14C_External_Validation_Summary.txt"
)


##############################################################
# Load Discovered Biomarkers
##############################################################

cat("Loading discovered AML biomarkers...\n")

biomarkers <- readr::read_csv(
  BIOMARKER_FILE,
  show_col_types = FALSE
)


cat(
  "Biomarkers Loaded :",
  nrow(biomarkers),
  "\n"
)


if (!"GeneID" %in% colnames(biomarkers)) {
  
  stop(
    "GeneID column was not found in Ranked_Biomarkers.csv."
  )
  
}


##############################################################
# Clean Ensembl IDs
##############################################################

biomarkers <- biomarkers |>
  
  dplyr::filter(
    !is.na(GeneID)
  ) |>
  
  dplyr::mutate(
    
    GeneID = sub(
      "\\..*$",
      "",
      GeneID
    )
    
  )


biomarkers <- biomarkers |>
  
  dplyr::distinct(
    GeneID,
    .keep_all = TRUE
  )


cat(
  "Unique Biomarkers :",
  nrow(biomarkers),
  "\n\n"
)


##############################################################
# Read GEO File
##############################################################

cat("Reading GSE13159 expression matrix...\n")


geo_lines <- readLines(
  gzfile(GEO_FILE),
  warn = FALSE
)


cat(
  "Total GEO lines :",
  length(geo_lines),
  "\n\n"
)


##############################################################
# Reconstruct GEO Metadata
##############################################################

cat("Reconstructing GEO sample metadata...\n")


accession_line <- geo_lines[
  grep(
    "^!Sample_geo_accession",
    geo_lines
  )[1]
]


title_line <- geo_lines[
  grep(
    "^!Sample_title",
    geo_lines
  )[1]
]


platform_line <- geo_lines[
  grep(
    "^!Sample_platform_id",
    geo_lines
  )[1]
]


characteristic_lines <- geo_lines[
  grep(
    "^!Sample_characteristics_ch1",
    geo_lines
  )
]


if (length(accession_line) == 0) {
  
  stop(
    "Could not locate GEO sample accession metadata."
  )
  
}


if (length(title_line) == 0) {
  
  stop(
    "Could not locate GEO sample title metadata."
  )
  
}


if (length(platform_line) == 0) {
  
  stop(
    "Could not locate GEO platform metadata."
  )
  
}


##############################################################
# Helper Function
##############################################################

parse_geo_line <- function(line) {
  
  fields <- strsplit(
    line,
    "\t",
    fixed = TRUE
  )[[1]]
  
  fields <- gsub(
    '^"|"$',
    "",
    fields
  )
  
  fields
  
}


##############################################################
# Parse Metadata
##############################################################

accessions <- parse_geo_line(
  accession_line
)


titles <- parse_geo_line(
  title_line
)


platforms <- parse_geo_line(
  platform_line
)


accessions <- accessions[
  -1
]


titles <- titles[
  -1
]


platforms <- platforms[
  -1
]


cat(
  "GEO Samples Detected :",
  length(accessions),
  "\n"
)


cat(
  "Sample Titles       :",
  length(titles),
  "\n"
)


cat(
  "Platform Entries    :",
  length(platforms),
  "\n"
)


cat(
  "Characteristic Rows :",
  length(characteristic_lines),
  "\n\n"
)


##############################################################
# Build Metadata Table
##############################################################

metadata <- data.frame(
  Sample = accessions,
  Title = titles,
  Platform = platforms,
  stringsAsFactors = FALSE
)


##############################################################
# Add Characteristics
##############################################################

if (length(characteristic_lines) > 0) {
  
  for (i in seq_along(characteristic_lines)) {
    
    fields <- parse_geo_line(
      characteristic_lines[i]
    )
    
    values <- fields[
      -1
    ]
    
    key <- paste0(
      "Characteristic_",
      i
    )
    
    metadata[[key]] <- values
    
  }
  
}


##############################################################
# Validate Metadata
##############################################################

metadata <- metadata |>
  
  dplyr::filter(
    !is.na(Sample),
    Sample != ""
  )


cat(
  "Metadata Samples :",
  nrow(metadata),
  "\n\n"
)


##############################################################
# Save Raw Metadata
##############################################################

readr::write_csv(
  metadata,
  METADATA_RAW_FILE
)


cat(
  "Metadata saved   :",
  METADATA_RAW_FILE,
  "\n\n"
)


##############################################################
# Locate Expression Matrix
##############################################################

cat("Locating expression matrix...\n")


header_candidates <- grep(
  "ID_REF",
  geo_lines,
  fixed = TRUE
)


header_index <- NA_integer_


for (idx in header_candidates) {
  
  current_line <- geo_lines[idx]
  
  fields <- strsplit(
    current_line,
    "\t",
    fixed = TRUE
  )[[1]]
  
  
  if (length(fields) >= 2) {
    
    first_field <- gsub(
      '"',
      "",
      trimws(fields[1])
    )
    
    
    if (
      identical(
        first_field,
        "ID_REF"
      )
    ) {
      
      header_index <- idx
      
      break
      
    }
    
  }
  
}


if (is.na(header_index)) {
  
  stop(
    paste0(
      "Expression matrix header ID_REF not found.\n",
      "GEO metadata was successfully reconstructed, ",
      "but the expression matrix could not be located."
    )
  )
  
}


cat(
  "Expression matrix header line :",
  header_index,
  "\n\n"
)


##############################################################
# Extract Expression Section
##############################################################

expression_lines <- geo_lines[
  header_index:length(geo_lines)
]


##############################################################
# Remove GEO Comment Lines
##############################################################

expression_lines <- expression_lines[
  !grepl(
    "^!",
    expression_lines
  )
]


##############################################################
# Remove Empty Lines
##############################################################

expression_lines <- expression_lines[
  nzchar(
    trimws(expression_lines)
  )
]


##############################################################
# Read Expression Matrix
##############################################################

expression_text <- paste(
  expression_lines,
  collapse = "\n"
)


expression_connection <- textConnection(
  expression_text
)


expression_data <- read.delim(
  expression_connection,
  header = TRUE,
  sep = "\t",
  quote = "\"",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)


close(
  expression_connection
)


##############################################################
# Validate Expression Matrix
##############################################################

probe_column <- colnames(
  expression_data
)[1]


cat(
  "Probe identifier column :",
  probe_column,
  "\n"
)


cat(
  "Probes detected  :",
  nrow(expression_data),
  "\n"
)


cat(
  "Samples detected :",
  ncol(expression_data) - 1,
  "\n\n"
)


if (
  probe_column != "ID_REF"
) {
  
  stop(
    paste0(
      "Unexpected expression matrix identifier column: ",
      probe_column
    )
  )
  
}


##############################################################
# Match Expression Samples to Metadata
##############################################################

cat(
  "Matching expression samples to GEO metadata...\n"
)


expression_samples <- colnames(
  expression_data
)[-1]


common_samples <- intersect(
  expression_samples,
  metadata$Sample
)


cat(
  "Common samples :",
  length(common_samples),
  "\n"
)


if (
  length(common_samples) < 100
) {
  
  stop(
    paste0(
      "Too few expression samples matched GEO metadata.\n",
      "Matched samples: ",
      length(common_samples)
    )
  )
  
}


##############################################################
# Reorder Expression Matrix
##############################################################

expression_data <- expression_data[
  ,
  c(
    "ID_REF",
    common_samples
  ),
  drop = FALSE
]


metadata <- metadata[
  match(
    common_samples,
    metadata$Sample
  ),
  ,
  drop = FALSE
]


rownames(metadata) <- metadata$Sample


cat(
  "Sample matching successful.\n\n"
)


##############################################################
# Identify Leukemia Classification
##############################################################

cat(
  "Identifying leukemia classification...\n\n"
)


characteristic_columns <- grep(
  "Characteristic_",
  colnames(metadata)
)


if (
  length(characteristic_columns) == 0
) {
  
  stop(
    "No GEO characteristic columns were found."
  )
  
}


leukemia_column <- NA_character_


for (col in colnames(metadata)[characteristic_columns]) {
  
  values <- metadata[[col]]
  
  if (
    any(
      grepl(
        "leukemia class",
        values,
        ignore.case = TRUE
      ),
      na.rm = TRUE
    )
  ) {
    
    leukemia_column <- col
    
    break
    
  }
  
}


if (is.na(leukemia_column)) {
  
  stop(
    "Could not identify the GEO leukemia class column."
  )
  
}


metadata$Leukemia_Class <- metadata[[
  leukemia_column
]]


##############################################################
# Clean Classification Labels
##############################################################

metadata$Leukemia_Class <- gsub(
  '^"|"$',
  "",
  metadata$Leukemia_Class
)


metadata$Leukemia_Class <- trimws(
  metadata$Leukemia_Class
)


##############################################################
# Identify AML Samples
##############################################################

metadata$AML_Status <- ifelse(
  
  grepl(
    "^leukemia class: AML",
    metadata$Leukemia_Class,
    ignore.case = TRUE
  ),
  
  "AML",
  
  "Non_AML"
  
)


##############################################################
# Classification Summary
##############################################################

classification_summary <- sort(
  table(
    metadata$Leukemia_Class
  ),
  decreasing = TRUE
)


cat(
  "Leukemia classes detected:\n\n"
)


print(
  classification_summary
)


cat("\n")


cat(
  "AML samples     :",
  sum(
    metadata$AML_Status == "AML"
  ),
  "\n"
)


cat(
  "Non AML samples :",
  sum(
    metadata$AML_Status == "Non_AML"
  ),
  "\n\n"
)


##############################################################
# Validate Class Counts
##############################################################

aml_count <- sum(
  metadata$AML_Status == "AML"
)


non_aml_count <- sum(
  metadata$AML_Status == "Non_AML"
)


if (
  aml_count < 10
) {
  
  stop(
    paste0(
      "Fewer than 10 AML samples were identified.\n",
      "AML samples: ",
      aml_count
    )
  )
  
}


if (
  non_aml_count < 10
) {
  
  stop(
    paste0(
      "Fewer than 10 non AML samples were identified.\n",
      "Non AML samples: ",
      non_aml_count
    )
  )
  
}


##############################################################
# Save Processed Metadata
##############################################################

readr::write_csv(
  metadata,
  METADATA_PROCESSED_FILE
)


cat(
  "Processed metadata saved :",
  METADATA_PROCESSED_FILE,
  "\n\n"
)


##############################################################
# Save Expression Matrix
##############################################################

readr::write_csv(
  expression_data,
  EXPRESSION_FILE
)


##############################################################
# Probe Annotation
##############################################################

cat(
  "Annotating GPL570 probes...\n"
)


probe_ids <- expression_data$ID_REF


probe_annotation <- AnnotationDbi::select(
  
  hgu133plus2.db,
  
  keys = probe_ids,
  
  keytype = "PROBEID",
  
  columns = c(
    "ENTREZID",
    "ENSEMBL",
    "SYMBOL"
  )
  
)


probe_annotation <- probe_annotation |>
  
  dplyr::filter(
    !is.na(ENSEMBL)
  ) |>
  
  dplyr::distinct(
    PROBEID,
    ENSEMBL,
    .keep_all = TRUE
  )


cat(
  "Annotated probes :",
  nrow(probe_annotation),
  "\n\n"
)


##############################################################
# Map Biomarkers to Probes
##############################################################

cat(
  "Mapping discovered biomarkers to GPL570 probes...\n"
)


biomarker_probe_map <- biomarkers |>
  
  dplyr::inner_join(
    probe_annotation,
    by = c(
      "GeneID" = "ENSEMBL"
    )
  )


biomarker_probe_map <- biomarker_probe_map |>
  
  dplyr::distinct(
    GeneID,
    PROBEID,
    .keep_all = TRUE
  )


cat(
  "Biomarkers mapped to probes :",
  nrow(
    biomarker_probe_map
  ),
  "\n"
)


cat(
  "Unique biomarkers mapped    :",
  dplyr::n_distinct(
    biomarker_probe_map$GeneID
  ),
  "\n"
)


cat(
  "Unique probes mapped        :",
  dplyr::n_distinct(
    biomarker_probe_map$PROBEID
  ),
  "\n\n"
)


if (
  nrow(biomarker_probe_map) < 10
) {
  
  stop(
    paste0(
      "Fewer than 10 discovered biomarkers could be mapped ",
      "to GPL570 probes."
    )
  )
  
}


##############################################################
# Save Probe Mapping
##############################################################

readr::write_csv(
  biomarker_probe_map,
  MAPPING_FILE
)


##############################################################
# Construct Biomarker Expression Matrix
##############################################################

cat(
  "Extracting biomarker expression profiles...\n"
)


mapped_probes <- unique(
  biomarker_probe_map$PROBEID
)


biomarker_expression <- expression_data |>
  
  dplyr::filter(
    ID_REF %in% mapped_probes
  )


##############################################################
# Convert to Numeric Matrix
##############################################################

expression_matrix <- as.matrix(
  biomarker_expression[
    ,
    -1,
    drop = FALSE
  ]
)


mode(
  expression_matrix
) <- "numeric"


rownames(
  expression_matrix
) <- biomarker_expression$ID_REF


##############################################################
# Z Score Normalization
##############################################################

cat(
  "Normalizing biomarker expression...\n"
)


expression_matrix_z <- t(
  scale(
    t(
      expression_matrix
    )
  )
)


##############################################################
# Handle Missing Values
##############################################################

expression_matrix_z[
  is.na(expression_matrix_z)
] <- 0


##############################################################
# Biomarker Direction
##############################################################

biomarker_direction <- biomarker_probe_map |>
  
  dplyr::select(
    GeneID,
    PROBEID,
    log2FoldChange
  ) |>
  
  dplyr::distinct(
    PROBEID,
    .keep_all = TRUE
  )


##############################################################
# Calculate Signature Score
##############################################################

cat(
  "Calculating biomarker signature scores...\n"
)


signature_scores <- numeric(
  ncol(expression_matrix_z)
)


for (
  i in seq_len(
    ncol(expression_matrix_z)
  )
) {
  
  sample_values <- expression_matrix_z[
    ,
    i
  ]
  
  
  probe_weights <- biomarker_direction$log2FoldChange
  
  
  weighted_values <-
    sample_values[
      biomarker_direction$PROBEID
    ] *
    sign(
      probe_weights
    )
  
  
  signature_scores[i] <- mean(
    weighted_values,
    na.rm = TRUE
  )
  
}


##############################################################
# Signature Data Frame
##############################################################

signature_data <- data.frame(
  
  Sample = colnames(
    expression_matrix_z
  ),
  
  Signature_Score = signature_scores,
  
  stringsAsFactors = FALSE
  
)


signature_data <- signature_data |>
  
  dplyr::left_join(
    metadata |>
      dplyr::select(
        Sample,
        AML_Status,
        Leukemia_Class
      ),
    by = "Sample"
  )


##############################################################
# Validate Signature
##############################################################

if (
  any(
    is.na(
      signature_data$AML_Status
    )
  )
) {
  
  stop(
    "Some expression samples do not have AML classification."
  )
  
}


##############################################################
# Save Signature Scores
##############################################################

readr::write_csv(
  signature_data,
  SIGNATURE_FILE
)


##############################################################
# ROC Analysis
##############################################################

cat("\n")
cat(
  "Running ROC analysis...\n"
)


roc_object <- pROC::roc(
  
  response = signature_data$AML_Status,
  
  predictor = signature_data$Signature_Score,
  
  levels = c(
    "Non_AML",
    "AML"
  ),
  
  direction = "<",
  
  quiet = TRUE
  
)


##############################################################
# AUC
##############################################################

auc_value <- as.numeric(
  pROC::auc(
    roc_object
  )
)


auc_ci <- pROC::ci.auc(
  roc_object,
  conf.level = 0.95
)


##############################################################
# Optimal Threshold
##############################################################

optimal_coords <- pROC::coords(
  
  roc_object,
  
  x = "best",
  
  best.method = "youden",
  
  ret = c(
    "threshold",
    "sensitivity",
    "specificity",
    "accuracy"
  ),
  
  transpose = FALSE
  
)


##############################################################
# ROC Results
##############################################################

roc_results <- data.frame(
  
  Metric = c(
    "AUC",
    "AUC_95CI_Lower",
    "AUC_95CI_Upper",
    "Optimal_Threshold",
    "Sensitivity",
    "Specificity",
    "Accuracy"
  ),
  
  Value = c(
    
    auc_value,
    
    as.numeric(
      auc_ci[1]
    ),
    
    as.numeric(
      auc_ci[3]
    ),
    
    as.numeric(
      optimal_coords$threshold
    ),
    
    as.numeric(
      optimal_coords$sensitivity
    ),
    
    as.numeric(
      optimal_coords$specificity
    ),
    
    as.numeric(
      optimal_coords$accuracy
    )
    
  )
  
)


##############################################################
# Save ROC Results
##############################################################

readr::write_csv(
  roc_results,
  ROC_RESULTS_FILE
)


##############################################################
# ROC Plot
##############################################################

png(
  
  ROC_PLOT_FILE,
  
  width = 2200,
  
  height = 1800,
  
  res = 300
  
)


plot(
  
  roc_object,
  
  col = "black",
  
  lwd = 3,
  
  main = paste0(
    "External Validation of AML Biomarker Signature\n",
    "GSE13159"
  ),
  
  legacy.axes = TRUE
  
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
      auc_value,
      3
    ),
    " (95% CI ",
    round(
      auc_ci[1],
      3
    ),
    "–",
    round(
      auc_ci[3],
      3
    ),
    ")"
  ),
  
  bty = "n"
  
)


dev.off()


##############################################################
# Validation Summary
##############################################################

summary_lines <- c(
  
  "=========================================",
  
  "Stage 14C External Validation Summary",
  
  "=========================================",
  
  "",
  
  paste(
    "Dataset : GSE13159"
  ),
  
  paste(
    "Total Samples :",
    nrow(
      signature_data
    )
  ),
  
  paste(
    "AML Samples :",
    aml_count
  ),
  
  paste(
    "Non AML Samples :",
    non_aml_count
  ),
  
  "",
  
  paste(
    "Discovered Biomarkers :",
    nrow(
      biomarkers
    )
  ),
  
  paste(
    "Mapped Biomarkers :",
    dplyr::n_distinct(
      biomarker_probe_map$GeneID
    )
  ),
  
  paste(
    "Mapped Probes :",
    dplyr::n_distinct(
      biomarker_probe_map$PROBEID
    )
  ),
  
  "",
  
  paste(
    "AUC :",
    round(
      auc_value,
      4
    )
  ),
  
  paste(
    "AUC 95% CI :",
    round(
      auc_ci[1],
      4
    ),
    "to",
    round(
      auc_ci[3],
      4
    )
  ),
  
  paste(
    "Optimal Threshold :",
    round(
      optimal_coords$threshold,
      4
    )
  ),
  
  paste(
    "Sensitivity :",
    round(
      optimal_coords$sensitivity,
      4
    )
  ),
  
  paste(
    "Specificity :",
    round(
      optimal_coords$specificity,
      4
    )
  ),
  
  paste(
    "Accuracy :",
    round(
      optimal_coords$accuracy,
      4
    )
  ),
  
  "",
  
  paste(
    "Signature Scores :",
    SIGNATURE_FILE
  ),
  
  paste(
    "ROC Results :",
    ROC_RESULTS_FILE
  ),
  
  paste(
    "ROC Plot :",
    ROC_PLOT_FILE
  ),
  
  "",
  
  paste(
    "Completed :",
    Sys.time()
  )
  
)


writeLines(
  summary_lines,
  SUMMARY_FILE
)


##############################################################
# Final Console Summary
##############################################################

cat("\n")
cat("=========================================\n")
cat("Stage 14C Completed Successfully\n")
cat("=========================================\n\n")


cat(
  "Dataset           : GSE13159\n"
)


cat(
  "Total Samples     :",
  nrow(signature_data),
  "\n"
)


cat(
  "AML Samples       :",
  aml_count,
  "\n"
)


cat(
  "Non AML Samples   :",
  non_aml_count,
  "\n"
)


cat(
  "Mapped Biomarkers :",
  dplyr::n_distinct(
    biomarker_probe_map$GeneID
  ),
  "\n"
)


cat(
  "Mapped Probes     :",
  dplyr::n_distinct(
    biomarker_probe_map$PROBEID
  ),
  "\n\n"
)


cat(
  "AUC               :",
  round(
    auc_value,
    4
  ),
  "\n"
)


cat(
  "AUC 95% CI        :",
  round(
    auc_ci[1],
    4
  ),
  "to",
  round(
    auc_ci[3],
    4
  ),
  "\n"
)


cat(
  "Sensitivity       :",
  round(
    optimal_coords$sensitivity,
    4
  ),
  "\n"
)


cat(
  "Specificity       :",
  round(
    optimal_coords$specificity,
    4
  ),
  "\n\n"
)


cat(
  "ROC Plot          :",
  ROC_PLOT_FILE,
  "\n"
)


cat(
  "Summary            :",
  SUMMARY_FILE,
  "\n\n"
)


cat(
  "Time Elapsed      :",
  round(
    Sys.time() - start_time,
    2
  ),
  "\n\n"
)

list.files(
  "C:/Users/USER/Bioinformatics/AML_Project/results/validation",
  recursive = TRUE,
  full.names = TRUE
)
