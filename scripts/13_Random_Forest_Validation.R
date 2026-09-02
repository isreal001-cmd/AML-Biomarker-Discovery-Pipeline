##############################################################
# AML Biomarker Discovery Pipeline
#
# Script:
# 13_Random_Forest_Validation.R
#
# Purpose:
# Validate the predictive performance of the Top 20
# biomarker candidates using a Random Forest classifier.
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
source("functions/plotting_functions.R")

cat("\n")
cat("=========================================\n")
cat("Stage 13 : Random Forest Validation\n")
cat("=========================================\n\n")

start_time <- Sys.time()

##############################################################
# Reproducibility
##############################################################

set.seed(12345)

##############################################################
# Directories
##############################################################

RF_DIR <- file.path(
  RESULTS_DIR,
  "random_forest"
)

dir.create(
  RF_DIR,
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

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  
  stop(
    paste(
      "Missing input file(s):",
      paste(missing_files, collapse = "\n")
    )
  )
  
}

##############################################################
# Load Data
##############################################################

cat("Loading input files...\n")

vsd <- readRDS(VSD_FILE)

metadata <- readr::read_csv(
  METADATA_FILE,
  show_col_types = FALSE
)

biomarkers <- readr::read_csv(
  BIOMARKER_FILE,
  show_col_types = FALSE
)

##############################################################
# Expression Matrix
##############################################################

expr <- SummarizedExperiment::assay(vsd)

##############################################################
# Standardise Ensembl IDs
##############################################################

clean_ensembl_id <- function(x) {
  
  x <- as.character(x)
  
  x <- sub(
    "\\..*$",
    "",
    x
  )
  
  x
  
}

##############################################################
# Normalise Expression Matrix Row IDs
##############################################################

expr_ids_original <- rownames(expr)

expr_ids_clean <- clean_ensembl_id(
  expr_ids_original
)

##############################################################
# Normalise Biomarker IDs
##############################################################

biomarker_ids_original <- biomarkers$GeneID

biomarker_ids_clean <- clean_ensembl_id(
  biomarker_ids_original
)

##############################################################
# Match Biomarkers to Expression Matrix
##############################################################

match_index <- match(
  biomarker_ids_clean,
  expr_ids_clean
)

matched <- !is.na(match_index)

cat("\n")
cat("Biomarker Mapping Summary\n")
cat("-------------------------\n")

cat(
  "Candidate biomarkers :",
  length(biomarker_ids_original),
  "\n"
)

cat(
  "Matched biomarkers   :",
  sum(matched),
  "\n"
)

cat(
  "Unmatched biomarkers :",
  sum(!matched),
  "\n\n"
)

##############################################################
# Report Unmatched Biomarkers
##############################################################

if (any(!matched)) {
  
  cat("Unmatched biomarker IDs:\n")
  
  print(
    biomarker_ids_original[!matched]
  )
  
  cat("\n")
  
}

##############################################################
# Stop if Nothing Matches
##############################################################

if (sum(matched) == 0) {
  
  stop(
    paste(
      "None of the biomarker genes matched the VST",
      "expression matrix after Ensembl ID normalisation."
    )
  )
  
}

##############################################################
# Extract Matched Expression Data
##############################################################

valid_genes <- expr_ids_original[
  match_index[matched]
]

expr <- expr[
  valid_genes,
  ,
  drop = FALSE
]

##############################################################
# Transpose for Random Forest
##############################################################

expr <- t(expr)

rf_data <- as.data.frame(expr)

##############################################################
# Add Phenotype
##############################################################

rf_data$Condition <- factor(
  metadata$Condition
)

cat(
  "Biomarkers Used :",
  ncol(rf_data) - 1,
  "\n\n"
)

##############################################################
# Train Random Forest
##############################################################

cat("Training Random Forest model...\n\n")

rf_model <- randomForest::randomForest(
  
  Condition ~ .,
  
  data = rf_data,
  
  importance = TRUE,
  
  ntree = 1000
  
)

cat("Random Forest training completed.\n\n")

##############################################################
# Variable Importance
##############################################################

importance_matrix <- randomForest::importance(rf_model)

importance_table <- data.frame(
  
  GeneID = rownames(importance_matrix),
  
  MeanDecreaseAccuracy =
    importance_matrix[, "MeanDecreaseAccuracy"],
  
  MeanDecreaseGini =
    importance_matrix[, "MeanDecreaseGini"]
  
)

importance_table <-
  
  importance_table |>
  
  dplyr::arrange(
    
    dplyr::desc(MeanDecreaseAccuracy)
    
  )

##############################################################
# Output Files
##############################################################

MODEL_FILE <- file.path(
  RF_DIR,
  "random_forest_model.rds"
)

IMPORTANCE_FILE <- file.path(
  RF_DIR,
  "Variable_Importance.csv"
)

PNG_FILE <- file.path(
  RF_DIR,
  "Variable_Importance.png"
)

PDF_FILE <- file.path(
  RF_DIR,
  "Variable_Importance.pdf"
)

##############################################################
# Save Outputs
##############################################################

saveRDS(
  
  rf_model,
  
  MODEL_FILE
  
)

readr::write_csv(
  
  importance_table,
  
  IMPORTANCE_FILE
  
)

##############################################################
# Plot Top 20
##############################################################

top20 <-
  
  dplyr::slice_head(
    
    importance_table,
    
    n = 20
    
  )

p <- ggplot2::ggplot(
  
  top20,
  
  ggplot2::aes(
    
    x = reorder(
      
      GeneID,
      
      MeanDecreaseAccuracy
      
    ),
    
    y = MeanDecreaseAccuracy
    
  )
  
) +
  
  ggplot2::geom_col(
    
    fill = AML_RED
    
  ) +
  
  ggplot2::coord_flip() +
  
  ggplot2::labs(
    
    title = "Random Forest Variable Importance",
    
    x = NULL,
    
    y = "Mean Decrease Accuracy"
    
  ) +
  
  publication_theme()

ggplot2::ggsave(
  
  filename = PNG_FILE,
  
  plot = p,
  
  width = FIG_WIDTH,
  
  height = FIG_HEIGHT,
  
  dpi = FIG_DPI
  
)

ggplot2::ggsave(
  
  filename = PDF_FILE,
  
  plot = p,
  
  width = FIG_WIDTH,
  
  height = FIG_HEIGHT
  
)

##############################################################
# Summary
##############################################################

end_time <- Sys.time()

cat("\n")
cat("=========================================\n")
cat("Stage 13 Completed Successfully\n")
cat("=========================================\n\n")

cat(
  "Top Biomarker :",
  top20$GeneID[1],
  "\n"
)

cat(
  "Mean Decrease Accuracy :",
  round(top20$MeanDecreaseAccuracy[1], 3),
  "\n\n"
)

cat(
  "Random Forest Model :",
  MODEL_FILE,
  "\n"
)

cat(
  "Importance Table    :",
  IMPORTANCE_FILE,
  "\n"
)

cat(
  "PNG Figure          :",
  PNG_FILE,
  "\n"
)

cat(
  "PDF Figure          :",
  PDF_FILE,
  "\n\n"
)

cat(
  "Time Elapsed :",
  round(end_time - start_time, 2),
  "\n\n"
)