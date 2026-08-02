##############################################################
# AML Biomarker Discovery Pipeline
# Global Configuration
#
# Author  : Isreal Oluwafemi Abiodun
# Version : 1.0.0
##############################################################

##############################################################
# PROJECT INFORMATION
##############################################################

PIPELINE_NAME <- "AML Biomarker Discovery Pipeline"

PIPELINE_VERSION <- "1.0.0"

AUTHOR <- "Isreal Oluwafemi Abiodun"

PIPELINE_DATE <- Sys.Date()

R_VERSION <- R.version.string

##############################################################
# PROJECT DIRECTORIES
##############################################################

PROJECT_DIR <- getwd()

DATA_DIR <- file.path(PROJECT_DIR, "data")

RAW_DATA_DIR <- file.path(DATA_DIR, "raw")

PROCESSED_DATA_DIR <- file.path(DATA_DIR, "processed")

VALIDATION_DIR <- file.path(DATA_DIR, "validation")

RESULTS_DIR <- file.path(PROJECT_DIR, "results")

##############################################################
# OUTPUT DIRECTORIES
##############################################################

FIGURE_DIR <- file.path(RESULTS_DIR, "figures")

TABLE_DIR <- file.path(RESULTS_DIR, "tables")

LOG_DIR <- file.path(RESULTS_DIR, "logs")

##############################################################
# FIGURE DIRECTORIES
##############################################################

PCA_DIR <- file.path(FIGURE_DIR, "PCA")

VOLCANO_DIR <- file.path(FIGURE_DIR, "Volcano")

HEATMAP_DIR <- file.path(FIGURE_DIR, "Heatmaps")

ROC_DIR <- file.path(FIGURE_DIR, "ROC")

ENRICHMENT_DIR <- file.path(FIGURE_DIR, "Enrichment")

RANDOM_FOREST_DIR <- file.path(FIGURE_DIR, "Random_Forest")

VALIDATION_FIGURE_DIR <- file.path(FIGURE_DIR, "Validation")

##############################################################
# VALIDATION DIRECTORIES
##############################################################

VALIDATION_GEO_DIR <- file.path(VALIDATION_DIR, "GEO")

VALIDATION_TCGA_DIR <- file.path(VALIDATION_DIR, "TCGA")

VALIDATION_GTEX_DIR <- file.path(VALIDATION_DIR, "GTEx")

VALIDATION_USER_DIR <- file.path(VALIDATION_DIR, "User")

##############################################################
# CREATE DIRECTORIES
##############################################################

dirs <- c(
  
  RESULTS_DIR,
  
  FIGURE_DIR,
  
  TABLE_DIR,
  
  LOG_DIR,
  
  PCA_DIR,
  
  VOLCANO_DIR,
  
  HEATMAP_DIR,
  
  ROC_DIR,
  
  ENRICHMENT_DIR,
  
  RANDOM_FOREST_DIR,
  
  VALIDATION_FIGURE_DIR,
  
  VALIDATION_GEO_DIR,
  
  VALIDATION_TCGA_DIR,
  
  VALIDATION_GTEX_DIR,
  
  VALIDATION_USER_DIR
  
)

for (d in dirs) {
  
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
}

##############################################################
# DESEQ2 SETTINGS
##############################################################

PADJ_THRESHOLD <- 0.05

LOG2FC_THRESHOLD <- 1

##############################################################
# HEATMAP SETTINGS
##############################################################

TOP_HEATMAP_GENES <- 50

##############################################################
# RANDOM FOREST SETTINGS
##############################################################

RF_TOP_GENES <- 20

RF_NTREES <- 500

##############################################################
# RANDOM SEED
##############################################################

RANDOM_SEED <- 123

##############################################################
# VALIDATION OPTIONS
##############################################################

USE_GEO <- TRUE

USE_TCGA <- TRUE

USE_GTEX <- FALSE

USE_USER_DATASET <- FALSE

##############################################################
# PARALLEL COMPUTING
##############################################################

N_THREADS <- max(1, parallel::detectCores() - 1)

##############################################################
# PUBLICATION FIGURE SETTINGS
##############################################################

FIG_WIDTH <- 8

FIG_HEIGHT <- 6

FIG_DPI <- 600

BASE_FONT_SIZE <- 14

##############################################################
# COLOUR PALETTE
##############################################################

AML_RED <- "#D73027"

HEALTHY_BLUE <- "#4575B4"

GREY <- "#808080"

##############################################################
# EXPORT SETTINGS
##############################################################

EXPORT_PNG <- TRUE

EXPORT_PDF <- TRUE

EXPORT_SVG <- FALSE

##############################################################
# SESSION SETTINGS
##############################################################

options(stringsAsFactors = FALSE)

set.seed(RANDOM_SEED)

##############################################################
# STARTUP MESSAGE
##############################################################

cat("\n")
cat("=========================================\n")
cat(PIPELINE_NAME, "\n")
cat("Version :", PIPELINE_VERSION, "\n")
cat("Author  :", AUTHOR, "\n")
cat("R        :", R_VERSION, "\n")
cat("Date     :", as.character(PIPELINE_DATE), "\n")
cat("=========================================\n\n")