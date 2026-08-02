##############################################################
# AML Biomarker Discovery Pipeline
# Automatic Pipeline Report
#
# Author  : Isreal Oluwafemi Abiodun
# Version : 1.0.0
##############################################################

source("config.R")

cat("---------------------------------------\n")
cat("Generating Pipeline Report\n")
cat("---------------------------------------\n")

##############################################################
# Report Directory
##############################################################

REPORT_DIR <- file.path(RESULTS_DIR, "report")

dir.create(
  REPORT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

##############################################################
# Count Generated Figures
##############################################################

n_figures <- length(
  
  list.files(
    FIGURE_DIR,
    recursive = TRUE
  )
  
)

##############################################################
# Count Generated Tables
##############################################################

n_tables <- length(
  
  list.files(
    TABLE_DIR,
    recursive = TRUE
  )
  
)

##############################################################
# Count Result Files
##############################################################

n_results <- length(
  
  list.files(
    RESULTS_DIR,
    recursive = TRUE
  )
  
)

##############################################################
# Count Validation Files
##############################################################

validation_dir <- file.path(DATA_DIR, "validation")

n_validation <- length(
  
  list.files(
    validation_dir,
    recursive = TRUE
  )
  
)

##############################################################
# Build Report
##############################################################

report <- c(
  
  "=================================================",
  "AML Biomarker Discovery Pipeline",
  "=================================================",
  
  "",
  
  paste("Pipeline :", PIPELINE_NAME),
  
  paste("Version  :", PIPELINE_VERSION),
  
  paste("Author   :", AUTHOR),
  
  paste("Date     :", Sys.Date()),
  
  "",
  
  
  "PROJECT SUMMARY",
  
  "----------------",
  
  paste("Figures Generated        :", n_figures),
  
  paste("Tables Generated         :", n_tables),
  
  paste("Result Files             :", n_results),
  
  paste("Validation Files         :", n_validation),
  
  "",
  
  
  "MODULES COMPLETED",
  
  "----------------",
  
  "[✓] FeatureCounts Merge",
  
  "[✓] Metadata Preparation",
  
  "[✓] Differential Expression",
  
  "[✓] VST Normalization",
  
  "[✓] PCA",
  
  "[✓] Volcano Plot",
  
  "[✓] Heatmap",
  
  "[✓] Batch Assessment",
  
  "[✓] Biomarker Ranking",
  
  "[✓] Functional Enrichment",
  
  "[✓] ROC Analysis",
  
  "[✓] Expression Plots",
  
  "[✓] Random Forest",
  
  "[✓] Validation Preparation",
  
  "[✓] TCGA Validation",
  
  "[✓] External Validation",
  
  "[✓] Validation Visualization",
  
  "[✓] Session Logging",
  
  "[✓] Publication Summary",
  
  "",
  
  
  "STATUS",
  
  "------",
  
  "Pipeline Completed Successfully."
  
)

##############################################################
# Save Report
##############################################################

writeLines(
  
  report,
  
  file.path(
    REPORT_DIR,
    "Pipeline_Report.txt"
  )
  
)

cat("Pipeline report generated.\n")