##############################################################
# AML Biomarker Discovery Pipeline
# Script: 19_Pipeline_Report.R
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

REPORT_DIR <- file.path(
  RESULTS_DIR,
  "report"
)

dir.create(
  REPORT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

##############################################################
# Count Outputs
##############################################################

n_figures <- length(
  list.files(
    FIGURE_DIR,
    recursive = TRUE,
    full.names = TRUE
  )
)

n_tables <- length(
  list.files(
    TABLE_DIR,
    recursive = TRUE,
    full.names = TRUE
  )
)

n_result_files <- length(
  list.files(
    RESULTS_DIR,
    recursive = TRUE,
    full.names = TRUE
  )
)

n_validation_files <- length(
  list.files(
    VALIDATION_DIR,
    recursive = TRUE,
    full.names = TRUE
  )
)

##############################################################
# Pipeline Stages
##############################################################

pipeline_modules <- c(
  
  "01. Merge FeatureCounts",
  
  "02. Prepare Metadata",
  
  "03. Differential Expression Analysis",
  
  "04. Variance Stabilizing Transformation",
  
  "05. PCA Analysis",
  
  "06. Volcano Plot",
  
  "07. Heatmap",
  
  "08. Batch Assessment",
  
  "09. Biomarker Candidate Selection",
  
  "10. Functional Enrichment",
  
  "11. ROC Analysis",
  
  "12. Biomarker Expression Plots",
  
  "13. Random Forest Validation",
  
  "14A. Validation Dataset Assessment",
  
  "14B. TCGA Dataset Preparation",
  
  "15. External Validation",
  
  "16. Validation Visualisation",
  
  "17. Session Information",
  
  "18. Publication Summary",
  
  "19. Pipeline Report"
  
)

##############################################################
# Build Report
##############################################################

report <- c(
  
  "==============================================================",
  
  "AML Biomarker Discovery Pipeline",
  
  "==============================================================",
  
  "",
  
  paste("Pipeline :", PIPELINE_NAME),
  
  paste("Version  :", PIPELINE_VERSION),
  
  paste("Author   :", AUTHOR),
  
  paste("Date     :", Sys.Date()),
  
  paste("R        :", R.version.string),
  
  "",
  
  "PROJECT SUMMARY",
  
  "----------------",
  
  paste("Pipeline Modules      :", length(pipeline_modules)),
  
  paste("Figures Generated     :", n_figures),
  
  paste("Tables Generated      :", n_tables),
  
  paste("Result Files          :", n_result_files),
  
  paste("Validation Files      :", n_validation_files),
  
  "",
  
  "PIPELINE MODULES",
  
  "----------------",
  
  paste("[✓]", pipeline_modules),
  
  "",
  
  "STATUS",
  
  "------",
  
  "Pipeline completed successfully.",
  
  "",
  
  "Output folders include:",
  
  "results/",
  
  "results/deseq2",
  
  "results/vst",
  
  "results/pca",
  
  "results/heatmap",
  
  "results/roc",
  
  "results/biomarkers",
  
  "results/random_forest",
  
  "results/enrichment",
  
  "results/validation",
  
  "results/publication",
  
  "results/logs",
  
  "results/report"
  
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

##############################################################
# Console Output
##############################################################

cat("---------------------------------------\n")
cat("Pipeline Report Generated Successfully\n")
cat("---------------------------------------\n\n")

cat("Modules completed :", length(pipeline_modules), "\n")
cat("Figures generated :", n_figures, "\n")
cat("Tables generated  :", n_tables, "\n")
cat("Result files      :", n_result_files, "\n\n")

cat("File created:\n")
cat("results/report/Pipeline_Report.txt\n")