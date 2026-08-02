##############################################################
# AML Biomarker Discovery Pipeline
# Master Pipeline Runner
#
# Author  : Isreal Oluwafemi Abiodun
# Version : 1.0.0
##############################################################

##############################################################
# Clean Environment
##############################################################

rm(list = ls())

gc()

##############################################################
# Load Global Configuration
##############################################################

source("config.R")

cat("Master Pipeline Started\n\n")

##############################################################
# Pipeline Start Time
##############################################################

start_time <- Sys.time()

##############################################################
# Pipeline Scripts
##############################################################

pipeline_scripts <- c(
  
  "scripts/00_check_packages.R",
  
  "scripts/01_merge_featureCounts.R",
  
  "scripts/02_prepare_metadata.R",
  
  "scripts/03_DESeq2_analysis.R",
  
  "scripts/04_VST_normalization.R",
  
  "scripts/05_PCA_analysis.R",
  
  "scripts/06_Volcano_plot.R",
  
  "scripts/07_Heatmap_top50.R",
  
  "scripts/08_Batch_Assessment.R",
  
  "scripts/09_Biomarker_Candidate_Selection.R",
  
  "scripts/10_Functional_Enrichment.R",
  
  "scripts/11_ROC_Analysis.R",
  
  "scripts/12_Biomarker_Expression_Plots.R",
  
  "scripts/13_Random_Forest_Validation.R",
  
  "scripts/14A_Prepare_Validation_Data.R",
  
  "scripts/14B_Download_TCGA_Data.R",
  
  "scripts/15_External_Validation.R",
  
  "scripts/16_Validation_Visualization.R",
  
  "scripts/17_Save_Session_Info.R",
  
  "scripts/18_Publication_Summary.R",
  
  "scripts/19_Generate_Report.R"
  
)

##############################################################
# Execute Pipeline
##############################################################

TOTAL_STEPS <- length(pipeline_scripts)

cat("Total Pipeline Steps :", TOTAL_STEPS, "\n\n")

for (i in seq_along(pipeline_scripts)) {
  
  script <- pipeline_scripts[i]
  
  cat("=========================================\n")
  cat("Step", i, "of", TOTAL_STEPS, "\n")
  cat("Running :", basename(script), "\n")
  cat("=========================================\n")
  
  script_start <- Sys.time()
  
  tryCatch(
    
    {
      
      source(script)
      
      script_end <- Sys.time()
      
      elapsed <- round(
        as.numeric(
          difftime(
            script_end,
            script_start,
            units = "secs"
          )
        ),
        2
      )
      
      cat("✓ Finished :", basename(script), "\n")
      cat("Elapsed    :", elapsed, "seconds\n\n")
      
    },
    
    error = function(e) {
      
      cat("\n")
      cat("=========================================\n")
      cat("PIPELINE FAILED\n")
      cat("=========================================\n")
      cat("Step       :", i, "of", TOTAL_STEPS, "\n")
      cat("Script     :", basename(script), "\n")
      cat("Reason     :", e$message, "\n")
      cat("=========================================\n\n")
      
      stop(e)
      
    }
    
  )
  
}

##############################################################
# Pipeline Finished
##############################################################

end_time <- Sys.time()

total_time <- round(
  as.numeric(
    difftime(
      end_time,
      start_time,
      units = "mins"
    )
  ),
  2
)

cat("=========================================\n")
cat("PIPELINE COMPLETED SUCCESSFULLY\n")
cat("=========================================\n")
cat("Pipeline :", PIPELINE_NAME, "\n")
cat("Version  :", PIPELINE_VERSION, "\n")
cat("Started  :", format(start_time), "\n")
cat("Finished :", format(end_time), "\n")
cat("Duration :", total_time, "minutes\n")
cat("=========================================\n")