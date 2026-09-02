##############################################################
# AML Biomarker Discovery Pipeline
# Script: 18_Publication_Summary.R
#
# Author  : Isreal Oluwafemi Abiodun
# Version : 1.0.0
##############################################################

source("config.R")

cat("---------------------------------------\n")
cat("Generating Publication Summary\n")
cat("---------------------------------------\n")

##############################################################
# Create Publication Directory
##############################################################

PUBLICATION_DIR <- file.path(
  RESULTS_DIR,
  "publication"
)

dir.create(
  PUBLICATION_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

##############################################################
# Pipeline Summary
##############################################################

summary_text <- c(
  
  "AML Biomarker Discovery Pipeline",
  
  "========================================",
  
  paste("Pipeline :", PIPELINE_NAME),
  
  paste("Version  :", PIPELINE_VERSION),
  
  paste("Author   :", AUTHOR),
  
  paste("Date     :", Sys.Date()),
  
  "",
  
  "Pipeline Modules",
  
  "----------------",
  
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
  
  "17. Save Session Information",
  
  "18. Publication Summary",
  
  "",
  
  "Pipeline completed successfully."
  
)

writeLines(
  summary_text,
  file.path(
    PUBLICATION_DIR,
    "Pipeline_Summary.txt"
  )
)

##############################################################
# Figure Index
##############################################################

figure_files <- list.files(
  
  FIGURE_DIR,
  
  recursive = TRUE,
  
  full.names = TRUE
  
)

if(length(figure_files) > 0){
  
  figure_index <- data.frame(
    
    Figure = basename(figure_files),
    
    Folder = dirname(figure_files),
    
    stringsAsFactors = FALSE
    
  )
  
  write.csv(
    
    figure_index,
    
    file.path(
      PUBLICATION_DIR,
      "Figure_Index.csv"
    ),
    
    row.names = FALSE
    
  )
  
}

##############################################################
# Important Results Manifest
##############################################################

important_files <- c(
  
  "results/deseq2/DESeq2_all_results.csv",
  
  "results/deseq2/DESeq2_significant_results.csv",
  
  "results/biomarkers/Ranked_Biomarkers.csv",
  
  "results/biomarkers/Top50_Biomarkers.csv",
  
  "results/biomarkers/Top20_Biomarkers.csv",
  
  "results/roc/ROC_AUC_Table.csv",
  
  "results/random_forest/Variable_Importance.csv",
  
  "results/validation/Validation_Summary.csv"
  
)

manifest <- data.frame(
  
  File = important_files,
  
  Exists = file.exists(important_files),
  
  stringsAsFactors = FALSE
  
)

write.csv(
  
  manifest,
  
  file.path(
    PUBLICATION_DIR,
    "Results_Manifest.csv"
  ),
  
  row.names = FALSE
  
)

##############################################################
# Pipeline Statistics
##############################################################

pipeline_stats <- data.frame(
  
  Metric = c(
    
    "Pipeline Version",
    
    "Author",
    
    "Analysis Date",
    
    "Total Figures",
    
    "Total Key Result Files"
    
  ),
  
  Value = c(
    
    PIPELINE_VERSION,
    
    AUTHOR,
    
    as.character(Sys.Date()),
    
    length(figure_files),
    
    sum(manifest$Exists)
    
  ),
  
  stringsAsFactors = FALSE
  
)

write.csv(
  
  pipeline_stats,
  
  file.path(
    PUBLICATION_DIR,
    "Pipeline_Statistics.csv"
  ),
  
  row.names = FALSE
  
)

##############################################################
# Console Output
##############################################################

cat("---------------------------------------\n")
cat("Publication Summary Generated\n")
cat("---------------------------------------\n\n")

cat("Files created:\n")

cat("results/publication/Pipeline_Summary.txt\n")
cat("results/publication/Figure_Index.csv\n")
cat("results/publication/Results_Manifest.csv\n")
cat("results/publication/Pipeline_Statistics.csv\n")