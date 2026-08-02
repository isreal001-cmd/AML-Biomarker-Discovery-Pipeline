##############################################################
# AML Biomarker Discovery Pipeline
# Publication Summary Export
#
# Author  : Isreal Oluwafemi Abiodun
# Version : 1.0.0
##############################################################

source("config.R")

cat("---------------------------------------\n")
cat("Generating Publication Summary\n")
cat("---------------------------------------\n")

##############################################################
# Create publication directory
##############################################################

PUBLICATION_DIR <- file.path(RESULTS_DIR, "publication")

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
  
  "----------------------------------------",
  
  paste("Pipeline :", PIPELINE_NAME),
  
  paste("Version  :", PIPELINE_VERSION),
  
  paste("Author   :", AUTHOR),
  
  paste("Date     :", Sys.Date()),
  
  "",
  
  "Analysis Completed Successfully.",
  
  "",
  
  "Pipeline Modules:",
  
  "1. FeatureCounts Merge",
  
  "2. Metadata Preparation",
  
  "3. Differential Expression Analysis",
  
  "4. VST Normalization",
  
  "5. PCA",
  
  "6. Volcano Plot",
  
  "7. Heatmap",
  
  "8. Batch Assessment",
  
  "9. Biomarker Ranking",
  
  "10. Functional Enrichment",
  
  "11. ROC Analysis",
  
  "12. Biomarker Expression",
  
  "13. Random Forest Validation",
  
  "14. Validation Dataset Preparation",
  
  "15. TCGA Data Processing",
  
  "16. External Validation",
  
  "17. Validation Visualization",
  
  "18. Session Information",
  
  "",
  
  "Framework completed successfully."
  
)

writeLines(
  
  summary_text,
  
  file.path(
    PUBLICATION_DIR,
    "Results_Summary.txt"
  )
  
)

##############################################################
# Figure Index
##############################################################

figures <- list.files(
  
  FIGURE_DIR,
  
  recursive = TRUE,
  
  full.names = TRUE
  
)

if(length(figures) > 0){
  
  figure_table <- data.frame(
    
    Figure = basename(figures),
    
    Folder = dirname(figures),
    
    stringsAsFactors = FALSE
    
  )
  
  write.csv(
    
    figure_table,
    
    file.path(
      PUBLICATION_DIR,
      "Figure_Index.csv"
    ),
    
    row.names = FALSE
    
  )
  
}

##############################################################
# Export Important Result Tables
##############################################################

important_files <- c(
  
  "results/DEGs/Significant_DEGs.csv",
  
  "results/biomarkers/Top_Biomarkers.csv",
  
  "results/ROC/ROC_Results.csv",
  
  "results/random_forest/Feature_Importance.csv"
  
)

manifest <- data.frame(
  
  File = important_files,
  
  Exists = file.exists(important_files)
  
)

write.csv(
  
  manifest,
  
  file.path(
    PUBLICATION_DIR,
    "Results_Manifest.csv"
  ),
  
  row.names = FALSE
  
)

cat("Publication summary created.\n")