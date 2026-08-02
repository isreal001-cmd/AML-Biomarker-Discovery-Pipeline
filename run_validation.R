##############################################################
# AML Biomarker Discovery Pipeline
# Validation Pipeline Runner
#
# Author  : Isreal Oluwafemi Abiodun
# Version : 1.0.0
##############################################################

rm(list = ls())

cat("\n=========================================\n")
cat(" External Validation Pipeline\n")
cat("=========================================\n\n")

validation_scripts <- c(
  
  "scripts/14A_Prepare_Validation_Data.R",
  
  "scripts/14B_Download_TCGA_Data.R",
  
  "scripts/15_External_Validation.R",
  
  "scripts/16_Validation_Visualization.R"
  
)

for(script in validation_scripts){
  
  cat("\nRunning:", basename(script), "\n")
  
  if(file.exists(script)){
    
    source(script)
    
  }else{
    
    stop(paste("Missing:", script))
    
  }
  
}

cat("\nValidation Pipeline Complete\n")