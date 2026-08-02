###############################################################
# AML Biomarker Discovery Pipeline
# Script: 04_VST_normalization.R
###############################################################

###############################################################
# 1. Clear Environment
###############################################################

# rm(list = ls())

###############################################################
# 2. Load Packages
###############################################################

packages <- c(
  "DESeq2"
)

for(pkg in packages){
  
  if(!require(pkg, character.only = TRUE)){
    
    if(pkg == "DESeq2"){
      
      if(!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      
      BiocManager::install("DESeq2")
      
    }
    
    library(pkg, character.only = TRUE)
    
  }
  
}

###############################################################
# 3. Create Output Folder
###############################################################

dir.create(
  "results/vst",
  recursive = TRUE,
  showWarnings = FALSE
)

###############################################################
# 4. Load DESeq2 Object
###############################################################

dds <- readRDS(
  "results/deseq2/dds.rds"
)

cat("---------------------------------------\n")
cat("DESeq2 Object Loaded Successfully\n")
cat("---------------------------------------\n")

###############################################################
# 5. Variance Stabilizing Transformation
###############################################################

cat("\nRunning VST...\n")

vsd <- vst(
  dds,
  blind = FALSE
)

cat("VST Completed.\n")

###############################################################
# 6. Extract Matrix
###############################################################

vst_matrix <- assay(vsd)

cat("\nDimensions:\n")

print(dim(vst_matrix))

###############################################################
# 7. Save Outputs
###############################################################

saveRDS(
  vsd,
  "results/vst/vsd.rds"
)

write.csv(
  vst_matrix,
  "results/vst/vst_expression_matrix.csv"
)

###############################################################
# 8. Summary
###############################################################

cat("\n---------------------------------------\n")

cat("VST Normalization Completed Successfully\n")

cat("---------------------------------------\n")

cat("\nFiles created:\n")

cat("results/vst/vsd.rds\n")

cat("results/vst/vst_expression_matrix.csv\n")