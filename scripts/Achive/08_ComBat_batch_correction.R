###############################################################
# AML Biomarker Discovery Pipeline
# Script: 08_ComBat_batch_correction.R
# Purpose:
# Correct sequencing platform batch effects using ComBat
###############################################################

#==============================================================
# 1. INITIALISE ENVIRONMENT
#==============================================================

rm(list = ls())

source("config.R")

packages <- c(
  "sva",
  "DESeq2",
  "readr"
)

for(pkg in packages){
  
  if(!require(pkg, character.only = TRUE)){
    
    if(pkg == "sva"){
      
      if(!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      
      BiocManager::install("sva")
      
    } else{
      
      install.packages(pkg)
      
    }
    
    library(pkg, character.only = TRUE)
    
  }
  
}

#==============================================================
# 2. LOAD INPUT FILES
#==============================================================

vsd <- readRDS(
  "results/vst/vsd.rds"
)

metadata <- read_csv(
  "data/metadata/sample_metadata.csv",
  show_col_types = FALSE
)

cat("---------------------------------------\n")
cat("Input Files Loaded Successfully\n")
cat("---------------------------------------\n\n")

#==============================================================
# 3. PREPARE EXPRESSION MATRIX
#==============================================================

expr <- assay(vsd)

cat("Expression Matrix Dimensions:\n")
print(dim(expr))

#==============================================================
# 4. PREPARE METADATA
#==============================================================

metadata <- as.data.frame(metadata)

rownames(metadata) <- metadata$Sample

metadata <- metadata[colnames(expr), ]

metadata$Condition <- factor(
  metadata$Condition,
  levels = c("Healthy","AML")
)

metadata$Platform <- factor(metadata$Platform)

cat("\nMetadata aligned successfully.\n")

#==============================================================
# 5. CREATE MODEL MATRIX
#==============================================================

mod <- model.matrix(
  ~ Condition,
  data = metadata
)

#==============================================================
# 6. RUN COMBAT
#==============================================================

cat("\nRunning ComBat Batch Correction...\n\n")

combat_expr <- ComBat(
  
  dat = expr,
  
  batch = metadata$Platform,
  
  mod = mod,
  
  par.prior = TRUE,
  
  prior.plots = FALSE
  
)

cat("ComBat completed successfully.\n")

#==============================================================
# 7. SAVE OUTPUT
#==============================================================

dir.create(
  "results/combat",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  combat_expr,
  "results/combat/combat_expression.rds"
)

write.csv(
  combat_expr,
  "results/combat/combat_expression.csv"
)

#==============================================================
# 8. SUMMARY
#==============================================================

cat("\n---------------------------------------\n")
cat("ComBat Batch Correction Completed\n")
cat("---------------------------------------\n\n")

cat("Files created:\n")

cat("results/combat/combat_expression.rds\n")

cat("results/combat/combat_expression.csv\n")