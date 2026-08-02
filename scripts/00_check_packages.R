##############################################################
# AML Biomarker Discovery Pipeline
#
# Script:
# 00_check_packages.R
#
# Purpose:
# Check, install and load all required packages.
#
# Author:
# Isreal Oluwafemi Abiodun
##############################################################

source("config.R")

cat("=========================================\n")
cat("Checking Required Packages\n")
cat("=========================================\n\n")

##############################################################
# CRAN Packages
##############################################################

cran_packages <- c(
  
  "readr",
  "dplyr",
  "tidyr",
  "purrr",
  "ggplot2",
  "stringr",
  "pROC",
  "randomForest",
  "pheatmap",
  "RColorBrewer",
  "circlize"
  
)

##############################################################
# Bioconductor Packages
##############################################################

bioc_packages <- c(
  
  "DESeq2",
  "SummarizedExperiment",
  "clusterProfiler",
  "org.Hs.eg.db",
  "ComplexHeatmap",
  "EnhancedVolcano",
  "TCGAbiolinks"
  
)

##############################################################
# Install CRAN Packages
##############################################################

for(pkg in cran_packages){
  
  if(!requireNamespace(pkg, quietly = TRUE)){
    
    cat("Installing CRAN package:", pkg, "\n")
    
    install.packages(pkg)
    
  }
  
}

##############################################################
# Install Bioconductor Packages
##############################################################

if(!requireNamespace("BiocManager", quietly = TRUE)){
  
  install.packages("BiocManager")
  
}

for(pkg in bioc_packages){
  
  if(!requireNamespace(pkg, quietly = TRUE)){
    
    cat("Installing Bioconductor package:", pkg, "\n")
    
    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )
    
  }
  
}

##############################################################
# Load Packages
##############################################################

all_packages <- c(
  cran_packages,
  bioc_packages
)

for(pkg in all_packages){
  
  suppressPackageStartupMessages(
    
    library(
      pkg,
      character.only = TRUE
    )
    
  )
  
  cat("✓", pkg, "\n")
  
}

cat("\n")
cat("=========================================\n")
cat("All Required Packages Loaded Successfully\n")
cat("=========================================\n\n")