###############################################################
# AML Biomarker Discovery Pipeline
# Script: 10_Functional_Enrichment.R
# Purpose:
# Functional enrichment analysis using GO and KEGG
###############################################################

# rm(list = ls())

###############################################################
# Load configuration
###############################################################

source("config.R")
#source("functions/enrichment_functions.R")

###############################################################
# Load packages
###############################################################

packages <- c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot",
  "ggplot2",
  "readr",
  "dplyr"
)

for(pkg in packages){
  
  if(!require(pkg, character.only = TRUE)){
    
    if(pkg %in% c(
      "clusterProfiler",
      "org.Hs.eg.db",
      "enrichplot"
    )){
      
      if(!requireNamespace("BiocManager",
                           quietly = TRUE))
        install.packages("BiocManager")
      
      BiocManager::install(
        pkg,
        ask = FALSE,
        update = FALSE
      )
      
    }else{
      
      install.packages(pkg)
      
    }
    
    library(pkg,
            character.only = TRUE)
    
  }
  
}

###############################################################
# Create output folder
###############################################################

dir.create(
  "results/enrichment",
  recursive = TRUE,
  showWarnings = FALSE
)

###############################################################
# Load biomarkers
###############################################################

biomarkers <- read_csv(
  "results/biomarkers/Ranked_Biomarkers.csv",
  show_col_types = FALSE
)

cat("---------------------------------------\n")
cat("Biomarkers Loaded Successfully\n")
cat("---------------------------------------\n\n")

cat("Genes:", nrow(biomarkers), "\n")