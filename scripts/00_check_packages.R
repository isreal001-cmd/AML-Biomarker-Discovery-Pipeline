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
#
# Version:
# 1.0.1
##############################################################

##############################################################
# Load Global Configuration
##############################################################

source("config.R")

cat("\n")
cat("=========================================\n")
cat("AML Biomarker Discovery Pipeline\n")
cat("Package Check and Installation\n")
cat("=========================================\n\n")

start_time <- Sys.time()

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
# Install BiocManager
##############################################################

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  
  cat("Installing BiocManager...\n")
  
  install.packages(
    "BiocManager",
    repos = "https://cran.rstudio.com/"
  )
  
}

##############################################################
# Check Bioconductor Version
##############################################################

cat(
  "Bioconductor Version :",
  as.character(BiocManager::version()),
  "\n\n"
)

##############################################################
# Install Missing CRAN Packages
##############################################################

cat("Checking CRAN packages...\n\n")

for (pkg in cran_packages) {
  
  if (requireNamespace(pkg, quietly = TRUE)) {
    
    cat(
      "Already installed :",
      pkg,
      "\n"
    )
    
  } else {
    
    cat(
      "Installing CRAN package :",
      pkg,
      "\n"
    )
    
    tryCatch(
      
      {
        
        install.packages(
          pkg,
          repos = "https://cran.rstudio.com/"
        )
        
      },
      
      error = function(e) {
        
        stop(
          paste0(
            "\nFailed to install CRAN package: ",
            pkg,
            "\nReason: ",
            conditionMessage(e)
          )
        )
        
      }
      
    )
    
    if (!requireNamespace(pkg, quietly = TRUE)) {
      
      stop(
        paste0(
          "\nPackage installation failed: ",
          pkg
        )
      )
      
    }
    
  }
  
}

cat("\n")

##############################################################
# Install Missing Bioconductor Packages
##############################################################

cat("Checking Bioconductor packages...\n\n")

for (pkg in bioc_packages) {
  
  if (requireNamespace(pkg, quietly = TRUE)) {
    
    cat(
      "Already installed :",
      pkg,
      "\n"
    )
    
  } else {
    
    cat(
      "Installing Bioconductor package :",
      pkg,
      "\n"
    )
    
    tryCatch(
      
      {
        
        BiocManager::install(
          pkg,
          ask = FALSE,
          update = FALSE
        )
        
      },
      
      error = function(e) {
        
        stop(
          paste0(
            "\nFailed to install Bioconductor package: ",
            pkg,
            "\nReason: ",
            conditionMessage(e)
          )
        )
        
      }
      
    )
    
    if (!requireNamespace(pkg, quietly = TRUE)) {
      
      stop(
        paste0(
          "\nPackage installation failed: ",
          pkg
        )
      )
      
    }
    
  }
  
}

cat("\n")

##############################################################
# Load and Validate All Packages
##############################################################

all_packages <- c(
  
  cran_packages,
  bioc_packages
  
)

cat("Loading and validating packages...\n\n")

loaded_packages <- character()
failed_packages <- character()

for (pkg in all_packages) {
  
  load_result <- tryCatch(
    
    {
      
      suppressPackageStartupMessages(
        
        library(
          pkg,
          character.only = TRUE
        )
        
      )
      
      TRUE
      
    },
    
    error = function(e) {
      
      cat(
        "FAILED :",
        pkg,
        "\n"
      )
      
      cat(
        "Reason :",
        conditionMessage(e),
        "\n\n"
      )
      
      FALSE
      
    }
    
  )
  
  if (load_result) {
    
    loaded_packages <- c(
      loaded_packages,
      pkg
    )
    
    cat(
      "[",
      length(loaded_packages),
      "/",
      length(all_packages),
      "] ",
      pkg,
      " loaded successfully\n",
      sep = ""
    )
    
  } else {
    
    failed_packages <- c(
      failed_packages,
      pkg
    )
    
  }
  
}

##############################################################
# Final Validation
##############################################################

if (length(failed_packages) > 0) {
  
  cat("\n")
  cat("=========================================\n")
  cat("PACKAGE VALIDATION FAILED\n")
  cat("=========================================\n\n")
  
  cat(
    "Failed Packages :",
    paste(failed_packages, collapse = ", "),
    "\n\n"
  )
  
  stop(
    "One or more required packages could not be loaded."
  )
  
}

##############################################################
# Summary
##############################################################

end_time <- Sys.time()

cat("\n")
cat("=========================================\n")
cat("All Required Packages Loaded Successfully\n")
cat("=========================================\n\n")

cat(
  "CRAN Packages        :",
  length(cran_packages),
  "\n"
)

cat(
  "Bioconductor Packages:",
  length(bioc_packages),
  "\n"
)

cat(
  "Total Packages       :",
  length(all_packages),
  "\n"
)

cat(
  "Time Elapsed         :",
  round(
    as.numeric(
      difftime(
        end_time,
        start_time,
        units = "secs"
      )
    ),
    2
  ),
  "seconds\n\n"
)

cat("Package check complete.\n\n")