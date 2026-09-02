##############################################################
# AML Biomarker Discovery Pipeline
#
# Script:
# 05_PCA_analysis.R
#
# Purpose:
# Perform Principal Component Analysis using the
# variance stabilized expression matrix and generate
# publication quality PCA figures.
#
# Author:
# Isreal Oluwafemi Abiodun
#
# Version:
# 1.0.0
##############################################################

##############################################################
# Load Global Configuration
##############################################################

source("config.R")
source("functions/plotting_functions.R")

cat("\n")
cat("=========================================\n")
cat("Stage 5 : Principal Component Analysis\n")
cat("=========================================\n\n")

start_time <- Sys.time()

##############################################################
# Output Directory
##############################################################

PCA_RESULTS_DIR <- file.path(
  RESULTS_DIR,
  "pca"
)

dir.create(
  PCA_RESULTS_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

##############################################################
# Input File
##############################################################

VSD_FILE <- file.path(
  RESULTS_DIR,
  "vst",
  "vsd.rds"
)

##############################################################
# Output Files
##############################################################

PCA_COORD_FILE <- file.path(
  PCA_RESULTS_DIR,
  "PCA_coordinates.csv"
)

PCA_PNG_FILE <- file.path(
  PCA_RESULTS_DIR,
  "PCA_before_ComBat.png"
)

PCA_PDF_FILE <- file.path(
  PCA_RESULTS_DIR,
  "PCA_before_ComBat.pdf"
)

##############################################################
# Validate Input
##############################################################

if (!file.exists(VSD_FILE)) {
  
  stop(
    "VST object not found.\nRun 04_VST_normalization.R first."
  )
  
}

##############################################################
# Load VST Object
##############################################################

cat("Loading VST object...\n")

vsd <- readRDS(
  VSD_FILE
)

cat("VST Object Successfully Loaded\n\n")

##############################################################
# Principal Component Analysis
##############################################################

cat("Calculating Principal Components...\n")

pcaData <- DESeq2::plotPCA(
  
  vsd,
  
  intgroup = c(
    "Condition",
    "Platform"
  ),
  
  returnData = TRUE
  
)

percentVar <- round(
  
  100 *
    
    attr(
      pcaData,
      "percentVar"
    )
  
)

cat("Variance Explained\n")

print(percentVar)

##############################################################
# Save PCA Coordinates
##############################################################

readr::write_csv(
  
  pcaData,
  
  PCA_COORD_FILE
  
)

##############################################################
# Generate Publication Figure
##############################################################

p <- ggplot2::ggplot(
  
  pcaData,
  
  ggplot2::aes(
    
    PC1,
    
    PC2,
    
    colour = Condition,
    
    shape = Platform
    
  )
  
) +
  
  ggplot2::geom_point(
    size = 4
  ) +
  
  ggplot2::labs(
    
    title = "Principal Component Analysis",
    
    subtitle = "Before Batch Correction Assessment",
    
    x = paste0(
      "PC1 (",
      percentVar[1],
      "%)"
    ),
    
    y = paste0(
      "PC2 (",
      percentVar[2],
      "%)"
    )
    
  ) +
  
  ggplot2::scale_colour_manual(
    
    values = c(
      
      Healthy = HEALTHY_BLUE,
      
      AML = AML_RED
      
    )
    
  ) +
  
  publication_theme()

##############################################################
# Export Figures
##############################################################

ggplot2::ggsave(
  
  filename = PCA_PNG_FILE,
  
  plot = p,
  
  width = FIG_WIDTH,
  
  height = FIG_HEIGHT,
  
  dpi = FIG_DPI
  
)

ggplot2::ggsave(
  
  filename = PCA_PDF_FILE,
  
  plot = p,
  
  width = FIG_WIDTH,
  
  height = FIG_HEIGHT
  
)

##############################################################
# Summary
##############################################################

end_time <- Sys.time()

cat("\n")
cat("=========================================\n")
cat("Stage 5 Completed Successfully\n")
cat("=========================================\n\n")

cat(
  "Coordinates      :",
  PCA_COORD_FILE,
  "\n"
)

cat(
  "PNG Figure       :",
  PCA_PNG_FILE,
  "\n"
)

cat(
  "PDF Figure       :",
  PCA_PDF_FILE,
  "\n\n"
)

cat(
  "Variance PC1     :",
  percentVar[1],
  "%\n"
)

cat(
  "Variance PC2     :",
  percentVar[2],
  "%\n\n"
)

cat(
  "Time Elapsed     :",
  round(end_time - start_time, 2),
  "\n\n"
)