###############################################################
# AML Biomarker Discovery Pipeline
# Script: 05_PCA_analysis.R
###############################################################

##############################
# 1. Clear Environment
##############################

# rm(list = ls())

##############################
# 2. Load Configuration
##############################

source("config.R")
source("functions/plotting_functions.R")

##############################
# 3. Load Packages
##############################

packages <- c(
  "DESeq2",
  "ggplot2"
)

for(pkg in packages){
  
  if(!require(pkg, character.only = TRUE)){
    
    if(pkg == "DESeq2"){
      
      if(!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      
      BiocManager::install("DESeq2")
      
    } else{
      
      install.packages(pkg)
      
    }
    
    library(pkg, character.only = TRUE)
    
  }
  
}

###############################################################
# 4. Load VST Object
###############################################################

vsd <- readRDS(
  "results/vst/vsd.rds"
)

cat("---------------------------------------\n")
cat("VST Object Loaded Successfully\n")
cat("---------------------------------------\n")

###############################################################
# 5. PCA Calculation
###############################################################

pcaData <- plotPCA(
  vsd,
  intgroup = c("Condition", "Platform"),
  returnData = TRUE
)

percentVar <- round(
  100 * attr(pcaData, "percentVar")
)

cat("\nVariance Explained\n")

print(percentVar)

###############################################################
# 6. Save PCA Coordinates
###############################################################

dir.create(
  "results/pca",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  pcaData,
  "results/pca/PCA_coordinates.csv",
  row.names = FALSE
)

###############################################################
# 7. Publication PCA Plot
###############################################################

p <- ggplot(
  pcaData,
  aes(
    PC1,
    PC2,
    colour = Condition,
    shape = Platform
  )
) +
  
  geom_point(size = 4) +
  
  xlab(
    paste0(
      "PC1 (",
      percentVar[1],
      "%)"
    )
  ) +
  
  ylab(
    paste0(
      "PC2 (",
      percentVar[2],
      "%)"
    )
  ) +
  
  ggtitle("PCA Before ComBat Correction") +
  
  scale_colour_manual(
    values = c(
      Healthy = HEALTHY_BLUE,
      AML = AML_RED
    )
  ) +
  
  publication_theme()

###############################################################
# 8. Save Figure
###############################################################

ggsave(
  filename = file.path(
    "results/pca",
    "PCA_before_ComBat.png"
  ),
  plot = p,
  width = FIG_WIDTH,
  height = FIG_HEIGHT,
  dpi = FIG_DPI
)

ggsave(
  filename = file.path(
    "results/pca",
    "PCA_before_ComBat.pdf"
  ),
  plot = p,
  width = FIG_WIDTH,
  height = FIG_HEIGHT
)

cat("\n---------------------------------------\n")
cat("PCA Analysis Completed Successfully\n")
cat("---------------------------------------\n")

cat("\nFiles created:\n")

cat("results/pca/PCA_before_ComBat.png\n")
cat("results/pca/PCA_before_ComBat.pdf\n")
cat("results/pca/PCA_coordinates.csv\n")
