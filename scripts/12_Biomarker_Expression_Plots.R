###############################################################
# AML Biomarker Discovery Pipeline
# Script: 12_Biomarker_Expression_Plots.R
# Purpose:
# Generate publication-quality boxplots for Top 20 biomarkers
###############################################################

# rm(list = ls())

###############################################################
# Load configuration
###############################################################

source("config.R")
source("functions/plotting_functions.R")

###############################################################
# Load packages
###############################################################

packages <- c(
  "DESeq2",
  "readr",
  "dplyr",
  "ggplot2"
)

for(pkg in packages){
  
  if(!require(pkg, character.only = TRUE)){
    
    if(pkg == "DESeq2"){
      
      if(!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      
      BiocManager::install("DESeq2")
      
    }else{
      
      install.packages(pkg)
      
    }
    
    library(pkg, character.only = TRUE)
    
  }
  
}

###############################################################
# Load input files
###############################################################

vsd <- readRDS(
  "results/vst/vsd.rds"
)

metadata <- read_csv(
  "data/metadata/sample_metadata.csv",
  show_col_types = FALSE
)

biomarkers <- read_csv(
  "results/biomarkers/Top20_Biomarkers.csv",
  show_col_types = FALSE
)

cat("---------------------------------------\n")
cat("Input Files Loaded Successfully\n")
cat("---------------------------------------\n\n")

###############################################################
# Prepare output folder
###############################################################

dir.create(
  "results/expression_plots",
  recursive = TRUE,
  showWarnings = FALSE
)

###############################################################
# Extract expression matrix
###############################################################

expr <- assay(vsd)

###############################################################
# Generate plots
###############################################################

for(gene in biomarkers$GeneID){
  
  values <- expr[gene, ]
  
  plot_data <- data.frame(
    
    Sample = metadata$Sample,
    
    Condition = metadata$Condition,
    
    Expression = values
    
  )
  
  p <- ggplot(
    plot_data,
    aes(
      x = Condition,
      y = Expression,
      fill = Condition
    )
  ) +
    
    geom_violin(
      alpha = 0.35,
      trim = FALSE
    ) +
    
    geom_boxplot(
      width = 0.18,
      outlier.shape = NA,
      colour = "black"
    ) +
    
    geom_jitter(
      width = 0.08,
      size = 2,
      alpha = 0.8
    ) +
    
    scale_fill_manual(
      values = c(
        Healthy = HEALTHY_BLUE,
        AML = AML_RED
      )
    ) +
    
    labs(
      title = gene,
      x = "",
      y = "VST Expression"
    ) +
    
    publication_theme()
  
  ggsave(
    filename = file.path(
      "results/expression_plots",
      paste0(gene, ".png")
    ),
    plot = p,
    width = FIG_WIDTH,
    height = FIG_HEIGHT,
    dpi = FIG_DPI
  )
  
  ggsave(
    filename = file.path(
      "results/expression_plots",
      paste0(gene, ".pdf")
    ),
    plot = p,
    width = FIG_WIDTH,
    height = FIG_HEIGHT
  )
  
}

###############################################################
# Summary
###############################################################

cat("---------------------------------------\n")
cat("Expression Plots Generated Successfully\n")
cat("---------------------------------------\n\n")

cat("Genes plotted :", nrow(biomarkers), "\n")

cat("\nOutput folder:\n")
cat("results/expression_plots/\n")