##############################################################
# AML Biomarker Discovery Pipeline
#
# Script:
# 07_Heatmap_top50.R
#
# Purpose:
# Generate a publication quality heatmap of the
# Top 25 Upregulated and Top 25 Downregulated genes.
#
# Author:
# Isreal Oluwafemi Abiodun
#
# Version:
# 1.0.0
##############################################################

##############################################################
# Load Configuration
##############################################################

source("config.R")
source("functions/plotting_functions.R")

cat("\n")
cat("=========================================\n")
cat("Stage 7 : Top 50 Heatmap\n")
cat("=========================================\n\n")

start_time <- Sys.time()

##############################################################
# Directories
##############################################################

HEATMAP_RESULTS_DIR <- file.path(
  RESULTS_DIR,
  "heatmap"
)

dir.create(
  HEATMAP_RESULTS_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

##############################################################
# Input Files
##############################################################

RESULT_FILE <- file.path(
  RESULTS_DIR,
  "deseq2",
  "DESeq2_all_results.csv"
)

VSD_FILE <- file.path(
  RESULTS_DIR,
  "vst",
  "vsd.rds"
)

METADATA_FILE <- file.path(
  "data",
  "metadata",
  "sample_metadata.csv"
)

##############################################################
# Output Files
##############################################################

PNG_FILE <- file.path(
  HEATMAP_RESULTS_DIR,
  "Top50_DEGs_heatmap.png"
)

PDF_FILE <- file.path(
  HEATMAP_RESULTS_DIR,
  "Top50_DEGs_heatmap.pdf"
)

##############################################################
# Validate Inputs
##############################################################

required_files <- c(
  RESULT_FILE,
  VSD_FILE,
  METADATA_FILE
)

missing_files <- required_files[!file.exists(required_files)]

if(length(missing_files) > 0){
  
  stop(
    paste(
      "Missing input files:\n",
      paste(missing_files, collapse = "\n")
    )
  )
  
}

##############################################################
# Load Data
##############################################################

cat("Loading input files...\n")

results <- readr::read_csv(
  RESULT_FILE,
  show_col_types = FALSE
)

vsd <- readRDS(
  VSD_FILE
)

metadata <- readr::read_csv(
  METADATA_FILE,
  show_col_types = FALSE
)

cat("Input files successfully loaded.\n\n")

##############################################################
# Select Top Differentially Expressed Genes
##############################################################

results <- dplyr::filter(
  
  results,
  
  !is.na(padj),
  
  padj < PADJ_THRESHOLD,
  
  abs(log2FoldChange) >= LOG2FC_THRESHOLD
  
)

top_up <- results |>
  
  dplyr::filter(
    log2FoldChange > LOG2FC_THRESHOLD
  ) |>
  
  dplyr::arrange(padj) |>
  
  dplyr::slice_head(n = 25)

top_down <- results |>
  
  dplyr::filter(
    log2FoldChange < -LOG2FC_THRESHOLD
  ) |>
  
  dplyr::arrange(padj) |>
  
  dplyr::slice_head(n = 25)

top50 <- dplyr::bind_rows(
  top_up,
  top_down
)

cat(
  "Selected genes :",
  nrow(top50),
  "\n\n"
)

##############################################################
# Extract Expression Matrix
##############################################################

expr <- assay(vsd)

available_genes <- intersect(
  top50$GeneID,
  rownames(expr)
)

expr <- expr[
  available_genes,
]

##############################################################
# Z Score Transformation
##############################################################

expr_scaled <- t(
  
  scale(
    
    t(expr)
    
  )
  
)

##############################################################
# Sample Annotation
##############################################################

annotation <- ComplexHeatmap::HeatmapAnnotation(
  
  Condition = metadata$Condition,
  
  col = list(
    
    Condition = c(
      
      Healthy = HEALTHY_BLUE,
      
      AML = AML_RED
      
    )
    
  )
  
)

##############################################################
# Build Heatmap
##############################################################

heatmap_plot <-
  
  ComplexHeatmap::Heatmap(
    
    expr_scaled,
    
    name = "Z score",
    
    top_annotation = annotation,
    
    cluster_rows = TRUE,
    
    cluster_columns = TRUE,
    
    show_row_names = TRUE,
    
    show_column_names = TRUE,
    
    row_names_gp = grid::gpar(
      fontsize = 8
    ),
    
    column_names_gp = grid::gpar(
      fontsize = 9
    ),
    
    column_title =
      
      "AML versus Healthy",
    
    row_title =
      
      "Top 50 Differentially Expressed Genes"
    
  )

##############################################################
# Export PNG
##############################################################

png(
  
  PNG_FILE,
  
  width = 2200,
  
  height = 1800,
  
  res = 300
  
)

ComplexHeatmap::draw(
  heatmap_plot
)

dev.off()

##############################################################
# Export PDF
##############################################################

pdf(
  
  PDF_FILE,
  
  width = 10,
  
  height = 8
  
)

ComplexHeatmap::draw(
  heatmap_plot
)

dev.off()

##############################################################
# Summary
##############################################################

end_time <- Sys.time()

cat("\n")
cat("=========================================\n")
cat("Stage 7 Completed Successfully\n")
cat("=========================================\n\n")

cat(
  "Genes plotted :",
  nrow(expr_scaled),
  "\n"
)

cat(
  "PNG Figure    :",
  PNG_FILE,
  "\n"
)

cat(
  "PDF Figure    :",
  PDF_FILE,
  "\n\n"
)

cat(
  "Time Elapsed  :",
  round(end_time - start_time, 2),
  "\n\n"
)