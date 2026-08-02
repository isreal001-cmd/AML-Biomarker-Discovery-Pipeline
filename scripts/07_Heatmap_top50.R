###############################################################
# AML Biomarker Discovery Pipeline
# Script: 07_Heatmap_top50.R
# Purpose:
# Generate a publication-quality heatmap of the
# Top 25 Upregulated + Top 25 Downregulated genes
###############################################################

#==============================================================
# 1. INITIALISE ENVIRONMENT
#==============================================================

# rm(list = ls())

source("config.R")
source("functions/plotting_functions.R")

library(ComplexHeatmap)
library(circlize)
library(DESeq2)
library(readr)
library(dplyr)

#==============================================================
# 2. LOAD INPUT FILES
#==============================================================

results <- read_csv(
  "results/deseq2/DESeq2_all_results.csv",
  show_col_types = FALSE
)

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
# 3. SELECT TOP DIFFERENTIALLY EXPRESSED GENES
#==============================================================

results <- results %>%
  filter(
    !is.na(padj),
    abs(log2FoldChange) >= 1
  )

top_up <- results %>%
  filter(log2FoldChange > 1) %>%
  arrange(padj) %>%
  slice_head(n = 25)

top_down <- results %>%
  filter(log2FoldChange < -1) %>%
  arrange(padj) %>%
  slice_head(n = 25)

top50 <- bind_rows(top_up, top_down)

cat("Top genes selected:", nrow(top50), "\n\n")

#==============================================================
# 4. EXTRACT VST EXPRESSION MATRIX
#==============================================================

expr <- assay(vsd)

expr <- expr[top50$GeneID, ]

#==============================================================
# 5. Z-SCORE NORMALISATION
#==============================================================

expr_scaled <- t(scale(t(expr)))

#==============================================================
# 6. SAMPLE ANNOTATION
#==============================================================

annotation <- HeatmapAnnotation(
  
  Condition = metadata$Condition,
  
  col = list(
    
    Condition = c(
      
      Healthy = HEALTHY_BLUE,
      
      AML = AML_RED
      
    )
    
  )
  
)

#==============================================================
# 7. CREATE OUTPUT DIRECTORY
#==============================================================

dir.create(
  "results/heatmap",
  recursive = TRUE,
  showWarnings = FALSE
)

#==============================================================
# 8. BUILD HEATMAP OBJECT
#==============================================================

heatmap_plot <- Heatmap(
  
  expr_scaled,
  
  name = "Z-score",
  
  top_annotation = annotation,
  
  cluster_rows = TRUE,
  
  cluster_columns = TRUE,
  
  show_row_names = TRUE,
  
  show_column_names = TRUE,
  
  row_names_gp = grid::gpar(fontsize = 8),
  
  column_names_gp = grid::gpar(fontsize = 9),
  
  column_title = "AML vs Healthy",
  
  row_title = "Top 50 Differentially Expressed Genes"
  
)

#==============================================================
# 9. SAVE PNG
#==============================================================

png(
  "results/heatmap/Top50_DEGs_heatmap.png",
  width = 2200,
  height = 1800,
  res = 300
)

draw(heatmap_plot)

dev.off()

#==============================================================
# 10. SAVE PDF
#==============================================================

pdf(
  "results/heatmap/Top50_DEGs_heatmap.pdf",
  width = 10,
  height = 8
)

draw(heatmap_plot)

dev.off()

#==============================================================
# 11. SUMMARY
#==============================================================

cat("---------------------------------------\n")
cat("Heatmap Generated Successfully\n")
cat("---------------------------------------\n\n")

cat("Files created:\n")

cat("results/heatmap/Top50_DEGs_heatmap.png\n")

cat("results/heatmap/Top50_DEGs_heatmap.pdf\n")