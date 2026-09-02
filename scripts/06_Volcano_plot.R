##############################################################
# AML Biomarker Discovery Pipeline
#
# Script:
# 06_Volcano_plot.R
#
# Purpose:
# Generate publication quality volcano plot showing
# differentially expressed genes.
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
cat("Stage 6 : Volcano Plot\n")
cat("=========================================\n\n")

start_time <- Sys.time()

##############################################################
# Output Directory
##############################################################

dir.create(
  VOLCANO_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

##############################################################
# Input File
##############################################################

RESULT_FILE <- file.path(
  RESULTS_DIR,
  "deseq2",
  "DESeq2_all_results.csv"
)

##############################################################
# Output Files
##############################################################

PNG_FILE <- file.path(
  VOLCANO_DIR,
  "Volcano_plot.png"
)

PDF_FILE <- file.path(
  VOLCANO_DIR,
  "Volcano_plot.pdf"
)

##############################################################
# Validate Input
##############################################################

if (!file.exists(RESULT_FILE)) {
  
  stop(
    "DESeq2 results not found.\nRun 03_DESeq2_analysis.R first."
  )
  
}

##############################################################
# Load Results
##############################################################

cat("Loading DESeq2 results...\n")

results <- readr::read_csv(
  RESULT_FILE,
  show_col_types = FALSE
)

cat("Genes Loaded :", nrow(results), "\n\n")

##############################################################
# Prepare Volcano Data
##############################################################

results <- dplyr::filter(
  results,
  !is.na(padj)
)

results$Significance <- dplyr::case_when(
  
  results$padj < PADJ_THRESHOLD &
    results$log2FoldChange >= LOG2FC_THRESHOLD ~
    
    "Up",
  
  results$padj < PADJ_THRESHOLD &
    results$log2FoldChange <= -LOG2FC_THRESHOLD ~
    
    "Down",
  
  TRUE ~
    
    "Not Significant"
  
)

results$minusLog10Padj <-
  
  -log10(results$padj)

##############################################################
# Candidate Biomarkers
##############################################################

genes_to_label <- c(
  
  "SNHG15",
  
  "MMP13",
  
  "XPO1"
  
)

label_data <- dplyr::filter(
  
  results,
  
  GeneID %in% genes_to_label
  
)

##############################################################
# Volcano Plot
##############################################################

p <- ggplot2::ggplot(
  
  results,
  
  ggplot2::aes(
    
    log2FoldChange,
    
    minusLog10Padj,
    
    colour = Significance
    
  )
  
) +
  
  ggplot2::geom_point(
    
    alpha = 0.70,
    
    size = 2
    
  ) +
  
  ggplot2::geom_vline(
    
    xintercept = c(
      
      -LOG2FC_THRESHOLD,
      
      LOG2FC_THRESHOLD
      
    ),
    
    linetype = 2
    
  ) +
  
  ggplot2::geom_hline(
    
    yintercept =
      
      -log10(PADJ_THRESHOLD),
    
    linetype = 2
    
  ) +
  
  ggrepel::geom_text_repel(
    
    data = label_data,
    
    ggplot2::aes(
      
      label = GeneID
      
    ),
    
    size = 4,
    
    max.overlaps = Inf
    
  ) +
  
  ggplot2::scale_colour_manual(
    
    values = c(
      
      Up = AML_RED,
      
      Down = HEALTHY_BLUE,
      
      "Not Significant" = GREY
      
    )
    
  ) +
  
  ggplot2::labs(
    
    title =
      
      "Differential Gene Expression",
    
    subtitle =
      
      "AML versus Healthy",
    
    x =
      
      "Log2 Fold Change",
    
    y =
      
      expression(-log[10]("Adjusted P value"))
    
  ) +
  
  publication_theme()

##############################################################
# Export Figures
##############################################################

ggplot2::ggsave(
  
  PNG_FILE,
  
  p,
  
  width = FIG_WIDTH,
  
  height = FIG_HEIGHT,
  
  dpi = FIG_DPI
  
)

ggplot2::ggsave(
  
  PDF_FILE,
  
  p,
  
  width = FIG_WIDTH,
  
  height = FIG_HEIGHT
  
)

##############################################################
# Summary
##############################################################

end_time <- Sys.time()

cat("\n")
cat("=========================================\n")
cat("Stage 6 Completed Successfully\n")
cat("=========================================\n\n")

print(table(results$Significance))

cat("\n")

cat("PNG Figure :", PNG_FILE, "\n")

cat("PDF Figure :", PDF_FILE, "\n\n")

cat(
  
  "Time Elapsed :",
  
  round(end_time - start_time, 2),
  
  "\n\n"
)