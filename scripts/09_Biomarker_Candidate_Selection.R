##############################################################
# AML Biomarker Discovery Pipeline
#
# Script:
# 09_Biomarker_Candidate_Selection.R
#
# Purpose:
# Rank significant differentially expressed genes and
# identify the strongest AML biomarker candidates.
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

cat("\n")
cat("=========================================\n")
cat("Stage 9 : Biomarker Candidate Selection\n")
cat("=========================================\n\n")

start_time <- Sys.time()

##############################################################
# Directories
##############################################################

BIOMARKER_DIR <- file.path(
  RESULTS_DIR,
  "biomarkers"
)

dir.create(
  BIOMARKER_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

##############################################################
# Input File
##############################################################

INPUT_FILE <- file.path(
  RESULTS_DIR,
  "deseq2",
  "DESeq2_significant_results.csv"
)

##############################################################
# Validate Input
##############################################################

if (!file.exists(INPUT_FILE)) {
  
  stop(
    "Significant DESeq2 results not found.\nRun 03_DESeq2_analysis.R first."
  )
  
}

##############################################################
# Load Data
##############################################################

cat("Loading significant biomarkers...\n")

results <- readr::read_csv(
  INPUT_FILE,
  show_col_types = FALSE
)

cat(
  "Significant Genes :",
  nrow(results),
  "\n\n"
)

##############################################################
# Remove Missing Values
##############################################################

results <- dplyr::filter(
  
  results,
  
  !is.na(padj)
  
)

##############################################################
# Add Gene Symbols
##############################################################

cat("Mapping Ensembl IDs to Gene Symbols...\n")

results$GeneID <- sub("\\..*$", "", results$GeneID)

gene_map <- AnnotationDbi::select(
  
  org.Hs.eg.db,
  
  keys = unique(results$GeneID),
  
  keytype = "ENSEMBL",
  
  columns = c("SYMBOL")
  
)

gene_map <- gene_map |>
  
  dplyr::filter(!is.na(SYMBOL)) |>
  
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

results <- dplyr::left_join(
  
  results,
  
  gene_map,
  
  by = c("GeneID" = "ENSEMBL")
  
)

results <- results |>
  
  dplyr::rename(
    
    GeneSymbol = SYMBOL
    
  )

cat(
  
  "Gene Symbols mapped :",
  
  sum(!is.na(results$GeneSymbol)),
  
  "\n\n"
  
)

##############################################################
# Rank Biomarkers
##############################################################

biomarkers <-
  
  results |>
  
  dplyr::arrange(
    
    padj,
    
    dplyr::desc(abs(log2FoldChange))
    
  )

##############################################################
# Select Top Biomarkers
##############################################################

top50 <-
  
  dplyr::slice_head(
    
    biomarkers,
    
    n = 50
    
  )

top20 <-
  
  dplyr::slice_head(
    
    biomarkers,
    
    n = 20
    
  )

top_up <-
  
  biomarkers |>
  
  dplyr::filter(
    
    log2FoldChange > 0
    
  ) |>
  
  dplyr::slice_head(
    
    n = 25
    
  )

top_down <-
  
  biomarkers |>
  
  dplyr::filter(
    
    log2FoldChange < 0
    
  ) |>
  
  dplyr::slice_head(
    
    n = 25
    
  )

##############################################################
# Output Files
##############################################################

RANKED_FILE <- file.path(
  BIOMARKER_DIR,
  "Ranked_Biomarkers.csv"
)

TOP50_FILE <- file.path(
  BIOMARKER_DIR,
  "Top50_Biomarkers.csv"
)

TOP20_FILE <- file.path(
  BIOMARKER_DIR,
  "Top20_Biomarkers.csv"
)

UP_FILE <- file.path(
  BIOMARKER_DIR,
  "Top25_Upregulated.csv"
)

DOWN_FILE <- file.path(
  BIOMARKER_DIR,
  "Top25_Downregulated.csv"
)

##############################################################
# Save Results
##############################################################

readr::write_csv(
  biomarkers,
  RANKED_FILE
)

readr::write_csv(
  top50,
  TOP50_FILE
)

readr::write_csv(
  top20,
  TOP20_FILE
)

readr::write_csv(
  top_up,
  UP_FILE
)

readr::write_csv(
  top_down,
  DOWN_FILE
)

##############################################################
# Summary
##############################################################

end_time <- Sys.time()

cat("\n")
cat("=========================================\n")
cat("Stage 9 Completed Successfully\n")
cat("=========================================\n\n")

cat(
  "Significant Genes :",
  nrow(results),
  "\n"
)

cat(
  "Top 50 Biomarkers :",
  nrow(top50),
  "\n"
)

cat(
  "Top Upregulated   :",
  nrow(top_up),
  "\n"
)

cat(
  "Top Downregulated :",
  nrow(top_down),
  "\n\n"
)

cat(
  "Ranked Biomarkers :",
  RANKED_FILE,
  "\n"
)

cat(
  "Top 50            :",
  TOP50_FILE,
  "\n"
)

cat(
  "Top 20            :",
  TOP20_FILE,
  "\n"
)

cat(
  "Top Upregulated   :",
  UP_FILE,
  "\n"
)

cat(
  "Top Downregulated :",
  DOWN_FILE,
  "\n\n"
)

cat(
  "Time Elapsed :",
  round(end_time - start_time, 2),
  "\n\n"
)