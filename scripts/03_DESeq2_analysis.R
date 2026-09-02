##############################################################
# AML Biomarker Discovery Pipeline
#
# Script:
# 03_DESeq2_analysis.R
#
# Purpose:
# Perform differential expression analysis using DESeq2
# to identify differentially expressed genes between
# AML and Healthy samples.
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

cat("\n")
cat("=========================================\n")
cat("Stage 3 : Differential Expression Analysis\n")
cat("=========================================\n\n")

start_time <- Sys.time()

##############################################################
# Output Directory
##############################################################

DESEQ_DIR <- file.path(
  RESULTS_DIR,
  "deseq2"
)

dir.create(
  DESEQ_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

##############################################################
# Output Files
##############################################################

ALL_RESULTS_FILE <- file.path(
  DESEQ_DIR,
  "DESeq2_all_results.csv"
)

SIG_RESULTS_FILE <- file.path(
  DESEQ_DIR,
  "DESeq2_significant_results.csv"
)

DDS_FILE <- file.path(
  DESEQ_DIR,
  "dds.rds"
)

##############################################################
# Input Files
##############################################################

COUNT_FILE <- file.path(
  DATA_DIR,
  "processed",
  "merged_counts.csv"
)

METADATA_FILE <- file.path(
  DATA_DIR,
  "metadata",
  "sample_metadata.csv"
)

##############################################################
# Validate Input Files
##############################################################

if (!file.exists(COUNT_FILE)) {
  
  stop(
    "Merged count matrix not found.\nRun 01_merge_featureCounts.R first."
  )
  
}

if (!file.exists(METADATA_FILE)) {
  
  stop(
    "Sample metadata not found.\nRun 02_prepare_metadata.R first."
  )
  
}

##############################################################
# Read Input Files
##############################################################

cat("Loading count matrix...\n")

counts <- readr::read_csv(
  COUNT_FILE,
  show_col_types = FALSE
)

cat("Loading metadata...\n")

metadata <- readr::read_csv(
  METADATA_FILE,
  show_col_types = FALSE
)

cat("\n")
cat("-----------------------------------------\n")
cat("Input Files Loaded Successfully\n")
cat("-----------------------------------------\n")

cat("Genes   :", nrow(counts), "\n")
cat("Samples :", ncol(counts) - 1, "\n\n")

##############################################################
# Prepare Count Matrix
##############################################################

count_matrix <- as.data.frame(counts)

rownames(count_matrix) <- count_matrix$Geneid

count_matrix <- count_matrix[, -1]

count_matrix <- as.matrix(count_matrix)

storage.mode(count_matrix) <- "integer"

cat("Count Matrix Dimensions\n")

print(dim(count_matrix))

##############################################################
# Prepare Metadata
##############################################################

metadata <- as.data.frame(metadata)

rownames(metadata) <- metadata$Sample

metadata <- metadata[, c(
  "Condition",
  "Platform"
)]

metadata$Condition <- factor(
  metadata$Condition,
  levels = c(
    "Healthy",
    "AML"
  )
)

metadata$Platform <- factor(
  metadata$Platform
)

##############################################################
# Match Sample Order
##############################################################

metadata <- metadata[
  colnames(count_matrix),
]

##############################################################
# Validate Sample Matching
##############################################################

if (!all(
  colnames(count_matrix) ==
  rownames(metadata)
)) {
  
  stop(
    "Sample names do not match between count matrix and metadata."
  )
  
}

cat("\n")
cat("Metadata Successfully Verified\n")

##############################################################
# Filter Low Count Genes
##############################################################

keep <- rowSums(
  count_matrix >= 10
) >= 3

count_matrix_filtered <- count_matrix[
  keep,
]

cat("\n")
cat("-----------------------------------------\n")
cat("Gene Filtering Summary\n")
cat("-----------------------------------------\n")

cat(
  "Genes Before Filtering :",
  nrow(count_matrix),
  "\n"
)

cat(
  "Genes After Filtering  :",
  nrow(count_matrix_filtered),
  "\n"
)

cat(
  "Genes Removed          :",
  nrow(count_matrix) -
    nrow(count_matrix_filtered),
  "\n\n"
)

##############################################################
# Construct DESeq Dataset
##############################################################

dds <- DESeq2::DESeqDataSetFromMatrix(
  
  countData = count_matrix_filtered,
  
  colData = metadata,
  
  design = ~ Condition
  
)

cat("DESeqDataSet Successfully Created\n")

##############################################################
# Differential Expression Analysis
##############################################################

cat("\n")
cat("Running DESeq2 Analysis...\n\n")

dds <- DESeq2::DESeq(dds)

cat("DESeq2 Analysis Completed Successfully\n")

##############################################################
# Extract Results
##############################################################

res <- DESeq2::results(dds)

cat("\n")
cat("-----------------------------------------\n")
cat("Differential Expression Summary\n")
cat("-----------------------------------------\n\n")

print(summary(res))

##############################################################
# Convert Results to Data Frame
##############################################################

res <- as.data.frame(res)

res$GeneID <- rownames(res)

res <- res[, c(
  
  "GeneID",
  
  "baseMean",
  
  "log2FoldChange",
  
  "lfcSE",
  
  "stat",
  
  "pvalue",
  
  "padj"
  
)]

##############################################################
# Significant Differentially Expressed Genes
##############################################################

sig_res <- subset(
  
  res,
  
  !is.na(padj) &
    
    padj < PADJ_THRESHOLD &
    
    abs(log2FoldChange) >= LOG2FC_THRESHOLD
  
)

##############################################################
# Upregulated Genes
##############################################################

upregulated <- subset(
  
  sig_res,
  
  log2FoldChange > LOG2FC_THRESHOLD
  
)

##############################################################
# Downregulated Genes
##############################################################

downregulated <- subset(
  
  sig_res,
  
  log2FoldChange < -LOG2FC_THRESHOLD
  
)

##############################################################
# Summary Statistics
##############################################################

cat("\n")
cat("-----------------------------------------\n")
cat("Significant DEG Summary\n")
cat("-----------------------------------------\n")

cat(
  
  "Total Significant :",
  
  nrow(sig_res),
  
  "\n"
  
)

cat(
  
  "Upregulated       :",
  
  nrow(upregulated),
  
  "\n"
  
)

cat(
  
  "Downregulated     :",
  
  nrow(downregulated),
  
  "\n"
  
)

##############################################################
# Save Results
##############################################################

readr::write_csv(
  
  res,
  
  ALL_RESULTS_FILE
  
)

readr::write_csv(
  
  sig_res,
  
  SIG_RESULTS_FILE
  
)

saveRDS(
  
  dds,
  
  DDS_FILE
  
)

##############################################################
# Final Summary
##############################################################

end_time <- Sys.time()

cat("\n")
cat("=========================================\n")
cat("Stage 3 Completed Successfully\n")
cat("=========================================\n\n")

cat(
  
  "Genes Tested        :",
  
  nrow(res),
  
  "\n"
  
)

cat(
  
  "Significant Genes   :",
  
  nrow(sig_res),
  
  "\n"
  
)

cat(
  
  "Upregulated Genes   :",
  
  nrow(upregulated),
  
  "\n"
  
)

cat(
  
  "Downregulated Genes :",
  
  nrow(downregulated),
  
  "\n\n"
  
)

cat(
  
  "All Results         :",
  
  ALL_RESULTS_FILE,
  
  "\n"
  
)

cat(
  
  "Significant Results :",
  
  SIG_RESULTS_FILE,
  
  "\n"
  
)

cat(
  
  "DESeq2 Object       :",
  
  DDS_FILE,
  
  "\n\n"
  
)

cat(
  
  "Time Elapsed        :",
  
  round(end_time - start_time, 2),
  
  "\n\n"
)