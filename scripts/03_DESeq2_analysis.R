###############################################################
# AML Biomarker Discovery Pipeline
# Script: 03_DESeq2_analysis.R
# Purpose: Differential Expression Analysis using DESeq2
# Author: Isreal Oluwafemi Abiodun
###############################################################

###############################################################
# 1. Clear Environment
###############################################################

# rm(list = ls())

###############################################################
# 2. Load Required Packages
###############################################################

packages <- c(
  "DESeq2",
  "readr",
  "dplyr",
  "tibble"
)

for(pkg in packages){
  
  if(!require(pkg, character.only = TRUE)){
    
    if(pkg == "DESeq2"){
      
      if(!requireNamespace("BiocManager", quietly =TRUE))
        install.packages("BiocManager")
      
      BiocManager::install("DESeq2")
      
    } else{
      
      install.packages(pkg)
      
    }
    
    library(pkg, character.only = TRUE)
    
  }
  
}

###############################################################
# 3. Create Output Folder
###############################################################

dir.create(
  "results/deseq2",
  recursive = TRUE,
  showWarnings = FALSE
)

###############################################################
# 4. Load Input Files
###############################################################

counts <- read_csv(
  "data/processed/merged_counts.csv",
  show_col_types = FALSE
)

metadata <- read_csv(
  "data/metadata/sample_metadata.csv",
  show_col_types = FALSE
)

cat("---------------------------------------\n")
cat("Input Files Loaded Successfully\n")
cat("---------------------------------------\n")

cat("Genes   :", nrow(counts), "\n")
cat("Samples :", ncol(counts)-1, "\n\n")

###############################################################
# 5. Prepare Count Matrix
###############################################################

count_matrix <- as.data.frame(counts)

rownames(count_matrix) <- count_matrix$Geneid

count_matrix <- count_matrix[, -1]

count_matrix <- as.matrix(count_matrix)

storage.mode(count_matrix) <- "integer"

cat("Count Matrix Dimensions:\n")

print(dim(count_matrix))

###############################################################
# 6. Prepare Metadata
###############################################################

metadata <- as.data.frame(metadata)

rownames(metadata) <- metadata$Sample

metadata <- metadata[, c("Condition","Platform")]

metadata$Condition <- factor(
  metadata$Condition,
  levels = c("Healthy","AML")
)

metadata$Platform <- factor(metadata$Platform)

metadata <- metadata[colnames(count_matrix),]

if(!all(colnames(count_matrix)==rownames(metadata))){
  
  stop("Sample names do not match!")
  
}

cat("\nMetadata verified successfully.\n")

###############################################################
# 7. Filter Low Count Genes
###############################################################

keep <- rowSums(count_matrix >=10) >=3

count_matrix_filtered <- count_matrix[keep,]

cat("\n---------------------------------------\n")

cat("Gene Filtering Summary\n")

cat("---------------------------------------\n")

cat("Before :", nrow(count_matrix), "\n")

cat("After  :", nrow(count_matrix_filtered), "\n")

cat("Removed:", nrow(count_matrix)-nrow(count_matrix_filtered), "\n")

###############################################################
# 8. Create DESeq2 Dataset
###############################################################

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_filtered,
  colData = metadata,
  design = ~ Condition
)

cat("\n---------------------------------------\n")
cat("DESeqDataSet Created Successfully\n")
cat("---------------------------------------\n")

###############################################################
# 9. Run Differential Expression Analysis
###############################################################

cat("\nRunning DESeq2...\n\n")

dds <- DESeq(dds)

cat("\n---------------------------------------\n")
cat("DESeq2 Analysis Completed Successfully\n")
cat("---------------------------------------\n")

###############################################################
# 10. Extract Results
###############################################################

res <- results(dds)

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

cat("\nSummary of Differential Expression\n\n")

print(summary(results(dds)))

###############################################################
# 11. Significant Genes
###############################################################

sig_res <- subset(
  res,
  padj < 0.05 &
    abs(log2FoldChange) >= 1
)

upregulated <- subset(
  sig_res,
  log2FoldChange > 1
)

downregulated <- subset(
  sig_res,
  log2FoldChange < -1
)

cat("\n---------------------------------------\n")
cat("Significant DEG Summary\n")
cat("---------------------------------------\n")

cat("Total Significant :", nrow(sig_res), "\n")
cat("Upregulated       :", nrow(upregulated), "\n")
cat("Downregulated     :", nrow(downregulated), "\n")

###############################################################
# 12. Save Results
###############################################################

write.csv(
  res,
  "results/deseq2/DESeq2_all_results.csv",
  row.names = FALSE
)

write.csv(
  sig_res,
  "results/deseq2/DESeq2_significant_results.csv",
  row.names = FALSE
)

saveRDS(
  dds,
  "results/deseq2/dds.rds"
)

cat("\n---------------------------------------\n")
cat("Files Successfully Saved\n")
cat("---------------------------------------\n")

cat("\nCreated:\n")

cat("results/deseq2/DESeq2_all_results.csv\n")
cat("results/deseq2/DESeq2_significant_results.csv\n")
cat("results/deseq2/dds.rds\n")

