###############################################################
# AML Biomarker Discovery Pipeline
# Script: 09_Biomarker_Candidate_Selection.R
# Purpose:
# Identify the strongest AML biomarker candidates
###############################################################

# rm(list = ls())

source("config.R")

packages <- c(
  "readr",
  "dplyr"
)

for(pkg in packages){
  
  if(!require(pkg, character.only = TRUE)){
    
    install.packages(pkg)
    
    library(pkg, character.only = TRUE)
    
  }
  
}

###############################################################
# Load DESeq2 Results
###############################################################

results <- read_csv(
  "results/deseq2/DESeq2_significant_results.csv",
  show_col_types = FALSE
)

cat("---------------------------------------\n")
cat("Significant DEGs Loaded\n")
cat("---------------------------------------\n")

cat("Genes:", nrow(results), "\n\n")

###############################################################
# Remove Missing Values
###############################################################

results <- results %>%
  
  filter(!is.na(padj))

###############################################################
# Rank Biomarkers
###############################################################

biomarkers <- results %>%
  
  arrange(
    padj,
    desc(abs(log2FoldChange))
  )

###############################################################
# Select Top 50
###############################################################

top50 <- biomarkers %>%
  
  dplyr::slice(1:50)

###############################################################
# Top 20
###############################################################

top20 <- biomarkers %>%
  
  dplyr::slice(1:20)

###############################################################
# Upregulated
###############################################################

top_up <- biomarkers %>%
  
  filter(log2FoldChange > 0) %>%
  
  dplyr::slice(1:25)

###############################################################
# Downregulated
###############################################################

top_down <- biomarkers %>%
  
  filter(log2FoldChange < 0) %>%
  
  dplyr::slice(1:25)

###############################################################
# Create Output Folder
###############################################################

dir.create(
  "results/biomarkers",
  recursive = TRUE,
  showWarnings = FALSE
)

###############################################################
# Save Results
###############################################################

write.csv(
  biomarkers,
  "results/biomarkers/Ranked_Biomarkers.csv",
  row.names = FALSE
)

write.csv(
  top50,
  "results/biomarkers/Top50_Biomarkers.csv",
  row.names = FALSE
)

write.csv(
  top20,
  "results/biomarkers/Top20_Biomarkers.csv",
  row.names = FALSE
)

write.csv(
  top_up,
  "results/biomarkers/Top25_Upregulated.csv",
  row.names = FALSE
)

write.csv(
  top_down,
  "results/biomarkers/Top25_Downregulated.csv",
  row.names = FALSE
)

###############################################################
# Summary
###############################################################

cat("---------------------------------------\n")
cat("Biomarker Candidates Selected\n")
cat("---------------------------------------\n\n")

cat("Total Significant Genes :", nrow(results), "\n")
cat("Top Biomarkers Saved    :", nrow(top50), "\n")
cat("Top Upregulated         :", nrow(top_up), "\n")
cat("Top Downregulated       :", nrow(top_down), "\n\n")

cat("Files created:\n")
cat("results/biomarkers/Ranked_Biomarkers.csv\n")
cat("results/biomarkers/Top50_Biomarkers.csv\n")
cat("results/biomarkers/Top20_Biomarkers.csv\n")
cat("results/biomarkers/Top25_Upregulated.csv\n")
cat("results/biomarkers/Top25_Downregulated.csv\n")
