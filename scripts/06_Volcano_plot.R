###############################################################
# AML Biomarker Discovery Pipeline
# Script: 06_Volcano_plot.R
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
  "ggplot2",
  "ggrepel",
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
# 4. Load Results
###############################################################

results <- read_csv(
  "results/deseq2/DESeq2_all_results.csv",
  show_col_types = FALSE
)

cat("---------------------------------------\n")
cat("DESeq2 Results Loaded\n")
cat("---------------------------------------\n")

cat("Genes:", nrow(results), "\n")

###############################################################
# 5. Prepare Volcano Data
###############################################################

results <- results %>%
  
  filter(!is.na(padj))

results$Significance <- "Not Significant"

results$Significance[
  results$padj < 0.05 &
    results$log2FoldChange >= 1
] <- "Up"

results$Significance[
  results$padj < 0.05 &
    results$log2FoldChange <= -1
] <- "Down"

results$minusLog10Padj <- -log10(results$padj)

table(results$Significance)

###############################################################
# 6. Label Biomarkers
###############################################################

genes_to_label <- c(
  
  "SNHG15",
  
  "MMP13",
  
  "XPO1"
  
)

label_data <- results %>%
  
  filter(GeneID %in% genes_to_label)

###############################################################
# 7. Volcano Plot
###############################################################

p <- ggplot(
  
  results,
  
  aes(
    
    x = log2FoldChange,
    
    y = minusLog10Padj,
    
    colour = Significance
    
  )
  
) +
  
  geom_point(
    
    alpha = 0.7,
    
    size = 2
    
  ) +
  
  scale_colour_manual(
    
    values = c(
      
      "Up" = AML_RED,
      
      "Down" = HEALTHY_BLUE,
      
      "Not Significant" = GREY
      
    )
    
  ) +
  
  geom_vline(
    
    xintercept = c(-1,1),
    
    linetype = 2
    
  ) +
  
  geom_hline(
    
    yintercept = -log10(0.05),
    
    linetype = 2
    
  ) +
  
  geom_text_repel(
    
    data = label_data,
    
    aes(label = GeneID),
    
    size = 4,
    
    max.overlaps = Inf
    
  ) +
  
  ggtitle("Differential Gene Expression in AML") +
  
  xlab("Log2 Fold Change") +
  
  ylab("-Log10 Adjusted P-value") +
  
  publication_theme()

###############################################################
# 8. Save Figure
###############################################################

ggsave(
  
  filename = file.path(
    
    VOLCANO_DIR,
    
    "Volcano_plot.png"
    
  ),
  
  plot = p,
  
  width = FIG_WIDTH,
  
  height = FIG_HEIGHT,
  
  dpi = FIG_DPI
  
)

ggsave(
  
  filename = file.path(
    
    VOLCANO_DIR,
    
    "Volcano_plot.pdf"
    
  ),
  
  plot = p,
  
  width = FIG_WIDTH,
  
  height = FIG_HEIGHT
  
)

cat("\n---------------------------------------\n")

cat("Volcano Plot Completed Successfully\n")

cat("---------------------------------------\n")