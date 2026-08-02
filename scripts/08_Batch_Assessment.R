###############################################################
# AML Biomarker Discovery Pipeline
# Script: 08_Batch_Assessment.R
# Purpose:
# Assess sequencing platform confounding and document why
# ComBat batch correction cannot be applied.
###############################################################

#==============================================================
# 1. INITIALISE ENVIRONMENT
#==============================================================

# rm(list = ls())

source("config.R")
source("functions/plotting_functions.R")

packages <- c(
  "readr",
  "dplyr",
  "ggplot2"
)

for(pkg in packages){
  
  if(!require(pkg, character.only = TRUE)){
    
    install.packages(pkg)
    
    library(pkg, character.only = TRUE)
    
  }
  
}

#==============================================================
# 2. LOAD DATA
#==============================================================

metadata <- read_csv(
  "data/metadata/sample_metadata.csv",
  show_col_types = FALSE
)

cat("---------------------------------------\n")
cat("Metadata Loaded Successfully\n")
cat("---------------------------------------\n\n")

cat("Samples:", nrow(metadata), "\n\n")

#==============================================================
# 3. SUMMARISE SAMPLE DISTRIBUTION
#==============================================================

distribution <- table(
  metadata$Condition,
  metadata$Platform
)

print(distribution)

#==============================================================
# 4. CHECK FOR COMPLETE CONFOUNDING
#==============================================================

confounded <- all(rowSums(distribution > 0) == 1) &
  all(colSums(distribution > 0) == 1)

cat("\n")

if(confounded){
  
  cat("Complete confounding detected.\n")
  
}else{
  
  cat("No complete confounding detected.\n")
  
}

#==============================================================
# 5. CREATE OUTPUT DIRECTORY
#==============================================================

dir.create(
  "results/batch_assessment",
  recursive = TRUE,
  showWarnings = FALSE
)

#==============================================================
# 6. SAVE CONTINGENCY TABLE
#==============================================================

write.csv(
  as.data.frame.matrix(distribution),
  "results/batch_assessment/Condition_Platform_Table.csv"
)

#==============================================================
# 7. CREATE BARPLOT
#==============================================================

plot_data <- metadata %>%
  dplyr::count(Condition, Platform)

p <- ggplot(
  plot_data,
  aes(
    x = Condition,
    y = n,
    fill = Platform
  )
) +
  
  geom_col(
    width = 0.7
  ) +
  
  geom_text(
    aes(label = n),
    vjust = -0.4,
    size = 5
  ) +
  
  scale_fill_manual(
    values = c(
      HiSeq1000 = "#4F81BD",
      NextSeq550 = "#C0504D"
    )
  ) +
  
  labs(
    title = "Distribution of Samples by Condition and Sequencing Platform",
    x = "Condition",
    y = "Number of Samples"
  ) +
  
  publication_theme()

ggsave(
  filename = "results/batch_assessment/Condition_vs_Platform.png",
  plot = p,
  width = FIG_WIDTH,
  height = FIG_HEIGHT,
  dpi = FIG_DPI
)

ggsave(
  filename = "results/batch_assessment/Condition_vs_Platform.pdf",
  plot = p,
  width = FIG_WIDTH,
  height = FIG_HEIGHT
)

#==============================================================
# 8. WRITE REPORT
#==============================================================

report <- c(
  
  "AML Biomarker Discovery Pipeline",
  
  "", 
  
  paste("Date:", Sys.time()),
  
  "",
  
  "Batch Assessment Summary",
  
  "",
  
  paste("Total samples:", nrow(metadata)),
  
  paste("AML samples:", sum(metadata$Condition == "AML")),
  
  paste("Healthy samples:", sum(metadata$Condition == "Healthy")),
  
  "",
  
  "Condition vs Platform:",
  
  capture.output(distribution),
  
  "",
  
  if(confounded){
    
    c(
      "Result: COMPLETE CONFOUNDING DETECTED.",
      "",
      "Interpretation:",
      "Sequencing platform is perfectly associated with biological condition.",
      "AML samples originate exclusively from NextSeq550.",
      "Healthy samples originate exclusively from HiSeq1000.",
      "",
      "ComBat empirical Bayes correction cannot separate technical effects from biological effects.",
      "",
      "Therefore ComBat was intentionally NOT applied.",
      "",
      "Independent external validation will instead be used to demonstrate robustness of identified biomarkers."
    )
    
  } else {
    
    "Result: No complete confounding detected."
    
  }
  
)

writeLines(
  report,
  "results/batch_assessment/Batch_Assessment_Report.txt"
)

#==============================================================
# 9. SUMMARY
#==============================================================

cat("\n---------------------------------------\n")
cat("Batch Assessment Completed Successfully\n")
cat("---------------------------------------\n\n")

cat("Files created:\n")

cat("results/batch_assessment/Condition_vs_Platform.png\n")
cat("results/batch_assessment/Condition_vs_Platform.pdf\n")
cat("results/batch_assessment/Condition_Platform_Table.csv\n")
cat("results/batch_assessment/Batch_Assessment_Report.txt\n")
