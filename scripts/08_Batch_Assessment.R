##############################################################
# AML Biomarker Discovery Pipeline
#
# Script:
# 08_Batch_Assessment.R
#
# Purpose:
# Assess sequencing platform confounding and document
# why ComBat batch correction cannot be applied.
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
cat("Stage 8 : Batch Assessment\n")
cat("=========================================\n\n")

start_time <- Sys.time()

##############################################################
# Directories
##############################################################

BATCH_DIR <- file.path(
  RESULTS_DIR,
  "batch_assessment"
)

dir.create(
  BATCH_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

##############################################################
# Input File
##############################################################

METADATA_FILE <- file.path(
  "data",
  "metadata",
  "sample_metadata.csv"
)

##############################################################
# Output Files
##############################################################

TABLE_FILE <- file.path(
  BATCH_DIR,
  "Condition_Platform_Table.csv"
)

PNG_FILE <- file.path(
  BATCH_DIR,
  "Condition_vs_Platform.png"
)

PDF_FILE <- file.path(
  BATCH_DIR,
  "Condition_vs_Platform.pdf"
)

REPORT_FILE <- file.path(
  BATCH_DIR,
  "Batch_Assessment_Report.txt"
)

##############################################################
# Validate Input
##############################################################

if (!file.exists(METADATA_FILE)) {
  
  stop(
    "Sample metadata not found.\nRun 02_prepare_metadata.R first."
  )
  
}

##############################################################
# Load Metadata
##############################################################

cat("Loading metadata...\n")

metadata <- readr::read_csv(
  METADATA_FILE,
  show_col_types = FALSE
)

cat(
  "Samples Loaded :",
  nrow(metadata),
  "\n\n"
)

##############################################################
# Contingency Table
##############################################################

distribution <- table(
  
  metadata$Condition,
  
  metadata$Platform
  
)

cat("-----------------------------------------\n")
cat("Condition × Platform\n")
cat("-----------------------------------------\n\n")

print(distribution)

##############################################################
# Confounding Assessment
##############################################################

confounded <-
  
  all(rowSums(distribution > 0) == 1) &
  
  all(colSums(distribution > 0) == 1)

cat("\n")

if (confounded) {
  
  cat("Complete confounding detected.\n\n")
  
} else {
  
  cat("No complete confounding detected.\n\n")
  
}

##############################################################
# Save Contingency Table
##############################################################

readr::write_csv(
  
  as.data.frame.matrix(distribution),
  
  TABLE_FILE
  
)

##############################################################
# Plot Sample Distribution
##############################################################

plot_data <-
  
  dplyr::count(
    
    metadata,
    
    Condition,
    
    Platform
    
  )

p <- ggplot2::ggplot(
  
  plot_data,
  
  ggplot2::aes(
    
    x = Condition,
    
    y = n,
    
    fill = Platform
    
  )
  
) +
  
  ggplot2::geom_col(
    
    width = 0.70
    
  ) +
  
  ggplot2::geom_text(
    
    ggplot2::aes(label = n),
    
    vjust = -0.4,
    
    size = 5
    
  ) +
  
  ggplot2::scale_fill_manual(
    
    values = c(
      
      HiSeq1000 = HEALTHY_BLUE,
      
      NextSeq550 = AML_RED
      
    )
    
  ) +
  
  ggplot2::labs(
    
    title = "Distribution of Samples by Sequencing Platform",
    
    subtitle = "Assessment of Platform Confounding",
    
    x = "Condition",
    
    y = "Number of Samples"
    
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
# Report
##############################################################

report <- c(
  
  "AML Biomarker Discovery Pipeline",
  
  "",
  
  paste("Date:", Sys.time()),
  
  "",
  
  "Batch Assessment Summary",
  
  "",
  
  paste("Total Samples :", nrow(metadata)),
  
  paste("AML Samples :", sum(metadata$Condition == "AML")),
  
  paste("Healthy Samples :", sum(metadata$Condition == "Healthy")),
  
  "",
  
  "Condition versus Platform",
  
  capture.output(distribution),
  
  "",
  
  if (confounded)
    
    c(
      
      "RESULT: COMPLETE CONFOUNDING DETECTED.",
      
      "",
      
      "Sequencing platform is perfectly associated",
      
      "with biological condition.",
      
      "",
      
      "AML samples were generated on NextSeq550.",
      
      "Healthy samples were generated on HiSeq1000.",
      
      "",
      
      "ComBat empirical Bayes correction cannot",
      
      "distinguish technical effects from",
      
      "true biological effects.",
      
      "",
      
      "Batch correction was intentionally omitted.",
      
      "",
      
      "External validation will be used instead."
      
    )
  
  else
    
    "RESULT: No complete confounding detected."
  
)

writeLines(
  
  report,
  
  REPORT_FILE
  
)

##############################################################
# Summary
##############################################################

end_time <- Sys.time()

cat("\n")
cat("=========================================\n")
cat("Stage 8 Completed Successfully\n")
cat("=========================================\n\n")

cat(
  "Contingency Table :",
  TABLE_FILE,
  "\n"
)

cat(
  "PNG Figure        :",
  PNG_FILE,
  "\n"
)

cat(
  "PDF Figure        :",
  PDF_FILE,
  "\n"
)

cat(
  "Report            :",
  REPORT_FILE,
  "\n\n"
)

cat(
  "Complete Confounding :",
  ifelse(confounded, "YES", "NO"),
  "\n\n"
)

cat(
  "Time Elapsed :",
  round(end_time - start_time, 2),
  "\n\n"
)