##############################################################
# AML Biomarker Discovery Pipeline
#
# Script:
# 12_Biomarker_Expression_Plots.R
#
# Purpose:
# Generate publication quality expression plots for
# the Top 20 biomarker candidates.
#
# Author:
# Isreal Oluwafemi Abiodun
#
# Version:
# 1.0.1  (fixed Ensembl ID version-mismatch bug)
##############################################################

##############################################################
# Load Configuration
##############################################################

source("config.R")
source("functions/plotting_functions.R")

cat("\n")
cat("=========================================\n")
cat("Stage 12 : Biomarker Expression Plots\n")
cat("=========================================\n\n")

start_time <- Sys.time()

##############################################################
# Directories
##############################################################

EXPRESSION_DIR <- file.path(
  RESULTS_DIR,
  "expression_plots"
)

dir.create(
  EXPRESSION_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

##############################################################
# Input Files
##############################################################

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

BIOMARKER_FILE <- file.path(
  RESULTS_DIR,
  "biomarkers",
  "Top20_Biomarkers.csv"
)

##############################################################
# Validate Inputs
##############################################################

required_files <- c(
  
  VSD_FILE,
  
  METADATA_FILE,
  
  BIOMARKER_FILE
  
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  
  stop(
    paste(
      "Missing input file(s):",
      paste(missing_files, collapse = "\n")
    )
  )
  
}

##############################################################
# Load Data
##############################################################

cat("Loading input files...\n")

vsd <- readRDS(VSD_FILE)

metadata <- readr::read_csv(
  METADATA_FILE,
  show_col_types = FALSE
)

biomarkers <- readr::read_csv(
  BIOMARKER_FILE,
  show_col_types = FALSE
)

cat(
  "Biomarkers Loaded :",
  nrow(biomarkers),
  "\n\n"
)

##############################################################
# Expression Matrix
##############################################################

expr <- SummarizedExperiment::assay(vsd)

##############################################################
# Standardize Ensembl Gene Identifiers
#
# The VST expression matrix stores versioned Ensembl IDs
# (e.g. "ENSG00000153815.14"), while the biomarker candidate
# table stores bare IDs (e.g. "ENSG00000153815"). Matching
# must be done on the stripped (unversioned) ID so that
# candidates are correctly located in the expression matrix.
##############################################################

strip_version <- function(x) {
  sub("\\..*$", "", x)
}

expr_ids_stripped <- strip_version(rownames(expr))

biomarkers$GeneID_stripped <- strip_version(biomarkers$GeneID)

# Map each stripped biomarker ID to its corresponding row name
# in the expression matrix (first match if duplicates exist)
id_lookup <- rownames(expr)[
  match(biomarkers$GeneID_stripped, expr_ids_stripped)
]

biomarkers$Matrix_RowID <- id_lookup

matched_n <- sum(!is.na(biomarkers$Matrix_RowID))
unmatched_n <- sum(is.na(biomarkers$Matrix_RowID))

cat(
  "Candidate biomarkers :",
  nrow(biomarkers),
  "\n"
)

cat(
  "Matched biomarkers   :",
  matched_n,
  "\n"
)

cat(
  "Unmatched biomarkers :",
  unmatched_n,
  "\n\n"
)

if (unmatched_n > 0) {
  
  cat(
    "Unmatched Gene IDs:\n"
  )
  
  print(
    biomarkers$GeneID[is.na(biomarkers$Matrix_RowID)]
  )
  
  cat("\n")
  
}

##############################################################
# Generate Expression Plots
##############################################################

plots_created <- 0

for (i in seq_len(nrow(biomarkers))) {
  
  gene_label <- biomarkers$GeneID[i]
  
  matrix_row <- biomarkers$Matrix_RowID[i]
  
  if (is.na(matrix_row)) {
    
    warning(
      paste(
        gene_label,
        "not found in expression matrix."
      )
    )
    
    next
    
  }
  
  values <- expr[matrix_row, ]
  
  plot_data <- data.frame(
    
    Sample = metadata$Sample,
    
    Condition = metadata$Condition,
    
    Expression = values
    
  )
  
  p <- ggplot2::ggplot(
    
    plot_data,
    
    ggplot2::aes(
      
      x = Condition,
      
      y = Expression,
      
      fill = Condition
      
    )
    
  ) +
    
    ggplot2::geom_violin(
      
      alpha = 0.35,
      
      trim = FALSE
      
    ) +
    
    ggplot2::geom_boxplot(
      
      width = 0.18,
      
      outlier.shape = NA,
      
      colour = "black"
      
    ) +
    
    ggplot2::geom_jitter(
      
      width = 0.08,
      
      size = 2,
      
      alpha = 0.80
      
    ) +
    
    ggplot2::scale_fill_manual(
      
      values = c(
        
        Healthy = HEALTHY_BLUE,
        
        AML = AML_RED
        
      )
      
    ) +
    
    ggplot2::labs(
      
      title = gene_label,
      
      x = NULL,
      
      y = "VST Expression"
      
    ) +
    
    publication_theme()
  
  # Use the bare (unversioned) gene ID for output filenames
  # so file names stay consistent with the biomarker tables
  # used elsewhere in the pipeline.
  
  png_file <- file.path(
    EXPRESSION_DIR,
    paste0(gene_label, ".png")
  )
  
  pdf_file <- file.path(
    EXPRESSION_DIR,
    paste0(gene_label, ".pdf")
  )
  
  ggplot2::ggsave(
    
    filename = png_file,
    
    plot = p,
    
    width = FIG_WIDTH,
    
    height = FIG_HEIGHT,
    
    dpi = FIG_DPI
    
  )
  
  ggplot2::ggsave(
    
    filename = pdf_file,
    
    plot = p,
    
    width = FIG_WIDTH,
    
    height = FIG_HEIGHT
    
  )
  
  plots_created <- plots_created + 1
  
}

##############################################################
# Summary
##############################################################

end_time <- Sys.time()

cat("\n")
cat("=========================================\n")
cat("Stage 12 Completed Successfully\n")
cat("=========================================\n\n")

cat(
  "Expression Plots Generated :",
  plots_created,
  "\n\n"
)

cat(
  "Output Directory :",
  EXPRESSION_DIR,
  "\n\n"
)

cat(
  "Time Elapsed :",
  round(end_time - start_time, 2),
  "\n\n"
)