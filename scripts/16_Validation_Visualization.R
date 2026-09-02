###############################################################
# AML Biomarker Discovery Pipeline
#
# Script:
# 16_Validation_Visualization.R
#
# Purpose:
# Publication quality visualization of external validation
#
# Author:
# Isreal Oluwafemi Abiodun
###############################################################


###############################################################
# Load Configuration
###############################################################

source("config.R")
source("functions/plotting_functions.R")


###############################################################
# Required Packages
###############################################################

packages <- c(
  "ggplot2",
  "readr",
  "dplyr",
  "ComplexHeatmap",
  "circlize"
)

install_and_load_packages(packages)


###############################################################
# Pipeline Header
###############################################################

cat("\n")
cat("=========================================\n")
cat("Stage 16 : Validation Visualization\n")
cat("=========================================\n\n")


start_time <- Sys.time()


###############################################################
# Load Validation Dataset
###############################################################

validation_dataset <- "TCGA"

cat(
  "Validation Dataset :",
  validation_dataset,
  "\n\n"
)

dataset <- load_validation_dataset(
  validation_dataset
)

expr <- dataset$expression

metadata <- dataset$metadata


###############################################################
# Validate Expression Dataset
###############################################################

if (!"GeneID" %in% names(expr)) {
  
  stop(
    "Expression dataset must contain a GeneID column."
  )
  
}


###############################################################
# Validate Metadata
###############################################################

if (nrow(metadata) == 0) {
  
  stop(
    "Validation metadata contains no samples."
  )
  
}


###############################################################
# Build Expression Matrix
###############################################################

expr_matrix <- as.matrix(
  expr[
    ,
    setdiff(
      names(expr),
      "GeneID"
    ),
    drop = FALSE
  ]
)


###############################################################
# Assign Gene IDs
###############################################################

rownames(expr_matrix) <- expr$GeneID


###############################################################
# Validate Sample Dimensions
###############################################################

if (ncol(expr_matrix) != nrow(metadata)) {
  
  stop(
    paste0(
      "Expression and metadata sample counts do not match.\n",
      "Expression samples : ",
      ncol(expr_matrix),
      "\n",
      "Metadata samples   : ",
      nrow(metadata)
    )
  )
  
}


###############################################################
# Load Validation Results
###############################################################

detection <- readr::read_csv(
  "results/validation/Biomarker_Detection.csv",
  show_col_types = FALSE
)

statistics <- readr::read_csv(
  "results/validation/Expression_Statistics.csv",
  show_col_types = FALSE
)


###############################################################
# Validate Detection Table
###############################################################

if (!"GeneID" %in% names(detection)) {
  
  stop(
    "Biomarker_Detection.csv must contain a GeneID column."
  )
  
}


###############################################################
# Validate Statistics Table
###############################################################

if (!"GeneID" %in% names(statistics)) {
  
  stop(
    "Expression_Statistics.csv must contain a GeneID column."
  )
  
}


###############################################################
# Output Folder
###############################################################

output_dir <- file.path(
  RESULTS_DIR,
  "validation",
  "plots"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


###############################################################
# Standardise Gene IDs
###############################################################

clean_gene_id <- function(x) {
  
  x <- as.character(x)
  
  sub(
    "\\..*$",
    "",
    x
  )
  
}


###############################################################
# Clean Expression Gene IDs
###############################################################

expr_ids_clean <- clean_gene_id(
  rownames(expr_matrix)
)


###############################################################
# Clean Detection Gene IDs
###############################################################

detection_ids_clean <- clean_gene_id(
  detection$GeneID
)


###############################################################
# Match Detection Biomarkers
###############################################################

detection_match <- match(
  detection_ids_clean,
  expr_ids_clean
)


###############################################################
# Report Matching
###############################################################

cat(
  "Candidate biomarkers :",
  length(detection_ids_clean),
  "\n"
)

cat(
  "Matched biomarkers   :",
  sum(!is.na(detection_match)),
  "\n"
)

cat(
  "Unmatched biomarkers :",
  sum(is.na(detection_match)),
  "\n\n"
)


###############################################################
# Stop If Nothing Matches
###############################################################

if (all(is.na(detection_match))) {
  
  stop(
    paste(
      "No biomarkers in Biomarker_Detection.csv",
      "could be matched to the validation",
      "expression matrix."
    )
  )
  
}


###############################################################
# Keep Matched Biomarkers
###############################################################

matched_detection <- detection[
  !is.na(detection_match),
  ,
  drop = FALSE
]

matched_indices <- detection_match[
  !is.na(detection_match)
]


###############################################################
# Expression Matrix For Validated Biomarkers
###############################################################

validated_expr_matrix <- expr_matrix[
  matched_indices,
  ,
  drop = FALSE
]


###############################################################
# Use Biomarker IDs As Row Names
###############################################################

rownames(validated_expr_matrix) <-
  matched_detection$GeneID


###############################################################
# Remove Zero Variance Genes
###############################################################

row_sd <- apply(
  validated_expr_matrix,
  1,
  sd,
  na.rm = TRUE
)


validated_expr_matrix <- validated_expr_matrix[
  is.finite(row_sd) &
    row_sd > 0,
  ,
  drop = FALSE
]


###############################################################
# Detection Rate Plot
###############################################################

p1 <-
  
  detection %>%
  
  arrange(
    desc(DetectionRate)
  ) %>%
  
  ggplot(
    
    aes(
      reorder(
        GeneID,
        DetectionRate
      ),
      DetectionRate
    )
    
  ) +
  
  geom_col(
    fill = AML_RED
  ) +
  
  coord_flip() +
  
  labs(
    
    title = "Biomarker Detection Rate",
    
    x = "",
    
    y = "Detection Rate"
    
  ) +
  
  publication_theme()


ggsave(
  
  file.path(
    output_dir,
    "Detection_Rate.png"
  ),
  
  p1,
  
  width = FIG_WIDTH,
  
  height = FIG_HEIGHT,
  
  dpi = FIG_DPI
  
)


ggsave(
  
  file.path(
    output_dir,
    "Detection_Rate.pdf"
  ),
  
  p1,
  
  width = FIG_WIDTH,
  
  height = FIG_HEIGHT
  
)


###############################################################
# Mean Expression Plot
###############################################################

if ("Mean" %in% names(statistics)) {
  
  p2 <-
    
    statistics %>%
    
    arrange(
      desc(Mean)
    ) %>%
    
    ggplot(
      
      aes(
        reorder(
          GeneID,
          Mean
        ),
        Mean
      )
      
    ) +
    
    geom_col(
      fill = HEALTHY_BLUE
    ) +
    
    coord_flip() +
    
    labs(
      
      title = "Mean Expression",
      
      x = "",
      
      y = "Mean Expression"
      
    ) +
    
    publication_theme()
  
  
  ggsave(
    
    file.path(
      output_dir,
      "Mean_Expression.png"
    ),
    
    p2,
    
    width = FIG_WIDTH,
    
    height = FIG_HEIGHT,
    
    dpi = FIG_DPI
    
  )
  
  
  ggsave(
    
    file.path(
      output_dir,
      "Mean_Expression.pdf"
    ),
    
    p2,
    
    width = FIG_WIDTH,
    
    height = FIG_HEIGHT
    
  )
  
}


###############################################################
# Top 20 Detection Plot
###############################################################

top20 <-
  
  detection %>%
  
  arrange(
    desc(DetectionRate)
  ) %>%
  
  slice_head(
    n = 20
  )


p3 <-
  
  ggplot(
    
    top20,
    
    aes(
      
      reorder(
        GeneID,
        DetectionRate
      ),
      
      DetectionRate
      
    )
    
  ) +
  
  geom_col(
    fill = "#2E8B57"
  ) +
  
  coord_flip() +
  
  labs(
    
    title = "Top 20 Validated Biomarkers",
    
    x = "",
    
    y = "Detection Rate"
    
  ) +
  
  publication_theme()


ggsave(
  
  file.path(
    output_dir,
    "Top20_Detection.png"
  ),
  
  p3,
  
  width = FIG_WIDTH,
  
  height = FIG_HEIGHT,
  
  dpi = FIG_DPI
  
)


ggsave(
  
  file.path(
    output_dir,
    "Top20_Detection.pdf"
  ),
  
  p3,
  
  width = FIG_WIDTH,
  
  height = FIG_HEIGHT
  
)


###############################################################
# Heatmap
###############################################################

cat(
  "Generating validation heatmap...\n"
)


###############################################################
# Scale Expression
###############################################################

expr_scaled <- t(
  scale(
    t(validated_expr_matrix)
  )
)


###############################################################
# Replace Non Finite Values
###############################################################

expr_scaled[
  !is.finite(expr_scaled)
] <- 0


###############################################################
# Generate Heatmap
###############################################################

heatmap_plot <-
  
  ComplexHeatmap::Heatmap(
    
    expr_scaled,
    
    name = "Z-score",
    
    cluster_rows = TRUE,
    
    cluster_columns = TRUE,
    
    show_row_names = TRUE,
    
    show_column_names = FALSE,
    
    column_title = paste(
      validation_dataset,
      "Validation"
    )
    
  )


###############################################################
# Save Heatmap PNG
###############################################################

png(
  
  file.path(
    output_dir,
    "Validation_Heatmap.png"
  ),
  
  width = 2200,
  
  height = 1800,
  
  res = 300
  
)

ComplexHeatmap::draw(
  heatmap_plot
)

dev.off()


###############################################################
# Save Heatmap PDF
###############################################################

pdf(
  
  file.path(
    output_dir,
    "Validation_Heatmap.pdf"
  ),
  
  width = 10,
  
  height = 8
  
)

ComplexHeatmap::draw(
  heatmap_plot
)

dev.off()


###############################################################
# Boxplots
###############################################################

if ("Condition" %in% names(metadata)) {
  
  cat(
    "Generating biomarker boxplots...\n"
  )
  
  
  selected <- matched_detection$GeneID[
    seq_len(
      min(
        10,
        nrow(matched_detection)
      )
    )
  ]
  
  
  for (gene in selected) {
    
    
    ###########################################################
    # Match Gene
    ###########################################################
    
    gene_clean <- clean_gene_id(
      gene
    )
    
    
    gene_index <- match(
      gene_clean,
      clean_gene_id(
        rownames(expr_matrix)
      )
    )
    
    
    if (is.na(gene_index)) {
      
      next
      
    }
    
    
    ###########################################################
    # Expression Values
    ###########################################################
    
    gene_values <- as.numeric(
      expr_matrix[
        gene_index,
      ]
    )
    
    
    ###########################################################
    # Build Plot Data
    ###########################################################
    
    plot_data <- data.frame(
      
      Condition = metadata$Condition,
      
      Expression = gene_values,
      
      stringsAsFactors = FALSE
      
    )
    
    
    ###########################################################
    # Generate Plot
    ###########################################################
    
    p <-
      
      ggplot(
        
        plot_data,
        
        aes(
          
          Condition,
          
          Expression,
          
          fill = Condition
          
        )
        
      ) +
      
      geom_violin(
        alpha = 0.3
      ) +
      
      geom_boxplot(
        width = 0.15
      ) +
      
      geom_jitter(
        width = 0.08
      ) +
      
      scale_fill_manual(
        
        values = c(
          
          Healthy = HEALTHY_BLUE,
          
          AML = AML_RED
          
        ),
        
        na.value = "grey70"
        
      ) +
      
      labs(
        
        title = gene,
        
        x = "",
        
        y = "Expression"
        
      ) +
      
      publication_theme()
    
    
    ###########################################################
    # Save Plot
    ###########################################################
    
    ggsave(
      
      file.path(
        
        output_dir,
        
        paste0(
          gene,
          "_boxplot.png"
        )
        
      ),
      
      p,
      
      width = FIG_WIDTH,
      
      height = FIG_HEIGHT,
      
      dpi = FIG_DPI
      
    )
    
  }
  
}


###############################################################
# Summary
###############################################################

end_time <- Sys.time()


cat("\n")
cat("---------------------------------------\n")
cat("Validation Visualisation Completed\n")
cat("---------------------------------------\n\n")

cat(
  "Validation Dataset :",
  validation_dataset,
  "\n"
)

cat(
  "Expression Genes   :",
  nrow(expr_matrix),
  "\n"
)

cat(
  "Validation Samples :",
  ncol(expr_matrix),
  "\n"
)

cat(
  "Candidate Biomarkers :",
  length(detection_ids_clean),
  "\n"
)

cat(
  "Matched Biomarkers   :",
  sum(!is.na(detection_match)),
  "\n"
)

cat(
  "Output Folder:\n"
)

cat(
  output_dir,
  "\n"
)

cat(
  "Time Elapsed :",
  round(
    end_time - start_time,
    2
  ),
  "\n\n"
)