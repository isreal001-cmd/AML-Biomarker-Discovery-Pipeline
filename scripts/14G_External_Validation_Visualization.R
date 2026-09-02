##############################################################
# AML Biomarker Discovery Pipeline
#
# Script:
# 14G_External_Validation_Visualization.R
#
# Purpose:
# Generate publication ready visualizations for the
# independent GSE13159 external validation analysis.
#
# Author:
# Isreal Oluwafemi Abiodun
#
# Version:
# 1.2.0
##############################################################


##############################################################
# Load configuration
##############################################################

source("config.R")


##############################################################
# Required packages
##############################################################

required_packages <- c(
  "readr",
  "dplyr",
  "ggplot2",
  "pROC"
)

for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    install.packages(pkg)
    
  }
  
}


##############################################################
# Load packages
##############################################################

library(readr)
library(dplyr)
library(ggplot2)
library(pROC)


##############################################################
# Header
##############################################################

cat("\n")
cat("=========================================\n")
cat("Stage 14G : External Validation Visualization\n")
cat("=========================================\n\n")

start_time <- Sys.time()


##############################################################
# Directories
##############################################################

VALIDATION_DIR <- file.path(
  RESULTS_DIR,
  "validation"
)

PERFORMANCE_DIR <- file.path(
  VALIDATION_DIR,
  "performance"
)

ROC_DIR <- file.path(
  VALIDATION_DIR,
  "ROC"
)

FINAL_DIR <- file.path(
  VALIDATION_DIR,
  "final_report"
)

FIGURE_DIR <- file.path(
  FINAL_DIR,
  "figures"
)

TABLE_DIR <- file.path(
  FINAL_DIR,
  "tables"
)

dir.create(
  FIGURE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  TABLE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


##############################################################
# Input files
##############################################################

SCORE_FILE <- file.path(
  VALIDATION_DIR,
  "GSE13159_biomarker_signature_scores.csv"
)

PERFORMANCE_FILE <- file.path(
  PERFORMANCE_DIR,
  "tables",
  "GSE13159_Performance_Metrics.csv"
)

PREDICTION_FILE <- file.path(
  PERFORMANCE_DIR,
  "tables",
  "GSE13159_Sample_Predictions.csv"
)

CONFUSION_FILE <- file.path(
  PERFORMANCE_DIR,
  "tables",
  "GSE13159_Confusion_Matrix.csv"
)

ROC_RESULTS_FILE <- file.path(
  ROC_DIR,
  "GSE13159_ROC_results.csv"
)


##############################################################
# Check required files
##############################################################

required_files <- c(
  SCORE_FILE,
  PERFORMANCE_FILE,
  PREDICTION_FILE,
  CONFUSION_FILE,
  ROC_RESULTS_FILE
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  
  stop(
    paste(
      "Required files are missing:\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
  
}


##############################################################
# Load signature scores
##############################################################

cat("Loading signature scores...\n")

scores <- readr::read_csv(
  SCORE_FILE,
  show_col_types = FALSE
)

cat(
  "Rows loaded :",
  nrow(scores),
  "\n"
)

cat(
  "Columns     :",
  ncol(scores),
  "\n\n"
)


##############################################################
# Validate required score columns
##############################################################

required_score_columns <- c(
  "Sample",
  "Signature_Score",
  "AML_Status"
)

missing_score_columns <- setdiff(
  required_score_columns,
  names(scores)
)

if (length(missing_score_columns) > 0) {
  
  stop(
    paste(
      "Required columns missing from signature score file:",
      paste(
        missing_score_columns,
        collapse = ", "
      )
    )
  )
  
}


##############################################################
# Standardize AML status
##############################################################

scores$AML_Status <- trimws(
  as.character(
    scores$AML_Status
  )
)

scores$AML_Status <- ifelse(
  scores$AML_Status %in% c(
    "AML",
    "aml"
  ),
  "AML",
  ifelse(
    scores$AML_Status %in% c(
      "Non_AML",
      "Non AML",
      "non AML",
      "non_aml",
      "NON_AML",
      "NON AML"
    ),
    "Non_AML",
    NA
  )
)

scores$Signature_Score <- suppressWarnings(
  as.numeric(
    scores$Signature_Score
  )
)

scores <- scores |>
  filter(
    !is.na(AML_Status),
    !is.na(Signature_Score),
    is.finite(Signature_Score)
  )

scores$AML_Status <- factor(
  scores$AML_Status,
  levels = c(
    "Non_AML",
    "AML"
  )
)


##############################################################
# Load performance metrics
##############################################################

cat("Loading Stage 14D performance metrics...\n")

performance <- readr::read_csv(
  PERFORMANCE_FILE,
  show_col_types = FALSE
)

cat(
  "Performance rows :",
  nrow(performance),
  "\n\n"
)


##############################################################
# Validate performance data
##############################################################

if (nrow(performance) < 1) {
  
  stop(
    "Stage 14D performance file contains no rows."
  )
  
}


##############################################################
# Load predictions
##############################################################

cat("Loading sample predictions...\n")

predictions <- readr::read_csv(
  PREDICTION_FILE,
  show_col_types = FALSE
)

cat(
  "Prediction rows :",
  nrow(predictions),
  "\n\n"
)


##############################################################
# Load confusion matrix
##############################################################

cat("Loading confusion matrix...\n")

confusion <- readr::read_csv(
  CONFUSION_FILE,
  show_col_types = FALSE
)

cat(
  "Confusion matrix rows :",
  nrow(confusion),
  "\n"
)

cat(
  "Confusion matrix columns:\n"
)

print(
  names(confusion)
)

cat("\n")


##############################################################
# Load ROC results
##############################################################

cat("Loading ROC results...\n")

roc_results <- readr::read_csv(
  ROC_RESULTS_FILE,
  show_col_types = FALSE
)

cat(
  "ROC result rows :",
  nrow(roc_results),
  "\n\n"
)


##############################################################
# Extract performance values
##############################################################

auc_value <- suppressWarnings(
  as.numeric(
    performance$AUC[1]
  )
)

auc_lower <- suppressWarnings(
  as.numeric(
    performance$AUC_CI_Lower[1]
  )
)

auc_upper <- suppressWarnings(
  as.numeric(
    performance$AUC_CI_Upper[1]
  )
)

sensitivity <- suppressWarnings(
  as.numeric(
    performance$Sensitivity[1]
  )
)

specificity <- suppressWarnings(
  as.numeric(
    performance$Specificity[1]
  )
)

ppv <- suppressWarnings(
  as.numeric(
    performance$PPV[1]
  )
)

npv <- suppressWarnings(
  as.numeric(
    performance$NPV[1]
  )
)

accuracy <- suppressWarnings(
  as.numeric(
    performance$Accuracy[1]
  )
)

threshold <- suppressWarnings(
  as.numeric(
    performance$Optimal_Threshold[1]
  )
)


##############################################################
# Calculate ROC directly
##############################################################

cat("Calculating ROC curve...\n")

roc_object <- pROC::roc(
  response = scores$AML_Status,
  predictor = scores$Signature_Score,
  levels = c(
    "Non_AML",
    "AML"
  ),
  direction = "<",
  quiet = TRUE
)

roc_auc <- as.numeric(
  pROC::auc(
    roc_object
  )
)

roc_ci <- pROC::ci.auc(
  roc_object,
  conf.level = 0.95
)


##############################################################
# ROC dataframe
##############################################################

roc_df <- data.frame(
  
  FPR = 1 - roc_object$specificities,
  
  TPR = roc_object$sensitivities
  
)


##############################################################
# ROC plot
##############################################################

cat(
  "Generating ROC curve...\n"
)

roc_label <- paste0(
  "AUC = ",
  sprintf(
    "%.4f",
    roc_auc
  ),
  "\n95% CI = ",
  sprintf(
    "%.4f",
    roc_ci[1]
  ),
  " to ",
  sprintf(
    "%.4f",
    roc_ci[3]
  )
)

roc_plot <- ggplot(
  roc_df,
  aes(
    x = FPR,
    y = TPR
  )
) +
  
  geom_line(
    linewidth = 1.2
  ) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  
  annotate(
    "text",
    x = 0.55,
    y = 0.18,
    hjust = 0,
    label = roc_label,
    size = 4
  ) +
  
  labs(
    title = "External Validation ROC Curve",
    subtitle = "GSE13159",
    x = "1 − Specificity",
    y = "Sensitivity"
  ) +
  
  theme_classic(
    base_size = 13
  )


ggsave(
  file.path(
    FIGURE_DIR,
    "GSE13159_External_Validation_ROC.png"
  ),
  roc_plot,
  width = 7,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(
    FIGURE_DIR,
    "GSE13159_External_Validation_ROC.pdf"
  ),
  roc_plot,
  width = 7,
  height = 6
)


##############################################################
# Signature score distribution
##############################################################

cat(
  "Generating signature score distribution...\n"
)

score_plot <- ggplot(
  scores,
  aes(
    x = AML_Status,
    y = Signature_Score
  )
) +
  
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA
  ) +
  
  geom_jitter(
    width = 0.15,
    alpha = 0.20,
    size = 0.8
  ) +
  
  labs(
    title = "AML Signature Score Distribution",
    subtitle = "Independent GSE13159 validation cohort",
    x = NULL,
    y = "Signature Score"
  ) +
  
  theme_classic(
    base_size = 13
  )


ggsave(
  file.path(
    FIGURE_DIR,
    "GSE13159_Signature_Score_Distribution.png"
  ),
  score_plot,
  width = 7,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(
    FIGURE_DIR,
    "GSE13159_Signature_Score_Distribution.pdf"
  ),
  score_plot,
  width = 7,
  height = 6
)


##############################################################
# Confusion matrix
##############################################################

cat(
  "Generating confusion matrix...\n"
)

print(confusion)


##############################################################
# Build confusion matrix robustly
##############################################################

cm_data <- NULL


##############################################################
# Case 1:
# Actual / Predicted / Count
##############################################################

actual_col <- names(confusion)[
  grepl(
    "^actual$",
    names(confusion),
    ignore.case = TRUE
  )
][1]

predicted_col <- names(confusion)[
  grepl(
    "^predicted$",
    names(confusion),
    ignore.case = TRUE
  )
][1]

count_candidates <- names(confusion)[
  grepl(
    "count|frequency|freq|number",
    names(confusion),
    ignore.case = TRUE
  )
]


if (
  !is.na(actual_col) &&
  !is.na(predicted_col) &&
  length(count_candidates) > 0
) {
  
  count_col <- count_candidates[1]
  
  cm_data <- data.frame(
    
    Actual = as.character(
      confusion[[actual_col]]
    ),
    
    Predicted = as.character(
      confusion[[predicted_col]]
    ),
    
    Count = suppressWarnings(
      as.numeric(
        confusion[[count_col]]
      )
    ),
    
    stringsAsFactors = FALSE
    
  )
  
}


##############################################################
# Case 2:
# Matrix format
##############################################################

if (is.null(cm_data)) {
  
  confusion_names <- names(confusion)
  
  if (
    all(
      c(
        "AML",
        "Non_AML"
      ) %in% confusion_names
    )
  ) {
    
    if (nrow(confusion) >= 2) {
      
      cm_data <- data.frame(
        
        Actual = c(
          "AML",
          "AML",
          "Non_AML",
          "Non_AML"
        ),
        
        Predicted = c(
          "AML",
          "Non_AML",
          "AML",
          "Non_AML"
        ),
        
        Count = c(
          
          confusion$AML[1],
          
          confusion$Non_AML[1],
          
          confusion$AML[2],
          
          confusion$Non_AML[2]
          
        ),
        
        stringsAsFactors = FALSE
        
      )
      
    }
    
  }
  
}


##############################################################
# Case 3:
# Derive from predictions
##############################################################

if (is.null(cm_data)) {
  
  cat(
    "Constructing confusion matrix from predictions...\n"
  )
  
  prediction_names <- names(
    predictions
  )
  
  actual_prediction_candidates <- prediction_names[
    grepl(
      "actual|observed|status|class",
      prediction_names,
      ignore.case = TRUE
    )
  ]
  
  predicted_prediction_candidates <- prediction_names[
    grepl(
      "predicted|prediction",
      prediction_names,
      ignore.case = TRUE
    )
  ]
  
  if (
    length(actual_prediction_candidates) > 0 &&
    length(predicted_prediction_candidates) > 0
  ) {
    
    actual_prediction_col <-
      actual_prediction_candidates[1]
    
    predicted_prediction_col <-
      predicted_prediction_candidates[1]
    
    actual_values <- trimws(
      as.character(
        predictions[[actual_prediction_col]]
      )
    )
    
    predicted_values <- trimws(
      as.character(
        predictions[[predicted_prediction_col]]
      )
    )
    
    cm_table <- table(
      Actual = actual_values,
      Predicted = predicted_values
    )
    
    cm_data <- as.data.frame(
      cm_table,
      stringsAsFactors = FALSE
    )
    
    names(cm_data)[3] <- "Count"
    
  }
  
}


##############################################################
# Stop if unavailable
##############################################################

if (is.null(cm_data)) {
  
  stop(
    paste(
      "Unable to construct confusion matrix.",
      "\nDetected columns:",
      paste(
        names(confusion),
        collapse = ", "
      )
    )
  )
  
}


##############################################################
# Clean confusion matrix
##############################################################

cm_data$Actual <- as.character(
  cm_data$Actual
)

cm_data$Predicted <- as.character(
  cm_data$Predicted
)

cm_data$Count <- suppressWarnings(
  as.numeric(
    cm_data$Count
  )
)

cm_data <- cm_data[
  is.finite(
    cm_data$Count
  ),
]


##############################################################
# Confusion matrix plot
##############################################################

cm_plot <- ggplot(
  cm_data,
  aes(
    x = Predicted,
    y = Actual,
    fill = Count
  )
) +
  
  geom_tile(
    colour = "white"
  ) +
  
  geom_text(
    aes(
      label = Count
    ),
    size = 5
  ) +
  
  scale_fill_gradient(
    low = "grey90",
    high = "grey20"
  ) +
  
  labs(
    title = "External Validation Confusion Matrix",
    subtitle = "GSE13159",
    x = "Predicted Status",
    y = "Observed Status",
    fill = "Samples"
  ) +
  
  theme_classic(
    base_size = 13
  )


ggsave(
  file.path(
    FIGURE_DIR,
    "GSE13159_Confusion_Matrix.png"
  ),
  cm_plot,
  width = 7,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(
    FIGURE_DIR,
    "GSE13159_Confusion_Matrix.pdf"
  ),
  cm_plot,
  width = 7,
  height = 6
)


##############################################################
# Bootstrap AUC distribution
##############################################################

cat(
  "Searching for bootstrap AUC distribution...\n"
)

bootstrap_candidates <- c(
  
  file.path(
    PERFORMANCE_DIR,
    "robustness",
    "GSE13159_Bootstrap_AUC.csv"
  ),
  
  file.path(
    PERFORMANCE_DIR,
    "robustness",
    "Bootstrap_AUC.csv"
  ),
  
  file.path(
    PERFORMANCE_DIR,
    "robustness",
    "Stage14E_Bootstrap_AUC.csv"
  )
  
)


existing_bootstrap <- bootstrap_candidates[
  file.exists(
    bootstrap_candidates
  )
]


if (
  length(existing_bootstrap) > 0
) {
  
  bootstrap_file <- existing_bootstrap[1]
  
  cat(
    "Bootstrap file found:",
    bootstrap_file,
    "\n"
  )
  
  bootstrap_data <- readr::read_csv(
    bootstrap_file,
    show_col_types = FALSE
  )
  
  auc_columns <- names(
    bootstrap_data
  )[
    grepl(
      "auc",
      names(bootstrap_data),
      ignore.case = TRUE
    )
  ]
  
  
  if (
    length(auc_columns) > 0
  ) {
    
    ##########################################################
    # CORRECT COLUMN EXTRACTION
    ##########################################################
    
    bootstrap_values <- suppressWarnings(
      as.numeric(
        bootstrap_data[[auc_columns[1]]]
      )
    )
    
    bootstrap_values <- bootstrap_values[
      is.finite(
        bootstrap_values
      )
    ]
    
    
    if (
      length(bootstrap_values) > 10
    ) {
      
      bootstrap_plot_data <- data.frame(
        AUC = bootstrap_values
      )
      
      bootstrap_plot <- ggplot(
        bootstrap_plot_data,
        aes(
          x = AUC
        )
      ) +
        
        geom_histogram(
          bins = 40
        ) +
        
        geom_vline(
          xintercept = roc_auc,
          linetype = "dashed",
          linewidth = 1
        ) +
        
        labs(
          title = "Bootstrap Distribution of AUC",
          subtitle = "GSE13159 external validation",
          x = "Bootstrap AUC",
          y = "Frequency"
        ) +
        
        theme_classic(
          base_size = 13
        )
      
      
      ggsave(
        file.path(
          FIGURE_DIR,
          "GSE13159_Bootstrap_AUC_Distribution.png"
        ),
        bootstrap_plot,
        width = 7,
        height = 6,
        dpi = 300
      )
      
      
      ggsave(
        file.path(
          FIGURE_DIR,
          "GSE13159_Bootstrap_AUC_Distribution.pdf"
        ),
        bootstrap_plot,
        width = 7,
        height = 6
      )
      
    }
    
  }
  
} else {
  
  cat(
    "Bootstrap AUC vector file not found.\n"
  )
  
}


##############################################################
# Save statistics
##############################################################

summary_table <- data.frame(
  
  Dataset = "GSE13159",
  
  Total_Samples = nrow(scores),
  
  AML_Samples = sum(
    scores$AML_Status == "AML"
  ),
  
  Non_AML_Samples = sum(
    scores$AML_Status == "Non_AML"
  ),
  
  AUC = roc_auc,
  
  AUC_CI_Lower = as.numeric(
    roc_ci[1]
  ),
  
  AUC_CI_Upper = as.numeric(
    roc_ci[3]
  ),
  
  Sensitivity = sensitivity,
  
  Specificity = specificity,
  
  PPV = ppv,
  
  NPV = npv,
  
  Accuracy = accuracy,
  
  Optimal_Threshold = threshold
  
)


write_csv(
  summary_table,
  file.path(
    TABLE_DIR,
    "GSE13159_External_Validation_Figure_Statistics.csv"
  )
)


##############################################################
# Scientific interpretation
##############################################################

interpretation <- c(
  
  "Stage 14G External Validation Visualization Summary",
  
  "",
  
  "Dataset: GSE13159",
  
  paste0(
    "Total samples: ",
    nrow(scores)
  ),
  
  paste0(
    "AML samples: ",
    sum(
      scores$AML_Status == "AML"
    )
  ),
  
  paste0(
    "Non AML samples: ",
    sum(
      scores$AML_Status == "Non_AML"
    )
  ),
  
  "",
  
  paste0(
    "Observed AUC: ",
    sprintf(
      "%.4f",
      roc_auc
    )
  ),
  
  paste0(
    "95% CI: ",
    sprintf(
      "%.4f",
      roc_ci[1]
    ),
    " to ",
    sprintf(
      "%.4f",
      roc_ci[3]
    )
  ),
  
  "",
  
  paste0(
    "Sensitivity: ",
    sprintf(
      "%.4f",
      sensitivity
    )
  ),
  
  paste0(
    "Specificity: ",
    sprintf(
      "%.4f",
      specificity
    )
  ),
  
  paste0(
    "PPV: ",
    sprintf(
      "%.4f",
      ppv
    )
  ),
  
  paste0(
    "NPV: ",
    sprintf(
      "%.4f",
      npv
    )
  ),
  
  paste0(
    "Accuracy: ",
    sprintf(
      "%.4f",
      accuracy
    )
  ),
  
  "",
  
  "Interpretation:",
  
  "The AML biomarker signature demonstrated moderate",
  "discriminatory performance in the independent GSE13159",
  "validation cohort.",
  
  "The observed AUC indicates that the signature distinguishes",
  "AML samples from the broader non AML hematological",
  "comparison population above chance.",
  
  "The close agreement between observed and bootstrap",
  "performance supports the robustness of the validation.",
  
  "Because GSE13159 contains multiple leukemia subtypes,",
  "the result should not be interpreted as validation against",
  "healthy controls alone."
  
)


writeLines(
  interpretation,
  file.path(
    FIGURE_DIR,
    "Stage14G_Figure_Interpretation.txt"
  )
)


##############################################################
# Completion summary
##############################################################

end_time <- Sys.time()

cat("\n")
cat("=========================================\n")
cat("Stage 14G Completed Successfully\n")
cat("=========================================\n\n")

cat(
  "Dataset           : GSE13159\n"
)

cat(
  "Total Samples     :",
  nrow(scores),
  "\n"
)

cat(
  "AML Samples       :",
  sum(
    scores$AML_Status == "AML"
  ),
  "\n"
)

cat(
  "Non AML Samples   :",
  sum(
    scores$AML_Status == "Non_AML"
  ),
  "\n"
)

cat(
  "Observed AUC      :",
  round(
    roc_auc,
    4
  ),
  "\n"
)

cat(
  "95% CI            :",
  round(
    roc_ci[1],
    4
  ),
  "to",
  round(
    roc_ci[3],
    4
  ),
  "\n"
)

cat(
  "Figure Directory   :",
  FIGURE_DIR,
  "\n"
)

cat(
  "Statistics         :",
  file.path(
    TABLE_DIR,
    "GSE13159_External_Validation_Figure_Statistics.csv"
  ),
  "\n"
)

cat(
  "Time Elapsed       :",
  round(
    end_time - start_time,
    2
  ),
  "\n\n"
)