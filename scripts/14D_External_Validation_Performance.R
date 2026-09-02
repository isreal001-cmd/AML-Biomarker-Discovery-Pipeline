##############################################################
# AML Biomarker Discovery Pipeline
#
# Script:
# 14D_External_Validation_Performance.R
#
# Purpose:
# Evaluate external validation performance of the AML
# biomarker signature using the GSE13159 dataset.
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
cat("Stage 14D : External Validation Performance\n")
cat("=========================================\n\n")

start_time <- Sys.time()

##############################################################
# Required Packages
##############################################################

if (!requireNamespace("pROC", quietly = TRUE)) {
  
  stop(
    "Package 'pROC' is required.\n",
    "Run 00_check_packages.R first."
  )
  
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  
  stop(
    "Package 'ggplot2' is required.\n",
    "Run 00_check_packages.R first."
  )
  
}

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

PLOT_DIR <- file.path(
  PERFORMANCE_DIR,
  "plots"
)

TABLE_DIR <- file.path(
  PERFORMANCE_DIR,
  "tables"
)

dir.create(
  PERFORMANCE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  PLOT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  TABLE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

##############################################################
# Locate Stage 14C Score File
##############################################################

SCORE_FILE <- file.path(
  VALIDATION_DIR,
  "GSE13159_biomarker_signature_scores.csv"
)

if (!file.exists(SCORE_FILE)) {
  
  stop(
    "Stage 14C signature score file was not found:\n",
    SCORE_FILE
  )
  
}

##############################################################
# Load Stage 14C Scores
##############################################################

cat("Loading Stage 14C signature scores...\n")

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

cat("Columns detected:\n")

print(names(scores))

##############################################################
# Validate Required Columns
##############################################################

required_columns <- c(
  "Sample",
  "Signature_Score",
  "AML_Status"
)

missing_columns <- setdiff(
  required_columns,
  names(scores)
)

if (length(missing_columns) > 0) {
  
  stop(
    "Required columns missing:\n",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
  
}

##############################################################
# Clean Data
##############################################################

scores <- scores |>
  
  dplyr::filter(
    !is.na(Signature_Score),
    !is.na(AML_Status)
  )

cat(
  "\nSamples after filtering :",
  nrow(scores),
  "\n\n"
)

##############################################################
# Inspect Original AML Status
##############################################################

cat("Original AML status distribution:\n\n")

print(
  table(
    scores$AML_Status,
    useNA = "ifany"
  )
)

##############################################################
# Standardize AML Status
##############################################################

scores$AML_Status_Clean <- dplyr::case_when(
  
  scores$AML_Status %in% c(
    "AML",
    "aml",
    "1",
    1,
    TRUE
  ) ~ "AML",
  
  scores$AML_Status %in% c(
    "Non_AML",
    "Non AML",
    "Non-AML",
    "non_aml",
    "non aml",
    "Healthy",
    "Control",
    "0",
    0,
    FALSE
  ) ~ "Non_AML",
  
  TRUE ~ NA_character_
  
)

##############################################################
# Check Standardized Status
##############################################################

cat("\nStandardized AML status distribution:\n\n")

print(
  table(
    scores$AML_Status_Clean,
    useNA = "ifany"
  )
)

##############################################################
# Stop if Classes Were Not Properly Identified
##############################################################

if (
  sum(scores$AML_Status_Clean == "AML", na.rm = TRUE) < 10
) {
  
  stop(
    "Fewer than 10 AML samples were identified."
  )
  
}

if (
  sum(scores$AML_Status_Clean == "Non_AML", na.rm = TRUE) < 10
) {
  
  stop(
    "Fewer than 10 non AML samples were identified."
  )
  
}

##############################################################
# Keep Valid Classes
##############################################################

scores <- scores |>
  
  dplyr::filter(
    AML_Status_Clean %in% c(
      "AML",
      "Non_AML"
    )
  )

##############################################################
# Create Binary Outcome
##############################################################

scores$AML_Binary <- ifelse(
  
  scores$AML_Status_Clean == "AML",
  
  1,
  
  0
  
)

##############################################################
# Final Sample Counts
##############################################################

AML_N <- sum(
  scores$AML_Binary == 1
)

NON_AML_N <- sum(
  scores$AML_Binary == 0
)

cat("\n")
cat("=========================================\n")
cat("Validation Cohort\n")
cat("=========================================\n\n")

cat(
  "Total Samples :",
  nrow(scores),
  "\n"
)

cat(
  "AML Samples   :",
  AML_N,
  "\n"
)

cat(
  "Non AML       :",
  NON_AML_N,
  "\n\n"
)

##############################################################
# ROC Analysis
##############################################################

cat("Running ROC analysis...\n")

roc_obj <- pROC::roc(
  
  response = scores$AML_Binary,
  
  predictor = scores$Signature_Score,
  
  direction = "auto",
  
  quiet = TRUE
  
)

##############################################################
# AUC
##############################################################

auc_value <- as.numeric(
  pROC::auc(roc_obj)
)

auc_ci <- pROC::ci.auc(
  roc_obj,
  conf.level = 0.95
)

##############################################################
# Optimal Threshold
##############################################################

optimal_coords <- pROC::coords(
  
  roc_obj,
  
  x = "best",
  
  best.method = "youden",
  
  ret = c(
    "threshold",
    "sensitivity",
    "specificity",
    "ppv",
    "npv",
    "accuracy"
  ),
  
  transpose = FALSE
  
)

##############################################################
# Extract Performance
##############################################################

optimal_threshold <- as.numeric(
  optimal_coords["threshold"]
)

sensitivity <- as.numeric(
  optimal_coords["sensitivity"]
)

specificity <- as.numeric(
  optimal_coords["specificity"]
)

ppv <- as.numeric(
  optimal_coords["ppv"]
)

npv <- as.numeric(
  optimal_coords["npv"]
)

accuracy <- as.numeric(
  optimal_coords["accuracy"]
)

##############################################################
# Classification Using Optimal Threshold
##############################################################

scores$Predicted_Status <- ifelse(
  
  scores$Signature_Score >= optimal_threshold,
  
  "AML",
  
  "Non_AML"
  
)

##############################################################
# Confusion Matrix
##############################################################

confusion_matrix <- table(
  
  Actual = scores$AML_Status_Clean,
  
  Predicted = scores$Predicted_Status
  
)

##############################################################
# Calculate Counts Safely
##############################################################

TP <- ifelse(
  "AML" %in% rownames(confusion_matrix) &&
    "AML" %in% colnames(confusion_matrix),
  confusion_matrix["AML", "AML"],
  0
)

TN <- ifelse(
  "Non_AML" %in% rownames(confusion_matrix) &&
    "Non_AML" %in% colnames(confusion_matrix),
  confusion_matrix["Non_AML", "Non_AML"],
  0
)

FP <- ifelse(
  "Non_AML" %in% rownames(confusion_matrix) &&
    "AML" %in% colnames(confusion_matrix),
  confusion_matrix["Non_AML", "AML"],
  0
)

FN <- ifelse(
  "AML" %in% rownames(confusion_matrix) &&
    "Non_AML" %in% colnames(confusion_matrix),
  confusion_matrix["AML", "Non_AML"],
  0
)

##############################################################
# Recalculate Performance
##############################################################

accuracy_calc <- (
  TP + TN
) / (
  TP + TN + FP + FN
)

sensitivity_calc <- (
  TP
) / (
  TP + FN
)

specificity_calc <- (
  TN
) / (
  TN + FP
)

##############################################################
# Save ROC Results
##############################################################

roc_results <- data.frame(
  
  Dataset = "GSE13159",
  
  Total_Samples = nrow(scores),
  
  AML_Samples = AML_N,
  
  Non_AML_Samples = NON_AML_N,
  
  AUC = auc_value,
  
  AUC_CI_Lower = as.numeric(auc_ci[1]),
  
  AUC_CI_Upper = as.numeric(auc_ci[3]),
  
  Optimal_Threshold = optimal_threshold,
  
  Sensitivity = sensitivity_calc,
  
  Specificity = specificity_calc,
  
  PPV = ppv,
  
  NPV = npv,
  
  Accuracy = accuracy_calc,
  
  stringsAsFactors = FALSE
  
)

ROC_RESULTS_FILE <- file.path(
  TABLE_DIR,
  "GSE13159_Performance_Metrics.csv"
)

readr::write_csv(
  roc_results,
  ROC_RESULTS_FILE
)

##############################################################
# Save Confusion Matrix
##############################################################

CONFUSION_FILE <- file.path(
  TABLE_DIR,
  "GSE13159_Confusion_Matrix.csv"
)

readr::write_csv(
  
  as.data.frame(confusion_matrix),
  
  CONFUSION_FILE
  
)

##############################################################
# Save Sample Level Predictions
##############################################################

PREDICTIONS_FILE <- file.path(
  TABLE_DIR,
  "GSE13159_Sample_Predictions.csv"
)

readr::write_csv(
  
  scores,
  
  PREDICTIONS_FILE
  
)

##############################################################
# ROC Plot
##############################################################

ROC_PLOT <- file.path(
  PLOT_DIR,
  "GSE13159_External_Validation_ROC.png"
)

png(
  ROC_PLOT,
  width = 2400,
  height = 2000,
  res = 300
)

plot(
  roc_obj,
  col = "black",
  lwd = 3,
  main = "External Validation of AML Biomarker Signature",
  xlab = "1 - Specificity",
  ylab = "Sensitivity"
)

abline(
  a = 0,
  b = 1,
  lty = 2
)

legend(
  "bottomright",
  legend = paste0(
    "AUC = ",
    round(auc_value, 4),
    "\n95% CI = ",
    round(auc_ci[1], 4),
    " to ",
    round(auc_ci[3], 4)
  ),
  bty = "n"
)

dev.off()

##############################################################
# Performance Summary
##############################################################

SUMMARY_FILE <- file.path(
  PERFORMANCE_DIR,
  "Stage14D_External_Validation_Performance_Summary.txt"
)

summary_lines <- c(
  
  "AML EXTERNAL VALIDATION PERFORMANCE",
  
  "=========================================",
  
  "",
  
  "Dataset : GSE13159",
  
  paste(
    "Total Samples :",
    nrow(scores)
  ),
  
  paste(
    "AML Samples :",
    AML_N
  ),
  
  paste(
    "Non AML Samples :",
    NON_AML_N
  ),
  
  "",
  
  paste(
    "AUC :",
    round(auc_value, 4)
  ),
  
  paste(
    "AUC 95% CI :",
    round(auc_ci[1], 4),
    "to",
    round(auc_ci[3], 4)
  ),
  
  paste(
    "Optimal Threshold :",
    round(optimal_threshold, 6)
  ),
  
  paste(
    "Sensitivity :",
    round(sensitivity_calc, 4)
  ),
  
  paste(
    "Specificity :",
    round(specificity_calc, 4)
  ),
  
  paste(
    "PPV :",
    round(ppv, 4)
  ),
  
  paste(
    "NPV :",
    round(npv, 4)
  ),
  
  paste(
    "Accuracy :",
    round(accuracy_calc, 4)
  ),
  
  "",
  
  "Confusion Matrix",
  
  "",
  
  paste(
    "True Positives :",
    TP
  ),
  
  paste(
    "True Negatives :",
    TN
  ),
  
  paste(
    "False Positives :",
    FP
  ),
  
  paste(
    "False Negatives :",
    FN
  ),
  
  "",
  
  paste(
    "Performance Metrics :",
    ROC_RESULTS_FILE
  ),
  
  paste(
    "Confusion Matrix :",
    CONFUSION_FILE
  ),
  
  paste(
    "Predictions :",
    PREDICTIONS_FILE
  ),
  
  paste(
    "ROC Plot :",
    ROC_PLOT
  )
  
)

writeLines(
  summary_lines,
  SUMMARY_FILE
)

##############################################################
# Console Summary
##############################################################

cat("\n")
cat("=========================================\n")
cat("Stage 14D Completed Successfully\n")
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
  AML_N,
  "\n"
)

cat(
  "Non AML Samples   :",
  NON_AML_N,
  "\n\n"
)

cat(
  "AUC               :",
  round(auc_value, 4),
  "\n"
)

cat(
  "AUC 95% CI        :",
  round(auc_ci[1], 4),
  "to",
  round(auc_ci[3], 4),
  "\n"
)

cat(
  "Optimal Threshold :",
  round(optimal_threshold, 6),
  "\n"
)

cat(
  "Sensitivity        :",
  round(sensitivity_calc, 4),
  "\n"
)

cat(
  "Specificity        :",
  round(specificity_calc, 4),
  "\n"
)

cat(
  "PPV                :",
  round(ppv, 4),
  "\n"
)

cat(
  "NPV                :",
  round(npv, 4),
  "\n"
)

cat(
  "Accuracy           :",
  round(accuracy_calc, 4),
  "\n\n"
)

cat(
  "Performance Table  :",
  ROC_RESULTS_FILE,
  "\n"
)

cat(
  "Confusion Matrix   :",
  CONFUSION_FILE,
  "\n"
)

cat(
  "Predictions        :",
  PREDICTIONS_FILE,
  "\n"
)

cat(
  "ROC Plot           :",
  ROC_PLOT,
  "\n"
)

cat(
  "Summary            :",
  SUMMARY_FILE,
  "\n\n"
)

cat(
  "Time Elapsed       :",
  round(
    Sys.time() - start_time,
    2
  ),
  "\n\n"
)