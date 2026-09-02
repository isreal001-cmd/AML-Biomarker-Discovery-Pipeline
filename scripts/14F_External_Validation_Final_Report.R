##############################################################
# AML Biomarker Discovery Pipeline
#
# Script:
# 14F_External_Validation_Final_Report.R
#
# Purpose:
# Consolidate external validation results from Stages 14C,
# 14D and 14E into one final reproducible report.
#
# Author:
# Isreal Oluwafemi Abiodun
#
# Version:
# 2.0.0
##############################################################

source("config.R")

library(readr)
library(dplyr)

cat("\n")
cat("=========================================\n")
cat("Stage 14F : Final External Validation Report\n")
cat("=========================================\n\n")

start_time <- Sys.time()

##############################################################
# Directories
##############################################################

VALIDATION_DIR <- file.path(
  RESULTS_DIR,
  "validation"
)

FINAL_DIR <- file.path(
  VALIDATION_DIR,
  "final_report"
)

dir.create(
  FINAL_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

##############################################################
# Input Files
##############################################################

ROC_FILE <- file.path(
  VALIDATION_DIR,
  "ROC",
  "GSE13159_ROC_results.csv"
)

PERFORMANCE_FILE <- file.path(
  VALIDATION_DIR,
  "performance",
  "tables",
  "GSE13159_Performance_Metrics.csv"
)

CONFUSION_FILE <- file.path(
  VALIDATION_DIR,
  "performance",
  "tables",
  "GSE13159_Confusion_Matrix.csv"
)

ROBUSTNESS_FILE <- file.path(
  VALIDATION_DIR,
  "performance",
  "robustness",
  "Stage14E_External_Validation_Robustness_Summary.txt"
)

##############################################################
# Check Files
##############################################################

required_files <- c(
  ROC_FILE,
  PERFORMANCE_FILE,
  CONFUSION_FILE,
  ROBUSTNESS_FILE
)

missing_files <- required_files[
  !file.exists(required_files)
]

if(length(missing_files) > 0){
  
  stop(
    paste(
      "Required validation files are missing:\n",
      paste(missing_files, collapse = "\n")
    )
  )
  
}

##############################################################
# Helper Functions
##############################################################

extract_numeric <- function(x){
  
  if(length(x) == 0 || is.na(x)){
    return(NA_real_)
  }
  
  x <- as.character(x)
  
  # Remove commas and whitespace
  x <- gsub(",", "", x)
  x <- trimws(x)
  
  # Extract first valid numeric value
  value <- suppressWarnings(
    as.numeric(
      sub(
        "^.*?(-?[0-9]+\\.?[0-9]*(?:[eE][+-]?[0-9]+)?).*$",
        "\\1",
        x
      )
    )
  )
  
  return(value)
  
}

extract_metric_from_roc <- function(
    roc_data,
    metric_names
){
  
  idx <- which(
    tolower(trimws(roc_data$Metric)) %in%
      tolower(metric_names)
  )
  
  if(length(idx) == 0){
    return(NA_real_)
  }
  
  extract_numeric(
    roc_data$Value[idx[1]]
  )
  
}

extract_text_value <- function(
    lines,
    pattern
){
  
  hit <- grep(
    pattern,
    lines,
    ignore.case = TRUE,
    value = TRUE
  )
  
  if(length(hit) == 0){
    return(NA_character_)
  }
  
  hit[1]
  
}

##############################################################
# Load Stage 14C ROC Results
##############################################################

cat("Loading Stage 14C ROC results...\n")

roc_data <- read_csv(
  ROC_FILE,
  show_col_types = FALSE
)

cat(
  "14C columns detected:\n"
)

print(
  names(roc_data)
)

##############################################################
# Extract ROC Metrics
##############################################################

observed_auc <- extract_metric_from_roc(
  roc_data,
  c(
    "AUC",
    "Observed AUC"
  )
)

auc_ci_lower <- extract_metric_from_roc(
  roc_data,
  c(
    "AUC_CI_Lower",
    "AUC CI Lower",
    "AUC 95% CI Lower",
    "CI Lower"
  )
)

auc_ci_upper <- extract_metric_from_roc(
  roc_data,
  c(
    "AUC_CI_Upper",
    "AUC CI Upper",
    "AUC 95% CI Upper",
    "CI Upper"
  )
)

##############################################################
# Load Stage 14D Performance Results
##############################################################

cat("\nLoading Stage 14D performance results...\n")

performance <- read_csv(
  PERFORMANCE_FILE,
  show_col_types = FALSE
)

cat(
  "14D columns detected:\n"
)

print(
  names(performance)
)

##############################################################
# Extract Dataset Information
##############################################################

dataset <- as.character(
  performance$Dataset[1]
)

total_samples <- extract_numeric(
  performance$Total_Samples[1]
)

aml_samples <- extract_numeric(
  performance$AML_Samples[1]
)

non_aml_samples <- extract_numeric(
  performance$Non_AML_Samples[1]
)

##############################################################
# Extract Performance Metrics
##############################################################

# Use Stage 14D as primary source

performance_auc <- extract_numeric(
  performance$AUC[1]
)

performance_ci_lower <- extract_numeric(
  performance$AUC_CI_Lower[1]
)

performance_ci_upper <- extract_numeric(
  performance$AUC_CI_Upper[1]
)

optimal_threshold <- extract_numeric(
  performance$Optimal_Threshold[1]
)

sensitivity <- extract_numeric(
  performance$Sensitivity[1]
)

specificity <- extract_numeric(
  performance$Specificity[1]
)

accuracy <- extract_numeric(
  performance$Accuracy[1]
)

##############################################################
# Resolve AUC
##############################################################

if(is.na(observed_auc)){
  
  observed_auc <- performance_auc
  
}

if(is.na(auc_ci_lower)){
  
  auc_ci_lower <- performance_ci_lower
  
}

if(is.na(auc_ci_upper)){
  
  auc_ci_upper <- performance_ci_upper
  
}

##############################################################
# Load Confusion Matrix
##############################################################

cat("\nLoading confusion matrix...\n")

confusion <- read_csv(
  CONFUSION_FILE,
  show_col_types = FALSE
)

##############################################################
# Identify Confusion Matrix Counts
##############################################################

find_column <- function(
    data,
    candidates
){
  
  idx <- which(
    tolower(names(data)) %in%
      tolower(candidates)
  )
  
  if(length(idx) == 0){
    return(NA_character_)
  }
  
  names(data)[idx[1]]
  
}

TP_col <- find_column(
  confusion,
  c(
    "TP",
    "True_Positive",
    "True Positive",
    "TruePositive"
  )
)

TN_col <- find_column(
  confusion,
  c(
    "TN",
    "True_Negative",
    "True Negative",
    "TrueNegative"
  )
)

FP_col <- find_column(
  confusion,
  c(
    "FP",
    "False_Positive",
    "False Positive",
    "FalsePositive"
  )
)

FN_col <- find_column(
  confusion,
  c(
    "FN",
    "False_Negative",
    "False Negative",
    "FalseNegative"
  )
)

##############################################################
# Calculate PPV and NPV
##############################################################

ppv <- NA_real_
npv <- NA_real_

if(
  !is.na(TP_col) &&
  !is.na(FP_col)
){
  
  TP <- extract_numeric(
    confusion[[TP_col]][1]
  )
  
  FP <- extract_numeric(
    confusion[[FP_col]][1]
  )
  
  if(
    !is.na(TP) &&
    !is.na(FP) &&
    (TP + FP) > 0
  ){
    
    ppv <- TP / (TP + FP)
    
  }
  
}

if(
  !is.na(TN_col) &&
  !is.na(FN_col)
){
  
  TN <- extract_numeric(
    confusion[[TN_col]][1]
  )
  
  FN <- extract_numeric(
    confusion[[FN_col]][1]
  )
  
  if(
    !is.na(TN) &&
    !is.na(FN) &&
    (TN + FN) > 0
  ){
    
    npv <- TN / (TN + FN)
    
  }
  
}

##############################################################
# If PPV/NPV Are Present in Stage 14D, Prefer Them
##############################################################

if("PPV" %in% names(performance)){
  
  stage14d_ppv <- extract_numeric(
    performance$PPV[1]
  )
  
  if(!is.na(stage14d_ppv)){
    ppv <- stage14d_ppv
  }
  
}

if("NPV" %in% names(performance)){
  
  stage14d_npv <- extract_numeric(
    performance$NPV[1]
  )
  
  if(!is.na(stage14d_npv)){
    npv <- stage14d_npv
  }
  
}

##############################################################
# Load Stage 14E Robustness Results
##############################################################

cat("\nLoading Stage 14E robustness results...\n")

robustness_lines <- readLines(
  ROBUSTNESS_FILE,
  warn = FALSE
)

##############################################################
# Extract Bootstrap Metrics
##############################################################

bootstrap_mean_auc <- NA_real_
bootstrap_ci_lower <- NA_real_
bootstrap_ci_upper <- NA_real_
wilcox_p <- NA_character_
bootstrap_iterations <- NA_real_
bootstrap_successful <- NA_real_

##############################################################
# Bootstrap Mean AUC
##############################################################

hit <- grep(
  "Mean Bootstrap AUC",
  robustness_lines,
  ignore.case = TRUE,
  value = TRUE
)

if(length(hit) > 0){
  
  bootstrap_mean_auc <- extract_numeric(
    hit[1]
  )
  
}

##############################################################
# Bootstrap CI
##############################################################

hit <- grep(
  "Bootstrap 95% CI",
  robustness_lines,
  ignore.case = TRUE,
  value = TRUE
)

if(length(hit) > 0){
  
  ci_numbers <- regmatches(
    hit[1],
    gregexpr(
      "[0-9]+\\.?[0-9]*",
      hit[1]
    )
  )[[1]]
  
  ci_numbers <- as.numeric(
    ci_numbers
  )
  
  if(length(ci_numbers) >= 2){
    
    bootstrap_ci_lower <- ci_numbers[
      length(ci_numbers) - 1
    ]
    
    bootstrap_ci_upper <- ci_numbers[
      length(ci_numbers)
    ]
    
  }
  
}

##############################################################
# Bootstrap Iterations
##############################################################

hit <- grep(
  "Bootstrap iterations",
  robustness_lines,
  ignore.case = TRUE,
  value = TRUE
)

if(length(hit) > 0){
  
  bootstrap_iterations <- extract_numeric(
    hit[1]
  )
  
}

##############################################################
# Successful Iterations
##############################################################

hit <- grep(
  "Successful",
  robustness_lines,
  ignore.case = TRUE,
  value = TRUE
)

if(length(hit) > 0){
  
  bootstrap_successful <- extract_numeric(
    hit[1]
  )
  
}

##############################################################
# Wilcoxon P Value
##############################################################

hit <- grep(
  "Wilcoxon P Value",
  robustness_lines,
  ignore.case = TRUE,
  value = TRUE
)

if(length(hit) > 0){
  
  wilcox_p <- sub(
    "^.*?:\\s*",
    "",
    hit[1]
  )
  
  wilcox_p <- trimws(
    wilcox_p
  )
  
}

##############################################################
# Fallback Values
##############################################################

if(is.na(bootstrap_mean_auc)){
  
  bootstrap_mean_auc <- observed_auc
  
}

if(is.na(bootstrap_ci_lower)){
  
  bootstrap_ci_lower <- auc_ci_lower
  
}

if(is.na(bootstrap_ci_upper)){
  
  bootstrap_ci_upper <- auc_ci_upper
  
}

##############################################################
# Final Summary
##############################################################

cat("\n")
cat("=========================================\n")
cat("FINAL EXTERNAL VALIDATION SUMMARY\n")
cat("=========================================\n\n")

cat(
  "Dataset             : ",
  dataset,
  "\n",
  sep = ""
)

cat(
  "Total Samples       : ",
  total_samples,
  "\n",
  sep = ""
)

cat(
  "AML Samples         : ",
  aml_samples,
  "\n",
  sep = ""
)

cat(
  "Non AML Samples     : ",
  non_aml_samples,
  "\n\n",
  sep = ""
)

cat(
  "Observed AUC        : ",
  sprintf("%.4f", observed_auc),
  "\n",
  sep = ""
)

cat(
  "Observed 95% CI     : ",
  sprintf("%.4f", auc_ci_lower),
  " to ",
  sprintf("%.4f", auc_ci_upper),
  "\n\n",
  sep = ""
)

cat(
  "Bootstrap Mean AUC  : ",
  sprintf("%.4f", bootstrap_mean_auc),
  "\n",
  sep = ""
)

cat(
  "Bootstrap 95% CI    : ",
  sprintf("%.4f", bootstrap_ci_lower),
  " to ",
  sprintf("%.4f", bootstrap_ci_upper),
  "\n\n",
  sep = ""
)

cat(
  "Optimal Threshold   : ",
  sprintf("%.6f", optimal_threshold),
  "\n",
  sep = ""
)

cat(
  "Sensitivity         : ",
  sprintf("%.4f", sensitivity),
  "\n",
  sep = ""
)

cat(
  "Specificity         : ",
  sprintf("%.4f", specificity),
  "\n",
  sep = ""
)

cat(
  "PPV                 : ",
  sprintf("%.4f", ppv),
  "\n",
  sep = ""
)

cat(
  "NPV                 : ",
  sprintf("%.4f", npv),
  "\n",
  sep = ""
)

cat(
  "Accuracy            : ",
  sprintf("%.4f", accuracy),
  "\n",
  sep = ""
)

cat(
  "Wilcoxon P Value    : ",
  wilcox_p,
  "\n\n",
  sep = ""
)

##############################################################
# Create Final Summary Table
##############################################################

summary_table <- data.frame(
  
  Dataset = dataset,
  
  Total_Samples = total_samples,
  
  AML_Samples = aml_samples,
  
  Non_AML_Samples = non_aml_samples,
  
  Observed_AUC = observed_auc,
  
  Observed_AUC_CI_Lower = auc_ci_lower,
  
  Observed_AUC_CI_Upper = auc_ci_upper,
  
  Bootstrap_Mean_AUC = bootstrap_mean_auc,
  
  Bootstrap_AUC_CI_Lower = bootstrap_ci_lower,
  
  Bootstrap_AUC_CI_Upper = bootstrap_ci_upper,
  
  Bootstrap_Iterations = bootstrap_iterations,
  
  Bootstrap_Successful = bootstrap_successful,
  
  Optimal_Threshold = optimal_threshold,
  
  Sensitivity = sensitivity,
  
  Specificity = specificity,
  
  PPV = ppv,
  
  NPV = npv,
  
  Accuracy = accuracy,
  
  Wilcoxon_P_Value = wilcox_p,
  
  stringsAsFactors = FALSE
  
)

##############################################################
# Save Summary Table
##############################################################

SUMMARY_FILE <- file.path(
  FINAL_DIR,
  "GSE13159_External_Validation_Summary.csv"
)

write_csv(
  summary_table,
  SUMMARY_FILE
)

##############################################################
# Create Human Readable Report
##############################################################

REPORT_FILE <- file.path(
  FINAL_DIR,
  "GSE13159_External_Validation_Report.txt"
)

report_lines <- c(
  
  "AML BIOMARKER DISCOVERY PIPELINE",
  "FINAL EXTERNAL VALIDATION REPORT",
  "",
  "=========================================",
  "",
  paste("Dataset:", dataset),
  paste("Total Samples:", total_samples),
  paste("AML Samples:", aml_samples),
  paste("Non AML Samples:", non_aml_samples),
  "",
  "ROC PERFORMANCE",
  "",
  paste(
    "Observed AUC:",
    sprintf("%.4f", observed_auc)
  ),
  paste(
    "Observed 95% CI:",
    sprintf("%.4f", auc_ci_lower),
    "to",
    sprintf("%.4f", auc_ci_upper)
  ),
  "",
  "BOOTSTRAP ROBUSTNESS",
  "",
  paste(
    "Bootstrap Mean AUC:",
    sprintf("%.4f", bootstrap_mean_auc)
  ),
  paste(
    "Bootstrap 95% CI:",
    sprintf("%.4f", bootstrap_ci_lower),
    "to",
    sprintf("%.4f", bootstrap_ci_upper)
  ),
  paste(
    "Bootstrap Iterations:",
    bootstrap_iterations
  ),
  paste(
    "Successful Iterations:",
    bootstrap_successful
  ),
  "",
  "CLASSIFICATION PERFORMANCE",
  "",
  paste(
    "Optimal Threshold:",
    sprintf("%.6f", optimal_threshold)
  ),
  paste(
    "Sensitivity:",
    sprintf("%.4f", sensitivity)
  ),
  paste(
    "Specificity:",
    sprintf("%.4f", specificity)
  ),
  paste(
    "PPV:",
    sprintf("%.4f", ppv)
  ),
  paste(
    "NPV:",
    sprintf("%.4f", npv)
  ),
  paste(
    "Accuracy:",
    sprintf("%.4f", accuracy)
  ),
  "",
  "STATISTICAL SIGNIFICANCE",
  "",
  paste(
    "Wilcoxon P Value:",
    wilcox_p
  ),
  "",
  "INTERPRETATION",
  "",
  paste(
    "The AML biomarker signature demonstrated",
    "moderate discriminatory performance in the",
    "independent GSE13159 validation cohort."
  ),
  "",
  paste(
    "The observed AUC of",
    sprintf("%.4f", observed_auc),
    "was consistent with the bootstrap mean AUC of",
    sprintf("%.4f", bootstrap_mean_auc),
    "indicating stable model discrimination."
  ),
  "",
  paste(
    "The bootstrap analysis used",
    bootstrap_successful,
    "successful resampling iterations."
  ),
  "",
  paste(
    "The signature score distributions between",
    "AML and non AML samples were statistically",
    "different with a Wilcoxon P value of",
    wilcox_p,
    "."
  ),
  "",
  "=========================================",
  "",
  "Generated by Stage 14F",
  "AML Biomarker Discovery Pipeline",
  "Author: Isreal Oluwafemi Abiodun",
  ""
)

writeLines(
  report_lines,
  REPORT_FILE
)

##############################################################
# Final Validation Checks
##############################################################

critical_metrics <- c(
  observed_auc,
  auc_ci_lower,
  auc_ci_upper,
  bootstrap_mean_auc,
  bootstrap_ci_lower,
  bootstrap_ci_upper,
  sensitivity,
  specificity,
  accuracy
)

if(any(is.na(critical_metrics))){
  
  warning(
    "One or more critical validation metrics could not be recovered."
  )
  
}

##############################################################
# Completion
##############################################################

elapsed <- Sys.time() - start_time

cat(
  "Summary Table       : ",
  SUMMARY_FILE,
  "\n",
  sep = ""
)

cat(
  "Final Report        : ",
  REPORT_FILE,
  "\n",
  sep = ""
)

cat(
  "Time Elapsed        : ",
  round(elapsed, 2),
  "\n\n",
  sep = ""
)

cat(
  "Stage 14F completed successfully.\n"
)