##############################################################
# AML Biomarker Discovery Pipeline
#
# Script:
# 14E_External_Validation_Robustness.R
#
# Purpose:
# Assess the robustness and stability of the AML biomarker
# signature using bootstrap validation in the independent
# GSE13159 cohort.
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
cat("Stage 14E : External Validation Robustness\n")
cat("=========================================\n\n")

start_time <- Sys.time()


##############################################################
# Required Packages
##############################################################

required_packages <- c(
  "readr",
  "dplyr",
  "ggplot2",
  "pROC"
)

for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    stop(
      "Required package not installed: ",
      pkg,
      "\nRun 00_check_packages.R first."
    )
    
  }
  
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

ROBUSTNESS_DIR <- file.path(
  PERFORMANCE_DIR,
  "robustness"
)

PLOT_DIR <- file.path(
  ROBUSTNESS_DIR,
  "plots"
)

TABLE_DIR <- file.path(
  ROBUSTNESS_DIR,
  "tables"
)

dir.create(
  ROBUSTNESS_DIR,
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
# Input File
##############################################################

SCORE_FILE <- file.path(
  VALIDATION_DIR,
  "GSE13159_biomarker_signature_scores.csv"
)

if (!file.exists(SCORE_FILE)) {
  
  stop(
    "Stage 14C signature score file not found:\n",
    SCORE_FILE
  )
  
}


##############################################################
# Load Signature Scores
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
# Validate Columns
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
# Remove Unrecognized Classes
##############################################################

scores <- scores |>
  
  dplyr::filter(
    AML_Status_Clean %in% c(
      "AML",
      "Non_AML"
    )
  )


##############################################################
# Binary Outcome
##############################################################

scores$AML_Binary <- ifelse(
  scores$AML_Status_Clean == "AML",
  1,
  0
)


##############################################################
# Sample Summary
##############################################################

TOTAL_N <- nrow(scores)

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
  TOTAL_N,
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
# Original ROC
##############################################################

cat("Calculating observed ROC performance...\n")

observed_roc <- pROC::roc(
  
  response = scores$AML_Binary,
  
  predictor = scores$Signature_Score,
  
  direction = "auto",
  
  quiet = TRUE
  
)

observed_auc <- as.numeric(
  pROC::auc(observed_roc)
)

observed_auc_ci <- pROC::ci.auc(
  observed_roc,
  conf.level = 0.95
)


##############################################################
# Bootstrap Settings
##############################################################

set.seed(20260809)

N_BOOTSTRAP <- 2000

cat("\n")
cat(
  "Bootstrap iterations :",
  N_BOOTSTRAP,
  "\n"
)

cat(
  "Random seed          : 20260809\n\n"
)


##############################################################
# Bootstrap Function
##############################################################

bootstrap_results <- vector(
  "list",
  N_BOOTSTRAP
)


##############################################################
# Run Bootstrap
##############################################################

cat("Running bootstrap validation...\n\n")

progress_points <- unique(
  round(
    seq(
      1,
      N_BOOTSTRAP,
      length.out = 10
    )
  )
)

for (i in seq_len(N_BOOTSTRAP)) {
  
  ############################################################
  # Stratified Bootstrap
  ############################################################
  
  aml_indices <- which(
    scores$AML_Binary == 1
  )
  
  non_aml_indices <- which(
    scores$AML_Binary == 0
  )
  
  sampled_aml <- sample(
    aml_indices,
    size = length(aml_indices),
    replace = TRUE
  )
  
  sampled_non_aml <- sample(
    non_aml_indices,
    size = length(non_aml_indices),
    replace = TRUE
  )
  
  bootstrap_indices <- c(
    sampled_aml,
    sampled_non_aml
  )
  
  boot_data <- scores[
    bootstrap_indices,
  ]
  
  ############################################################
  # ROC
  ############################################################
  
  boot_roc <- tryCatch(
    
    pROC::roc(
      
      response = boot_data$AML_Binary,
      
      predictor = boot_data$Signature_Score,
      
      direction = observed_roc$direction,
      
      quiet = TRUE
      
    ),
    
    error = function(e) {
      
      NULL
      
    }
    
  )
  
  
  ############################################################
  # Store Metrics
  ############################################################
  
  if (!is.null(boot_roc)) {
    
    boot_auc <- as.numeric(
      pROC::auc(boot_roc)
    )
    
    boot_coords <- tryCatch(
      
      pROC::coords(
        
        boot_roc,
        
        x = "best",
        
        best.method = "youden",
        
        ret = c(
          "threshold",
          "sensitivity",
          "specificity"
        ),
        
        transpose = FALSE
        
      ),
      
      error = function(e) {
        
        NULL
        
      }
      
    )
    
    if (!is.null(boot_coords)) {
      
      bootstrap_results[[i]] <- data.frame(
        
        Iteration = i,
        
        AUC = boot_auc,
        
        Threshold =
          as.numeric(
            boot_coords["threshold"]
          ),
        
        Sensitivity =
          as.numeric(
            boot_coords["sensitivity"]
          ),
        
        Specificity =
          as.numeric(
            boot_coords["specificity"]
          )
        
      )
      
    } else {
      
      bootstrap_results[[i]] <- data.frame(
        
        Iteration = i,
        
        AUC = boot_auc,
        
        Threshold = NA_real_,
        
        Sensitivity = NA_real_,
        
        Specificity = NA_real_
        
      )
      
    }
    
  }
  
  ############################################################
  # Progress
  ############################################################
  
  if (i %in% progress_points) {
    
    cat(
      "Bootstrap progress :",
      i,
      "/",
      N_BOOTSTRAP,
      "\n"
    )
    
  }
  
}


##############################################################
# Combine Bootstrap Results
##############################################################

bootstrap_results <- dplyr::bind_rows(
  bootstrap_results
)


##############################################################
# Remove Failed Iterations
##############################################################

bootstrap_results <- bootstrap_results |>
  
  dplyr::filter(
    !is.na(AUC)
  )


##############################################################
# Bootstrap Summary
##############################################################

bootstrap_auc <- bootstrap_results$AUC

bootstrap_sensitivity <-
  bootstrap_results$Sensitivity[
    !is.na(
      bootstrap_results$Sensitivity
    )
  ]

bootstrap_specificity <-
  bootstrap_results$Specificity[
    !is.na(
      bootstrap_results$Specificity
    )
  ]


##############################################################
# Bootstrap Confidence Intervals
##############################################################

auc_boot_ci <- quantile(
  
  bootstrap_auc,
  
  probs = c(
    0.025,
    0.975
  ),
  
  na.rm = TRUE
  
)


sensitivity_boot_ci <- quantile(
  
  bootstrap_sensitivity,
  
  probs = c(
    0.025,
    0.975
  ),
  
  na.rm = TRUE
  
)


specificity_boot_ci <- quantile(
  
  bootstrap_specificity,
  
  probs = c(
    0.025,
    0.975
  ),
  
  na.rm = TRUE
  
)


##############################################################
# Bootstrap Statistics
##############################################################

auc_mean <- mean(
  bootstrap_auc,
  na.rm = TRUE
)

auc_sd <- sd(
  bootstrap_auc,
  na.rm = TRUE
)

sensitivity_mean <- mean(
  bootstrap_sensitivity,
  na.rm = TRUE
)

specificity_mean <- mean(
  bootstrap_specificity,
  na.rm = TRUE
)


##############################################################
# Signature Score Comparison
##############################################################

cat("\n")
cat("Comparing signature score distributions...\n")

aml_scores <- scores |>
  
  dplyr::filter(
    AML_Status_Clean == "AML"
  ) |>
  
  dplyr::pull(
    Signature_Score
  )


non_aml_scores <- scores |>
  
  dplyr::filter(
    AML_Status_Clean == "Non_AML"
  ) |>
  
  dplyr::pull(
    Signature_Score
  )


##############################################################
# Wilcoxon Test
##############################################################

wilcox_test <- wilcox.test(
  
  aml_scores,
  
  non_aml_scores,
  
  exact = FALSE
  
)


##############################################################
# Score Summary
##############################################################

score_summary <- scores |>
  
  dplyr::group_by(
    AML_Status_Clean
  ) |>
  
  dplyr::summarise(
    
    N = n(),
    
    Mean = mean(
      Signature_Score,
      na.rm = TRUE
    ),
    
    SD = sd(
      Signature_Score,
      na.rm = TRUE
    ),
    
    Median = median(
      Signature_Score,
      na.rm = TRUE
    ),
    
    Q1 = quantile(
      Signature_Score,
      0.25,
      na.rm = TRUE
    ),
    
    Q3 = quantile(
      Signature_Score,
      0.75,
      na.rm = TRUE
    ),
    
    Min = min(
      Signature_Score,
      na.rm = TRUE
    ),
    
    Max = max(
      Signature_Score,
      na.rm = TRUE
    ),
    
    .groups = "drop"
    
  )


##############################################################
# Save Bootstrap Results
##############################################################

BOOTSTRAP_FILE <- file.path(
  
  TABLE_DIR,
  
  "GSE13159_Bootstrap_Results.csv"
  
)

readr::write_csv(
  
  bootstrap_results,
  
  BOOTSTRAP_FILE
  
)


##############################################################
# Save Score Summary
##############################################################

SCORE_SUMMARY_FILE <- file.path(
  
  TABLE_DIR,
  
  "GSE13159_Signature_Score_Summary.csv"
  
)

readr::write_csv(
  
  score_summary,
  
  SCORE_SUMMARY_FILE
  
)


##############################################################
# Save Statistical Test
##############################################################

WILCOX_FILE <- file.path(
  
  TABLE_DIR,
  
  "GSE13159_Signature_Score_Wilcoxon.csv"
  
)

wilcox_results <- data.frame(
  
  Test = "Wilcoxon rank sum test",
  
  Statistic =
    as.numeric(
      wilcox_test$statistic
    ),
  
  P_Value =
    wilcox_test$p.value,
  
  stringsAsFactors = FALSE
  
)

readr::write_csv(
  
  wilcox_results,
  
  WILCOX_FILE
  
)


##############################################################
# Bootstrap AUC Plot
##############################################################

BOOTSTRAP_PLOT <- file.path(
  
  PLOT_DIR,
  
  "GSE13159_Bootstrap_AUC_Distribution.png"
  
)

bootstrap_plot <- ggplot2::ggplot(
  
  data.frame(
    AUC = bootstrap_auc
  ),
  
  ggplot2::aes(
    x = AUC
  )
  
) +
  
  ggplot2::geom_histogram(
    
    bins = 40,
    
    color = "black",
    
    fill = "grey80"
    
  ) +
  
  ggplot2::geom_vline(
    
    xintercept = observed_auc,
    
    linetype = "dashed",
    
    linewidth = 1
    
  ) +
  
  ggplot2::geom_vline(
    
    xintercept = auc_boot_ci[1],
    
    linetype = "dotted"
    
  ) +
  
  ggplot2::geom_vline(
    
    xintercept = auc_boot_ci[2],
    
    linetype = "dotted"
    
  ) +
  
  ggplot2::labs(
    
    title =
      "Bootstrap Distribution of AUC",
    
    subtitle = paste0(
      
      "Observed AUC = ",
      
      round(
        observed_auc,
        4
      ),
      
      " | Bootstrap 95% CI = ",
      
      round(
        auc_boot_ci[1],
        4
      ),
      
      " to ",
      
      round(
        auc_boot_ci[2],
        4
      )
      
    ),
    
    x = "Bootstrap AUC",
    
    y = "Frequency"
    
  ) +
  
  ggplot2::theme_classic()


ggplot2::ggsave(
  
  BOOTSTRAP_PLOT,
  
  bootstrap_plot,
  
  width = 8,
  
  height = 6,
  
  dpi = 300
  
)


##############################################################
# Signature Score Boxplot
##############################################################

BOXPLOT_FILE <- file.path(
  
  PLOT_DIR,
  
  "GSE13159_Signature_Score_Distribution.png"
  
)

boxplot <- ggplot2::ggplot(
  
  scores,
  
  ggplot2::aes(
    
    x = AML_Status_Clean,
    
    y = Signature_Score
    
  )
  
) +
  
  ggplot2::geom_boxplot(
    
    outlier.shape = NA
    
  ) +
  
  ggplot2::geom_jitter(
    
    width = 0.15,
    
    alpha = 0.15,
    
    size = 1
    
  ) +
  
  ggplot2::labs(
    
    title =
      "AML Biomarker Signature Score Distribution",
    
    subtitle = paste0(
      
      "Wilcoxon P = ",
      
      format.pval(
        wilcox_test$p.value,
        digits = 3
      )
      
    ),
    
    x = "Disease Status",
    
    y = "Signature Score"
    
  ) +
  
  ggplot2::theme_classic()


ggplot2::ggsave(
  
  BOXPLOT_FILE,
  
  boxplot,
  
  width = 8,
  
  height = 6,
  
  dpi = 300
  
)


##############################################################
# Density Plot
##############################################################

DENSITY_FILE <- file.path(
  
  PLOT_DIR,
  
  "GSE13159_Signature_Score_Density.png"
  
)

density_plot <- ggplot2::ggplot(
  
  scores,
  
  ggplot2::aes(
    
    x = Signature_Score,
    
    fill = AML_Status_Clean
    
  )
  
) +
  
  ggplot2::geom_density(
    
    alpha = 0.35
    
  ) +
  
  ggplot2::labs(
    
    title =
      "AML Signature Score Density",
    
    x = "Signature Score",
    
    y = "Density",
    
    fill = "Status"
    
  ) +
  
  ggplot2::theme_classic()


ggplot2::ggsave(
  
  DENSITY_FILE,
  
  density_plot,
  
  width = 8,
  
  height = 6,
  
  dpi = 300
  
)


##############################################################
# Robustness Interpretation
##############################################################

if (
  
  auc_boot_ci[1] > 0.70
  
) {
  
  robustness_statement <-
    "The bootstrap confidence interval remained above 0.70, supporting stable discriminatory performance."
  
} else if (
  
  auc_boot_ci[1] > 0.50
  
) {
  
  robustness_statement <-
    "The bootstrap confidence interval remained above random discrimination but indicates moderate robustness."
  
} else {
  
  robustness_statement <-
    "The bootstrap confidence interval approached or crossed random discrimination and warrants caution."
  
}


##############################################################
# Summary Report
##############################################################

SUMMARY_FILE <- file.path(
  
  ROBUSTNESS_DIR,
  
  "Stage14E_External_Validation_Robustness_Summary.txt"
  
)


summary_lines <- c(
  
  "AML EXTERNAL VALIDATION ROBUSTNESS ANALYSIS",
  
  "=========================================",
  
  "",
  
  "Dataset : GSE13159",
  
  paste(
    "Total Samples :",
    TOTAL_N
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
  
  "OBSERVED ROC PERFORMANCE",
  
  paste(
    "Observed AUC :",
    round(
      observed_auc,
      4
    )
  ),
  
  paste(
    "Observed AUC 95% CI :",
    round(
      observed_auc_ci[1],
      4
    ),
    "to",
    round(
      observed_auc_ci[3],
      4
    )
  ),
  
  "",
  
  "BOOTSTRAP VALIDATION",
  
  paste(
    "Bootstrap Iterations :",
    N_BOOTSTRAP
  ),
  
  paste(
    "Successful Iterations :",
    nrow(
      bootstrap_results
    )
  ),
  
  paste(
    "Mean Bootstrap AUC :",
    round(
      auc_mean,
      4
    )
  ),
  
  paste(
    "Bootstrap AUC SD :",
    round(
      auc_sd,
      4
    )
  ),
  
  paste(
    "Bootstrap AUC 95% CI :",
    round(
      auc_boot_ci[1],
      4
    ),
    "to",
    round(
      auc_boot_ci[2],
      4
    )
  ),
  
  "",
  
  "BOOTSTRAP CLASSIFICATION PERFORMANCE",
  
  paste(
    "Mean Sensitivity :",
    round(
      sensitivity_mean,
      4
    )
  ),
  
  paste(
    "Sensitivity 95% CI :",
    round(
      sensitivity_boot_ci[1],
      4
    ),
    "to",
    round(
      sensitivity_boot_ci[2],
      4
    )
  ),
  
  paste(
    "Mean Specificity :",
    round(
      specificity_mean,
      4
    )
  ),
  
  paste(
    "Specificity 95% CI :",
    round(
      specificity_boot_ci[1],
      4
    ),
    "to",
    round(
      specificity_boot_ci[2],
      4
    )
  ),
  
  "",
  
  "SIGNATURE SCORE COMPARISON",
  
  paste(
    "Wilcoxon P Value :",
    format.pval(
      wilcox_test$p.value,
      digits = 4
    )
  ),
  
  "",
  
  "INTERPRETATION",
  
  robustness_statement,
  
  "",
  
  "OUTPUTS",
  
  paste(
    "Bootstrap Results :",
    BOOTSTRAP_FILE
  ),
  
  paste(
    "Score Summary :",
    SCORE_SUMMARY_FILE
  ),
  
  paste(
    "Wilcoxon Test :",
    WILCOX_FILE
  ),
  
  paste(
    "Bootstrap AUC Plot :",
    BOOTSTRAP_PLOT
  ),
  
  paste(
    "Score Distribution Plot :",
    BOXPLOT_FILE
  ),
  
  paste(
    "Score Density Plot :",
    DENSITY_FILE
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
cat("Stage 14E Completed Successfully\n")
cat("=========================================\n\n")

cat(
  "Dataset             : GSE13159\n"
)

cat(
  "Total Samples       :",
  TOTAL_N,
  "\n"
)

cat(
  "AML Samples         :",
  AML_N,
  "\n"
)

cat(
  "Non AML Samples     :",
  NON_AML_N,
  "\n\n"
)

cat(
  "Observed AUC        :",
  round(
    observed_auc,
    4
  ),
  "\n"
)

cat(
  "Observed 95% CI     :",
  round(
    observed_auc_ci[1],
    4
  ),
  "to",
  round(
    observed_auc_ci[3],
    4
  ),
  "\n\n"
)

cat(
  "Bootstrap Iterations:",
  N_BOOTSTRAP,
  "\n"
)

cat(
  "Successful          :",
  nrow(
    bootstrap_results
  ),
  "\n"
)

cat(
  "Mean Bootstrap AUC  :",
  round(
    auc_mean,
    4
  ),
  "\n"
)

cat(
  "Bootstrap 95% CI    :",
  round(
    auc_boot_ci[1],
    4
  ),
  "to",
  round(
    auc_boot_ci[2],
    4
  ),
  "\n\n"
)

cat(
  "Mean Sensitivity    :",
  round(
    sensitivity_mean,
    4
  ),
  "\n"
)

cat(
  "Mean Specificity    :",
  round(
    specificity_mean,
    4
  ),
  "\n\n"
)

cat(
  "Wilcoxon P Value    :",
  format.pval(
    wilcox_test$p.value,
    digits = 4
  ),
  "\n\n"
)

cat(
  "Robustness Summary  :",
  SUMMARY_FILE,
  "\n\n"
)

cat(
  "Time Elapsed        :",
  round(
    Sys.time() - start_time,
    2
  ),
  "\n\n"
)