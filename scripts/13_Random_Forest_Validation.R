###############################################################
# AML Biomarker Discovery Pipeline
# Script: 13_Random_Forest_Validation.R
# Purpose:
# Random Forest validation of Top 20 biomarker candidates
###############################################################

# rm(list = ls())

###############################################################
# Load configuration
###############################################################

source("config.R")
source("functions/plotting_functions.R")

###############################################################
# Load packages
###############################################################

packages <- c(
  "randomForest",
  "DESeq2",
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

###############################################################
# Reproducibility
###############################################################

set.seed(12345)

###############################################################
# Load input files
###############################################################

vsd <- readRDS(
  "results/vst/vsd.rds"
)

metadata <- read_csv(
  "data/metadata/sample_metadata.csv",
  show_col_types = FALSE
)

biomarkers <- read_csv(
  "results/biomarkers/Top20_Biomarkers.csv",
  show_col_types = FALSE
)

cat("---------------------------------------\n")
cat("Input Files Loaded Successfully\n")
cat("---------------------------------------\n\n")

###############################################################
# Prepare expression matrix
###############################################################

expr <- assay(vsd)

expr <- expr[
  biomarkers$GeneID,
]

expr <- t(expr)

rf_data <- as.data.frame(expr)

rf_data$Condition <- factor(
  metadata$Condition
)

###############################################################
# Train Random Forest
###############################################################

cat("Training Random Forest...\n\n")

rf_model <- randomForest(
  
  Condition ~ .,
  
  data = rf_data,
  
  importance = TRUE,
  
  ntree = 1000
  
)

cat("Training completed.\n\n")

###############################################################
# Variable importance
###############################################################

importance_table <- data.frame(
  
  GeneID = rownames(
    importance(rf_model)
  ),
  
  MeanDecreaseAccuracy =
    importance(rf_model)[,"MeanDecreaseAccuracy"],
  
  MeanDecreaseGini =
    importance(rf_model)[,"MeanDecreaseGini"]
  
)

importance_table <- importance_table %>%
  
  arrange(desc(MeanDecreaseAccuracy))

###############################################################
# Output folder
###############################################################

dir.create(
  
  "results/random_forest",
  
  recursive = TRUE,
  
  showWarnings = FALSE
  
)

###############################################################
# Save outputs
###############################################################

saveRDS(
  
  rf_model,
  
  "results/random_forest/random_forest_model.rds"
  
)

write.csv(
  
  importance_table,
  
  "results/random_forest/Variable_Importance.csv",
  
  row.names = FALSE
  
)

###############################################################
# Plot Top 20
###############################################################

top20 <- importance_table %>%
  
  dplyr::slice(1:20)

p <- ggplot(
  
  top20,
  
  aes(
    
    x = reorder(
      GeneID,
      MeanDecreaseAccuracy
    ),
    
    y = MeanDecreaseAccuracy
    
  )
  
) +
  
  geom_col(
    
    fill = AML_RED
    
  ) +
  
  coord_flip() +
  
  labs(
    
    title = "Random Forest Variable Importance",
    
    x = "",
    
    y = "Mean Decrease Accuracy"
    
  ) +
  
  publication_theme()

ggsave(
  
  filename = "results/random_forest/Variable_Importance.png",
  
  plot = p,
  
  width = FIG_WIDTH,
  
  height = FIG_HEIGHT,
  
  dpi = FIG_DPI
  
)

ggsave(
  
  filename = "results/random_forest/Variable_Importance.pdf",
  
  plot = p,
  
  width = FIG_WIDTH,
  
  height = FIG_HEIGHT
  
)

###############################################################
# Console summary
###############################################################

cat("---------------------------------------\n")
cat("Random Forest Validation Completed\n")
cat("---------------------------------------\n\n")

cat("Top Biomarker:\n")

print(top20$GeneID[1])

cat("\nTop Mean Decrease Accuracy:\n")

print(top20$MeanDecreaseAccuracy[1])

cat("\nFiles created:\n")

cat("results/random_forest/random_forest_model.rds\n")
cat("results/random_forest/Variable_Importance.csv\n")
cat("results/random_forest/Variable_Importance.png\n")
cat("results/random_forest/Variable_Importance.pdf\n")