###############################################################
# AML Biomarker Discovery Pipeline
# Script: 16_Validation_Visualization.R
#
# Purpose:
# Generate publication-quality visualizations for
# external validation results.
###############################################################

# rm(list = ls())

###############################################################
# Load Configuration
###############################################################

source("config.R")
source("functions/plotting_functions.R")

###############################################################
# Required Packages
###############################################################

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

###############################################################
# Create Output Folder
###############################################################

dir.create(
  "results/validation/plots",
  recursive = TRUE,
  showWarnings = FALSE
)

###############################################################
# Load Results
###############################################################

detection <- read_csv(
  "results/validation/Biomarker_Detection.csv",
  show_col_types = FALSE
)

statistics <- read_csv(
  "results/validation/Expression_Statistics.csv",
  show_col_types = FALSE
)

###############################################################
# Detection Rate Plot
###############################################################

p1 <- detection %>%
  
  arrange(desc(DetectionRate)) %>%
  
  ggplot(
    aes(
      x = reorder(GeneID, DetectionRate),
      y = DetectionRate
    )
  ) +
  
  geom_col(fill = "#2E86AB") +
  
  coord_flip() +
  
  labs(
    title = "Biomarker Detection Rate",
    x = "",
    y = "Detection Rate"
  ) +
  
  publication_theme()

ggsave(
  "results/validation/plots/Detection_Rate.png",
  p1,
  width = FIG_WIDTH,
  height = FIG_HEIGHT,
  dpi = FIG_DPI
)

ggsave(
  "results/validation/plots/Detection_Rate.pdf",
  p1,
  width = FIG_WIDTH,
  height = FIG_HEIGHT
)

###############################################################
# Mean Expression Plot
###############################################################

p2 <- statistics %>%
  
  arrange(desc(Mean)) %>%
  
  ggplot(
    aes(
      x = reorder(GeneID, Mean),
      y = Mean
    )
  ) +
  
  geom_col(fill = "#D35400") +
  
  coord_flip() +
  
  labs(
    title = "Mean Expression of Validated Biomarkers",
    x = "",
    y = "Mean Counts"
  ) +
  
  publication_theme()

ggsave(
  "results/validation/plots/Mean_Expression.png",
  p2,
  width = FIG_WIDTH,
  height = FIG_HEIGHT,
  dpi = FIG_DPI
)

ggsave(
  "results/validation/plots/Mean_Expression.pdf",
  p2,
  width = FIG_WIDTH,
  height = FIG_HEIGHT
)

###############################################################
# Top 20 Detection Plot
###############################################################

top20 <- detection %>%
  
  arrange(desc(DetectionRate)) %>%
  
  slice(1:20)

p3 <- ggplot(
  top20,
  aes(
    x = reorder(GeneID, DetectionRate),
    y = DetectionRate
  )
) +
  
  geom_col(fill = "#27AE60") +
  
  coord_flip() +
  
  labs(
    title = "Top 20 Validated Biomarkers",
    x = "",
    y = "Detection Rate"
  ) +
  
  publication_theme()

ggsave(
  "results/validation/plots/Top20_Validated_Biomarkers.png",
  p3,
  width = FIG_WIDTH,
  height = FIG_HEIGHT,
  dpi = FIG_DPI
)

ggsave(
  "results/validation/plots/Top20_Validated_Biomarkers.pdf",
  p3,
  width = FIG_WIDTH,
  height = FIG_HEIGHT
)

###############################################################
# Summary
###############################################################

cat("---------------------------------------\n")
cat("Validation Visualization Completed\n")
cat("---------------------------------------\n\n")

cat("Files created:\n")

cat("results/validation/plots/Detection_Rate.png\n")
cat("results/validation/plots/Detection_Rate.pdf\n")
cat("results/validation/plots/Mean_Expression.png\n")
cat("results/validation/plots/Mean_Expression.pdf\n")
cat("results/validation/plots/Top20_Validated_Biomarkers.png\n")
cat("results/validation/plots/Top20_Validated_Biomarkers.pdf\n")