###############################################################
# AML Biomarker Discovery Pipeline
# File: validation_functions.R
#
# Reusable functions for external validation
###############################################################

###############################################################
# Check Required Files
###############################################################

check_validation_dataset <- function(dataset){
  
  dataset_path <- file.path(
    "data",
    "validation",
    dataset
  )
  
  expr_file <- file.path(
    dataset_path,
    "expression.csv"
  )
  
  meta_file <- file.path(
    dataset_path,
    "metadata.csv"
  )
  
  if(!file.exists(expr_file))
    stop("expression.csv not found.")
  
  if(!file.exists(meta_file))
    stop("metadata.csv not found.")
  
  invisible(TRUE)
  
}

###############################################################
# Load Validation Dataset
###############################################################

load_validation_dataset <- function(dataset){
  
  check_validation_dataset(dataset)
  
  expr <- readr::read_csv(
    file.path(
      "data",
      "validation",
      dataset,
      "expression.csv"
    ),
    show_col_types = FALSE
  )
  
  metadata <- readr::read_csv(
    file.path(
      "data",
      "validation",
      dataset,
      "metadata.csv"
    ),
    show_col_types = FALSE
  )
  
  list(
    
    expression = expr,
    
    metadata = metadata
    
  )
  
}

###############################################################
# Load Biomarkers
###############################################################

load_biomarkers <- function(){
  
  readr::read_csv(
    
    "results/biomarkers/Top50_Biomarkers.csv",
    
    show_col_types = FALSE
    
  )
  
}

###############################################################
# Match Biomarkers
###############################################################

match_biomarkers <- function(expr, biomarkers){
  
  expr %>%
    
    dplyr::filter(
      
      GeneID %in% biomarkers$GeneID
      
    )
  
}

###############################################################
# Detect Dataset Type
###############################################################

detect_dataset_type <- function(metadata){
  
  if(!"Condition" %in% names(metadata))
    return("single_group")
  
  groups <- unique(metadata$Condition)
  
  if(length(groups) < 2)
    return("single_group")
  
  "two_group"
  
}

###############################################################
# Detection Rate
###############################################################

calculate_detection_rate <- function(expr){
  
  detected <- rowSums(
    
    expr[,-1] > 0
    
  )
  
  data.frame(
    
    GeneID = expr$GeneID,
    
    SamplesDetected = detected,
    
    DetectionRate = detected /
      (ncol(expr)-1)
    
  )
  
}

###############################################################
# Basic Expression Statistics
###############################################################

calculate_expression_statistics <- function(expr){
  
  values <- expr[,-1]
  
  data.frame(
    
    GeneID = expr$GeneID,
    
    Mean = apply(values,1,mean),
    
    Median = apply(values,1,median),
    
    SD = apply(values,1,sd),
    
    Max = apply(values,1,max),
    
    Min = apply(values,1,min)
    
  )
  
}