###############################################################
# AML Biomarker Discovery Pipeline
# Script: 17_Save_Session_Information.R
#
# Purpose:
# Save session information, installed packages,
# pipeline version and reproducibility information.
#
# Author: Isreal Oluwafemi Abiodun
###############################################################

# rm(list = ls())

###############################################################
# Load Configuration
###############################################################

source("config.R")

###############################################################
# Create Log Directory
###############################################################

dir.create(
  LOG_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("---------------------------------------\n")
cat("Saving Session Information\n")
cat("---------------------------------------\n\n")

###############################################################
# Save sessionInfo()
###############################################################

sink(
  file.path(
    LOG_DIR,
    "SessionInfo.txt"
  )
)

sessionInfo()

sink()

###############################################################
# Save Installed Packages
###############################################################

pkgs <- installed.packages()

package_table <- data.frame(
  
  Package = pkgs[, "Package"],
  
  Version = pkgs[, "Version"],
  
  stringsAsFactors = FALSE
  
)

package_table <- package_table[
  order(package_table$Package),
]

write.csv(
  
  package_table,
  
  file.path(
    LOG_DIR,
    "Installed_Packages.csv"
  ),
  
  row.names = FALSE
  
)

###############################################################
# Save Pipeline Information
###############################################################

pipeline_info <- c(
  
  paste("Pipeline :", PIPELINE_NAME),
  
  paste("Version :", PIPELINE_VERSION),
  
  paste("Author :", AUTHOR),
  
  paste("Date :", Sys.Date()),
  
  paste("R Version :", R.version.string),
  
  paste("Working Directory :", getwd()),
  
  paste("Random Seed :", RANDOM_SEED),
  
  paste("Number of Threads :", N_THREADS)
  
)

writeLines(
  
  pipeline_info,
  
  file.path(
    LOG_DIR,
    "Pipeline_Version.txt"
  )
  
)

###############################################################
# Save Timestamp
###############################################################

writeLines(
  
  paste(
    "Pipeline completed on:",
    Sys.time()
  ),
  
  file.path(
    LOG_DIR,
    "Completion_Time.txt"
  )
  
)

###############################################################
# Summary
###############################################################

cat("---------------------------------------\n")
cat("Session Information Saved Successfully\n")
cat("---------------------------------------\n\n")

cat("Files created:\n")

cat("results/logs/SessionInfo.txt\n")
cat("results/logs/Installed_Packages.csv\n")
cat("results/logs/Pipeline_Version.txt\n")
cat("results/logs/Completion_Time.txt\n")