##############################################################
# AML Biomarker Discovery Pipeline
# Save Session Information
#
# Author  : Isreal Oluwafemi Abiodun
# Version : 1.0.0
##############################################################

source("config.R")

cat("---------------------------------------\n")
cat("Saving Session Information\n")
cat("---------------------------------------\n")

##############################################################
# Create Log Directory
##############################################################

dir.create(
  LOG_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

##############################################################
# Save sessionInfo()
##############################################################

sink(file.path(LOG_DIR, "SessionInfo.txt"))

sessionInfo()

sink()

##############################################################
# Save Installed Packages
##############################################################

pkgs <- installed.packages()

package_table <- data.frame(
  
  Package = pkgs[, "Package"],
  
  Version = pkgs[, "Version"],
  
  stringsAsFactors = FALSE
  
)

package_table <- package_table[order(package_table$Package), ]

write.csv(
  
  package_table,
  
  file.path(LOG_DIR, "Installed_Packages.csv"),
  
  row.names = FALSE
  
)

##############################################################
# Save Pipeline Version
##############################################################

writeLines(
  
  c(
    
    paste("Pipeline :", PIPELINE_NAME),
    
    paste("Version  :", PIPELINE_VERSION),
    
    paste("Author   :", AUTHOR),
    
    paste("Date     :", Sys.Date()),
    
    paste("R        :", R.version.string)
    
  ),
  
  file.path(LOG_DIR, "Pipeline_Version.txt")
  
)

cat("Session information saved.\n")