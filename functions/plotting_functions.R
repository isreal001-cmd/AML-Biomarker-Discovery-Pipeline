###############################################################
# AML Biomarker Discovery Pipeline
# plotting_functions.R
###############################################################

library(ggplot2)

publication_theme <- function(){
  
  theme_bw(base_size = 14) +
    
    theme(
      
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 16
      ),
      
      axis.title = element_text(
        face = "bold"
      ),
      
      axis.text = element_text(
        colour = "black"
      ),
      
      legend.title = element_text(
        face = "bold"
      ),
      
      panel.grid.major = element_blank(),
      
      panel.grid.minor = element_blank(),
      
      panel.border = element_rect(
        colour = "black",
        fill = NA
      )
      
    )
  
}