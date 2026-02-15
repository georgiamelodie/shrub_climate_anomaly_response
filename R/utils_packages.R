## Alnus and Salix shrub anatomical response to climate
## Author: Georgia Hole
## scripts/paths_setup.R
## R/utils_packages.R


load_project_packages <- function() {
  
  pkgs <- c(
    # base utilities used directly
    "utils", "here",
    
    # data handling 
    "dplyr", "tidyr", "readr", "tibble", "stringr", 
    
    # dates
    "lubridate",
    
    # plotting
    "ggplot2", "patchwork",
    
    # dendro / climate-growth analysis
    "dplR", "treeclim",
    
    # diagnostics / correlation plots
    "corrplot", "DHARMa",
    
    # modelling and summaries
    "lme4", "MASS", "MuMIn", "broom.mixed"
  )
  
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      message("Installing package: ", p)
      install.packages(p, dependencies = TRUE)
    }
    suppressPackageStartupMessages(
      library(p, character.only = TRUE)
    )
  }
  
  # Type-III sums-of-squares by default (if you use car::Anova etc.)
  options(contrasts = c("contr.sum", "contr.poly"))
}

