## wood anatomical response to climate in Alnus and Salix
##Author: Georgia Hole
##scripts/00_paths_setup.R

source("R/utils_packages.R")
load_project_packages()

# Root is the project directory (assuming you open the .Rproj here)
root <- here::here()

# Define subfolders
data_dir    <- file.path(root, "Data")
results_dir <- file.path(root, "results")
figs_dir     <- file.path(root, "figures")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figs_dir,     showWarnings = FALSE, recursive = TRUE)

message("Project root: ", root)
message("Data dir: ", data_dir)
message("Results dir: ", results_dir)
message("Figures dir: ", figs_dir)
