##############################################################################
# Alnus and Salix shrub anatomical response to climate
##scripts/00_paths_setup.R

# Root project directory
root <- here::here()

# Define subfolders 
data_dir    <- file.path(root, "data")
output_dir  <- file.path(root, "output")
figures_dir <- file.path(root, "figures")

dir.create(output_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

message("Project root: ", root)
message("Data dir: ", data_dir)
message("Output dir: ", output_dir)
message("Figures dir: ", figures_dir)
