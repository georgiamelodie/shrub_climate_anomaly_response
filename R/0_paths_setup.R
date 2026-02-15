## wood anatomical response to climate in Alnus and Salix
##Author: Georgia Hole
##scripts/0_paths_setup.R

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










list_files <- function(path = ".", depth = 3, prefix = "") {
  if (depth == 0) return(NULL)
  
  items <- list.files(path, all.files = FALSE)
  items <- items[!items %in% c(".git", ".Rhistory", ".RData")]
  
  out <- character()
  
  for (i in seq_along(items)) {
    item <- items[i]
    full_path <- file.path(path, item)
    connector <- if (i == length(items)) "└── " else "├── "
    
    out <- c(out, paste0(prefix, connector, item))
    
    if (file.info(full_path)$isdir) {
      extension <- if (i == length(items)) "    " else "│   "
      out <- c(out, list_files(full_path, depth - 1, paste0(prefix, extension)))
    }
  }
  
  out
}

cat(paste(list_files(".", depth = 3), collapse = "\n"))
