################################################################################
# IO Utilities for scRNA-seq Pipeline
################################################################################

#' @title Load and Validate PMH Manifest
#' @description Reads the input CSV and verifies that all specific data_paths exist.
#' @param path String. Path to the manifest CSV file.
#' @return A data frame containing the manifest data.
#' @export
load_manifest <- function(path) {
  if(!file.exists(path)) {
    stop(paste("Execution Halted: Manifest file not found at", path))
  }

  manifest <- read.csv(path, stringsAsFactors=FALSE) 

  # check whether required columns exist
  required_cols <- c("sample_name", "ident1", "ident2", "data_path")
  if (!all(required_cols %in% colnames(manifest))) {
    missing <- setdiff(required_cols, colnames(manifest))
    stop(paste("CRITICAL ERROR: Manifest is missing columns:", paste(missing, collapse=", ")))
  }

  # verify each directory path
  for (i in 1:nrow(manifest)) {
    if (!dir.exists(manifest$data_path[i])) {
      warning(paste("PATH WARNING: Folder not found for", 
                    manifest$sample_name[i], "at", manifest$data_path[i]))
    }
  }
  
  message(paste("Successfully loaded manifest with", nrow(manifest), "samples."))
  return(manifest) 

}

#' @title Load 10x Data with Metadata
#' @description Imports raw counts and injects ident1/ident2 metadata immediately.
#' @param data_path String. Path to the 10x directory.
#' @param sample_name String. Project name for the Seurat object.
#' @param ident1 String. Condition metadata (e.g., Healthy_Fresh).
#' @param ident2 String. Patient/Sample ID (e.g., HTY131).
#' @return A Seurat object with metadata attached.
#' @export
load_10x_data <- function(data_path, sample_name, ident1, ident2) {
  message(paste("Loading 10x data for:", sample_name))
  
  counts <- Seurat::Read10X(data.dir = data_path)
  obj <- Seurat::CreateSeuratObject(
    counts = counts, 
    project = sample_name,
    min.cells = 3, 
    min.features = 200
  )
  
  # inject CSV metadata directly into the object
  obj$condition <- ident1
  obj$patient_id <- ident2
  
  return(obj)
}

#' @title Save Seurat Object
#' @description Saves object with a timestamp to prevent accidental overwrites.
#' @export
save_seurat_object <- function(obj, output_dir, step_name) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
  file_path <- file.path(output_dir, paste0(step_name, "_", timestamp, ".rds"))
  
  message(paste("Saving result to:", file_path))
  saveRDS(obj, file = file_path)
}