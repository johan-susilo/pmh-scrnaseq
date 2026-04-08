#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(harmony)
  library(clustree)
  library(ggplot2)
  library(future)
})

options(future.globals.maxSize = 100 * 1024^3)
# -----------------------------------------------------

# =========================================================================
# Define Command-Line Arguments
# =========================================================================
option_list <- list(
  make_option(c("--project_dir"), type="character", default=NULL, help="Path to Nextflow project"),
  make_option(c("--output_dir"), type="character", default="./output", help="Output directory"),
  make_option(c("--resolution"), type="numeric", default=0.4, help="Default clustering resolution"),
  make_option(c("--use_sct"), type="logical", default=TRUE, help="Use SCTransform"),
  make_option(c("--batch_var"), type="character", default="orig.ident2", help="Metadata column for Harmony batch correction")
)

opt_parser <- OptionParser(option_list = option_list)
parsed_args <- parse_args(opt_parser, positional_arguments = TRUE) 

# separate the named options from the positional file arguments
opt <- parsed_args$options
file_paths <- parsed_args$args

if (!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)

source(file.path(opt$project_dir, "utils", "plot.R"))

message("\n==================== Starting Harmony Integration ====================")
message("Loading ", length(file_paths), " clean samples for integration...")

sample_list <- lapply(file_paths, readRDS)

# prefix barcodes to prevent Seurat from crashing if two samples have the exact same 10x barcode
sample_list <- lapply(seq_along(sample_list), function(i) {
  obj <- sample_list[[i]]
  
  # safely look for identifiers without triggering Seurat v5 aborts
  tag <- NA
  meta <- obj@meta.data
  if ("orig.ident2" %in% colnames(meta)) {
    tag <- as.character(unique(meta$orig.ident2)[1])
  } else if ("orig.ident1" %in% colnames(meta)) {
    tag <- as.character(unique(meta$orig.ident1)[1])
  } else {
    tag <- as.character(unique(meta$orig.ident)[1])
  }
  
  if (is.null(tag) || is.na(tag)) tag <- paste0("sample", i)
  
  # Ensure the tag is safe for column names
  tag_clean <- make.names(tag)
  colnames(obj) <- paste0(tag_clean, "_", colnames(obj))
  
  DefaultAssay(obj) <- "RNA"
  if ("SCT" %in% names(obj)) obj[["SCT"]] <- NULL
  return(obj)
})

# =========================================================================
# Merge & Cell Cycle Scoring
# =========================================================================
message("Merging objects...")

# merge multiple seurat_obj if it contains more than 1 object
if (length(sample_list) > 1) {
  TN.combined <- merge(x = sample_list[[1]], y = sample_list[2:length(sample_list)])
} else {
  TN.combined <- sample_list[[1]]
}

# use trycatch for join layers so it can capture errors
TN.combined <- tryCatch({
  JoinLayers(TN.combined)
}, error = function(e) {
  message("JoinLayers skipped: ", e$message)
  TN.combined  # return the original object unchanged
})

message("Calculating Cell Cycle Scores...")
TN.combined <- NormalizeData(TN.combined, assay = "RNA", verbose = FALSE)

# assess each cell's position in the cell cycle (G1, S, G2/M)
# cell cycle effects can confound downstream analyses like clustering or integration
if (!exists("cc.genes")) data("cc.genes", package = "Seurat")
if (length(cc.genes$s.genes) > 0 && length(cc.genes$g2m.genes) > 0){
  TN.combined <- CellCycleScoring(
    TN.combined, 
    s.features = cc.genes$s.genes, 
    g2m.features = cc.genes$g2m.genes
  )
}

# =========================================================================
# Global Normalization & PCA
# =========================================================================

meta <- TN.combined@meta.data
if (opt$batch_var %in% colnames(meta)) {
  TN.combined@meta.data$batch <- meta[[opt$batch_var]]
} else {
  message("Warning: ", opt$batch_var, " not found. Falling back to orig.ident")
  TN.combined@meta.data$batch <- meta$orig.ident
}

if (isTRUE(opt$use_sct)) {
  message("Running Global SCTransform v2...")

  vars_regress <- c("percent.mt") # always regress percent mitochondrial genes
  if ("S.Score" %in% colnames(TN.combined@meta.data))
    vars_regress <- c(vars_regress, "S.Score", "G2M.Score") # regress cell cycle effects if available

  # run SCTransform normalization with glmGamPoi method and v2 flavor
  TN.combined <- SCTransform(
    TN.combined, 
    assay = "RNA", 
    vars.to.regress = vars_regress, 
    method = "glmGamPoi", 
    vst.flavor = "v2", 
    verbose = FALSE
  )

  # prepare for downstream marker detection
  TN.combined <- PrepSCTFindMarkers(TN.combined,  assay = "SCT", verbose = FALSE)

  # run pca on normalized data
  TN.combined <- RunPCA(TN.combined, assay = "SCT", npcs = 50, verbose = FALSE)
  # note which assay to use for Harmony batch correction
  harmony_assay <- "SCT"

} else {
  message("Running standard LogNormalize...")
  TN.combined <- NormalizeData(TN.combined) %>% 
    FindVariableFeatures() %>% 
    ScaleData(vars.to.regress = "percent.mt") %>% 
    RunPCA(npcs = 50, verbose = FALSE)

  harmony_assay <- "RNA"
}

# =========================================================================
# Harmony Batch Correction & Clustering
# =========================================================================

message("Running Harmony Integration...")
TN.combined <- RunHarmony(TN.combined, 
  group.by.vars = "batch", 
  assay.use = harmony_assay, 
  verbose = FALSE
)

message("Running UMAP and Multi-Resolution Clustering...")
TN.combined <- RunUMAP(TN.combined, 
                       reduction = "harmony", 
                       dims = 1:30, 
                       verbose = FALSE, 
                       umap.method = "uwot", 
                       metric = "cosine"
                       # metric cosine preserves similarity of expression patterns 
                       # across high-dimensional PCA space.
                       )


TN.combined <- FindNeighbors(TN.combined, reduction = "harmony", dims = 1:30, verbose = FALSE)

# precompute multiple resolutions so we don't need to rerun just to get different resolutions
resolutions <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8)
TN.combined <- FindClusters(TN.combined, resolution = resolutions, verbose = FALSE)

# =========================================================================
# Diagnostics & Final Save
# =========================================================================
cluster_prefix <- paste0(harmony_assay, "_snn_res.")

message("Generating Clustree...")
p_tree <- clustree(TN.combined, prefix = cluster_prefix) +
  ggtitle("Global Integration Clustree")
safe_save_plot(p_tree, file.path(opt$output_dir, "Clustree"), w = 12, h = 10)

# set default resolution
Idents(TN.combined) <- paste0(cluster_prefix, opt$resolution)

output_rds <- file.path(opt$output_dir, "TN.combined_dim30.rds")
saveRDS(TN.combined, output_rds)

message(paste("\nSuccessfully integrated data and saved to:", output_rds))