#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(clustree)
  library(harmony)
  library(ggrepel)
})

source("workflow/scripts/00_utils.R")

# ==============================================================================
# 1. COMMAND-LINE INTERFACE
# ==============================================================================
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Path to globally annotated RDS (e.g., TN.combined_annotated.rds)"),
  make_option(c("-c", "--celltype"), type = "character", help = "Target cell type string to subset (e.g., 'fibroblast')"),
  make_option(c("-o", "--outdir"), type = "character", help = "Base output directory for subset clusters"),
  make_option(c("-r", "--resolution"), type = "numeric", default = 0.2, help = "Default clustering resolution for the subset"),
  make_option(c("--pcs"), type = "integer", default = 20, help = "Number of PCs for subsetting"),
  make_option(c("--resolutions"), type = "character", default = "0.2,0.4,0.6,0.8,1.0,1.2", help = "Comma-separated resolutions to calculate"),
  make_option(c("--seed"), type = "integer", default = 42, help = "Global random seed for reproducibility")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input) || is.null(opt$celltype) || is.null(opt$outdir)) {
  stop("Missing required arguments: --input, --celltype, or --outdir")
}

PROOF_GENES <- c("LUM", "DCN", "COL1A1", "KRT14", "KRT1", "PTPRC", "CD68", "PECAM1", "VWF")

# Parse --resolutions ("0.2,0.4,0.6,0.8,1.0,1.2") into a numeric vector once.
# This is what actually gets passed down into FindClusters().
resolutions_vec <- as.numeric(trimws(strsplit(opt$resolutions, ",")[[1]]))
if (any(is.na(resolutions_vec))) {
  stop("Could not parse --resolutions into numeric values: ", opt$resolutions)
}

# ==============================================================================
# 2. HELPER FUNCTIONS & CORE PROCESSING
# ==============================================================================
create_proportion_barplot <- function(seurat_obj, group_col, title) {
  prop_df <- as.data.frame(prop.table(table(Idents(seurat_obj), seurat_obj[[group_col]][[1]]), margin = 2))
  colnames(prop_df) <- c("Cluster", "Group", "Proportion")
  ggplot(prop_df %>% filter(Proportion > 0), aes(x = Group, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Cluster, size = Proportion), position = position_stack(vjust = 0.5)) +
    scale_size_continuous(range = c(2, 6), guide = "none") + 
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_y_continuous(labels = scales::percent) +
    ggtitle(title)
}

run_dim_reduction <- function(sub_obj, resolutions) {
  DefaultAssay(sub_obj) <- "RNA"
  sub_obj <- SCTransform(sub_obj, method = "glmGamPoi", vst.flavor = "v2", verbose = FALSE)
  sub_obj <- RunPCA(sub_obj, assay = "SCT", seed.use = opt$seed, verbose = FALSE)
  sub_obj <- RunHarmony(sub_obj, group.by.vars = "orig.ident2", assay.use = "SCT", verbose = FALSE)

  # Clamp requested PCs to what's actually available, and use the SAME pcs
  # for both UMAP and the neighbor graph (previously RunUMAP used opt$pcs
  # unclamped while FindNeighbors used a separately-clamped pcs_to_use).
  pcs_to_use <- min(opt$pcs, ncol(sub_obj) - 1)

  sub_obj <- RunUMAP(sub_obj, reduction = "harmony", dims = 1:pcs_to_use, seed.use = opt$seed, verbose = FALSE)

  # FindNeighbors MUST run before FindClusters — FindClusters consumes the
  # SNN graph FindNeighbors builds. The original order (FindClusters before
  # FindNeighbors) had no graph to cluster on.
  sub_obj <- FindNeighbors(sub_obj, reduction = "harmony", dims = 1:pcs_to_use, verbose = FALSE)

  # `resolutions` here is the function parameter (parsed to numeric by the
  # caller), not the raw opt$resolutions CLI string — FindClusters needs a
  # numeric vector, not "0.2,0.4,0.6,0.8,1.0,1.2" as literal text.
  sub_obj <- FindClusters(sub_obj, resolution = resolutions, random.seed = opt$seed, verbose = FALSE)

  return(sub_obj)
}

assign_metadata <- function(sub_obj) {
  sub_obj$Detailed_Condition <- case_when(
    grepl("HTY|UA", sub_obj$orig.ident2, ignore.case = TRUE) ~ "Healthy",
    grepl("AC", sub_obj$orig.ident2, ignore.case = TRUE) ~ "Acute",
    grepl("CH", sub_obj$orig.ident2, ignore.case = TRUE) ~ "Chronic",
    TRUE ~ "Unknown"
  )
  sub_obj$Condition <- ifelse(sub_obj$Detailed_Condition == "Healthy", "Healthy", "PMH")
  return(sub_obj)
}

process_subset <- function(seurat_obj, subset_clusters, prefix, out_base_dir,
                           resolutions = resolutions_vec, default_res = 0.2) {
  message(sprintf("\n=== Starting Pipeline for: %s ===", toupper(prefix)))
  
  out_dirs <- list(
    processed = file.path(out_base_dir, prefix, "processed"),
    plots     = file.path(out_base_dir, prefix, "plots"),
    tables    = file.path(out_base_dir, prefix, "tables"),
    dge       = file.path(out_base_dir, prefix, "dge")
  )
  lapply(out_dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
  
  if (!"global_cluster" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$global_cluster <- Idents(seurat_obj)
  }
  
  # Subset based on the matched dynamic labels
  sub_obj <- subset(seurat_obj, idents = subset_clusters)
  
  if (ncol(sub_obj) < 10) {
     stop("Too few cells found for this cell type. Aborting.")
  }

  sub_obj <- run_dim_reduction(sub_obj, resolutions)
  
  cluster_prefix <- ifelse(paste0("SCT_snn_res.", default_res) %in% colnames(sub_obj@meta.data), 
                           "SCT_snn_res.", "RNA_snn_res.")
  res_col <- paste0(cluster_prefix, default_res)
  Idents(sub_obj) <- res_col
  sub_obj$seurat_clusters <- sub_obj[[res_col]]
  
p_tree <- clustree(sub_obj, prefix = cluster_prefix) +
    ggtitle(paste(prefix, "- Clustree Resolution Tracker")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  save_plot(p_tree, file.path(out_dirs$plots, paste0(prefix, "_clustree")), w = 12, h = 10)

  valid_features <- intersect(PROOF_GENES, rownames(sub_obj))
  p_val <- DotPlot(sub_obj, features = valid_features) + RotatedAxis() + 
    ggtitle(paste(prefix, "Lineage Validation (Pre-Cleaning)"))
  save_plot(p_val, file.path(out_dirs$plots, paste0(prefix, "_validation_dotplot")), w = 10, h = 5)

  sub_obj <- assign_metadata(sub_obj)

  p_cond <- DimPlot(sub_obj, split.by = "Condition", label = TRUE, label.size = 5) +
    ggtitle(paste(prefix, "- Split by Condition (Res:", default_res, ")"))
  save_plot(p_cond, file.path(out_dirs$plots, paste0(prefix, "_umap_condition")), w = 12, h = 6)
  
  p_id2 <- DimPlot(sub_obj, split.by = "orig.ident2", label = TRUE, ncol = 3)
  save_plot(p_id2, file.path(out_dirs$plots, paste0(prefix, "_umap_origident2")), w = 15, h = 10)

  p_bar_stage <- create_proportion_barplot(sub_obj, "Detailed_Condition", paste(prefix, "Composition by Disease Stage"))
  save_plot(p_bar_stage, file.path(out_dirs$plots, paste0(prefix, "_stage_proportions")), w = 8, h = 7)

  p_bar_sample <- create_proportion_barplot(sub_obj, "orig.ident2", paste(tools::toTitleCase(prefix), "Composition by Sample"))
  save_plot(p_bar_sample, file.path(out_dirs$plots, paste0(prefix, "_sample_proportions")), w = 10, h = 7)

  DefaultAssay(sub_obj) <- "RNA"
  try({ sub_obj <- JoinLayers(sub_obj) }, silent = TRUE)
  sub_obj <- NormalizeData(sub_obj, verbose = FALSE)

  all_markers <- FindAllMarkers(sub_obj, assay = "RNA", only.pos = TRUE, 
                                min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)

  if (nrow(all_markers) > 0) {
    write.csv(all_markers, file.path(out_dirs$dge, paste0(prefix, "_markers.csv")), row.names = FALSE)
    top10 <- all_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    sub_obj <- ScaleData(sub_obj, features = top10$gene, verbose = FALSE)
    p_heat <- DoHeatmap(sub_obj, features = top10$gene, assay = "RNA") + NoLegend() + 
      ggtitle(paste(prefix, "Top 10 Markers (Res:", default_res, ")"))
    save_plot(p_heat, file.path(out_dirs$plots, paste0(prefix, "_heatmap_top10")), w = 15, h = 15)
  }

  cluster_compare <- table(sub_obj$global_cluster, Idents(sub_obj))
  write.csv(cluster_compare, file.path(out_dirs$tables, paste0(prefix, "_global_vs_new_clusters.csv")))
  saveRDS(sub_obj, file.path(out_dirs$processed, paste0(prefix, "_subset_processed.rds")))
  
  message(sprintf("=== Finished Pipeline for: %s ===\n", toupper(prefix)))
  return(sub_obj)
}

# ==============================================================================
# 3. DYNAMIC EXECUTION BLOCK
# ==============================================================================
message("Loading master annotated object: ", opt$input)
pmh_obj <- readRDS(opt$input)

# Apply QC filtering as in original
pmh_obj <- subset(pmh_obj, subset = orig.ident2 != "HTY244")

# --- FIX: Switch active identity to 'cell_type_full' to catch semicolon-separated ties ---
if ("cell_type_full" %in% colnames(pmh_obj@meta.data)) {
  Idents(pmh_obj) <- "cell_type_full"
} else if ("Cell_Type" %in% colnames(pmh_obj@meta.data)) {
  Idents(pmh_obj) <- "Cell_Type"
}

# Dynamically search the Idents levels for the target cell type (case-insensitive partial match)
available_idents <- levels(Idents(pmh_obj))
target_idents <- grep(opt$celltype, available_idents, ignore.case = TRUE, value = TRUE)

if (length(target_idents) == 0) {
  stop(sprintf("No clusters found matching '%s'.\nAvailable full clusters:\n%s", 
               opt$celltype, paste(available_idents, collapse="\n")))
}

message("Subsetting the following dynamically matched clusters:\n - ", paste(target_idents, collapse="\n - "))

# Run the pipeline module
subset_obj <- process_subset(
  seurat_obj = pmh_obj, 
  subset_clusters = target_idents, 
  prefix = opt$celltype, 
  out_base_dir = opt$outdir,
  resolutions = resolutions_vec,
  default_res = opt$resolution
)