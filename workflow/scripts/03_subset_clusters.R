#!/usr/bin/env Rscript
# Usage:
#   Rscript 03_subset_clusters.R -i results/02_annotation/res_0.2/TN.combined_annotated.rds \
#     -c fibroblast -o results/03_subsets
#
# CONTRACT WITH 02: reads `cell_type_full` (exact strings produced by
# 02_global_annotation.R's apply_labels step) as the primary key for
# subsetting. This replaces the old bare grep() substring match, which could
# silently pull in unrelated clusters whose label happened to contain the
# search string, or silently miss clusters labeled "Low-confidence: X" /
# "X (mixed)". See match_target_clusters() below.
#
# CONTRACT WITH 04: writes
#   <outdir>/<celltype>/processed/<celltype>_subset_processed.rds
# with `Condition` metadata set, matching what 04_detail_annotation.R expects.

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
  make_option(c("-i", "--input"),      type = "character", help = "Path to TN.combined_annotated.rds (from 02_global_annotation.R apply_labels step)"),
  make_option(c("-c", "--celltype"),   type = "character", help = "Target cell type to subset (e.g. 'Fibroblasts'); matched against cell_type_full"),
  make_option(c("-o", "--outdir"),     type = "character", default = "results/03_subsets", help = "Base output directory for subset clusters"),
  make_option(c("--config"),           type = "character", default = "config/config.yaml", help = "Path to config.yaml"),
  make_option(c("-r", "--resolution"), type = "numeric",   default = NULL, help = "Default clustering resolution for the subset [config: subsetting.pcs / resolution]"),
  make_option(c("--pcs"),              type = "integer",   default = NULL, help = "Number of PCs for subsetting [config: subsetting.pcs]"),
  make_option(c("--resolutions"),      type = "character", default = NULL, help = "Comma-separated resolutions to sweep [config: subsetting.resolutions]"),
  make_option(c("--exclude_samples"),  type = "character", default = "", help = "Comma-separated orig.ident2 values to drop before subsetting (e.g. known-bad samples)"),
  make_option(c("--seed"),             type = "integer",   default = NULL, help = "Global random seed [config: reproducibility.random_seed]")
)
opt <- parse_args(OptionParser(option_list = option_list))
cfg <- get_config(opt$config)

`%||%` <- function(a, b) if (is.null(a)) b else a
opt$resolution <- opt$resolution %||% cfg_get(cfg, "preprocessing", "resolution", default = 0.2)
opt$pcs        <- opt$pcs        %||% cfg_get(cfg, "subsetting", "pcs", default = 20)
opt$seed       <- opt$seed       %||% cfg_get(cfg, "reproducibility", "random_seed", default = 42)
if (is.null(opt$resolutions)) {
  res_cfg <- cfg_get(cfg, "subsetting", "resolutions", default = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2))
  opt$resolutions <- paste(res_cfg, collapse = ",")
}

if (is.null(opt$input) || is.null(opt$celltype) || is.null(opt$outdir)) {
  stop("Missing required arguments: --input, --celltype, or --outdir")
}

set.seed(opt$seed)

PROOF_GENES <- c("LUM", "DCN", "COL1A1", "KRT14", "KRT1", "PTPRC", "CD68", "PECAM1", "VWF")

resolutions_vec <- as.numeric(trimws(strsplit(opt$resolutions, ",")[[1]]))
if (any(is.na(resolutions_vec))) {
  stop("Could not parse --resolutions into numeric values: ", opt$resolutions)
}
if (!(opt$resolution %in% resolutions_vec)) resolutions_vec <- sort(c(resolutions_vec, opt$resolution))

# ==============================================================================
# 2. CLUSTER MATCHING — the fixed contract described at the top of this file
# ==============================================================================

# Ambiguous/low-confidence consensus labels that should never be silently
# swept into a subset via substring matching (see 02_global_annotation.R).
AMBIGUOUS_LABEL_PATTERN <- "Ambiguous|\\(mixed\\)|^Low-confidence:|^Unknown$"

# Decide which cluster labels in `available_idents` belong to the requested
# celltype. Strategy, in order:
#   1. Exact match against `requested` (case-sensitive) — always preferred.
#   2. Case-insensitive exact match.
#   3. Substring match, but ONLY among labels that are NOT ambiguous/
#      low-confidence, and every match is logged explicitly so it's auditable
#      which labels were pulled in and why.
# Ambiguous-pattern labels are never silently substring-matched; if the user
# wants those included they must pass the exact label string.
match_target_clusters <- function(requested, available_idents) {
  exact <- available_idents[available_idents == requested]
  if (length(exact) > 0) {
    message("Exact match for '", requested, "': ", paste(exact, collapse = ", "))
    return(exact)
  }

  exact_ci <- available_idents[tolower(available_idents) == tolower(requested)]
  if (length(exact_ci) > 0) {
    message("Case-insensitive exact match for '", requested, "': ", paste(exact_ci, collapse = ", "))
    return(exact_ci)
  }

  safe_pool <- available_idents[!grepl(AMBIGUOUS_LABEL_PATTERN, available_idents)]
  candidates <- grep(requested, safe_pool, ignore.case = TRUE, value = TRUE)

  if (length(candidates) == 0) {
    ambiguous_hits <- grep(requested, available_idents, ignore.case = TRUE, value = TRUE)
    msg <- sprintf(
      "No confident cluster labels matched '%s'.\nAll available labels:\n%s",
      requested, paste(available_idents, collapse = "\n")
    )
    if (length(ambiguous_hits) > 0) {
      msg <- paste0(msg, "\n\nNOTE: these labels contain '", requested,
                     "' but were excluded as ambiguous/low-confidence:\n  - ",
                     paste(ambiguous_hits, collapse = "\n  - "),
                     "\nPass the exact label string via --celltype if you want to include them.")
    }
    stop(msg)
  }

  message("Substring-matched labels for '", requested, "' (excluding ambiguous/low-confidence labels): ",
          paste(candidates, collapse = ", "))
  candidates
}

# ==============================================================================
# 3. HELPER FUNCTIONS & CORE PROCESSING
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
  # Purge global embeddings to prevent cross-contamination
  sub_obj@reductions <- list()
  DefaultAssay(sub_obj) <- "RNA"
  sub_obj <- SCTransform(sub_obj, method = "glmGamPoi", vst.flavor = "v2", verbose = FALSE)
  sub_obj <- RunPCA(sub_obj, assay = "SCT", seed.use = opt$seed, verbose = FALSE)
  sub_obj <- RunHarmony(sub_obj, group.by.vars = "orig.ident2", assay.use = "SCT", verbose = FALSE)

  # Clamp requested PCs to what's actually available; use the SAME pcs for
  # both UMAP and the neighbor graph.
  pcs_to_use <- min(opt$pcs, ncol(sub_obj) - 1)
  if (pcs_to_use < opt$pcs)
    message("Requested pcs=", opt$pcs, " exceeds ncol-1=", ncol(sub_obj) - 1, "; clamping to ", pcs_to_use)

  sub_obj <- RunUMAP(sub_obj, reduction = "harmony", dims = 1:pcs_to_use, seed.use = opt$seed, verbose = FALSE)
  sub_obj <- FindNeighbors(sub_obj, reduction = "harmony", dims = 1:pcs_to_use, verbose = FALSE)
  sub_obj <- FindClusters(sub_obj, resolution = resolutions, random.seed = opt$seed, verbose = FALSE)

  sub_obj
}

process_subset <- function(seurat_obj, subset_clusters, prefix, out_base_dir,
                           resolutions, default_res) {
  message(sprintf("\n=== Starting Pipeline for: %s ===", toupper(prefix)))

  out_dirs <- make_stage_dirs(file.path(out_base_dir, prefix),
                               c("00_data", "01_subset_QC/plots", "01_subset_QC/tables", "01_subset_QC/markers"))
  names(out_dirs) <- c("processed", "plots", "tables", "dge")

  if (!"global_cluster" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$global_cluster <- Idents(seurat_obj)
  }

  sub_obj <- subset(seurat_obj, idents = subset_clusters)
  message("Subsetted ", ncol(sub_obj), " cells from clusters: ", paste(subset_clusters, collapse = ", "))

  if (ncol(sub_obj) < 10) {
     stop("Too few cells found for this cell type (", ncol(sub_obj), " < 10). Aborting.")
  }

  sub_obj <- run_dim_reduction(sub_obj, resolutions)

  cluster_prefix <- ifelse(paste0("SCT_snn_res.", default_res) %in% colnames(sub_obj@meta.data),
                           "SCT_snn_res.", "RNA_snn_res.")
  res_col <- paste0(cluster_prefix, default_res)
  if (!res_col %in% colnames(sub_obj@meta.data))
    stop("Resolution column '", res_col, "' not produced by FindClusters. Requested resolutions: ",
         paste(resolutions, collapse = ", "))
  Idents(sub_obj) <- res_col
  sub_obj$seurat_clusters <- sub_obj[[res_col]]

  p_tree <- clustree(sub_obj, prefix = cluster_prefix) +
    ggtitle(paste(prefix, "- Clustree Resolution Tracker")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  save_plot(p_tree, file.path(out_dirs$plots, paste0(prefix, "_clustree")), w = 12, h = 10)

  valid_features <- intersect(PROOF_GENES, rownames(sub_obj))
  if (length(valid_features) > 0) {
    p_val <- DotPlot(sub_obj, features = valid_features) + RotatedAxis() +
      ggtitle(paste(prefix, "Lineage Validation (Pre-Cleaning)"))
    save_plot(p_val, file.path(out_dirs$plots, paste0(prefix, "_validation_dotplot")), w = 10, h = 5)
  } else {
    message("No PROOF_GENES found in subset; skipping validation dotplot")
  }

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
  sub_obj <- try_join_layers(sub_obj)
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
  } else {
    message("No positive markers found for subset; skipping heatmap")
  }

  cluster_compare <- table(sub_obj$global_cluster, Idents(sub_obj))
  write.csv(cluster_compare, file.path(out_dirs$tables, paste0(prefix, "_global_vs_new_clusters.csv")))

  sub_obj@misc$pipeline_seed        <- opt$seed
  sub_obj@misc$subset_source_labels <- subset_clusters
  saveRDS(sub_obj, file.path(out_dirs$processed, paste0(prefix, "_subset_processed.rds")))

  message(sprintf("=== Finished Pipeline for: %s ===\n", toupper(prefix)))
  sub_obj
}

# ==============================================================================
# 4. DYNAMIC EXECUTION BLOCK
# ==============================================================================
message("Loading annotated object: ", opt$input)
if (!file.exists(opt$input))
  stop("Input not found at '", opt$input, "'. Run 02_global_annotation.R's 'apply_labels' step first.")
pmh_obj <- readRDS(opt$input)

exclude_samples <- trimws(strsplit(opt$exclude_samples, ",")[[1]])
exclude_samples <- exclude_samples[exclude_samples != ""]
if (length(exclude_samples) > 0 && "orig.ident2" %in% colnames(pmh_obj@meta.data)) {
  message("Excluding samples: ", paste(exclude_samples, collapse = ", "))
  pmh_obj <- subset(pmh_obj, subset = !(orig.ident2 %in% exclude_samples))
}

# Use cell_type_full (the exact, canonical label from 02) as the identity to
# subset on. The old "Cell_Type" fallback is removed: 02_global_annotation.R
# never produces a column by that name, so that branch was dead code.
if ("cell_type_full" %in% colnames(pmh_obj@meta.data)) {
  Idents(pmh_obj) <- "cell_type_full"
} else {
  stop("Input object is missing 'cell_type_full'. Re-run 02_global_annotation.R's ",
       "'apply_labels' step, which is the only step that creates this column.")
}

available_idents <- levels(Idents(pmh_obj))
search_term <- gsub("_", " ", opt$celltype)
target_idents <- match_target_clusters(search_term, available_idents)

message("Subsetting the following matched cluster label(s):\n - ", paste(target_idents, collapse = "\n - "))

subset_obj <- process_subset(
  seurat_obj      = pmh_obj,
  subset_clusters = target_idents,
  prefix          = opt$celltype,
  out_base_dir    = opt$outdir,
  resolutions     = resolutions_vec,
  default_res     = opt$resolution
)

message("Step complete at ", Sys.time())
