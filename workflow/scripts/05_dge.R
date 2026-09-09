#!/usr/bin/env Rscript
# Unified Pseudobulk DGE Factory
#
# Runs DESeq2 (macro "all combined" + per-subcluster) for every celltype
# folder found under --indir. Fully generic: derives celltype names from
# the directory structure produced by 03_subset_clusters.R / 04_detail_annotation.R,
# never assumes fibroblast/macrophage specifically.
#
# CONTRACT WITH 03/04: --indir must be the SAME base directory both of those
# scripts wrote to (default "results/03_subsets" for all three) -- this
# script looks for <indir>/<celltype>/processed/<celltype>_detailed_annotated.rds
# first, falling back to <celltype>_subset_processed.rds if 04 hasn't run yet.
#
# Usage:
#   Rscript 05_dge.R --indir results/03_subsets
#   Rscript 05_dge.R --indir results/03_subsets --celltypes fibroblast,macrophage

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(DESeq2)
  library(tidyr)
  library(ggrepel)
  library(stringr)
})

source("workflow/scripts/00_utils.R")

# ==============================================================================
# 1. COMMAND-LINE INTERFACE
# ==============================================================================
option_list <- list(
  make_option(c("-i", "--indir"), type = "character", default = "results/03_subsets",
              help = "Base directory containing one folder per celltype (SAME dir as 03/04 -o) [default: results/03_subsets]"),
  make_option(c("-c", "--celltypes"), type = "character", default = NULL,
              help = "Optional comma-separated whitelist of celltype folder names to process (default: all subfolders of --indir)"),
  make_option(c("--config"), type = "character", default = "config/config.yaml", help = "Path to config.yaml"),
  make_option(c("--compare_by"), type = "character", default = NULL, help = "Metadata column to compare on [config: dge.compare_by]"),
  make_option(c("--ident_1"), type = "character", default = NULL, help = "Primary group [config: dge.ident_1]"),
  make_option(c("--ident_2"), type = "character", default = NULL, help = "Control group [config: dge.ident_2]"),
  make_option(c("--seed"), type = "integer", default = NULL, help = "Global random seed [config: reproducibility.random_seed]")
)
opt <- parse_args(OptionParser(option_list = option_list))
cfg <- get_config(opt$config)

`%||%` <- function(a, b) if (is.null(a)) b else a
opt$compare_by <- opt$compare_by %||% cfg_get(cfg, "dge", "compare_by", default = "Condition")
opt$ident_1    <- opt$ident_1    %||% cfg_get(cfg, "dge", "ident_1",    default = "PMH")
opt$ident_2    <- opt$ident_2    %||% cfg_get(cfg, "dge", "ident_2",    default = "Healthy")
opt$seed       <- opt$seed       %||% cfg_get(cfg, "reproducibility", "random_seed", default = 42)
set.seed(opt$seed)

if (is.null(opt$indir)) {
  stop("Missing required argument: --indir (base 03_subsets directory)")
}
if (!dir.exists(opt$indir)) {
  stop(sprintf("--indir does not exist: %s", opt$indir))
}

message("\n==================================================================")
message("=== Starting Unified Pseudobulk DGE Factory ===")
message("=== Comparing '", opt$compare_by, "': ", opt$ident_1, " vs ", opt$ident_2, " ===")
message("==================================================================")

base_dir <- opt$indir

all_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)
# never treat marker/log/scratch folders as celltypes
all_dirs <- all_dirs[!all_dirs %in% c("logs") & !startsWith(all_dirs, "_")]

if (!is.null(opt$celltypes)) {
  requested <- trimws(strsplit(opt$celltypes, ",")[[1]])
  missing <- setdiff(requested, all_dirs)
  if (length(missing) > 0) {
    warning(sprintf("Requested celltype(s) not found under %s and will be skipped: %s",
                     base_dir, paste(missing, collapse = ", ")))
  }
  cell_types <- intersect(requested, all_dirs)
} else {
  cell_types <- all_dirs
}

if (length(cell_types) == 0) {
  stop("No celltype folders found to process. Check --indir / --celltypes.")
}

message("Celltypes to process: ", paste(cell_types, collapse = ", "))

# ------------------------------------------------------------------------------
# Helper Function: Run DESeq2, Save CSV, and Plot Volcano
# ------------------------------------------------------------------------------
run_and_save_deseq2 <- function(counts_matrix, meta, comparison_name, out_dir, title_prefix,
                                 ident_1 = opt$ident_1, ident_2 = opt$ident_2) {

  valid_idx <- !is.na(meta$condition)
  meta <- meta[valid_idx, , drop = FALSE]
  counts_matrix <- counts_matrix[, valid_idx, drop = FALSE]

  ident1_count <- sum(meta$condition == ident_1)
  ident2_count <- sum(meta$condition == ident_2)

  ident1_count <- sum(meta$condition == ident_1)
  ident2_count <- sum(meta$condition == ident_2)

  if (ident1_count < 2 | ident2_count < 2) {
    message(sprintf("   -> SKIPPING %s: Insufficient replicates (%s: %d, %s: %d)",
                    comparison_name, ident_1, ident1_count, ident_2, ident2_count))
    return(NULL)
  }

  message(sprintf("   -> Running %s (%s: %d, %s: %d)...", comparison_name, ident_1, ident1_count, ident_2, ident2_count))

  dds <- DESeqDataSetFromMatrix(countData = counts_matrix, colData = meta, design = ~ condition)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", ident_1, ident_2))

  up_label   <- paste0("Upregulated in ", ident_1)
  down_label <- paste0("Downregulated in ", ident_1)

  res_df <- as.data.frame(res) %>%
    filter(!is.na(padj)) %>%
    mutate(gene = rownames(.)) %>%
    mutate(
      significance = case_when(
        padj < 0.05 & log2FoldChange > 1 ~ up_label,
        padj < 0.05 & log2FoldChange < -1 ~ down_label,
        TRUE ~ "Not Significant"
      )
    ) %>%
    dplyr::select(gene, everything()) %>%
    arrange(padj)

  # Save CSV
  write.csv(res_df, file.path(out_dir, paste0("deseq2_", comparison_name, ".csv")), row.names = FALSE)

  # Save Volcano Plot
  top_genes <- res_df %>% filter(significance != "Not Significant") %>%
    group_by(significance) %>% slice_head(n = 20) %>% ungroup()

  volcano_pal <- setNames(c("red", "blue", "grey80"), c(up_label, down_label, "Not Significant"))
  p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_text_repel(data = top_genes, aes(label = gene), color = "black", box.padding = 0.5, max.overlaps = Inf) +
    scale_color_manual(values = volcano_pal) +
    theme_minimal() +
    labs(title = paste(title_prefix, "(", ident_1, "vs", ident_2, ")"))

  ggsave(filename = file.path(out_dir, paste0("volcano_", comparison_name, ".png")), plot = p_volcano, width = 7, height = 6)
}

# ------------------------------------------------------------------------------
# Main Execution Loop
# ------------------------------------------------------------------------------
for (cell_type in cell_types) {

  message(paste("\n========================================================"))
  message(paste("=== Processing Cell Type:", toupper(cell_type), "==="))
  message(paste("========================================================"))

  annotated_path <- file.path(base_dir, cell_type, "00_data", paste0(cell_type, "_detailed_annotated.rds"))
  processed_path <- file.path(base_dir, cell_type, "00_data", paste0(cell_type, "_subset_processed.rds"))

  if (file.exists(annotated_path)) {
    message(paste("   -> Loading annotated RDS for", cell_type))
    sub_obj <- readRDS(annotated_path)
    cluster_col <- "Detailed_Label"
  } else if (file.exists(processed_path)) {
    message(paste("   -> Loading basic processed RDS for", cell_type))
    sub_obj <- readRDS(processed_path)
    cluster_col <- "seurat_clusters"
  } else {
    message(paste("WARNING: Could not find any processed RDS for", cell_type, "- Skipping."))
    next
  }

  if (!opt$compare_by %in% colnames(sub_obj@meta.data)) {
    message(paste("WARNING:", cell_type, "object is missing '", opt$compare_by, "' metadata - Skipping."))
    next
  }

  sub_obj@meta.data[[opt$compare_by]] <- factor(sub_obj@meta.data[[opt$compare_by]],
                                                 levels = c(opt$ident_2, opt$ident_1))

  out_dir <- file.path(base_dir, cell_type, "03_dge_pseudobulk")
  if (!dir.exists(out_dir)) { dir.create(out_dir, recursive = TRUE) }


  # ==============================================================================
  # PART A: "All Combined" DGE (The Macro View)
  # ==============================================================================
  message("\n--- Aggregating ALL cells (Macro View) ---")
  pb_all <- AggregateExpression(sub_obj, assays = "RNA", slot = "counts",
                                group.by = c("orig.ident2", opt$compare_by), return.seurat = FALSE)$RNA

  ident_regex <- paste0("(", opt$ident_1, "|", opt$ident_2, ")$")
  meta_all <- data.frame(pseudobulk_id = colnames(pb_all)) %>%
    mutate(
      condition = str_extract(pseudobulk_id, ident_regex),
      sample_id = str_remove(pseudobulk_id, paste0("_?", ident_regex))
    )
  rownames(meta_all) <- meta_all$pseudobulk_id
  meta_all$condition <- factor(meta_all$condition, levels = c(opt$ident_2, opt$ident_1))

  run_and_save_deseq2(
    counts_matrix = pb_all,
    meta = meta_all,
    comparison_name = paste0("ALL_COMBINED_", cell_type),
    out_dir = out_dir,
    title_prefix = paste("All", tools::toTitleCase(cell_type), "Combined")
  )


  # ==============================================================================
  # PART B: Sub-cluster DGE (The Micro View)
  # ==============================================================================
  message("\n--- Aggregating by Sub-cluster (Micro View) ---")
  pb_sub <- AggregateExpression(sub_obj, assays = "RNA", slot = "counts",
                                group.by = c(cluster_col, "orig.ident2", opt$compare_by), return.seurat = FALSE)$RNA

  # Safely extract the cluster levels whether it's Detailed_Label or seurat_clusters
  clusters <- levels(as.factor(sub_obj@meta.data[[cluster_col]]))
  all_sub_cols <- colnames(pb_sub)

  for (cluster_id in clusters) {

    seurat_safe_id <- str_replace_all(cluster_id, "_", "-")

    # Isolate columns for this specific cluster using the dashed ID
    cluster_cols <- all_sub_cols[grepl(paste0("^g?", seurat_safe_id, "_"), all_sub_cols)]

    if (length(cluster_cols) == 0) { next }

    counts_sub <- pb_sub[, cluster_cols, drop = FALSE]

    meta_sub <- data.frame(pseudobulk_id = cluster_cols) %>%
      mutate(
        condition = str_extract(pseudobulk_id, ident_regex),
        # Use the dashed ID to correctly extract the sample name
        sample_id = str_remove(pseudobulk_id, paste0("_?", ident_regex)) %>%
                    str_remove(paste0("^g?", seurat_safe_id, "_"))
      )
    rownames(meta_sub) <- meta_sub$pseudobulk_id
    meta_sub$condition <- factor(meta_sub$condition, levels = c(opt$ident_2, opt$ident_1))

    run_and_save_deseq2(
      counts_matrix = counts_sub,
      meta = meta_sub,
      comparison_name = paste0(cell_type, "_subcluster_", cluster_id), # Keep the original name for the file!
      out_dir = out_dir,
      title_prefix = paste(tools::toTitleCase(cell_type), "Cluster", cluster_id)
    )
  }

  # Clean up memory before loading the next cell type
  rm(sub_obj, pb_all, pb_sub); gc()
}

message("\n==================================================================")
message("=== Unified DGE Factory Complete! ===")
message("==================================================================")