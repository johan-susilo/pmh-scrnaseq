# ==============================================================================
# 00_utils.R — Shared helpers for the entire scRNA-seq pipeline (01-04)
#
# CONTRACT: every other script sources this file first, and reads its
# defaults from config.yaml via get_config()/cfg_get() rather than hardcoding
# them a second time in make_option(). This is the fix for the
# config-vs-CLI-default drift problem: config.yaml is the single source of
# truth, and CLI flags are optional overrides on top of it.
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(RColorBrewer)
  library(yaml)
})

# ==============================================================================
# 1. CONFIG LOADING (single source of truth for all default parameters)
# ==============================================================================

# Reads config.yaml once per R session and caches it, so repeated calls in a
# loop don't re-hit disk. Every script calls this the same way:
#   cfg <- get_config(opt$config)
#   default_res <- cfg_get(cfg, "preprocessing", "resolution", default = 0.2)
.cfg_cache <- new.env(parent = emptyenv())

get_config <- function(path = "config/config.yaml") {
  key <- normalizePath(path, mustWork = FALSE)
  if (!is.null(.cfg_cache[[key]])) return(.cfg_cache[[key]])
  if (!file.exists(path)) {
    warning("config.yaml not found at '", path, "'; falling back to hardcoded defaults")
    .cfg_cache[[key]] <- list()
    return(.cfg_cache[[key]])
  }
  cfg <- yaml::read_yaml(path)
  .cfg_cache[[key]] <- cfg
  cfg
}

# Safe nested lookup: cfg_get(cfg, "preprocessing", "resolution", default = 0.2)
# Never errors on a missing key -- just falls back to `default`, so config.yaml
# can be extended without every script needing a matching update.
cfg_get <- function(cfg, ..., default = NULL) {
  keys <- list(...)
  node <- cfg
  for (k in keys) {
    if (is.null(node) || !is.list(node) || is.null(node[[k]])) return(default)
    node <- node[[k]]
  }
  node
}

# ==============================================================================
# 2. SHARED COLOR PALETTES
# ==============================================================================
mycolor <- c(
  "#CCCCCC", "#A6CEE3", "#FF7F00", "#09519C", "#FFFB33", "#556b2f",
  "#FF00FF", "#377EB8", "#bb8fce", "#666666", "#90EE90", "#ff4500",
  "#A6761D", "#E67E22", "#323695", "#E81E32", "#006837", "#CBC3E3",
  "#F1C40F", "#3498DB", "#34495E", "#FA9399", "#48C9B0", "#7D3C98",
  "#ff4500", "#8b4513", "#8a2be2", "#f0e68c", "#00ffff", "#32CD32", "#b03060"
)

doublet_color  <- c("Doublet" = "#e35473", "Singlet" = "#54cdeb")
volcano_colors <- c("Upregulated in PMH" = "red", "Downregulated in PMH" = "blue",
                     "Not Significant" = "grey80")

pick_colors <- function(n) {
  if (n <= length(mycolor)) mycolor[1:n] else colorRampPalette(mycolor)(n)
}

# ==============================================================================
# 3. UNIVERSAL GGPLOT THEMES
# ==============================================================================
theme_pipeline <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold", size = base_size + 2),
      axis.text.x      = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y      = element_text(color = "black"),
      strip.text       = element_text(size = base_size - 2, face = "bold"),
      strip.background = element_rect(fill = "lightgray"),
      legend.title     = element_blank()
    )
}

theme_trajectory <- function(base_size = 18) {
  theme(
    text = element_text(size = base_size),
    axis.title = element_text(size = base_size + 2, face = "bold"),
    axis.text = element_text(size = base_size - 2),
    legend.title = element_text(size = base_size, face = "bold"),
    legend.text = element_text(size = base_size - 2),
    plot.title = element_text(size = base_size + 4, face = "bold", hjust = 0.5),
    strip.text = element_text(size = base_size + 4, face = "bold")
  )
}

# ==============================================================================
# 4. UNIVERSAL PLOT SAVER (PDF + PNG together)
# ==============================================================================
save_plot <- function(plot_obj, base_filepath, w = 12, h = 8) {
  for (fmt in c("pdf", "png")) {
    out_path <- paste0(base_filepath, ".", fmt)
    tryCatch({
      if (fmt == "pdf") pdf(out_path, width = w, height = h)
      else png(out_path, width = w, height = h, units = "in", res = 300)
      print(plot_obj)
      dev.off()
      message("Saved: ", out_path)
    }, error = function(e) {
      message("Warning: could not save ", out_path, ": ", conditionMessage(e))
      if (length(dev.list()) > 0) dev.off()
    })
  }
}

# ==============================================================================
# 5. CLUSTER / RESOLUTION-COLUMN RESOLUTION (single canonical implementation;
#    01, 02, and 03 all previously reinvented a version of this)
# ==============================================================================

# Given a Seurat object and a numeric/character resolution, find the matching
# "<assay>_snn_res.<r>" column, preferring SCT over RNA. Returns NULL (not an
# error) if nothing matches, so callers can decide how to fall back.
resolve_res_col <- function(obj, resolution) {
  md_cols <- colnames(obj@meta.data)
  for (prefix in c("SCT_snn_res.", "RNA_snn_res.")) {
    col <- paste0(prefix, resolution)
    if (col %in% md_cols) return(col)
  }
  NULL
}

# Strip a leading "X" that data.frame/table column-name mangling adds to
# purely-numeric cluster IDs (e.g. table() colnames "0","1" -> "X0","X1").
# Used consistently everywhere a cluster ID crosses a data.frame boundary,
# instead of every script writing its own gsub("^X", ...).
strip_cluster_prefix <- function(x) {
  trimws(gsub("^X+", "", as.character(x)))
}

# Call JoinLayers when available (Seurat v5), silently skip otherwise.
try_join_layers <- function(obj) {
  if (exists("JoinLayers", where = asNamespace("Seurat"), mode = "function")) {
    tryCatch(JoinLayers(obj), error = function(e) {
      message("Warning: JoinLayers failed: ", conditionMessage(e)); obj
    })
  } else {
    message("JoinLayers() not available; skipping")
    obj
  }
}

# ==============================================================================
# 6. CANONICAL OUTPUT PATH LAYOUT (same shape used by 01/02/03/04)
#
#    <01-out>/{processed,plots,tables,logs}/
#    <01-out>/TN.combined_dim30.rds, TN.combined_dim30_full.rds
#      (<01-out> is whatever -o the 01_preprocessing.R rule is given, e.g.
#       results/01_preprocessing -- used AS-IS, no extra nesting)
#
#    <02-out>/res_<r>/{singleR,markers,celliD,scCATCH,consensus,
#                       annotation_plots,combined_plots,logs}/
#    <02-out>/res_<r>/TN.combined_annotated.rds
#      (<02-out> likewise used as-is, e.g. results/02_annotation)
#
#    results/03_subsets/<celltype>/{processed,plots,tables,dge}/
#    results/03_subsets/<celltype>/<celltype>_subset_processed.rds       (03)
#    results/03_subsets/<celltype>/<celltype>_detailed_annotated.rds     (04)
#    results/03_subsets/<celltype>/{1_UMAPs,2_Validation_Plots,
#                                    3_Summary_Plots,4_Mucin_ECM,
#                                    5_Proportions}/                     (04)
#      -- 03 and 04 MUST be given the same base outdir (default
#         "results/03_subsets" for both) so 05_dge.R / 06_go.R can find both
#         files for a celltype in one shared <celltype>/processed/ folder.
#
#    results/03_subsets/<celltype>/dge_pseudobulk/                       (05)
#    results/03_subsets/<celltype>/pathways/                             (06)
# ==============================================================================

make_stage_dirs <- function(base, subdirs) {
  dirs <- lapply(subdirs, function(s) file.path(base, s))
  names(dirs) <- subdirs
  lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
  dirs
}
