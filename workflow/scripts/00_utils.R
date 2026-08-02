# workflow/scripts/utils.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(RColorBrewer)
})

# ==============================================================================
# 1. SHARED COLOR PALETTES
# ==============================================================================
mycolor <- c(
  "#CCCCCC", "#A6CEE3", "#FF7F00", "#09519C", "#FFFB33", "#556b2f",
  "#FF00FF", "#377EB8", "#bb8fce", "#666666", "#90EE90", "#ff4500",
  "#A6761D", "#E67E22", "#323695", "#E81E32", "#006837", "#CBC3E3",
  "#F1C40F", "#3498DB", "#34495E", "#FA9399", "#48C9B0", "#7D3C98",
  "#ff4500", "#8b4513", "#8a2be2", "#f0e68c", "#00ffff", "#32CD32", "#b03060"
)

doublet_color <- c("Doublet" = "#e35473", "Singlet" = "#54cdeb")

volcano_colors <- c("Upregulated in PMH" = "red", "Downregulated in PMH" = "blue", "Not Significant" = "grey80")

# Helper to dynamically pick colors for clusters
pick_colors <- function(n) {
  if (n <= length(mycolor)) mycolor[1:n] else colorRampPalette(mycolor)(n)
}

# ==============================================================================
# 2. UNIVERSAL GGPLOT THEMES
# ==============================================================================

# Standard pipeline theme for UMAPs, Barplots, and DotPlots
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

# Enlarged theme specifically for Monocle3 Trajectory and Pseudotime plots
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
# 3. UNIVERSAL PLOT SAVER
# ==============================================================================
# Saves both PDF and PNG formats identically
save_plot <- function(plot_obj, base_filepath, w = 12, h = 8) {
  for (fmt in c("pdf", "png")) {
    out_path <- paste0(base_filepath, ".", fmt)
    tryCatch({
      if (fmt == "pdf") {
        pdf(out_path, width = w, height = h)
      } else {
        png(out_path, width = w, height = h, units = "in", res = 300)
      }
      print(plot_obj)
      dev.off()
      message("Saved: ", out_path)
    }, error = function(e) {
      message("Warning: could not save ", out_path, ": ", conditionMessage(e))
      if (length(dev.list()) > 0) dev.off()
    })
  }
}

REQUIRED_COLS_AFTER_02 <- c("cell_type_short", "cell_type_full", "cluster_label")

resolve_res_col <- function(obj, resolution) {
  md <- colnames(obj@meta.data)
  for (p in c("SCT_snn_res.", "RNA_snn_res.")) {
    col <- paste0(p, resolution)
    if (col %in% md) return(col)
  }
  NULL
}

read_config <- function(path = "config/config.yaml") yaml::read_yaml(path)