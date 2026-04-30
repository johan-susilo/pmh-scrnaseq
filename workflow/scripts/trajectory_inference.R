suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(monocle3)
  library(viridis)
})

# ==============================================================================
# 1. COMMAND LINE ARGUMENTS
# ==============================================================================
option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Path to the annotated RDS object"),
  make_option(c("-p", "--prefix"), type="character", help="Cell type prefix (e.g., 'fibroblast')"),
  make_option(c("-r", "--root_state"), type="character", default=NULL, help="Exact name of the cluster/annotation to start the trajectory (e.g., 'F2_Universal_Reticular')"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory for plots and the updated RDS")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (!is.null(opt$root_state)) {
  if (opt$root_state == "" || toupper(opt$root_state) == "NULL" || toupper(opt$root_state) == "NONE") {
    opt$root_state <- NULL
  }
}

message("\n==================================================================")
message(paste("=== Starting Monocle 3 Trajectory for:", toupper(opt$prefix), "==="))
message("==================================================================")

# Setup Output Directories
out_plots <- file.path(opt$outdir, "plots")
out_processed <- file.path(opt$outdir, "processed")
dir.create(out_plots, recursive = TRUE, showWarnings = FALSE)
dir.create(out_processed, recursive = TRUE, showWarnings = FALSE)

safe_save_plot <- function(plot_obj, base_filename, w = 8, h = 6) {
  pdf_path <- file.path(out_plots, paste0(base_filename, ".pdf"))
  png_path <- file.path(out_plots, paste0(base_filename, ".png"))
  tryCatch({ pdf(pdf_path, width = w, height = h); print(plot_obj); dev.off() }, error = function(e) {})
  tryCatch({ png(png_path, width = w, height = h, units = "in", res = 300); print(plot_obj); dev.off() }, error = function(e) {})
}

# ==============================================================================
# 2. LOAD & CONVERT SEURAT TO MONOCLE 3 CDS
# ==============================================================================
sub_obj <- readRDS(opt$input)

if (!"Detailed_Condition" %in% colnames(sub_obj@meta.data)) {
  stop("ERROR: 'Detailed_Condition' missing. Ensure the object was cleaned properly.")
}

message("Converting Seurat object to Monocle 3 CDS...")
expression_matrix <- GetAssayData(sub_obj, assay = "RNA", slot = "counts")
cell_metadata <- sub_obj@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix), row.names = rownames(expression_matrix))

cds <- new_cell_data_set(expression_matrix, cell_metadata = cell_metadata, gene_metadata = gene_annotation)

# CRITICAL: Transfer Harmony/Seurat UMAP coordinates directly into Monocle
reducedDims(cds)[["UMAP"]] <- sub_obj@reductions$umap@cell.embeddings

# ==============================================================================
# 3. GRAPH LEARNING & ROOT SELECTION
# ==============================================================================
message("Clustering cells and learning trajectory graph...")
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = TRUE, close_loop = FALSE)

# Determine the Root State
if (!is.null(opt$root_state)) {
  message(paste("Manual Root State requested:", opt$root_state))
  
  # Check SingleR_label first, then fall back to seurat_clusters
  if ("SingleR_label" %in% colnames(sub_obj@meta.data) && opt$root_state %in% sub_obj$SingleR_label) {
    root_barcodes <- rownames(sub_obj@meta.data[sub_obj$SingleR_label == opt$root_state, ])
  } else if (opt$root_state %in% sub_obj$seurat_clusters) {
    root_barcodes <- rownames(sub_obj@meta.data[sub_obj$seurat_clusters == opt$root_state, ])
  } else {
    stop(paste("ERROR: Root state", opt$root_state, "not found in SingleR_label or seurat_clusters!"))
  }
} else {
  message("Auto-detecting baseline root cluster...")
  prop_table <- prop.table(table(sub_obj$seurat_clusters, sub_obj$Detailed_Condition), margin = 1)
  
  if ("Healthy" %in% colnames(prop_table)) {
    root_cluster <- names(which.max(prop_table[, "Healthy"]))
    message(paste("Auto-detected Root:", root_cluster, "(Highest proportion of Healthy cells)"))
  } else {
    root_cluster <- levels(Idents(sub_obj))[1]
    message(paste("Warning: No 'Healthy' condition found. Defaulting to cluster:", root_cluster))
  }
  root_barcodes <- rownames(sub_obj@meta.data[sub_obj$seurat_clusters == root_cluster, ])
}

# ==============================================================================
# 4. CALCULATE PSEUDOTIME
# ==============================================================================
message("Ordering cells by Pseudotime...")
cds <- order_cells(cds, root_cells = root_barcodes)
sub_obj$Pseudotime <- pseudotime(cds)

# ==============================================================================
# 5. VISUALIZATIONS
# ==============================================================================
message("Generating trajectory plots...")

# Plot A: Monocle Graph
p_monocle <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 1.5) +
  ggtitle(paste(toupper(opt$prefix), "- Monocle 3 Trajectory Graph")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
safe_save_plot(p_monocle, paste0(opt$prefix, "_Monocle_Trajectory_Graph"))

# Plot B: Pseudotime FeaturePlot
p_traj <- FeaturePlot(sub_obj, features = "Pseudotime", pt.size = 0.6) +
  scale_color_viridis(option = "C", na.value = "grey90") +
  ggtitle(paste(toupper(opt$prefix), "- Pseudotime Heatmap")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
safe_save_plot(p_traj, paste0(opt$prefix, "_UMAP_Pseudotime_Heatmap"))

# Plot C: Density Progression
meta_df <- sub_obj@meta.data %>% filter(is.finite(Pseudotime))
meta_df$Detailed_Condition <- factor(meta_df$Detailed_Condition, levels = c("Healthy", "Acute", "Chronic"))

p_density <- ggplot(meta_df, aes(x = Pseudotime, fill = Detailed_Condition, color = Detailed_Condition)) +
  geom_density(alpha = 0.4, linewidth = 1) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("Healthy" = "#4daf4a", "Acute" = "#377eb8", "Chronic" = "#e41a1c")) +
  scale_color_manual(values = c("Healthy" = "#4daf4a", "Acute" = "#377eb8", "Chronic" = "#e41a1c")) +
  labs(
    title = paste(toupper(opt$prefix), "- Disease Progression Map"),
    subtitle = paste("Root State:", ifelse(is.null(opt$root_state), "Auto-Detected", opt$root_state)),
    x = "Pseudotime (0 = Origin, High = Terminal State)", y = "Density of Cells"
  )
safe_save_plot(p_density, paste0(opt$prefix, "_Density_Pseudotime_by_Stage"))

# ==============================================================================
# 6. SAVE OUTPUT
# ==============================================================================
message("Saving updated RDS with Pseudotime metadata...")
saveRDS(sub_obj, file.path(out_processed, paste0(opt$prefix, "_annotated_with_pseudotime.rds")))
message("=== Trajectory Analysis Complete ===")