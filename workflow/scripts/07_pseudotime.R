#!/usr/bin/env Rscript
# Pseudotime Analysis using Monocle3
# Usage: Rscript ./pseudotime.R -r TN.combined_dim30.rds -c 0,2,3,6,15 -o output/pseudotime
# Example: Rscript ./pseudotime.R -r /home/johan/output/skin_pmh/TN.combined_dim30.rds -c 0,1,2,3 -o /home/johan/output/skin_pmh/pseudotime
# Example: Rscript ./pseudotime.R -r TN.combined_dim30.rds --all_clusters -o output/pseudotime

#Rscript /home/johan/pipeline/scRNA/skin/pseudotime.R  -r "/home/johan/output/skin_pmh_harmony_sctransform2/subset_cluster/fibroblast/processed/fibroblast_detailed_annotated.rds"  --all_clusters  -o "/home/johan/output/skin_pmh_harmony_sctransform2/subset_cluster/fibroblast/pseudotime_global"

Sys.time()

suppressPackageStartupMessages({
  library(optparse)
  library(monocle3)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# Command-line Interface ----------------------------------------------------
option_list <- list(
  make_option(c("-r", "--rds"), type = "character", default = NULL,
              help = "Path to RDS file (TN.combined_dim30.rds)"),
  make_option(c("-c", "--clusters"), type = "character", default = NULL,
              help = "Comma-separated list of cluster IDs to include (e.g., '0,2,3,6,15')"),
  make_option(c("--all_clusters"), action = "store_true", default = FALSE,
              help = "Use all clusters for pseudotime analysis"),
  make_option(c("-o", "--output"), type = "character", default = "output/pseudotime",
              help = "Output directory [default: output/pseudotime]"),
  make_option(c("--method"), type = "character", default = "monocle3",
              help = "Trajectory inference method: monocle3 [default: monocle3]")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

# Create output directory
output_dir <- opt$output
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Set up logging
log_file <- file.path(output_dir, paste0("pseudotime_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_conn <- file(log_file, open = "wt")
sink(log_conn, type = "output", split = TRUE)
sink(log_conn, type = "message")

message("============================================================")
message("Pseudotime Analysis Pipeline")
message("============================================================")
message("Log file: ", log_file)
message("Started at: ", Sys.time())

# Helper function to safely save PDF plots
safe_save_pdf <- function(plot_obj, filepath, w = 15, h = 15) {
  tryCatch({
    pdf(filepath, width = w, height = h)
    on.exit(dev.off(), add = TRUE)
    print(plot_obj)
    message("Saved plot: ", filepath)
  }, error = function(e) {
    message("Warning: Failed to save plot ", filepath, ": ", conditionMessage(e))
    if (length(dev.list()) > 0) dev.off()
  })
}

# Load and subset data ------------------------------------------------------
message("\n============================================================")
message("Loading Seurat object")
message("============================================================")

if (is.null(opt$rds)) {
  stop("RDS file path must be specified with --rds")
}

TN.combined <- readRDS(file = opt$rds)
message("Loaded Seurat object: ", ncol(TN.combined), " cells, ", nrow(TN.combined), " features")
message("Total clusters: ", length(unique(Idents(TN.combined))))

# Determine which clusters to use
if (opt$all_clusters) {
  message("Using ALL clusters for pseudotime analysis")
  
  # Ensure we use detailed labels
  Idents(TN.combined) <- "Detailed_Label"
  
  TN.subset <- TN.combined
  selected_clusters <- as.character(unique(Idents(TN.combined)))
  
} else if (!is.null(opt$clusters)) {
  # Parse cluster IDs from comma-separated string (NO as.numeric)
  cluster_ids <- trimws(strsplit(opt$clusters, ",")[[1]])
  message("Subsetting to clusters: ", paste(cluster_ids, collapse = ", "))

  # Ensure we are using your detailed dictionary labels
  Idents(TN.combined) <- "Detailed_Label"
  
  TN.subset <- subset(TN.combined, idents = cluster_ids)
  selected_clusters <- cluster_ids

  if (ncol(TN.subset) == 0) {
    stop("No cells found in the specified clusters. Check cluster IDs.")
  }
  message("Subset contains ", ncol(TN.subset), " cells")
} else {
  stop("Must specify either --clusters or --all_clusters")
}


TN.subset$clusters <- as.character(Idents(TN.subset))

# Plot the subset
message("\nGenerating UMAP plot of selected clusters...")
p_umap <- DimPlot(TN.subset, reduction = "umap", label = TRUE, pt.size = 0.8) +
  ggtitle("Selected Clusters for Pseudotime")
safe_save_pdf(p_umap, file.path(output_dir, "subset_umap.pdf"))

# Rename clusters for clarity
message("\nRenaming cluster identities...")
cluster_names <- paste0("C", selected_clusters)
names(cluster_names) <- levels(TN.subset)
TN.subset <- RenameIdents(TN.subset, cluster_names)

# Update metadata cluster labels
if ("clusters" %in% colnames(TN.subset@meta.data)) {
  original_levels <- levels(TN.subset@meta.data$clusters)
  new_levels <- paste0("C", selected_clusters)
  # Create mapping for all possible levels
  full_mapping <- setNames(paste0("C", 0:(length(original_levels)-1)), original_levels)
  levels(TN.subset@meta.data$clusters) <- full_mapping[original_levels]
}

message("Cluster renaming complete")

# Convert to Monocle3 format ------------------------------------------------
message("\n============================================================")
message("Converting to Monocle3 CDS object")
message("============================================================")

# Set default assay to RNA
DefaultAssay(TN.subset) <- "RNA"

# Join layers to get proper data matrix
message("Joining RNA layers...")
Joined_subset <- JoinLayers(TN.subset)

# Extract expression data
message("Extracting expression data...")
data <- as(as.matrix(Joined_subset@assays$RNA$data), 'sparseMatrix')

# Create AnnotatedDataFrame for phenotype data (cell metadata)
message("Creating phenotype data...")
pd <- new('AnnotatedDataFrame', data = Joined_subset@meta.data)

# Create AnnotatedDataFrame for feature data (gene metadata)
message("Creating feature data...")
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

# Create CellDataSet object for Monocle3
message("Creating Monocle3 CellDataSet...")
cds <- new_cell_data_set(data,
                         cell_metadata = Joined_subset@meta.data,
                         gene_metadata = fData)

message("CDS object created successfully")
message("Cells: ", ncol(cds))
message("Genes: ", nrow(cds))

# ==============================================================================
# Porting Seurat Embeddings to Monocle3
# ==============================================================================
message("\nImporting Harmony-integrated UMAP from Seurat...")

# 1. Force the exact Seurat UMAP coordinates into Monocle3
reducedDims(cds)[["UMAP"]] <- TN.subset[["umap"]]@cell.embeddings

# 2. Force the exact detailed text labels as the Monocle3 clusters
cds@clusters$UMAP$clusters <- as.factor(TN.subset$Detailed_Label)

# 3. Force a single partition (This prevents broken/disconnected graphs)
recreate_partition <- c(rep(1, length(cds@clusters$UMAP$clusters)))
names(recreate_partition) <- names(cds@clusters$UMAP$clusters)
recreate_partition <- as.factor(recreate_partition)
cds@clusters$UMAP$partitions <- recreate_partition

# ==============================================================================
# Learning the Trajectory Graph
# ==============================================================================
message("\nLearning trajectory graph on imported UMAP...")

# Run learn_graph directly on the imported Seurat layout
cds <- learn_graph(cds, use_partition = FALSE)

message("Plotting trajectory...")
p_trajectory <- plot_cells(cds,
                          color_cells_by = "cluster",  # Uses the imported Detailed_Labels
                          label_groups_by_cluster = FALSE,
                          label_leaves = TRUE,
                          label_branch_points = TRUE,
                          graph_label_size = 6,
                          group_label_size = 6) +
  ggtitle("Learned Trajectory (Seurat UMAP)") +
  theme(
    text = element_text(size = 18),                    # <-- Global text size
    axis.title = element_text(size = 20, face = "bold"), # <-- Axis titles (UMAP 1 / UMAP 2)
    axis.text = element_text(size = 16),               # <-- Axis tick numbers
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5)
  )

safe_save_pdf(p_trajectory, file.path(output_dir, "trajectory_by_cluster.pdf"))



message("Ordering cells in pseudotime")

# Set the root to the FIRST cluster you listed in your -c argument
root_cluster <- selected_clusters[1]
message("Setting root cluster to: ", root_cluster)

# Extract the specific cell barcodes that belong to this root cluster
root_cells <- rownames(colData(cds)[as.character(colData(cds)$clusters) == as.character(root_cluster), ])

if (length(root_cells) == 0) {
  stop("Could not find any cells for the root cluster: ", root_cluster)
}

# Order cells non-interactively using the extracted root_cells
cds <- order_cells(cds, root_cells = root_cells)

message("Pseudotime ordering complete")

# Plot pseudotime
message("\nPlotting pseudotime...")
p_pseudotime <- plot_cells(cds,
                           color_cells_by = "pseudotime",
                           label_cell_groups = FALSE,
                           label_leaves = FALSE,
                           label_branch_points = FALSE,
                           label_roots = TRUE,
                           graph_label_size = 6,
group_label_size = 6) +
  ggtitle("Cells Ordered by Pseudotime") +
  theme(
    text = element_text(size = 18),                    # <-- Global text size
    axis.title = element_text(size = 20, face = "bold"), # <-- Axis titles (UMAP 1 / UMAP 2)
    axis.text = element_text(size = 16),               # <-- Axis tick numbers
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5)
  )

safe_save_pdf(p_pseudotime, file.path(output_dir, "pseudotime.pdf"))

# Combined trajectory plots
message("\nGenerating combined trajectory plots...")
p_combined <- plot_cells(cds,
                         color_cells_by = "clusters",
                         label_groups_by_cluster = TRUE,
                         label_leaves = TRUE,
                         label_branch_points = TRUE,
                         label_roots = TRUE,
                         graph_label_size = 6,
group_label_size = 6) +
  ggtitle("Trajectory with Cluster Labels") + 
  theme(
    text = element_text(size = 18),                    # <-- Global text size
    axis.title = element_text(size = 20, face = "bold"), # <-- Axis titles (UMAP 1 / UMAP 2)
    axis.text = element_text(size = 16),               # <-- Axis tick numbers
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5)
  )

safe_save_pdf(p_combined, file.path(output_dir, "trajectory_combined.pdf"))

# Split by sample if orig.ident1 exists
if ("orig.ident1" %in% colnames(colData(cds))) {
  message("\nPlotting trajectory split by sample group...")
  p_by_sample <- plot_cells(cds,
                            color_cells_by = "orig.ident1",
                            label_cell_groups = FALSE,
                            label_leaves = FALSE,
                            label_branch_points = TRUE,
                            graph_label_size = 6,
group_label_size = 6) +
    ggtitle("Trajectory Colored by Sample") +
    theme(
    text = element_text(size = 18),                    # <-- Global text size
    axis.title = element_text(size = 20, face = "bold"), # <-- Axis titles (UMAP 1 / UMAP 2)
    axis.text = element_text(size = 16),               # <-- Axis tick numbers
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5)
  )

  safe_save_pdf(p_by_sample, file.path(output_dir, "trajectory_by_sample.pdf"))
}

# Save results --------------------------------------------------------------
message("\n============================================================")
message("Saving results")
message("============================================================")

# Save CDS object
message("Saving Monocle3 CDS object...")
saveRDS(cds, file.path(output_dir, "monocle3_cds.rds"))

# Extract pseudotime values and save
message("Extracting pseudotime values...")

# Extract values and convert Inf to NA for unreachable cells
p_time <- pseudotime(cds)
p_time[is.infinite(p_time)] <- NA 

pseudotime_df <- data.frame(
  cell = colnames(cds),
  pseudotime = p_time,
  cluster = colData(cds)$seurat_clusters
)

# Add sample info if available
if ("orig.ident1" %in% colnames(colData(cds))) {
  pseudotime_df$sample <- colData(cds)$orig.ident1
}

write.csv(pseudotime_df, file.path(output_dir, "pseudotime_values.csv"), row.names = FALSE)
message("Pseudotime values saved")

# Summary statistics
message("\n============================================================")
message("Pseudotime Summary Statistics")
message("============================================================")
message("Min pseudotime: ", round(min(pseudotime_df$pseudotime, na.rm = TRUE), 3))
message("Max pseudotime: ", round(max(pseudotime_df$pseudotime, na.rm = TRUE), 3))
message("Mean pseudotime: ", round(mean(pseudotime_df$pseudotime, na.rm = TRUE), 3))
message("Median pseudotime: ", round(median(pseudotime_df$pseudotime, na.rm = TRUE), 3))

# Pseudotime by cluster
pseudotime_by_cluster <- pseudotime_df %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),
    mean_pseudotime = mean(pseudotime, na.rm = TRUE),
    median_pseudotime = median(pseudotime, na.rm = TRUE),
    min_pseudotime = min(pseudotime, na.rm = TRUE),
    max_pseudotime = max(pseudotime, na.rm = TRUE)
  )

write.csv(pseudotime_by_cluster, file.path(output_dir, "pseudotime_by_cluster.csv"), row.names = FALSE)
message("\nPseudotime statistics by cluster:")
print(pseudotime_by_cluster)

# ==============================================================================
# Create Violin Plot
# ==============================================================================
message("\nGenerating pseudotime violin plot...")

# Strictly filter out NA values so ggplot does not crash
plot_df <- pseudotime_df[!is.na(pseudotime_df$pseudotime), ]

p_violin <- ggplot(plot_df, aes(x = cluster, y = pseudotime, fill = cluster)) +
  geom_violin() +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "Pseudotime Distribution by Cluster",
       x = "Cluster",
       y = "Pseudotime") +
  coord_flip()

safe_save_pdf(p_violin, file.path(output_dir, "pseudotime_violin_by_cluster.pdf"))

# Find genes that change as a function of pseudotime ------------------------
message("\n============================================================")
message("Finding pseudotime-dependent genes")
message("============================================================")

message("This may take a while for large datasets...")
message("Testing genes for pseudotime dependence using graph_test...")

# Use graph_test to find genes that vary along the trajectory
tryCatch({
  gene_fits <- graph_test(cds, neighbor_graph = "principal_graph", cores = 1)

  # Filter significant genes
  sig_genes <- gene_fits %>%
    filter(q_value < 0.05) %>%
    arrange(q_value)

  message("Found ", nrow(sig_genes), " genes significantly changing along pseudotime")

  # Save results
  write.csv(gene_fits, file.path(output_dir, "pseudotime_gene_fits.csv"), row.names = FALSE)
  write.csv(sig_genes, file.path(output_dir, "pseudotime_significant_genes.csv"), row.names = FALSE)

  # Plot top genes
  if (nrow(sig_genes) > 0) {
    message("\nPlotting top pseudotime-dependent genes...")
    top_genes <- head(sig_genes$gene_short_name, 9)

    p_genes <- plot_genes_in_pseudotime(cds[top_genes,],
                                        color_cells_by = "clusters",
                                        min_expr = 0.5,
                                        ncol = 3) +
  theme(
    text = element_text(size = 20),                    # <-- Makes gene names (strip text) larger
    axis.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 16),
    legend.position = "bottom",                        # <-- Moves legend to bottom to save horizontal space
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 22, face = "bold") # <-- Makes the Gene Names at the top of each box much larger
  )

    safe_save_pdf(p_genes, file.path(output_dir, "top_pseudotime_genes.pdf"), w = 18, h = 18)
  }

}, error = function(e) {
  message("Warning: Could not compute pseudotime-dependent genes: ", e$message)
})

# Final summary -------------------------------------------------------------
message("\n============================================================")
message("Pseudotime Analysis Complete!")
message("============================================================")
message("Output directory: ", output_dir)
message("Files generated:")
message("  - subset_umap.pdf: UMAP of selected clusters")
message("  - monocle3_umap_clusters.pdf: Monocle3 UMAP colored by cluster")
message("  - trajectory_by_cluster.pdf: Learned trajectory")
message("  - pseudotime.pdf: Cells colored by pseudotime")
message("  - trajectory_combined.pdf: Combined trajectory visualization")
message("  - pseudotime_values.csv: Pseudotime values for each cell")
message("  - pseudotime_by_cluster.csv: Summary statistics by cluster")
message("  - pseudotime_violin_by_cluster.pdf: Distribution of pseudotime")
message("  - monocle3_cds.rds: Saved Monocle3 CDS object")

if (file.exists(file.path(output_dir, "pseudotime_significant_genes.csv"))) {
  message("  - pseudotime_gene_fits.csv: All gene test results")
  message("  - pseudotime_significant_genes.csv: Significant genes")
  message("  - top_pseudotime_genes.pdf: Expression of top genes")
}

message("\n============================================================")
message("Completed at: ", Sys.time())
message("============================================================")

# Close log file
sink(type = "message")
sink(type = "output")
close(log_conn)
