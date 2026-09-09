#!/usr/bin/env Rscript
# Global CellChat cross-talk analysis using the top-level annotated object.

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(CellChat)
  library(patchwork)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(circlize)
  library(RColorBrewer)
  library(ComplexHeatmap)
})
source("workflow/scripts/00_utils.R")

# ==============================================================================
# 0. DEFINE CUSTOM COLORS GLOBALLY
# ==============================================================================
GLOBAL_COLORS <- c(
  "Endothelial cells" = "#146B78", "Endothelial_cells" = "#146B78",
  "Epithelial cells"  = "#5FAF8D", "Epithelial_cells"  = "#5FAF8D",
  "Keratinocytes"     = "#70AAA9", 
  "Fibroblasts"       = "#F2A928", 
  "Macrophages"       = "#F37758", 
  "Mast cells"        = "#C43D73", "Mast_cells"        = "#C43D73",
  "T cells"           = "#A8D9D6", "T_cells"           = "#A8D9D6",
  "Melanocytes"       = "#1A1029"
)

# ==============================================================================
# 1. SETUP COMMAND LINE ARGUMENTS
# ==============================================================================
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Path to global TN.combined_annotated.rds"),
  make_option(c("-o", "--outdir"), type = "character", help = "Output directory"),
  make_option(c("--min_cells"), type = "integer", default = 10),
  make_option(c("--pval_thresh"), type = "numeric", default = 0.05),
  make_option(c("--seed"), type = "integer", default = 42)
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
message("Loading Global Annotated Object: ", opt$input)

combined <- readRDS(opt$input)

# Join layers before any splitting or extraction
combined[["RNA"]] <- JoinLayers(combined[["RNA"]])

# Set the active identity to the clean, top-level global cell types[cite: 1]
combined$Detailed_Label <- factor(as.character(combined$cell_type_short))
Idents(combined) <- "Detailed_Label"

# ==============================================================================
# 2. PLOTTING HELPER FUNCTIONS
# ==============================================================================
generate_chord_plot <- function(interaction_df, out_pdf, plot_title) {
  if (nrow(interaction_df) == 0) return(NULL)
  
  interaction_df <- interaction_df %>%
    mutate(
      source_group = as.character(source), 
      target_group = as.character(target)
    ) %>%
    filter(!is.na(prob), prob > 0, is.na(pval) | pval <= 0.05)[cite: 1]

  cell_network <- interaction_df %>%
    group_by(source_group, target_group) %>%
    summarise(interaction_strength = sum(prob, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(interaction_strength))[cite: 1]
    
  if(nrow(cell_network) == 0) return(NULL)
  
  cell_network_plot <- cell_network %>% slice_max(order_by = interaction_strength, n = 40, with_ties = FALSE)[cite: 1]
  cell_groups <- sort(unique(c(cell_network_plot$source_group, cell_network_plot$target_group)))[cite: 1]
  
  chord_matrix <- matrix(0, nrow = length(cell_groups), ncol = length(cell_groups), dimnames = list(cell_groups, cell_groups))[cite: 1]
  
  for (i in seq_len(nrow(cell_network_plot))) {
    src <- as.character(cell_network_plot$source_group[i])
    tgt <- as.character(cell_network_plot$target_group[i])
    chord_matrix[src, tgt] <- cell_network_plot$interaction_strength[i]
  }

  base_colors <- GLOBAL_COLORS
  missing_groups <- setdiff(cell_groups, names(base_colors))[cite: 1]
  if (length(missing_groups) > 0) {
    extra_colors <- rep(brewer.pal(max(3, min(8, length(missing_groups))), "Set3"), length.out = length(missing_groups))[cite: 1]
    names(extra_colors) <- missing_groups
    base_colors <- c(base_colors, extra_colors)[cite: 1]
  }
  grid_colors <- base_colors[cell_groups]
  
  min_strength <- min(cell_network_plot$interaction_strength)[cite: 1]
  max_strength <- max(cell_network_plot$interaction_strength)[cite: 1]
  if (min_strength == max_strength) max_strength <- min_strength + 0.01 
  
  link_color_func <- colorRamp2(breaks = c(min_strength, max_strength), colors = c("#398F8B", "#F37758"))[cite: 1]
  
  pdf(out_pdf, width = 12, height = 12)[cite: 1]
  circos.clear()
  circos.par(start.degree = 90, gap.after = rep(6, length(cell_groups)), track.margin = c(0.01, 0.01), cell.padding = c(0.01, 0.01, 0.01, 0.01))[cite: 1]
  
  chordDiagram(
    x = chord_matrix, grid.col = grid_colors, col = link_color_func, directional = 1,
    direction.type = c("arrows", "diffHeight"), link.arr.type = "big.arrow", link.sort = TRUE,
    link.largest.ontop = TRUE, transparency = 0.3, annotationTrack = c("grid", "name"),
    preAllocateTracks = list(track.height = 0.08)
  )[cite: 1]
  title(plot_title, cex.main = 1.6)[cite: 1]
  
  lgd <- Legend(col_fun = link_color_func, title = "Interaction\nStrength", direction = "horizontal")[cite: 1]
  draw(lgd, x = unit(0.5, "npc"), y = unit(0.05, "npc"), just = c("center", "bottom"))[cite: 1]
  circos.clear()
  dev.off()
}

# ==============================================================================
# 3. CORE CELLCHAT MATH ENGINE
# ==============================================================================
run_cellchat_workflow <- function(seu_obj, condition_name, outdir_base) {
  message(sprintf("\n=== Running CellChat Workflow for Condition: %s ===", condition_name))[cite: 1]
  
  seu_obj$Detailed_Label <- droplevels(factor(seu_obj$Detailed_Label))[cite: 1]
  
  if (ncol(seu_obj) < opt$min_cells) {
    message("Not enough cells to run CellChat for ", condition_name)[cite: 1]
    return(NULL)
  }

  set.seed(opt$seed)
  data.input <- GetAssayData(seu_obj, assay = "RNA", layer = "data")

  cc <- createCellChat(object = data.input, meta = seu_obj@meta.data, group.by = "Detailed_Label")[cite: 1]
  cc@options$color.use <- GLOBAL_COLORS[levels(cc@idents)]
  cc@DB <- subsetDB(CellChatDB.human, search = "Secreted Signaling")[cite: 1]

  message("Computing probabilities...")
  cc <- subsetData(cc)[cite: 1]
  cc <- identifyOverExpressedGenes(cc)[cite: 1]
  cc <- identifyOverExpressedInteractions(cc)[cite: 1]
  cc <- computeCommunProb(cc)[cite: 1]
  cc <- filterCommunication(cc, min.cells = opt$min_cells)[cite: 1]
  cc <- computeCommunProbPathway(cc)[cite: 1]
  cc <- aggregateNet(cc)[cite: 1]

  if (condition_name == "Global") {
    pfx <- "Global"
  } else if (condition_name == "PMH") {
    pfx <- "PMH_Only"
  } else if (condition_name == "Healthy") {
    pfx <- "Healthy_Only"
  } else {
    pfx <- condition_name
  }

  df <- subsetCommunication(cc, thresh = opt$pval_thresh)[cite: 1]
  write.csv(df, file.path(outdir_base, sprintf("Table_Interactions_%s.csv", pfx)), row.names = FALSE)[cite: 1]
  
  # 1. Generate Chord Plot
  generate_chord_plot(df, file.path(outdir_base, sprintf("Chord_Diagram_%s.pdf", pfx)), sprintf("Macro Cell-Cell Communication (%s)", condition_name))[cite: 1]

  # 2. Generate Circle Plot
  # 2. Generate Circle Plot
  message("Generating Network Circle Plot...")
  
  # Extract the exact cell types present in the network matrix
  active_celltypes <- rownames(cc@net$count)
  
  # Subset your global palette and explicitly name the vector
  my_colors <- GLOBAL_COLORS[active_celltypes]
  names(my_colors) <- active_celltypes
  
  # Catch any NA values (missing cell types) and assign a fallback color
  missing_colors <- is.na(my_colors)
  if (any(missing_colors)) {
    message("Warning: Missing predefined colors for: ", paste(active_celltypes[missing_colors], collapse = ", "))
    # Generate fallback colors dynamically
    fallback_palette <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(sum(missing_colors))
    my_colors[missing_colors] <- fallback_palette
  }

  pdf(file.path(outdir_base, sprintf("Network_Circle_%s.pdf", pfx)), width = 12, height = 12)
  netVisual_circle(
    cc@net$count, 
    weight.scale = TRUE, 
    label.edge = FALSE, 
    title.name = sprintf("Total Macro Interactions (%s)", condition_name),
    color.use = my_colors # Use the validated, gap-free color vector
  )
  dev.off()

  # 3. Generate Bubble Plot (NEW)
  message("Generating Bubble Plot...")
  bubble_plot <- netVisual_bubble(cc, remove.isolate = FALSE) + 
    ggtitle(sprintf("Significant Interactions (%s)", condition_name))
  ggsave(file.path(outdir_base, sprintf("Bubble_Plot_%s.pdf", pfx)), plot = bubble_plot, width = 12, height = 10, dpi = 300)

  # 4. Output the RDS
  rds_name <- ifelse(condition_name == "Global", "macro_global_cellchat.rds", sprintf("macro_global_cellchat_%s.rds", pfx))[cite: 1]
  saveRDS(cc, file.path(outdir_base, rds_name))[cite: 1]
  
  return(cc)
}

# ==============================================================================
# 4. EXECUTE THE INDIVIDUAL WORKFLOWS
# ==============================================================================
# 1. Run Global (All cells: Healthy + PMH)[cite: 1]
cc_global <- run_cellchat_workflow(combined, "Global", opt$outdir)[cite: 1]

if ("Condition" %in% colnames(combined@meta.data)) {
  
  # 2. Run PMH Cells Only[cite: 1]
  combined_pmh <- subset(combined, Condition == "PMH")[cite: 1]
  if (ncol(combined_pmh) > 0) {
    cc_pmh <- run_cellchat_workflow(combined_pmh, "PMH", opt$outdir)[cite: 1]
  } else {
    message("No PMH cells found to run condition-specific CellChat.")[cite: 1]
  }
  
  # 3. Run Healthy Cells Only[cite: 1]
  combined_healthy <- subset(combined, Condition == "Healthy")[cite: 1]
  if (ncol(combined_healthy) > 0) {
    cc_healthy <- run_cellchat_workflow(combined_healthy, "Healthy", opt$outdir)[cite: 1]
  } else {
    message("No Healthy cells found to run condition-specific CellChat.")[cite: 1]
  }
}

# ==============================================================================
# 5. MERGE AND COMPARE (NEW)
# ==============================================================================
# If both conditions ran successfully, merge them to generate a comparison plot
if (exists("cc_pmh") && exists("cc_healthy") && !is.null(cc_pmh) && !is.null(cc_healthy)) {
  message("\n=== Merging PMH and Healthy for Comparison ===")
  
  # Combine into a list
  object.list <- list(Healthy = cc_healthy, PMH = cc_pmh)
  
  # Lift up CellChat objects to ensure they have the same cell labels before merging
  cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))
  
  # Save the combined RDS
  saveRDS(cellchat_merged, file.path(opt$outdir, "macro_cellchat_merged_comparison.rds"))
  
  # Generate Comparative Bubble Plot
  message("Generating Comparative Bubble Plot (Healthy vs PMH)...")
  
  # This plots the probabilities of interactions in Healthy vs PMH
  bubble_comp <- netVisual_bubble(
    cellchat_merged, 
    comparison = c(1, 2), 
    angle.x = 45, 
    remove.isolate = FALSE, 
    title.name = "Increased/Decreased Signaling in PMH vs Healthy"
  )
  
  ggsave(file.path(opt$outdir, "Bubble_Plot_Comparison_Healthy_vs_PMH.pdf"), plot = bubble_comp, width = 14, height = 10, dpi = 300)
}

message("\nMulti-way CellChat Analysis Complete!")