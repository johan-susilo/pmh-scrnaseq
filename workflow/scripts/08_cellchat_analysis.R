#!/usr/bin/env Rscript
# Unified Dynamic Multi-way CellChat cross-talk analysis for Global and Sub-cluster conditions.

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
# 1. SETUP COMMAND LINE ARGUMENTS
# ==============================================================================
option_list <- list(
  make_option(c("--basedir"), type = "character", default = "results/03_subsets"),
  make_option(c("--celltypes"), type = "character", help = "Comma-separated list of cell types OR 'Global'"),
  make_option(c("--global_input"), type = "character", help = "Path to global annotated RDS (if celltypes == 'Global')"),
  make_option(c("--outdir"), type = "character", help = "Output directory"),
  make_option(c("--min_cells"), type = "integer", default = 10),
  make_option(c("--pval_thresh"), type = "numeric", default = 0.05),
  make_option(c("--seed"), type = "integer", default = 42)
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 2. CONDITIONAL LOAD: GLOBAL VS DETAILED SUBSETS
# ==============================================================================
if (opt$celltypes == "Global") {
  message("Loading Global Annotated Object: ", opt$global_input)
  combined <- readRDS(opt$global_input)
  
  # Set identity for global object
  combined$Detailed_Label <- factor(as.character(combined$cell_type_short))
  Idents(combined) <- "Detailed_Label"
  
} else {
  types <- trimws(unlist(strsplit(opt$celltypes, ",")))
  if (length(types) < 2) stop("You must specify at least 2 cell types for detailed CellChat.")
  message("Loading subset objects for: ", paste(types, collapse = ", "))
  
  obj_list <- list()
  for (i in seq_along(types)) {
    ctype <- types[i]
    rds_path <- file.path(opt$basedir, ctype, "00_data", paste0(ctype, "_detailed_annotated.rds"))
    if (!file.exists(rds_path)) stop("Missing RDS for ", ctype)
    
    obj <- readRDS(rds_path)
    assay_use <- ifelse("RNA" %in% Assays(obj), "RNA", DefaultAssay(obj))
    if (inherits(obj[[assay_use]], "Assay5")) try({ obj <- JoinLayers(obj) }, silent = TRUE)
    
    mat <- tryCatch(GetAssayData(obj, assay = assay_use, layer = "counts"), 
                    error = function(e) GetAssayData(obj, assay = assay_use, slot = "counts"))
    
    prefix <- paste0(toupper(substr(ctype, 1, 4)), i, "_")
    colnames(mat) <- paste0(prefix, colnames(mat))
    
    meta <- obj@meta.data[, intersect(c("orig.ident1", "orig.ident2", "Condition", "Detailed_Label"), colnames(obj@meta.data)), drop = FALSE]
    rownames(meta) <- paste0(prefix, rownames(meta))
    meta$source_celltype <- ctype
    
    obj_list[[ctype]] <- list(mat = mat, meta = meta)
  }

  common_genes <- Reduce(intersect, lapply(obj_list, function(x) rownames(x$mat)))
  combined_counts <- do.call(cbind, lapply(obj_list, function(x) x$mat[common_genes, ]))
  combined_meta <- do.call(rbind, lapply(obj_list, function(x) x$meta))

  message("Building new multi-way Seurat object...")
  combined <- CreateSeuratObject(counts = combined_counts, meta.data = combined_meta)
  combined <- NormalizeData(combined, verbose = FALSE)

  # Clean and collapse fragmented labels
  clean_labels <- function(x) {
    x <- str_replace_all(x, "^T_cells(_.*)?$|^T Cells(_.*)?$", "T cells")
    x <- str_replace_all(x, "^Mast_cells(_.*)?$|^Mast Cells(_.*)?$", "Mast cells")
    x <- str_replace_all(x, "^Melanocytes(_.*)?$", "Melanocytes")
    x <- str_replace_all(x, "^Keratinocytes(_.*)?$", "Keratinocytes")
    x <- str_replace_all(x, "^Endothelial_cells(_.*)?$", "Endothelial cells")
    return(x)
  }

  combined$Detailed_Label <- clean_labels(as.character(combined$Detailed_Label))
  combined$Detailed_Label <- factor(combined$Detailed_Label)
  Idents(combined) <- "Detailed_Label"
  
  rm(obj_list, combined_counts)
  gc()
}

# ==============================================================================
# 3. PLOTTING HELPER FUNCTIONS
# ==============================================================================
generate_chord_plot <- function(interaction_df, out_pdf, plot_title) {
  if (nrow(interaction_df) == 0) return(NULL)
  
  interaction_df <- interaction_df %>%
    mutate(source_group = as.character(source), target_group = as.character(target)) %>%
    filter(!is.na(prob), prob > 0, is.na(pval) | pval <= 0.05)
  
  cell_network <- interaction_df %>%
    group_by(source_group, target_group) %>%
    summarise(interaction_strength = sum(prob, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(interaction_strength))
    
  if(nrow(cell_network) == 0) return(NULL)
  
  cell_network_plot <- cell_network %>% slice_max(order_by = interaction_strength, n = 40, with_ties = FALSE)
  cell_groups <- sort(unique(c(cell_network_plot$source_group, cell_network_plot$target_group)))
  
  chord_matrix <- matrix(0, nrow = length(cell_groups), ncol = length(cell_groups), dimnames = list(cell_groups, cell_groups))
  for (i in seq_len(nrow(cell_network_plot))) {
    src <- as.character(cell_network_plot$source_group[i])
    tgt <- as.character(cell_network_plot$target_group[i])
    chord_matrix[src, tgt] <- cell_network_plot$interaction_strength[i]
  }
  
  # Merged color palette supporting both Global and Detailed clusters
  base_colors <- c(
    "T cells" = "#A8D9D6", "Mast cells" = "#C43D73", "Melanocytes" = "#1A1029",
    "Keratinocytes" = "#70AAA9", "Endothelial cells" = "#246A73", "Epithelial cells" = "#86af4d",
    "Fibroblasts" = "#F2542D", "Macrophages" = "#F6AE2D",
    "M_Homeostatic_Resident" = "#0072B2", "M_Chronic_GvHD_Inflammatory" = "#E69F00",
    "F1_Superficial" = "#146B78", "F3_FRC_like" = "#70AAA9", "F4_DS_DPEP1" = "#82B8B7",
    "F4_DP_HHIP" = "#A8D9D6", "F5_NGFR" = "#5FAF8D", "F6_Inflammatory_Myofibroblast" = "#E85D3F"
  )
  
  missing_groups <- setdiff(cell_groups, names(base_colors))
  if (length(missing_groups) > 0) {
    extra_colors <- rep(brewer.pal(max(3, min(8, length(missing_groups))), "Set3"), length.out = length(missing_groups))
    names(extra_colors) <- missing_groups
    base_colors <- c(base_colors, extra_colors)
  }
  grid_colors <- base_colors[cell_groups]
  
  min_strength <- min(cell_network_plot$interaction_strength)
  max_strength <- max(cell_network_plot$interaction_strength)
  if (min_strength == max_strength) max_strength <- min_strength + 0.01 
  
  link_color_func <- colorRamp2(breaks = c(min_strength, max_strength), colors = c("#398F8B", "#F37758"))
  
  pdf(out_pdf, width = 12, height = 12)
  circos.clear()
  circos.par(start.degree = 90, gap.after = rep(6, length(cell_groups)), track.margin = c(0.01, 0.01), cell.padding = c(0.01, 0.01, 0.01, 0.01))
  
  chordDiagram(
    x = chord_matrix, grid.col = grid_colors, col = link_color_func, directional = 1,
    direction.type = c("arrows", "diffHeight"), link.arr.type = "big.arrow", link.sort = TRUE,
    link.largest.ontop = TRUE, transparency = 0.3, annotationTrack = c("grid", "name"),
    preAllocateTracks = list(track.height = 0.08)
  )
  title(plot_title, cex.main = 1.6)
  
  lgd <- Legend(col_fun = link_color_func, title = "Interaction\nStrength", direction = "horizontal")
  draw(lgd, x = unit(0.5, "npc"), y = unit(0.05, "npc"), just = c("center", "bottom"))
  circos.clear()
  dev.off()
}

save_dynamic_bubble <- function(cellchat_obj, sources, targets, out_path, plot_title) {
  tryCatch({
    p <- netVisual_bubble(cellchat_obj, sources.use = sources, targets.use = targets, remove.isolate = TRUE, thresh = opt$pval_thresh)
    if (is.null(p)) return(NULL)
    
    p <- p + ggtitle(plot_title) + 
         theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
               plot.margin = margin(t = 20, r = 20, b = 60, l = 20))
    
    n_pairs <- length(unique(p$data$source.target))
    n_pathways <- length(unique(p$data$interaction_name_2))
    if (is.null(n_pairs) || n_pairs == 0) n_pairs <- 10
    if (is.null(n_pathways) || n_pathways == 0) n_pathways <- 20
    
    dyn_width <- max(10, n_pairs * 0.5 + 4)
    dyn_height <- max(8, n_pathways * 0.25 + 6) 
    
    pdf(out_path, width = dyn_width, height = dyn_height)
    print(p)
    dev.off()
  }, error = function(e) message("Skipping bubble plot ", basename(out_path), ": ", e$message))
}

# ==============================================================================
# 4. CORE CELLCHAT MATH & WORKFLOW ENGINE
# ==============================================================================
run_cellchat_workflow <- function(seu_obj, condition_name, outdir_base) {
  message(sprintf("\n=== Running CellChat Workflow for Condition: %s ===", condition_name))
  
  seu_obj$Detailed_Label <- droplevels(factor(seu_obj$Detailed_Label))
  
  if (ncol(seu_obj) < opt$min_cells) {
    message("Not enough cells to run CellChat for ", condition_name)
    return(NULL)
  }

  set.seed(opt$seed)
  data.input <- tryCatch(GetAssayData(seu_obj, assay = "RNA", layer = "data"), 
                         error = function(e) GetAssayData(seu_obj, assay = "RNA", slot = "data"))

  cc <- createCellChat(object = data.input, meta = seu_obj@meta.data, group.by = "Detailed_Label")
  cc@DB <- subsetDB(CellChatDB.human, search = "Secreted Signaling")

  message("Computing probabilities...")
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc)
  cc <- filterCommunication(cc, min.cells = opt$min_cells)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)

  # Map condition to filename prefix
  if (condition_name == "AllCells") {
    pfx <- "Global"
  } else if (condition_name == "PMH") {
    pfx <- "PMH_Only"
  } else if (condition_name == "Healthy") {
    pfx <- "Healthy_Only"
  } else {
    pfx <- condition_name
  }

  df <- subsetCommunication(cc, thresh = opt$pval_thresh)
  write.csv(df, file.path(outdir_base, sprintf("Table_Interactions_%s.csv", pfx)), row.names = FALSE)
  generate_chord_plot(df, file.path(outdir_base, sprintf("Chord_Diagram_%s.pdf", pfx)), sprintf("Cell-Cell Communication (%s)", condition_name))

  pdf(file.path(outdir_base, sprintf("Network_Circle_%s.pdf", pfx)), width = 12, height = 12)
  netVisual_circle(cc@net$count, weight.scale = TRUE, label.edge = FALSE, title.name = sprintf("Total Interactions (%s)", condition_name))
  dev.off()

  # Only generate dynamic sender/receiver bubble plots for the subset networks, not global macro
  if (opt$celltypes != "Global") {
    message("Generating dynamic bubble plots...")
    types <- trimws(unlist(strsplit(opt$celltypes, ",")))
    for (sender in types) {
      sender_label <- tools::toTitleCase(sender)
      sender_clusters <- as.character(unique(seu_obj$Detailed_Label[seu_obj$source_celltype == sender]))
      receiver_clusters <- as.character(unique(seu_obj$Detailed_Label[seu_obj$source_celltype != sender]))
      
      if (length(sender_clusters) > 0 && length(receiver_clusters) > 0) {
        save_dynamic_bubble(cc, sender_clusters, receiver_clusters, 
                            file.path(outdir_base, sprintf("Bubble_%s_%s_to_Others.pdf", pfx, sender_label)), 
                            sprintf("%s Signals: %s -> All Others", condition_name, sender_label))
      }
    }
  }

  rds_name <- ifelse(condition_name == "AllCells", "multiway_cellchat.rds", sprintf("multiway_cellchat_%s.rds", pfx))
  saveRDS(cc, file.path(outdir_base, rds_name))
  
  return(cc)
}

# ==============================================================================
# 5. EXECUTE THE WORKFLOW
# ==============================================================================
# 1. Run All Cells (Healthy + PMH)
cc_global <- run_cellchat_workflow(combined, "AllCells", opt$outdir)

if ("Condition" %in% colnames(combined@meta.data)) {
  
  # 2. Run PMH Cells Only
  combined_pmh <- subset(combined, Condition == "PMH")
  if (ncol(combined_pmh) > 0) {
    cc_pmh <- run_cellchat_workflow(combined_pmh, "PMH", opt$outdir)
  }
  
  # 3. Run Healthy Cells Only (Ensure "Healthy" matches your metadata exactly)
  combined_healthy <- subset(combined, Condition == "Healthy")
  if (ncol(combined_healthy) > 0) {
    cc_healthy <- run_cellchat_workflow(combined_healthy, "Healthy", opt$outdir)
  }
  
}

message("\nUnified CellChat Analysis Complete!")