suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(pheatmap)
})

source("workflow/scripts/00_utils.R")

#Rscript workflow/scripts/04_detail_annotation.R -c fibroblast  -i results/03_subsets/fibroblast/processed/fibroblast_subset_processed.rds -o results/03_subsets
# Rscript workflow/scripts/04_detail_annotation.R -c macrophage  -i results/03_subsets/macrophage/processed/macrophage_subset_processed.rds -o results/03_subsets
# Rscript workflow/scripts/04_detail_annotation.R -c Mast -i results/03_subsets/Mast/processed/Mast_subset_processed.rds -o results/03_subsets


# ==============================================================================
# 1. COMMAND LINE ARGUMENTS (SNAKEMAKE INTEGRATION)
# ==============================================================================
option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Path to the subset RDS object"),
  make_option(c("-c", "--celltype"), type="character", help="The cell type being processed (e.g., 'fibroblast' or 'macrophage')"),
  make_option(c("-o", "--outdir"), type="character", help="Base subset_cluster output directory"),
  make_option(c("--seed"), type = "integer", default = 42, help = "Global random seed for reproducibility"),
  make_option(c("--uncertain_margin"), type = "double", default = 0.10,
              help = "Relative margin between the best and second-best cluster module score below which a cluster is called 'Uncertain'. E.g. 0.10 means top2 must be more than 10%% below top1 (relative to top1) to avoid Uncertain. Only applied when top1 is comfortably above zero -- see --uncertain_abs_margin. [default: %default]"),
  make_option(c("--uncertain_abs_margin"), type = "double", default = 0.002,
              help = "Absolute margin between the best and second-best cluster module score below which a cluster is called 'Uncertain'. This is a floor that keeps the call reachable when top1 is small or near zero, where a purely relative margin would almost never be satisfied. The cluster is required to clear BOTH margins (blended as max(abs_margin, relative_margin * abs(top1))). [default: %default]")
)
opt <- parse_args(OptionParser(option_list=option_list))
set.seed(opt$seed)

# Create highly organized subdirectories dynamically
base_out_dir <- file.path(opt$outdir, opt$celltype, "haniffa_analysis")
dirs <- list(
  umaps_global  = file.path(base_out_dir, "1_UMAPs", "Global"),
  umaps_sample  = file.path(base_out_dir, "1_UMAPs", "Per_Sample"),
  val_vln       = file.path(base_out_dir, "2_Validation_Plots", "VlnPlots"),
  val_dot       = file.path(base_out_dir, "2_Validation_Plots", "DotPlots"),
  master_sum    = file.path(base_out_dir, "3_Summary_Plots"),
  mucin         = file.path(base_out_dir, "4_Mucin_ECM"),
  proportions   = file.path(base_out_dir, "5_Proportions")
)
lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)



# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

# Generate clean, uniform proportion barplots (Modified to accept your dictionary colors)
create_proportion_barplot <- function(seurat_obj, group_col, title, custom_colors = NULL) {
  prop_df <- as.data.frame(prop.table(table(Idents(seurat_obj), seurat_obj[[group_col]][[1]]), margin = 2))
  colnames(prop_df) <- c("Cluster", "Group", "Proportion")

  p <- ggplot(prop_df %>% filter(Proportion > 0), aes(x = Group, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
    # Removed the text overlay as it gets too cluttered with many detailed clusters
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
          axis.text.y = element_text(color = "black"),
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    scale_y_continuous(labels = scales::percent) +
    ggtitle(title)

  # Apply detailed annotation colors if provided
  if (!is.null(custom_colors)) {
    p <- p + scale_fill_manual(values = custom_colors)
  }

  return(p)
}

message(paste("Loading mathematically cleaned", toupper(opt$celltype), "object..."))
seu_obj <- readRDS(opt$input)
DefaultAssay(seu_obj) <- "RNA"

# ==============================================================================
# 2. DICTIONARIES & COLORS (FIBROBLAST & MACROPHAGE)
# ==============================================================================
biological_dictionary <- list(
  "F1_Superficial"                = c("APCDD1", "COL18A1", "COL23A1", "COL13A1", "NKD2", "RSPO1", "AXIN2", "WIF1"),
  "F2_Universal"                  = c("CD34", "PI16", "DPP4", "MFAP5", "PCOLCE2", "SLPI", "CD70", "LGR5"),
  "F2_F3_Perivascular"            = c("CXCL12", "APOE", "EFEMP1", "APOC1", "C7", "PLA2G2A", "PPARG", "MYOC", "GDF10"),
  "F3_FRC_like"                   = c("CCL19", "CD74", "CH25H", "TNFSF13B", "IL33", "IRF8", "IL15", "VCAM1", "HLA-DRA", "HLA-DRB1"),
  "F4_DS_DPEP1"                   = c("MEF2C", "DPEP1", "MYL4"),
  "F4_TNN_COCH"                   = c("TNN", "COCH", "CRABP1", "COL24A1", "RSPO4", "SLITRK6", "NRG3", "MKX", "TNMD"),
  "F4_DP_HHIP"                    = c("CORIN", "BMP7", "LEF1", "HHIP", "RSPO3", "PTCH1"),
  "F5_RAMP1"                      = c("RAMP1", "RELN", "PLEKHA6", "IGFBP2", "SFRP1"),
  "F5_NGFR"                       = c("EBF2", "NGFR", "ITGA6", "CDH19", "CLDN1"),
  "F6_Inflammatory_Myofibroblast" = c("IL11", "IL24", "CXCL5", "CXCL6", "CXCL8", "CXCL13", "MMP1", "MMP3", "IL7R", "CSF3", "TDO2", "WWC1", "CHI3L1", "CCL5", "CCL11", "FAM167A", "HIF1A"),
  "F7_Myofibroblast"              = c("ACTA2", "TAGLN", "COL3A1", "COL5A1", "COL8A1", "POSTN", "LRRC15", "RUNX2", "KIF26B", "ZNF469", "SULF1", "ADAM12", "ADAM19", "CREB3L1", "CCN4", "FABP5", "CDH2", "C1QTNF3", "CADM1", "LRRC17", "KCNMA1", "NRG1", "OGN", "WNT2", "COL10A1", "LAMP5"),
  "F8_Fascia_like_Myofibroblast"  = c("ACAN", "SCX", "THBS4", "ITGA10", "EVI2A", "FGF18", "PRG4", "CRTAC1")
)

# 1. Add "F8_..." to your factor levels
f1_f8_order <- c(
  "F1_Superficial", 
  "F2_Universal", 
  "F2_F3_Perivascular", 
  "F3_FRC_like",
  "F4_DS_DPEP1", 
  "F4_TNN_COCH", 
  "F4_DP_HHIP",
  "F5_RAMP1", 
  "F5_NGFR",
  "F6_Inflammatory_Myofibroblast", 
  "F7_Myofibroblast",
  "F8_Fascia_like_Myofibroblast"
)

f1_f8_colors <- c(
  "F1_Superficial"                = "#F5FB5C",
  "F2_Universal"                  = "#B1C8E1",
  "F2_F3_Perivascular"            = "#7F86BF",
  "F3_FRC_like"                   = "#E5C3D9",
  "F4_DS_DPEP1"                   = "#B5E595",
  "F4_TNN_COCH"                   = "#7FC97F",
  "F4_DP_HHIP"                    = "#278F48",
  "F5_RAMP1"                      = "#9E9AC8",
  "F5_NGFR"                       = "#796EB2",
  "F6_Inflammatory_Myofibroblast" = "#f22f15",
  "F7_Myofibroblast"              = "#5599FF",
  "F8_Fascia_like_Myofibroblast"  = "#F292BF"
)

# --- MACROPHAGE DICTIONARY (unchanged) ---
macrophage_dictionary <- list(
  "M0_Non_Polarized" = c("CD68","LYZ","CSF1R", "AIF1"),
  "M1_Inflammatory"  = c("CD80", "IL6", "CXCL9", "CXCL10"),
  "M2_Wound_Healing" = c("CD163", "MRC1", "FOLR2", "CD209", "IL10", "CCL18")
)
macrophage_order <- c("M0_Non_Polarized", "M1_Inflammatory", "M2_Wound_Healing")
macrophage_colors <- c("M0_Non_Polarized" = "#A6CEE3", "M1_Inflammatory"  = "#E31A1C", "M2_Wound_Healing" = "#33A02C")

# ==============================================================================
# 3. DYNAMIC AUTO-ANNOTATION ENGINE (AddModuleScore + Cluster Averaging)
# ==============================================================================
message("Running Cluster-Level Module Score Auto-Annotation...")

# Setup variables based on celltype
if (opt$celltype == "fibroblast") {
  active_dict <- biological_dictionary
  active_colors <- f1_f8_colors
  active_order <- f1_f8_order
} else if (opt$celltype == "macrophage") {
  active_dict <- macrophage_dictionary
  active_colors <- macrophage_colors
  active_order <- macrophage_order
} else {
  active_dict <- list()
}

if (opt$celltype %in% c("fibroblast", "macrophage")) {
  # 3a. Clean the dictionary and enforce quality control (Minimum 3 gene)
  clean_dictionary <- function(dict, obj, min_genes = 3) {
    cleaned_dict <- list()
    for (prog in names(dict)) {
      valid_genes <- intersect(dict[[prog]], rownames(obj))
      if (length(valid_genes) >= min_genes) {
        cleaned_dict[[prog]] <- valid_genes
        message(sprintf(" [KEEP] %s: %d valid genes", prog, length(valid_genes)))
      } else {
        message(sprintf(" [DROP] %s: Only %d valid genes detected (requires %d)", prog, length(valid_genes), min_genes))
      }
    }
    return(cleaned_dict)
  }

  message("\nEvaluating dictionary signatures against dataset...")
  active_dict_clean <- clean_dictionary(active_dict, seu_obj)

  # Safety check: Stop if all dictionaries were dropped
  if(length(active_dict_clean) == 0) {
    message("FATAL: No signatures had enough valid genes to proceed.")
    quit(save = "no", status = 1)
  }


  # 3b. Run AddModuleScore (computes continuous scores for each signature)
  seu_obj <- AddModuleScore(
    object = seu_obj,
    features = active_dict_clean,
    name = "ProgScore_",
    seed = opt$seed
  )

  # Seurat appends numbers to the name (ProgScore_1, ProgScore_2, etc.)
  prog_cols <- paste0("ProgScore_", 1:length(active_dict_clean))
  prog_names <- names(active_dict_clean)

  # Ensure seurat_clusters exists (fallback to current Idents if missing)
  if (!"seurat_clusters" %in% colnames(seu_obj@meta.data)) {
    seu_obj$seurat_clusters <- Idents(seu_obj)
  }

  # 3c. Average the module scores per structural cluster
  score_df <- seu_obj@meta.data %>%
    dplyr::select(seurat_clusters, dplyr::all_of(prog_cols)) %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(prog_cols), mean), .groups = "drop")

  # normalize by gene set size to penalize massive sets like F6_Myofibroblast
  gene_counts <- sapply(active_dict_clean, length)

  for (i in seq_along(prog_cols)) {
    score_df[[prog_cols[i]]] <- score_df[[prog_cols[i]]] / sqrt(gene_counts[i])
  }

  message("Generating QC Module Score Heatmap and UMAPs...")

  # A. Score Heatmap (Validates cluster specificity)
  score_matrix <- as.matrix(score_df[, prog_cols])
  rownames(score_matrix) <- score_df$seurat_clusters
  colnames(score_matrix) <- prog_names # Mapped from cleaned dictionary

  pdf(file.path(dirs$master_sum, "QC_ModuleScore_Heatmap.pdf"), width = 12, height = 8)
  pheatmap(score_matrix,
           scale = "column", # Z-scores each signature to highlight peak cluster
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           main = paste(toupper(opt$celltype), "Module Score Z-Scores per Cluster"))
  dev.off()

  # Raw Score FeaturePlots (Visualizes continuous gradients on UMAP)
  # Plots the first 4 successful signatures to prevent massive messy PDFs
  myo_idx <- which(names(active_dict_clean) == "F6_Myofibroblast")

  # Plot the first 3 programs PLUS the F6_Myofibroblast program (if present)
  plot_indices <- unique(c(1, 2, 3, myo_idx))
  cols_to_plot <- prog_cols[plot_indices]

  p_scores <- FeaturePlot(seu_obj, features = cols_to_plot, ncol = 2, order = TRUE)
  save_plot(p_scores, file.path(dirs$umaps_global, "UMAP_ModuleScores_Raw"), w = 12, h = 10)


  # 3d. Assign the highest-scoring program to each cluster.
  # Blends a RELATIVE margin (top2 must be more than `uncertain_margin`
  # fraction below top1) with an ABSOLUTE floor (`uncertain_abs_margin`).
  # A purely relative margin doesn't work once gene-set-size normalization
  # shrinks scores toward zero: dividing by a tiny/near-zero top1 makes the
  # required gap unreachable (and a hard "top1 <= 0 -> Uncertain" rule
  # discards clusters outright whenever the best program's average score
  # dips at or below its background-control baseline, which AddModuleScore
  # output does routinely). The blend takes whichever margin is larger, so:
  #   - when top1 is comfortably positive and large, the relative margin
  #     dominates (as originally intended), scaling with score magnitude.
  #   - when top1 is small or near zero, the absolute floor keeps the
  #     threshold reachable instead of silently failing every cluster.
  required_gap <- function(top1_score) {
    max(opt$uncertain_abs_margin, opt$uncertain_margin * abs(top1_score))
  }

  best_labels <- apply(score_df[ , prog_cols, drop = FALSE], 1, function(x) {
    ord <- order(x, decreasing = TRUE)
    top1 <- ord[1]
    top2 <- ord[2]

    gap <- x[top1] - x[top2]

    # If the gap between best and second-best doesn't clear the blended
    # margin, the call is too close to trust - mark as Uncertain.
    if (gap < required_gap(x[top1])) {
      return("Uncertain")
    }
    return(prog_names[top1])
  })

  # 3e. Map the cluster labels back to the individual cells in the Seurat object
  cluster_to_label <- setNames(best_labels, score_df$seurat_clusters)

  seu_obj$Detailed_Label <- unname(cluster_to_label[as.character(seu_obj$seurat_clusters)])

  # Convert to factor to maintain your custom ordering, then set as active identity
  active_levels <- intersect(active_order, unique(seu_obj$Detailed_Label))
  if("Uncertain" %in% unique(seu_obj$Detailed_Label)) {
    active_levels <- c(active_levels, "Uncertain")
    # Add a dedicated color for Uncertain (e.g., light grey)
    active_colors["Uncertain"] <- "#D3D3D3"
  }
  seu_obj$Detailed_Label <- factor(seu_obj$Detailed_Label, levels = active_levels)

  Idents(seu_obj) <- "Detailed_Label"
} else {
  message("Bypassing dictionary scoring. Numbering", opt$celltype, "clusters directly...")

  # Ensure seurat_clusters exists
  if (!"seurat_clusters" %in% colnames(seu_obj@meta.data)) {
    seu_obj$seurat_clusters <- Idents(seu_obj)
  }

  # Assign raw cluster numbers as the detailed label
  seu_obj$Detailed_Label <- factor(paste0("Cluster_", seu_obj$seurat_clusters))
  Idents(seu_obj) <- "Detailed_Label"

  # Generate dynamic colors and order so downstream UMAPs don't break
  active_order <- levels(seu_obj$Detailed_Label)
  active_colors <- setNames(scales::hue_pal()(length(active_order)), active_order)
}



# --- GLOBAL UMAP ---
p_global_umap <- DimPlot(seu_obj, group.by = "Detailed_Label", cols = active_colors, label = TRUE, repel = TRUE) +
  ggtitle(paste(toupper(opt$celltype), "- Global Annotated UMAP")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

save_plot(p_global_umap, file.path(dirs$umaps_global, "UMAP_Annotated_Global"), w = 10, h = 8)
# --- PER SAMPLE UMAP ---
p_per_sample <- DimPlot(seu_obj,
                        group.by = "Detailed_Label",
                        split.by = "orig.ident2",
                        cols     = active_colors,
                        pt.size  = 0.3,
                        ncol     = 3) +
  theme_void() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) +
  guides(color = guide_legend(override.aes = list(size = 4), nrow = 2)) +
  ggtitle(paste(toupper(opt$celltype), "- Per Sample"))

save_plot(p_per_sample, file.path(dirs$umaps_sample, "UMAP_PerSample_Annotated"), w = 18, h = ceiling(length(unique(seu_obj$orig.ident2)) / 3) * 5)

# ==============================================================================
# 4. UNIVERSAL VALIDATION PLOTS (VlnPlots & DotPlots)
# ==============================================================================
message("Generating Validation Plots...")

for (category_name in names(active_dict)) {
  genes_to_plot <- active_dict[[category_name]]
  valid_genes <- intersect(genes_to_plot, rownames(seu_obj))
  if (length(valid_genes) == 0) next

  exp_data <- GetAssayData(seu_obj, assay = "RNA", layer = "data")
  valid_genes <- valid_genes[rowSums(exp_data[valid_genes, , drop = FALSE]) > 0]
  if (length(valid_genes) == 0) next

  try({
    p_dot <- DotPlot(seu_obj, features = valid_genes) + RotatedAxis() + ggtitle(paste(toupper(opt$celltype), "-", category_name))
    ggsave(file.path(dirs$val_dot, paste0("DotPlot_", category_name, ".pdf")), plot = p_dot, width = max(6, length(valid_genes)*0.5 + 2), height = 6)
  }, silent = TRUE)

  try({
    p_vln <- VlnPlot(seu_obj, features = valid_genes, stack = TRUE, flip = TRUE, cols = active_colors) + theme(legend.position = "none") + ggtitle(paste(toupper(opt$celltype), "-", category_name))
    ggsave(file.path(dirs$val_vln, paste0("VlnPlot_", category_name, ".pdf")), plot = p_vln, width = 8, height = max(6, length(valid_genes)*1.5))
  }, silent = TRUE)
}

# ==============================================================================
# 5. FIBROBLAST-SPECIFIC DOWNSTREAM TASKS (Bypassed for Macrophages)
# ==============================================================================
if (opt$celltype == "fibroblast") {
  message("Running Fibroblast-specific downstream analysis (Mucin, Highlights, Lineages, Condition labels)...")

  # Publication Master DotPlot - updated to the new lineage marker set
  paper_signature_genes <- c(
    "APCDD1", "COL18A1", "WIF1",                       # F1_Superficial
    "PI16", "CD34", "MFAP5", "DPP4",                    # F2_Universal
    "PPARG", "CXCL12",                                  # F2_F3_Bridge
    "CCL19", "CD74",                                    # F3_FRC_like
    "DPEP1", "MYL4",                                    # F4_DS_DPEP1
    "TNN", "COCH",                                      # F4_TNN_COCH
    "CORIN", "HHIP",                                    # F4_DP_HHIP
    "RAMP1", "RELN",                                     # F5_RAMP1
    "NGFR", "ITGA6",                                     # F5_NGFR
    "IL11", "CXCL8", "MMP1",                             # F6_Inflammatory_Myofibroblast
    "ACTA2", "COL8A1", "LRRC15",                         # F6_Myofibroblast
    "ACAN", "ITGA10", "PRG4"                             # F7_Fascia
  )
  valid_paper_genes <- intersect(paper_signature_genes, rownames(seu_obj))
  master_dotplot <- DotPlot(seu_obj, features = valid_paper_genes, dot.scale = 6) + theme_minimal() + RotatedAxis() + scale_color_gradientn(colors = c("lightgrey", "blue", "darkred")) + labs(title = "Fibroblast Subpopulation Signatures in PMH", x = "Key Marker Genes", y = "Identified Subclusters") + theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "italic"), axis.text.y = element_text(size = 11, face = "bold"), legend.position = "right")
  save_plot(master_dotplot, file.path(dirs$master_sum, "Publication_Master_DotPlot"), w = 14, h = 8)



  # Macro-Lineages: collapse related sub-populations for a coarser overview UMAP
  seu_obj$Macro_Lineage <- dplyr::case_when(
    seu_obj$Detailed_Label == "F1_Superficial" ~ "F1_Superficial",
    seu_obj$Detailed_Label == "F2_Universal" ~ "F2_Universal",
    seu_obj$Detailed_Label == "F2_F3_Bridge" ~ "F2/F3_Bridge",
    seu_obj$Detailed_Label == "F3_FRC_like" ~ "F3_FRC_like",
    grepl("^F4", seu_obj$Detailed_Label) ~ "F4_HairFollicle",
    grepl("^F5", seu_obj$Detailed_Label) ~ "F5_Schwann",
    seu_obj$Detailed_Label == "F6_Inflammatory_Myofibroblast" ~ "F6_Inflammatory_Myo",
    seu_obj$Detailed_Label == "F6_Myofibroblast" ~ "F6_Myofibroblast",
    seu_obj$Detailed_Label == "F7_Fascia" ~ "F7_Fascia",
    TRUE ~ "Unknown"
  )

  macro_colors <- c("F1_Superficial" = "#1F77B4", "F2_Universal" = "#2CA02C", "F2/F3_Bridge" = "#FF7F00", "F3_FRC_like" = "#9467BD", "F4_HairFollicle" = "#E31A1C", "F5_Schwann" = "#E7298A", "F6_Inflammatory_Myo" = "#17BECF", "F6_Myofibroblast" = "#7F7F7F", "F7_Fascia" = "#8B0000", "Unknown" = "#D3D3D3")

  p_macro_umap <- DimPlot(seu_obj, group.by = "Macro_Lineage", label = TRUE, repel = TRUE, cols = macro_colors) + ggtitle("Global UMAP: Major Fibroblast Lineages") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  save_plot(p_macro_umap, file.path(dirs$umaps_global, "UMAP_MacroLineage_Global"), w = 9, h = 7)

  # Condition-qualified label, e.g. "F1_Superficial_Healthy" / "F1_Superficial_Disease".
  # Built directly from the ground-truth Condition column rather than a second,
  # overlapping gene signature - see the note above biological_dictionary.
  if ("Condition" %in% colnames(seu_obj@meta.data)) {
    seu_obj$Detailed_Label_Condition <- paste0(
      as.character(seu_obj$Detailed_Label), "_",
      ifelse(seu_obj$Condition == "PMH", "Disease", "Healthy")
    )
  }
}


# ==============================================================================
# 5.4 MUCIN
# ==============================================================================


if (opt$celltype == "fibroblast") {
  # Mucin & ECM Gene Analysis
  mucin_ecm_genes <- c("MUC1", "HAS1", "HAS2", "MMP1", "MUC12", "HAS3", "VCAN", "FN1", "CEMIP", "HYAL1", "HYAL2", "CTGF", "TGFBI", "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "SPARC", "POSTN", "ACTA2", "TAGLN", "LOX", "LOXL2")
  available_mucin <- intersect(mucin_ecm_genes, rownames(seu_obj))
  if(length(available_mucin) > 0) {
    mucin_dotplot <- DotPlot(seu_obj, features = available_mucin, dot.scale = 8) + theme_minimal() + RotatedAxis() + scale_color_gradientn(colors = c("lightgrey", "blue", "darkred")) + labs(title = "Mucin & ECM Production by Fibroblast Subtype", x = "Target Genes", y = "Fibroblast Subcluster") + theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5), axis.text.x = element_text(face = "italic", color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))
    save_plot(mucin_dotplot, file.path(dirs$mucin, "Mucin_DotPlot_Summary"), w = 10, h = 7)
  }
}


# ==============================================================================
# 5.5 PROPORTIONS: PER-SAMPLE COMPOSITION
# ==============================================================================
message("Generating Sample Proportion Barplots...")

p_bar_sample <- create_proportion_barplot(
  seurat_obj = seu_obj,
  group_col = "orig.ident2",
  title = paste(toupper(opt$celltype), "Detailed Composition by Sample"),
  custom_colors = active_colors
)

save_plot(p_bar_sample, file.path(dirs$proportions, paste0(opt$celltype, "_sample_proportions")), w = 12, h = 7)

# ==============================================================================
# 6. SAVE DETAILED ANNOTATED RDS
# ==============================================================================
message("Saving detailed annotated RDS...")

final_rds_path <- file.path(opt$outdir, opt$celltype, "processed", paste0(opt$celltype, "_detailed_annotated.rds"))

# This prevents the "cannot open the connection" crash!
dir.create(dirname(final_rds_path), recursive = TRUE, showWarnings = FALSE)

saveRDS(seu_obj, final_rds_path)
message("=== Detail Annotation Complete! ===")