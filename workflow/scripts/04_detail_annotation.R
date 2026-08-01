suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(pheatmap)
})

source("workflow/scripts/00_utils.R")

set.seed(42)

#Rscript /home/johan/pipeline/scRNA/skin/detail_annotation.R -c fibroblast  -i /home/johan/output/skin_pmh_harmony_sctransform2/subset_cluster/fibroblast/processed/fibroblast_subset_processed.rds -o /home/johan/output/skin_pmh_harmony_sctransform2/subset_cluster
# Rscript detail_annotation.R -c macrophage  -i /home/johan/output/skin_pmh_harmony_sctransform2/subset_cluster/macrophage/processed/macrophage_subset_processed.rds -o /home/johan/output/skin_pmh_harmony_sctransform2/subset_cluster
# Rscript detail_annotation.R -c mast_cells -i /home/johan/output/skin_pmh_harmony_sctransform2/subset_cluster/mast_cell/processed/mast_cell_subset_processed.rds -o /home/johan/output/skin_pmh_harmony_sctransform2/subset_cluster


# ==============================================================================
# 1. COMMAND LINE ARGUMENTS (SNAKEMAKE INTEGRATION)
# ==============================================================================
option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Path to the subset RDS object"),
  make_option(c("-c", "--celltype"), type="character", help="The cell type being processed (e.g., 'fibroblast' or 'macrophage')"),
  make_option(c("-o", "--outdir"), type="character", help="Base subset_cluster output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))

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
    # Removed the text overlay as it gets too cluttered with 22 detailed clusters
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
  "F1_Superficial_Healthy" = c("APCDD1", "COL18A1", "COL23A1", "COL13A1", "COMP", "NKD2", "RSPO1", "AXIN2"),
  "F1_Superficial_Disease" = c("CRABP1", "CYP26B1", "TNFRSF21", "CXCL1", "WNT5A", "COL18A1", "COL23A1", "COL13A1", "NKD2", "AXIN2", "RSPO1"),
  "F2_Reticular_Healthy" = c("CD34", "PI16", "MFAP5", "DPP4", "PCOLCE2", "LGR5", "SLPI", "CD70", "CTHRC1"),
  "F2_Reticular_Disease" = c("PI16", "DPP4", "PCOLCE2", "MFAP5", "CD70", "LGR5"),
  "F2_F3_Perivascular_Healthy" = c("CXCL12", "APOE", "EFEMP1", "APOC1", "C7", "PLA2G2A", "PPARG", "MYOC", "GDF10"),
  "F2_F3_Perivascular_Disease" = c("CXCL12", "APOE", "C7", "PLA2G2A", "EFEMP1", "GDF10", "MYOC"),
  "F3_FRC_like_Healthy" = c("CCL19", "CD74", "CH25H", "TNFSF13B", "IL33", "IRF8", "IL15", "VCAM1", "HLA-DRB1", "HLA-DRA"),
  "F3_FRC_like_Disease" = c("CCL19", "CD74", "CH25H", "TNFSF13B", "IL33", "HLA-DRA", "IRF8", "COX4I2", "RBP5", "ADAMDEC1", "CXCL9", "CXCL10", "APOE", "CXCL12"),
  "F4_HF_DPEP_Healthy" = c("DPEP1", "MYL4", "MEF2C", "COL11A1"),
  "F4_HF_DPEP_Disease" = c("MEF2C", "MYL4", "COL11A1", "POSTN", "DPEP1"),
  "F4_HF_TNN_Healthy" = c("TNN", "COCH", "TNMD", "MKX", "NRG3", "SLITRK6"),
  "F4_HF_TNN_Disease" = c("COCH", "CRABP1", "COL24A1", "RSPO4", "SLITRK6", "NRG3", "MKX", "TNMD"),
  "F4_HF_DP_Healthy" = c("CORIN", "HHIP", "BMP7", "WNT5A", "LEF1"),
  "F4_HF_DP_Disease" = c("CRABP1", "COL24A1", "RSPO4", "RSPO3", "BMP7", "WNT5A", "LEF1", "SOX18", "HHIP"),
  "F5_Schwann_NGFR_Healthy" = c("NGFR", "ITGA6", "SCN7A", "CDH19", "CLDN1", "SFRP4", "TENM2"),
  "F5_Schwann_NGFR_Disease" = c("NGFR", "TM4SF1", "SFRP4", "ANGPTL7", "ITGA6", "CDH19", "CLDN1", "EBF2", "OLFM2", "SCN7A"),
  "F5_Schwann_RAMP1_Healthy" = c("RAMP1", "RELN", "PLEKHA6", "IGFBP2", "FGFBP2", "SCN7A"),
  "F5_Schwann_RAMP1_Disease" = c("RAMP1", "IGFBP2", "RELN", "COL26A1", "PLEKHA6", "FMO2", "FGFBP2"),
  "F6_Inflammatory_Myo_Disease" = c("IL11", "IL24", "CXCL5", "CXCL6", "CXCL8", "MMP9", "WNT2", "COL10A1", "MMP1", "MMP3"),
  "F7_Myofibroblast" = c("ACTA2", "TAGLN", "CTHRC1", "RUNX2", "KIF26B", "SULF1", "ADAM12", "COL8A1", "LRRC15", "CCN4", "ASPN", "POSTN", "TNC", "COL3A1", "WNT2", "COL10A1"),
  "F8_Fascia_like_Myo_Disease" = c("ACAN", "ITGA10", "CDH2", "DPP4", "CCN3"),
  "F_Fascia_Disease" = c("ITGA10", "CCN3", "DPP4", "CDH13", "PRG4", "CRTAC1", "PCOLCE2", "LGR5")
)

f1_f8_order <- c("F1_Superficial_Healthy", "F1_Superficial_Disease", "F2_Reticular_Healthy", "F2_Reticular_Disease", "F2_F3_Perivascular_Healthy", "F2_F3_Perivascular_Disease", "F3_FRC_like_Healthy", "F3_FRC_like_Disease", "F4_HF_DPEP_Healthy", "F4_HF_DPEP_Disease", "F4_HF_TNN_Healthy", "F4_HF_TNN_Disease", "F4_HF_DP_Healthy", "F4_HF_DP_Disease", "F5_Schwann_NGFR_Healthy", "F5_Schwann_NGFR_Disease", "F5_Schwann_RAMP1_Healthy", "F5_Schwann_RAMP1_Disease", "F6_Inflammatory_Myo_Disease", "F7_Myofibroblast", "F8_Fascia_like_Myo_Disease", "F_Fascia_Disease")

f1_f8_colors <- c(
  # F1
  "F1_Superficial_Healthy"  = "#A6CEE3",
  "F1_Superficial_Disease"  = "#1F78B4",

  # F2 Reticular
  "F2_Reticular_Healthy"    = "#B2DF8A",
  "F2_Reticular_Disease"    = "#33A02C",

  # F2/F3 Perivascular
  "F2_F3_Perivascular_Healthy" = "#FDBF6F",
  "F2_F3_Perivascular_Disease" = "#FF7F00",

  # F3 FRC-like
  "F3_FRC_like_Healthy"     = "#CAB2D6",
  "F3_FRC_like_Disease"     = "#6A3D9A",

  # F4 HF DPEP
  "F4_HF_DPEP_Healthy"      = "#8DD3C7",
  "F4_HF_DPEP_Disease"      = "#1B9E77",

  # F4 HF TNN
  "F4_HF_TNN_Healthy"       = "#c59a82",
  "F4_HF_TNN_Disease"       = "#B15928",

  # F4 HF DP
  "F4_HF_DP_Healthy"        = "#f5fccd",
  "F4_HF_DP_Disease"        = "#c4c000",

  # F5 Schwann NGFR
  "F5_Schwann_NGFR_Healthy" = "#FBB4AE",
  "F5_Schwann_NGFR_Disease" = "#E7298A",

  # F5 Schwann RAMP1
  "F5_Schwann_RAMP1_Healthy" = "#D9D9D9",
  "F5_Schwann_RAMP1_Disease" = "#636363",

  # F6
  "F6_Inflammatory_Myo_Disease" = "#66A61E",

  # F7
  "F7_Myofibroblast" = "#E41A1C",

  # F8
  "F8_Fascia_like_Myo_Disease" = "#A65628",

  # Fascia
  "F_Fascia_Disease" = "#984EA3"
)

# --- MACROPHAGE DICTIONARY ---
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
    name = "ProgScore_"
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

  # normalize by gene set size to penalize massive sets like F7
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
  f7_idx <- which(names(active_dict_clean) == "F7_Myofibroblast")

  # Plot the first 3 programs PLUS the F7 program
  plot_indices <- unique(c(1, 2, 3, f7_idx))
  cols_to_plot <- prog_cols[plot_indices]

  p_scores <- FeaturePlot(seu_obj, features = cols_to_plot, ncol = 2, order = TRUE)
  save_plot(p_scores, file.path(dirs$umaps_global, "UMAP_ModuleScores_Raw"), w = 12, h = 10)


  # 3d. Assign the highest-scoring program to each cluster
  best_labels <- apply(score_df[ , prog_cols, drop = FALSE], 1, function(x) {
    ord <- order(x, decreasing = TRUE)
    top1 <- ord[1]
    top2 <- ord[2]

    #If the difference between the best and second-best score is tiny, mark as Uncertain
    if ((x[top1] - x[top2]) < 0.002) {
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
  message("Bypassing dictionary scoring. Numbering Mast Cell clusters directly...")
  
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
  message("Running Fibroblast-specific downstream analysis (Mucin, Highlights, Lineages)...")
  
  # Publication Master DotPlot
  paper_signature_genes <- c("APCDD1", "COL18A1", "CRABP1", "CYP26B1", "WNT5A", "PI16", "CD34", "MFAP5", "KLF5", "DPP4", "PPARG", "CXCL12", "CCL19", "CD74", "CXCL9", "DPEP1", "TNN", "CORIN", "RAMP1", "NGFR", "IL11", "CXCL8", "MMP1", "ACTA2", "COL11A1", "COMP", "LRRC15", "ACAN", "ITGA10", "THBS4")
  valid_paper_genes <- intersect(paper_signature_genes, rownames(seu_obj))
  master_dotplot <- DotPlot(seu_obj, features = valid_paper_genes, dot.scale = 6) + theme_minimal() + RotatedAxis() + scale_color_gradientn(colors = c("lightgrey", "blue", "darkred")) + labs(title = "Fibroblast Subpopulation Signatures in PMH", x = "Key Marker Genes", y = "Identified Subclusters") + theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "italic"), axis.text.y = element_text(size = 11, face = "bold"), legend.position = "right")
  save_plot(master_dotplot, file.path(dirs$master_sum, "Publication_Master_DotPlot"), w = 14, h = 8)



  # Macro-Lineages
  seu_obj$Macro_Lineage <- dplyr::case_when(
    grepl("^F1", seu_obj$Detailed_Label) ~ "F1_Superficial", 
    grepl("^F2_Reticular", seu_obj$Detailed_Label) ~ "F2_Reticular", 
    grepl("^F2_F3", seu_obj$Detailed_Label) ~ "F2/F3_Perivascular", 
    grepl("^F3", seu_obj$Detailed_Label) ~ "F3_FRC_like", 
    grepl("^F4", seu_obj$Detailed_Label) ~ "F4_HairFollicle", 
    grepl("^F5", seu_obj$Detailed_Label) ~ "F5_Schwann", 
    grepl("^F6", seu_obj$Detailed_Label) ~ "F6_Inflammatory_Myo", 
    grepl("^F7", seu_obj$Detailed_Label) ~ "F7_Myofibroblast", 
    grepl("^F8", seu_obj$Detailed_Label) ~ "F8_Fascia_like_Myo", 
    grepl("^F_Fascia", seu_obj$Detailed_Label) ~ "F_Fascia", 
    TRUE ~ "Unknown"
  )

  macro_colors <- c("F1_Superficial" = "#1F77B4", "F2_Reticular" = "#2CA02C", "F2/F3_Perivascular" = "#FF7F00", "F3_FRC_like" = "#9467BD", "F4_HairFollicle" = "#E31A1C", "F5_Schwann" = "#E7298A", "F6_Inflammatory_Myo" = "#17BECF", "F7_Myofibroblast" = "#7F7F7F", "F8_Fascia_like_Myo" = "#8C564B", "F_Fascia" = "#8B0000", "Unknown" = "#D3D3D3")

  p_macro_umap <- DimPlot(seu_obj, group.by = "Macro_Lineage", label = TRUE, repel = TRUE, cols = macro_colors) + ggtitle("Global UMAP: Major Fibroblast Lineages") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  save_plot(p_macro_umap, file.path(dirs$umaps_global, "UMAP_MacroLineage_Global"), w = 9, h = 7)
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