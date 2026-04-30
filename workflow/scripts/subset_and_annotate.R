suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(clustree)
  library(harmony)
  library(ggrepel)
  library(patchwork)
  library(SingleR)
  library(SingleCellExperiment)
})

# ==============================================================================
# 1. COMMAND LINE ARGUMENTS
# ==============================================================================
option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Path to the global integrated RDS"),
  make_option(c("-c", "--clusters"), type="character", help="Comma-separated clusters to subset (e.g., '1,3,8')"),
  make_option(c("-r", "--remove"), type="character", default=NULL, help="Comma-separated imposter clusters to remove"),
  make_option(c("-p", "--prefix"), type="character", help="Name of the cell type (e.g., 'fibroblast')"),
  make_option(c("-o", "--outdir"), type="character", help="Base output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))

target_clusters <- unlist(strsplit(opt$clusters, ","))
imposter_clusters <- if(!is.null(opt$remove)) unlist(strsplit(opt$remove, ",")) else NULL

# Setup Directories
out_base <- file.path(opt$outdir, opt$prefix)
dirs <- list(
  processed = file.path(out_base, "processed"),
  plots_raw = file.path(out_base, "plots", "1_raw_clustering"),
  plots_anno = file.path(out_base, "plots", "2_annotation"),
  tables = file.path(out_base, "tables"),
  dge = file.path(out_base, "dge")
)
lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

# Safe Plot Saver (PDF & High-Res PNG)
safe_save_plot <- function(plot_obj, base_filename, out_dir, w = 10, h = 6) {
  pdf_path <- file.path(out_dir, paste0(base_filename, ".pdf"))
  png_path <- file.path(out_dir, paste0(base_filename, ".png"))
  tryCatch({ pdf(pdf_path, width = w, height = h); print(plot_obj); dev.off() }, error = function(e) {})
  tryCatch({ png(png_path, width = w, height = h, units = "in", res = 300); print(plot_obj); dev.off() }, error = function(e) {})
}

# ==============================================================================
# 2. SUBSET & CLEAN
# ==============================================================================
message("\n=== Phase 1: Subsetting & Mathematical Cleaning ===")
seurat_obj <- readRDS(opt$input)

if (!"global_cluster" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj$global_cluster <- Idents(seurat_obj)
}

sub_obj <- subset(seurat_obj, idents = target_clusters)
DefaultAssay(sub_obj) <- "RNA"

message("Re-processing standard Seurat workflow...")
sub_obj <- NormalizeData(sub_obj, verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>% 
  RunPCA(verbose = FALSE)

sub_obj <- RunHarmony(sub_obj, group.by.vars = "orig.ident2", assay.use = "RNA", verbose = FALSE)

pcs_to_use <- min(20, ncol(Embeddings(sub_obj, "harmony")))
sub_obj <- RunUMAP(sub_obj, reduction = "harmony", dims = 1:pcs_to_use, verbose = FALSE)
sub_obj <- FindNeighbors(sub_obj, reduction = "harmony", dims = 1:pcs_to_use, verbose = FALSE)
sub_obj <- FindClusters(sub_obj, resolution = seq(0.2, 1.2, by = 0.2), verbose = FALSE)

# Clustree
p_tree <- clustree(sub_obj, prefix = "RNA_snn_res.") +
  ggtitle(paste(toupper(opt$prefix), "- Clustree Resolution Tracker")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
safe_save_plot(p_tree, paste0(opt$prefix, "_clustree"), dirs$plots_raw, w = 12, h = 10)

Idents(sub_obj) <- "RNA_snn_res.0.2"
sub_obj$seurat_clusters <- as.character(Idents(sub_obj))

# Imposter Cleaning
if (!is.null(imposter_clusters)) {
  message(paste("Cleaning out imposter clusters:", paste(imposter_clusters, collapse = ", ")))
  sub_obj <- subset(sub_obj, idents = imposter_clusters, invert = TRUE)
  sub_obj <- RunUMAP(sub_obj, reduction = "harmony", dims = 1:pcs_to_use, verbose = FALSE)
  sub_obj <- FindNeighbors(sub_obj, reduction = "harmony", dims = 1:pcs_to_use, verbose = FALSE)
  sub_obj <- FindClusters(sub_obj, resolution = 0.2, verbose = FALSE)
  Idents(sub_obj) <- "RNA_snn_res.0.2"
  sub_obj$seurat_clusters <- as.character(Idents(sub_obj))
}

# Standardize Conditions
sub_obj$Detailed_Condition <- case_when(
  grepl("HTY|UA", sub_obj$orig.ident2, ignore.case = TRUE) ~ "Healthy",
  grepl("AC", sub_obj$orig.ident2, ignore.case = TRUE) ~ "Acute",
  grepl("CH", sub_obj$orig.ident2, ignore.case = TRUE) ~ "Chronic",
  TRUE ~ "Unknown"
)
sub_obj$Condition <- ifelse(sub_obj$Detailed_Condition == "Healthy", "Healthy", "PMH")

# ==============================================================================
# 3. RAW VISUALIZATIONS & PROPORTIONS
# ==============================================================================
message("\n=== Phase 2: Raw Topology & Proportions ===")

p_umap <- DimPlot(sub_obj, group.by = "seurat_clusters", label = TRUE, label.size = 5) + 
  ggtitle(paste(toupper(opt$prefix), "- Raw Clusters"))
safe_save_plot(p_umap, "1_Raw_UMAP", dirs$plots_raw, w = 8, h = 6)

p_cond <- DimPlot(sub_obj, split.by = "Condition", label = TRUE) +
  ggtitle(paste(toupper(opt$prefix), "- Split by Condition"))
safe_save_plot(p_cond, "2_Raw_UMAP_Condition", dirs$plots_raw, w = 12, h = 6)

# Proportion Barplots
prop_sample <- as.data.frame(prop.table(table(Idents(sub_obj), sub_obj$orig.ident2), margin = 2))
colnames(prop_sample) <- c("Cluster", "Sample", "Proportion")
p_bar_sample <- ggplot(prop_sample %>% filter(Proportion > 0), aes(x = Sample, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = scales::percent) +
  ggtitle(paste(toupper(opt$prefix), "Composition by Sample"))
safe_save_plot(p_bar_sample, "3_Sample_Proportions", dirs$plots_raw, w = 10, h = 7)

# Raw DGE & Heatmap
sub_obj <- JoinLayers(sub_obj)
all_markers <- FindAllMarkers(sub_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
if (nrow(all_markers) > 0) {
  write.csv(all_markers, file.path(dirs$dge, paste0(opt$prefix, "_raw_markers.csv")), row.names = FALSE)
  top10 <- all_markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 10)
  p_heat <- DoHeatmap(sub_obj, features = unique(top10$gene)) + NoLegend() +
    ggtitle(paste(toupper(opt$prefix), "Top 10 Markers per Raw Cluster"))
  safe_save_plot(p_heat, "4_Heatmap_Top10", dirs$plots_raw, w = 15, h = 15)
}

# ==============================================================================
# 4. SINGLE-R BIOLOGICAL ANNOTATION
# ==============================================================================
if (opt$prefix == "fibroblast") {
  message("\n=== Phase 3: SingleR F1-F8 Auto-Annotation ===")
  
  biological_dictionary <- list(
    "F1_Superficial_Papillary" = c("COL13A1", "COL18A1", "COL23A1", "APCDD1", "WIF1", "NKD2", "CRABP1", "CYP26B1", "WNT5A", "COMP", "AXIN2", "RSPO1", "SFRP2", "IGFBP2", "TNFRSF21", "PDE4B", "GUCY1A1"),
    "F2_Universal_Reticular" = c("PI16", "CD34", "MFAP5", "KLF5", "DPP4", "PCOLCE2", "CTHRC1", "SLPI", "CD70", "LGR5"),
    "F2_F3_Perivascular" = c("PPARG", "CXCL12", "APOC1", "C7", "PLA2G2A", "APOE", "EFEMP1", "MYOC", "GDF10"),
    "F3_FRC_like" = c("CCL19", "CXCL12", "CH25H", "IL33", "IL15", "TNFSF13B", "VCAM1", "CD74", "HLA-DRA", "CXCL9", "ADAMDEC1", "IRF8", "HLA-DRB1"),
    "F4_HairFollicle_Associated" = c("ASPN", "COL11A1", "DPEP1", "MKX", "TNMD", "CORIN", "HHIP", "RSPO3", "LEF1", "MEF2C", "MYL4", "TNN", "COCH", "COL24A1", "RSPO4", "SLITRK6", "NRG3", "BMP7", "INHBA", "PTCH1"),
    "F5_Schwann_like" = c("SCN7A", "FMO2", "FGFBP2", "OLFML2A", "NGFR", "ITGA6", "EBF2", "RAMP1", "PEAR1", "RELN", "PLEKHA6", "IGFBP2", "SFRP1", "CDH19", "CLDN1"),
    "Myofibroblast_Disease_Signature" = c("ACTA2", "COL3A1", "COL5A1", "COL8A1", "POSTN", "CTHRC1", "LRRC15", "SFRP4", "ASPN", "RUNX2", "SCX"),
    "F6_Inflammatory_Myofibroblast" = c("IL11", "IL24", "IL7R", "CXCL5", "CXCL8", "CXCL13", "CCL11", "MMP1", "MMP3", "CSF3", "TDO2", "WWC1", "CHI3L1", "STAT4", "CCL5", "CCL3", "FAM167A"),
    "F7_Myofibroblast" = c("PIEZO2", "COL1A1", "COL3A1", "POSTN", "WNT2", "COL10A1", "LAMP5", "NRG1", "OGN", "TAGLN", "KIF26B", "ZNF469", "SULF1", "ADAM12", "CDH2", "LRRC17", "KCNMA1", "ADAM19", "CREB3L1", "CCN4", "FABP5", "C1QTNF3", "CADM1"),
    "F8_Fascia_like_Myofibroblast" = c("ACAN", "THBS4", "ITGA10", "FGF18", "PRG4", "CRTAC1"),
    "Prenatal_LTo_like" = c("CCL21", "CXCL13", "MADCAM1", "FDCSP", "TNFSF11")
  )
  
  validation_dictionary <- biological_dictionary
  validation_dictionary[["QC_Stress_Markers"]] <- c("MT2A", "MT1M", "MT1X", "HSP90AA1", "JUNB", "GADD45B", "IER3")

  ref_mat <- sapply(names(biological_dictionary), function(label) {
    genes <- intersect(biological_dictionary[[label]], rownames(sub_obj))
    vec <- setNames(rep(0, length(rownames(sub_obj))), rownames(sub_obj))
    vec[genes] <- 1
    vec
  })
  ref_sce <- SingleCellExperiment(assays = list(logcounts = as.matrix(ref_mat)))
  colData(ref_sce)$label <- colnames(ref_mat)
  
  singler_res <- SingleR(
    test = as.SingleCellExperiment(sub_obj, assay = "RNA"),
    ref = ref_sce,
    labels = ref_sce$label,
    clusters = sub_obj$seurat_clusters, 
    prune = FALSE 
  )
  
  label_map <- setNames(as.data.frame(singler_res)$labels, rownames(singler_res))
  sub_obj$SingleR_label <- unname(label_map[as.character(sub_obj$seurat_clusters)])
  Idents(sub_obj) <- "SingleR_label"
  
  p_anno_umap <- DimPlot(sub_obj, label = TRUE, repel = TRUE) + ggtitle("SingleR Auto-Annotated F1-F8 Fibroblasts")
  safe_save_plot(p_anno_umap, "1_Annotated_UMAP", dirs$plots_anno, w = 10, h = 7)

  # ==============================================================================
  # 5. HANIFFA VALIDATION LOOPS & MASTER DOTPLOT
  # ==============================================================================
  message("\n=== Phase 4: Crash-Proof Validation Plots ===")
  
  for (category_name in names(validation_dictionary)) {
    genes_to_plot <- validation_dictionary[[category_name]]
    valid_genes <- intersect(genes_to_plot, rownames(sub_obj))
    
    # Filter out 0-expression genes to prevent crashes
    exp_data <- GetAssayData(sub_obj, assay = "RNA", slot = "data")
    valid_genes <- valid_genes[rowSums(exp_data[valid_genes, , drop = FALSE]) > 0]
    
    if (length(valid_genes) == 0) next
    
    tryCatch({
      p_dot <- DotPlot(sub_obj, features = valid_genes) + RotatedAxis() + ggtitle(paste("Fibroblast -", category_name))
      safe_save_plot(p_dot, paste0("DotPlot_", category_name), dirs$plots_anno, width = max(6, length(valid_genes)*0.5 + 2), height = 6)
    }, error = function(e) {})
    
    tryCatch({
      p_vln <- VlnPlot(sub_obj, features = valid_genes, stack = TRUE, flip = TRUE) + theme(legend.position = "none") + ggtitle(paste("Fibroblast -", category_name))
      safe_save_plot(p_vln, paste0("VlnPlot_", category_name), dirs$plots_anno, width = 8, height = max(6, length(valid_genes)*1.5))
    }, error = function(e) {})
  }
  
  # Publication Master DotPlot
  paper_signature_genes <- c("APCDD1", "COL18A1", "CRABP1", "CYP26B1", "WNT5A", "PI16", "CD34", "MFAP5", "KLF5", "DPP4", "PPARG", "CXCL12", "CCL19", "CD74", "CXCL9", "DPEP1", "TNN", "CORIN", "RAMP1", "NGFR", "IL11", "CXCL8", "MMP1", "ACTA2", "COL11A1", "COMP", "LRRC15", "ACAN", "ITGA10", "THBS4")
  valid_paper_genes <- intersect(paper_signature_genes, rownames(sub_obj))
  
  expected_order <- c("F1_Superficial_Papillary", "F2_Universal_Reticular", "F2_F3_Perivascular", "F3_FRC_like", "F4_HairFollicle_Associated", "F5_Schwann_like", "F6_Inflammatory_Myofibroblast", "F7_Myofibroblast", "F8_Fascia_like_Myofibroblast", "Myofibroblast_Disease_Signature", "Prenatal_LTo_like", "QC_Stress_Markers", "Unassigned")
  ordered_levels <- intersect(expected_order, unique(as.character(Idents(sub_obj))))
  Idents(sub_obj) <- factor(Idents(sub_obj), levels = rev(ordered_levels))
  
  master_dotplot <- DotPlot(sub_obj, features = valid_paper_genes, dot.scale = 6) +
    theme_minimal() + RotatedAxis() + scale_color_gradientn(colors = c("lightgrey", "blue", "darkred")) +
    labs(title = "Fibroblast Subpopulation Signatures in PMH", x = "Key Marker Genes", y = "Identified Subclusters") +
    theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "italic"), axis.text.y = element_text(size = 11, face = "bold"), legend.position = "right")
  safe_save_plot(master_dotplot, "Publication_Master_DotPlot", dirs$plots_anno, w = 14, h = 8)

  # ==============================================================================
  # 6. TARGETED HIGHLIGHTS
  # ==============================================================================
  message("\n=== Phase 5: Specific Myofibroblast Highlights ===")
  sub_obj$Highlight_Group_1 <- case_when(
    sub_obj$SingleR_label %in% c("F7_Myofibroblast") ~ "F7_Myofibroblast",
    sub_obj$SingleR_label %in% c("F8_Fascia_like_Myofibroblast") ~ "F8_Fascia-like_Myofibroblast",
    TRUE ~ "F1-F6 (Healthy & Inflammatory)"
  )
  sub_obj$Highlight_Group_2 <- case_when(
    sub_obj$SingleR_label %in% c("F7_Myofibroblast", "F8_Fascia_like_Myofibroblast", "Myofibroblast_Disease_Signature") ~ "F7+F8 (All Terminal Myofibroblasts)",
    TRUE ~ "F1-F6 (Healthy & Inflammatory)"
  )
  
  col_1 <- c("F1-F6 (Healthy & Inflammatory)" = "#D9D9D9", "F7_Myofibroblast" = "#252525", "F8_Fascia-like_Myofibroblast" = "#810F7C")
  col_2 <- c("F1-F6 (Healthy & Inflammatory)" = "#D9D9D9", "F7+F8 (All Terminal Myofibroblasts)" = "#E31A1C")
  
  p_hi_1 <- DimPlot(sub_obj, group.by = "Highlight_Group_1", cols = col_1) + ggtitle("Myofibroblast Isolation (F7 & F8 Split)")
  safe_save_plot(p_hi_1, "UMAP_Highlight_F7_F8_Split", dirs$plots_anno, w = 9, h = 7)
  
  p_hi_2 <- DimPlot(sub_obj, group.by = "Highlight_Group_2", cols = col_2) + ggtitle("Myofibroblast Isolation (Combined)")
  safe_save_plot(p_hi_2, "UMAP_Highlight_Combined", dirs$plots_anno, w = 9, h = 7)
}

# ==============================================================================
# 7. SAVE FINAL OBJECT
# ==============================================================================
message("\n=== Phase 6: Saving Data ===")
saveRDS(sub_obj, file.path(dirs$processed, paste0(opt$prefix, "_annotated_final.rds")))
message("=== Pipeline Complete ===")