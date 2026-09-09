#!/usr/bin/env Rscript
# Usage:
#   Rscript 04_detail_annotation.R -c fibroblasts \
#     -i results/03_subsets/fibroblasts/processed/fibroblast_subset_processed.rds \
#     -o results/03_subsets
#
# CONTRACT WITH 03: expects `<celltype>_subset_processed.rds` with
# `seurat_clusters` and `Condition` metadata already set (both written by
# 03_subset_clusters.R's process_subset()). -o here MUST match 03's -o
# (default "results/03_subsets" for both) -- 05_dge.R and 06_go.R scan one
# shared <base>/<celltype>/processed/ folder, expecting BOTH
# <celltype>_subset_processed.rds (from 03) and
# <celltype>_detailed_annotated.rds (from 04) to live there together.
#
# CONTRACT WITH 05/06/07/08: writes
#   <outdir>/<celltype>/processed/<celltype>_detailed_annotated.rds
# with `Detailed_Label` as the active identity — this is what 05_dge.R,
# 07_pseudotime.R, and 08_cellchat_analysis.R all read.

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(pheatmap)
})

source("workflow/scripts/00_utils.R")

# ==============================================================================
# 1. COMMAND LINE ARGUMENTS
# ==============================================================================
option_list <- list(
  make_option(c("-i", "--input"),    type = "character", help = "Path to the subset RDS object (from 03_subset_clusters.R)"),
  make_option(c("-c", "--celltype"), type = "character", help = "The cell type being processed (e.g. 'fibroblasts' or 'macrophages')"),
  make_option(c("-o", "--outdir"),   type = "character", default = "results/03_subsets", help = "Base output directory. MUST be the SAME directory 03_subset_clusters.R wrote to -- 04 writes <outdir>/<celltype>/processed/<celltype>_detailed_annotated.rds alongside 03's <celltype>_subset_processed.rds in that same folder, and 05_dge.R/06_go.R scan that shared tree for both files."),
  make_option(c("--config"),         type = "character", default = "config/config.yaml", help = "Path to config.yaml"),
  make_option(c("--seed"),                type = "integer", default = NULL, help = "Global random seed [config: reproducibility.random_seed]"),
  make_option(c("--uncertain_margin"),    type = "double",  default = NULL, help = "Relative margin for calling a cluster Uncertain [config: subsetting.uncertain_margin]"),
  make_option(c("--uncertain_abs_margin"),type = "double",  default = NULL, help = "Absolute margin floor for calling a cluster Uncertain [config: subsetting.uncertain_abs_margin]")
)
opt <- parse_args(OptionParser(option_list = option_list))
cfg <- get_config(opt$config)

`%||%` <- function(a, b) if (is.null(a)) b else a
opt$seed                 <- opt$seed                 %||% cfg_get(cfg, "reproducibility", "random_seed", default = 42)
opt$uncertain_margin     <- opt$uncertain_margin     %||% cfg_get(cfg, "subsetting", "uncertain_margin", default = 0.10)
opt$uncertain_abs_margin <- opt$uncertain_abs_margin %||% cfg_get(cfg, "subsetting", "uncertain_abs_margin", default = 0.002)

if (is.null(opt$input) || is.null(opt$celltype))
  stop("Missing required arguments: --input, --celltype")

set.seed(opt$seed)

# Canonical output layout: results/04_detail/<celltype>/...
base_out_dir <- file.path(opt$outdir, opt$celltype)
dirs <- make_stage_dirs(base_out_dir, c(
  file.path("02_detailed_annotation", "1_UMAPs", "Global"),
  file.path("02_detailed_annotation", "1_UMAPs", "Per_Sample"),
  file.path("02_detailed_annotation", "2_Validation_Plots", "VlnPlots"),
  file.path("02_detailed_annotation", "2_Validation_Plots", "DotPlots"),
  file.path("02_detailed_annotation", "3_Summary_Plots"),
  file.path("02_detailed_annotation", "4_Mucin_ECM"),
  file.path("02_detailed_annotation", "5_Proportions"),
  "00_data",
  "logs"
))
names(dirs) <- c("umaps_global", "umaps_sample", "val_vln", "val_dot",
                  "master_sum", "mucin", "proportions", "processed", "logs")

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

create_proportion_barplot <- function(seurat_obj, group_col, title, custom_colors = NULL) {
  prop_df <- as.data.frame(prop.table(table(Idents(seurat_obj), seurat_obj[[group_col]][[1]]), margin = 2))
  colnames(prop_df) <- c("Cluster", "Group", "Proportion")

  p <- ggplot(prop_df %>% filter(Proportion > 0), aes(x = Group, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
          axis.text.y = element_text(color = "black"),
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    scale_y_continuous(labels = scales::percent) +
    ggtitle(title)

  if (!is.null(custom_colors)) p <- p + scale_fill_manual(values = custom_colors)
  p
}

message(paste("Loading subset-processed", toupper(opt$celltype), "object..."))
if (!file.exists(opt$input))
  stop("Input not found at '", opt$input, "'. Run 03_subset_clusters.R first.")
seu_obj <- readRDS(opt$input)
DefaultAssay(seu_obj) <- "RNA"

if (!"Condition" %in% colnames(seu_obj@meta.data))
  stop("Input object is missing 'Condition' metadata. It should have been set by ",
       "03_subset_clusters.R's assign_metadata().")

# ==============================================================================
# 3. DICTIONARIES & COLORS (FIBROBLAST & MACROPHAGE)
# ==============================================================================
fibroblast_dictionary <- list(
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

f1_f8_order <- c(
  "F1_Superficial", "F2_Universal", "F2_F3_Perivascular", "F3_FRC_like",
  "F4_DS_DPEP1", "F4_TNN_COCH", "F4_DP_HHIP", "F5_RAMP1", "F5_NGFR",
  "F6_Inflammatory_Myofibroblast", "F7_Myofibroblast", "F8_Fascia_like_Myofibroblast"
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

macrophage_dictionary <- list(
  # Steady-State / Healthy Skin Macrophages
  # Grounded in healthy skin quantitative proteomics and MERFISH resident profiles
  "M_Homeostatic" = c(
    "CSF1R", "MERTK", "F13A1",
    "C1QA", "C1QB", "C1QC",
    "FOLR2", "VSIG4", "C3AR1",
    "LILRB4", "MRC1", "CD163",
    "MS4A6A", "ABCA1", "GPX3" 
  ),

  # Acute GvHD / Tissue-Remodelling & Regulatory
  # Markers driving active wound healing, angiogenesis, and Treg signaling
  "M_Acute_Repair" = c(
    "CD163", "F13A1", "FOLR2",
    "MRC1", "MERTK", "VSIG4",
    "IL10", "TGFB1", "VEGFA",
    "TIMP1", "PLIN2", "CCL18",
    "CCL13", "MSR1", "MARCO",
    "LGALS9"                       # LGALS3 removed (associated with activated cDCs in human skin)
  ),

  # Chronic GvHD / Proinflammatory
  # Macrophages repolarized to a proinflammatory, interferon-responsive state
  "M_Chronic_Inflammatory" = c(
    "CCR7", "TREM1", "IL1B",
    "TNF", "CXCL8", "CXCL9",
    "CXCL10", "CXCL11", "CCL3",
    "CCL4", "IRF1", "GBP1", 
    "CD86", "HLA-DRA", "HLA-DRB1"   # GBP5 removed (strongly T-cell enriched in skin proteomics)
  ),

  # Proliferating Macrophages (Active in situ division)
  # Upregulated cell-cycle genes driving local macrophages proliferation in skin lesions
  "M_Proliferating" = c(
    "MKI67", "MYBL2", "CCND1",     # CCND1 and MYBL2 are the key paper-specific GvHD proliferation markers
    "TOP2A", "STMN1", "TYMS", 
    "PCNA", "MCM2", "MCM3", "MCM4", 
    "MCM5", "MCM6", "UBE2C", "BIRC5", 
    "CENPF", "HMGB2", "CCNB1", "CCNB2", 
    "NUSAP1"
  )
)

macrophage_order  <- c(
  "M_Homeostatic", 
  "M_Acute_Repair", 
  "M_Chronic_Inflammatory", 
  "M_Proliferating"
)

macrophage_colors <- c(
  "M_Homeostatic"      = "#A6CEE3", # Light Blue
  "M_Acute_Repair"         = "#33A02C", # Green (Repair/Resolving)
  "M_Chronic_Inflammatory" = "#E31A1C", # Red (Inflammation/Damage)
  "M_Proliferating"             = "#CAB2D6"  # Light Purple
)

t_cell_dictionary <- list(

  # Broad T-cell identity
  # Core genes for identifying conventional T cells in human skin scRNA-seq
  "T_Core" = c(
    "CD3D", "CD3E", "CD3G",
    "TRAC", "CD2", "CD247",
    "LCK", "MAL", "LTB",
    "IL32", "TRBC1", "TRBC2"
  ),


  # Naive / Central-Memory T Cells
  # Resting lymphoid-homing T cells; useful for identifying non-resident circulating-like cells
  "T_Naive_CentralMemory" = c(
    "CCR7", "SELL", "TCF7",
    "LEF1", "MAL", "LTB",
    "IL7R", "MALAT1",
    "NOSIP", "LTB", "TRBC2"
  ),


  # Skin-Resident Memory T Cells
  # Tissue-retention and skin-residency program reported in human skin and inflammatory dermatoses
  "T_Skin_TRM" = c(
    "CD69", "ITGAE", "ITGA1",
    "CXCR6", "ZNF683",
    "RUNX3", "CD44",
    "RGS1", "CCL5",
    "HOPX", "ITGB1"
  ),


  # CD4 Helper / Memory T Cells
  # Conventional CD4 helper and memory phenotype; interpret with subtype-specific modules
  "T_CD4_HelperMemory" = c(
    "CD4", "IL7R", "LTB",
    "MAL", "LTB", "CCR7",
    "LTB", "MALAT1",
    "IL32", "LTB"
  ),


  # Regulatory T Cells
  # Activated and tissue-associated Treg program
  "T_Treg_Regulatory" = c(
    "FOXP3", "IL2RA", "CTLA4",
    "TIGIT", "TNFRSF4",
    "IKZF2", "LAYN",
    "ICOS", "BATF",
    "IL7R", "TNFRSF18"
  ),


  # Th1-like / Type-1 Inflammatory T Cells
  # IFN-gamma and type-1 inflammatory program
  "T_Th1_Type1" = c(
    "TBX21", "CXCR3",
    "IFNG", "TNF",
    "CCL5", "PRDM1",
    "GZMK", "NKG7",
    "IL12RB2", "STAT4"
  ),


  # Th2 / Type-2 Skin-Inflammatory T Cells
  # Type-2 program particularly relevant to atopic dermatitis and allergic skin inflammation
  "T_Th2_Type2" = c(
    "GATA3", "CCR4",
    "IL4", "IL13",
    "IL7R", "KLRB1",
    "PTGDR2", "IL1RL1",
    "CCL17", "CCL22"
  ),


  # Th17 / Tc17 Inflammatory T Cells
  # IL-17-associated inflammatory program found in psoriasis and other skin diseases
  "T_Th17_Tc17" = c(
    "CCR6", "KLRB1",
    "RORA", "IL23R",
    "IL17A", "IL17F",
    "CCL20", "IL22",
    "CXCR6", "CXCL13",
    "IL7R", "AHR"
  ),


  # Cytotoxic CD8 T Cells
  # Effector cytotoxic program; distinguish from NK cells using CD3/TRAC expression
  "T_CD8_Cytotoxic" = c(
    "CD8A", "CD8B",
    "CCL5", "NKG7",
    "GZMK", "GZMB",
    "GZMH", "GNLY",
    "PRF1", "CTSW",
    "KLRD1", "FGFBP2"
  ),


  # GZMK-positive Memory / Effector T Cells
  # Less terminally differentiated cytotoxic-memory state
  "T_GZMK_MemoryEffector" = c(
    "GZMK", "CCL5",
    "IL7R", "LTB",
    "CD8A", "CD8B",
    "CCL4", "IL32",
    "NKG7", "CXCR6"
  ),


  # MAIT Cells
  # Require TRAV1-2 or SLC4A10 together with KLRB1; KLRB1 alone is not specific
  "T_MAIT" = c(
    "TRAV1-2", "SLC4A10",
    "KLRB1", "NCR3",
    "KLRD1", "IL7R",
    "GZMK", "CCL5",
    "NKG7", "TRBC1"
  ),


  # Gamma-Delta T Cells
  # TCR gamma-delta lineage; TRDC/TRGC genes are the most informative markers
  "T_GammaDelta" = c(
    "TRDC", "TRGC1", "TRGC2",
    "CD3D", "CD3E",
    "KLRB1", "RORA",
    "CCR6", "IL23R",
    "IL17A", "IL17F",
    "CCL5", "NKG7"
  ),


  # Activated T Cells
  # General activation and antigen-experience program
  "T_Activated" = c(
    "CD69", "IL2RA",
    "TNFRSF4", "TNFRSF9",
    "ICOS", "CD40LG",
    "HLA-DRA", "HLA-DRB1",
    "CD38", "MKI67"
  ),


  # Dysfunctional / Exhausted T Cells
  # Chronic stimulation-associated program; should not be assigned using PDCD1 alone
  "T_Exhausted_Dysfunctional" = c(
    "PDCD1", "TOX",
    "TOX2", "TIGIT",
    "LAG3", "CTLA4",
    "HAVCR2", "ENTPD1",
    "LAYN", "TNFRSF9",
    "BATF", "CXCL13"
  ),


  # Proliferating T Cells
  # Cell-cycle state rather than an independent T-cell lineage
  "T_Proliferating" = c(
    "MKI67", "TOP2A",
    "STMN1", "TYMS",
    "PCNA", "MCM2",
    "MCM3", "MCM4",
    "MCM5", "MCM6",
    "UBE2C", "BIRC5",
    "CENPF", "HMGB2",
    "CCNB1", "CCNB2",
    "NUSAP1", "TUBA1B"
  )
)

t_cell_order <- c(
  "T_Core",                     # <--- ADDED
  "T_Naive_CentralMemory",
  "T_Skin_TRM",
  "T_CD4_HelperMemory",         # <--- ADDED
  "T_Treg_Regulatory",
  "T_Th1_Type1",
  "T_Th2_Type2",
  "T_Th17_Tc17",
  "T_CD8_Cytotoxic",
  "T_GZMK_MemoryEffector",
  "T_MAIT",
  "T_GammaDelta",
  "T_Activated",
  "T_Exhausted_Dysfunctional",
  "T_Proliferating"
)

t_cell_colors <- c(
  "T_Core"                    = "#B2DF8A", # Light Green (ADDED)
  "T_Naive_CentralMemory"     = "#A6CEE3", # Light blue
  "T_Skin_TRM"                = "#1F78B4", # Blue
  "T_CD4_HelperMemory"        = "#FDBF6F", # Light orange (ADDED)
  "T_Treg_Regulatory"         = "#33A02C", # Green
  "T_Th1_Type1"               = "#E31A1C", # Red
  "T_Th2_Type2"               = "#FF7F00", # Orange
  "T_Th17_Tc17"               = "#6A3D9A", # Purple
  "T_CD8_Cytotoxic"           = "#B15928", # Brown
  "T_GZMK_MemoryEffector"     = "#FB9A99", # Salmon
  "T_MAIT"                    = "#CAB2D6", # Lavender
  "T_GammaDelta"              = "#F781BF", # Pink
  "T_Activated"               = "#FFD92F", # Yellow
  "T_Exhausted_Dysfunctional" = "#666666", # Dark gray
  "T_Proliferating"           = "#8DD3C7"  # Teal
)

# ==============================================================================
# 4. DYNAMIC AUTO-ANNOTATION ENGINE (AddModuleScore + Cluster Averaging)
# ==============================================================================
message("Running Cluster-Level Module Score Auto-Annotation...")

if (opt$celltype == "fibroblasts") {
  active_dict   <- fibroblast_dictionary
  active_colors <- f1_f8_colors
  active_order  <- f1_f8_order
} else if (opt$celltype == "macrophages") {
  active_dict   <- macrophage_dictionary
  active_colors <- macrophage_colors
  active_order  <- macrophage_order
} else if (opt$celltype %in% c("t_cells", "tcells", "T_cells")) {

  active_dict   <- t_cell_dictionary
  active_colors <- t_cell_colors
  active_order  <- t_cell_order

} else {
  active_dict <- list()
}

if (opt$celltype %in% c("fibroblasts", "macrophages", "t_cells")) {

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
    cleaned_dict
  }

  message("\nEvaluating dictionary signatures against dataset...")
  active_dict_clean <- clean_dictionary(active_dict, seu_obj)

  if (length(active_dict_clean) == 0) {
    stop("FATAL: No signatures had enough valid genes to proceed for celltype '", opt$celltype, "'.")
  }

  seu_obj <- AddModuleScore(
    object = seu_obj,
    features = active_dict_clean,
    name = "ProgScore_",
    seed = opt$seed
  )

  prog_cols  <- paste0("ProgScore_", 1:length(active_dict_clean))
  prog_names <- names(active_dict_clean)

  if (!"seurat_clusters" %in% colnames(seu_obj@meta.data)) {
    seu_obj$seurat_clusters <- Idents(seu_obj)
  }

  score_df <- seu_obj@meta.data %>%
    dplyr::select(seurat_clusters, dplyr::all_of(prog_cols)) %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(prog_cols), mean), .groups = "drop")

  message("Generating QC Module Score Heatmap and UMAPs...")

  score_matrix <- as.matrix(score_df[, prog_cols])
  rownames(score_matrix) <- score_df$seurat_clusters
  colnames(score_matrix) <- prog_names

  pdf(file.path(dirs$master_sum, "QC_ModuleScore_Heatmap.pdf"), width = 12, height = 8)
  pheatmap(score_matrix,
           scale = "column",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           main = paste(toupper(opt$celltype), "Module Score Z-Scores per Cluster"))
  dev.off()

  # Plot up to the first 4 signatures (safe against fewer than 3 programs,
  # and against a myofibroblast column that may not exist under this name).
  myo_idx <- which(names(active_dict_clean) %in% c("F6_Myofibroblast", "F7_Myofibroblast"))
  plot_indices <- unique(c(seq_len(min(3, length(prog_cols))), myo_idx))
  cols_to_plot <- prog_cols[plot_indices]

  p_scores <- FeaturePlot(seu_obj, features = cols_to_plot, ncol = 2, order = TRUE)
  save_plot(p_scores, file.path(dirs$umaps_global, "UMAP_ModuleScores_Raw"), w = 12, h = 10)

  # Blended margin: relative floor scales with score magnitude; absolute
  # floor keeps the threshold reachable when top1 is near/at zero (see
  # config.yaml comments for the full rationale).
  required_gap <- function(top1_score) {
    max(opt$uncertain_abs_margin, opt$uncertain_margin * abs(top1_score))
  }

  best_labels <- apply(score_df[ , prog_cols, drop = FALSE], 1, function(x) {
    ord  <- order(x, decreasing = TRUE)
    top1 <- ord[1]
    top2 <- ord[2]
    gap  <- x[top1] - x[top2]
    if (gap < required_gap(x[top1])) return("Uncertain")
    prog_names[top1]
  })

  cluster_to_label <- setNames(best_labels, score_df$seurat_clusters)
  seu_obj$Detailed_Label <- unname(cluster_to_label[as.character(seu_obj$seurat_clusters)])

  active_levels <- intersect(active_order, unique(seu_obj$Detailed_Label))
  if ("Uncertain" %in% unique(seu_obj$Detailed_Label)) {
    active_levels <- c(active_levels, "Uncertain")
    active_colors["Uncertain"] <- "#D3D3D3"
  }
  seu_obj$Detailed_Label <- factor(seu_obj$Detailed_Label, levels = active_levels)
  Idents(seu_obj) <- "Detailed_Label"

} else {
  message("Bypassing dictionary scoring. Numbering ", opt$celltype, " clusters directly...")

  if (!"seurat_clusters" %in% colnames(seu_obj@meta.data)) {
    seu_obj$seurat_clusters <- Idents(seu_obj)
  }

  seu_obj$Detailed_Label <- factor(paste0(tools::toTitleCase(opt$celltype), "_", seu_obj$seurat_clusters))
  Idents(seu_obj) <- "Detailed_Label"

  active_order  <- levels(seu_obj$Detailed_Label)
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

save_plot(p_per_sample, file.path(dirs$umaps_sample, "UMAP_PerSample_Annotated"),
          w = 18, h = ceiling(length(unique(seu_obj$orig.ident2)) / 3) * 5)

# ==============================================================================
# 5. UNIVERSAL VALIDATION PLOTS (VlnPlots & DotPlots)
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
    ggsave(file.path(dirs$val_dot, paste0("DotPlot_", category_name, ".pdf")), plot = p_dot,
           width = max(6, length(valid_genes) * 0.5 + 2), height = 6)
  }, silent = TRUE)

  try({
    p_vln <- VlnPlot(seu_obj, features = valid_genes, stack = TRUE, flip = TRUE, cols = active_colors) +
      theme(legend.position = "none") + ggtitle(paste(toupper(opt$celltype), "-", category_name))
    ggsave(file.path(dirs$val_vln, paste0("VlnPlot_", category_name, ".pdf")), plot = p_vln,
           width = 8, height = max(6, length(valid_genes) * 1.5))
  }, silent = TRUE)
}

# ==============================================================================
# 6. FIBROBLAST-SPECIFIC DOWNSTREAM TASKS (Bypassed for other celltypes)
# ==============================================================================
if (opt$celltype == "fibroblasts") {
  message("Running Fibroblast-specific downstream analysis (Mucin, Lineages, Condition labels)...")

  paper_signature_genes <- c(
    "APCDD1", "COL18A1", "WIF1", "PI16", "CD34", "MFAP5", "DPP4",
    "PPARG", "CXCL12", "CCL19", "CD74", "DPEP1", "MYL4", "TNN", "COCH",
    "CORIN", "HHIP", "RAMP1", "RELN", "NGFR", "ITGA6",
    "IL11", "CXCL8", "MMP1", "ACTA2", "COL8A1", "LRRC15",
    "ACAN", "ITGA10", "PRG4"
  )
  valid_paper_genes <- intersect(paper_signature_genes, rownames(seu_obj))
  if (length(valid_paper_genes) > 0) {
    master_dotplot <- DotPlot(seu_obj, features = valid_paper_genes, dot.scale = 6) +
      theme_minimal() + RotatedAxis() +
      scale_color_gradientn(colors = c("lightgrey", "blue", "darkred")) +
      labs(title = "Fibroblast Subpopulation Signatures in PMH",
           x = "Key Marker Genes", y = "Identified Subclusters") +
      theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "italic"),
            axis.text.y = element_text(size = 11, face = "bold"),
            legend.position = "right")
    save_plot(master_dotplot, file.path(dirs$master_sum, "Publication_Master_DotPlot"), w = 14, h = 8)
  }

  seu_obj$Macro_Lineage <- dplyr::case_when(
    seu_obj$Detailed_Label == "F1_Superficial" ~ "F1_Superficial",
    seu_obj$Detailed_Label == "F2_Universal" ~ "F2_Universal",
    seu_obj$Detailed_Label == "F2_F3_Perivascular" ~ "F2/F3_Bridge",
    seu_obj$Detailed_Label == "F3_FRC_like" ~ "F3_FRC_like",
    grepl("^F4", seu_obj$Detailed_Label) ~ "F4_HairFollicle",
    grepl("^F5", seu_obj$Detailed_Label) ~ "F5_Schwann",
    seu_obj$Detailed_Label == "F6_Inflammatory_Myofibroblast" ~ "F6_Inflammatory_Myo",
    seu_obj$Detailed_Label == "F7_Myofibroblast" ~ "F7_Myofibroblast",
    seu_obj$Detailed_Label == "F8_Fascia_like_Myofibroblast" ~ "F8_Fascia",
    TRUE ~ "Unknown"
  )

  macro_colors <- c("F1_Superficial" = "#1F77B4", "F2_Universal" = "#2CA02C",
                     "F2/F3_Bridge" = "#FF7F00", "F3_FRC_like" = "#9467BD",
                     "F4_HairFollicle" = "#E31A1C", "F5_Schwann" = "#E7298A",
                     "F6_Inflammatory_Myo" = "#17BECF", "F7_Myofibroblast" = "#7F7F7F",
                     "F8_Fascia" = "#8B0000", "Unknown" = "#D3D3D3")

  p_macro_umap <- DimPlot(seu_obj, group.by = "Macro_Lineage", label = TRUE, repel = TRUE, cols = macro_colors) +
    ggtitle("Global UMAP: Major Fibroblast Lineages") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  save_plot(p_macro_umap, file.path(dirs$umaps_global, "UMAP_MacroLineage_Global"), w = 9, h = 7)

  if ("Condition" %in% colnames(seu_obj@meta.data)) {
    seu_obj$Detailed_Label_Condition <- paste0(
      as.character(seu_obj$Detailed_Label), "_",
      ifelse(seu_obj$Condition == "PMH", "Disease", "Healthy")
    )
  }
}

# ==============================================================================
# 7. MUCIN & ECM ANALYSIS (fibroblasts only)
# ==============================================================================
if (opt$celltype == "fibroblasts") {
  mucin_ecm_genes <- c("MUC1", "HAS1", "HAS2", "MMP1", "MUC12", "HAS3", "VCAN", "FN1",
                        "CEMIP", "HYAL1", "HYAL2", "CTGF", "TGFBI", "COL1A1", "COL1A2",
                        "COL3A1", "COL5A1", "COL6A1", "SPARC", "POSTN", "ACTA2", "TAGLN",
                        "LOX", "LOXL2")
  available_mucin <- intersect(mucin_ecm_genes, rownames(seu_obj))
  if (length(available_mucin) > 0) {
    mucin_dotplot <- DotPlot(seu_obj, features = available_mucin, dot.scale = 8) +
      theme_minimal() + RotatedAxis() +
      scale_color_gradientn(colors = c("lightgrey", "blue", "darkred")) +
      labs(title = "Mucin & ECM Production by Fibroblast Subtype",
           x = "Target Genes", y = "Fibroblast Subcluster") +
      theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
            axis.text.x = element_text(face = "italic", color = "black", size = 12),
            axis.text.y = element_text(color = "black", size = 12))
    save_plot(mucin_dotplot, file.path(dirs$mucin, "Mucin_DotPlot_Summary"), w = 10, h = 7)
    p_mucin_feature <- FeaturePlot(seu_obj,
                features = c("HAS1", "HYAL2"),
                split.by = "Condition",
                pt.size = 0.8,
                order = TRUE,
                combine = TRUE
                ) +
      theme(
        legend.position = "right",
        legend.justification = "center"
      ) +
      labs(color = "Expression level")
    
    save_plot(p_mucin_feature, file.path(dirs$mucin, "Mucin_FeaturePlot_HAS1_HYAL2"), w = 10, h = 8)
  }
}

# ==============================================================================
# 8. PROPORTIONS: PER-SAMPLE COMPOSITION
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
# 9. SAVE DETAILED ANNOTATED RDS
# ==============================================================================
message("Saving detailed annotated RDS...")

final_rds_path <- file.path(dirs$processed, paste0(opt$celltype, "_detailed_annotated.rds"))
seu_obj@misc$pipeline_seed         <- opt$seed
seu_obj@misc$uncertain_margin      <- opt$uncertain_margin
seu_obj@misc$uncertain_abs_margin  <- opt$uncertain_abs_margin
saveRDS(seu_obj, final_rds_path)
message("Saved: ", final_rds_path)
message("=== Detail Annotation Complete! ===")
