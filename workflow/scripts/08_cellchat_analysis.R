#!/usr/bin/env Rscript
# Pairwise CellChat cross-talk analysis between any two detailed-annotated celltypes.
#
# FIX vs original cellchat_analysis.R: group membership for the "who signals
# to whom" plots used to be inferred by grepping label prefixes ("^F" for
# fibroblast, "^M" for macrophage), which only worked because those two
# specific dictionaries happen to produce F*/M*-prefixed Detailed_Label
# values. Any other celltype pair (e.g. Mast <-> keratinocyte) uses
# "Cluster_0", "Cluster_1", ... labels that match neither prefix, so the old
# script would silently produce empty group vectors and blank/erroring plots.
#
# Now group membership is taken directly from which input object a cell came
# from (recorded before the merge), and the "disease-focused" plots use the
# generic `Condition` metadata column (== "PMH") instead of grepping the
# word "Disease" out of fibroblast-specific label text.
#
# Requires the `cellchat` mamba/conda environment (kept separate from
# seurat_env due to CellChat's Matrix/igraph version pinning).
#
# Usage:
#   Rscript 08_cellchat_analysis.R \
#     --obj1 results/03_subsets/fibroblast/processed/fibroblast_detailed_annotated.rds --type1 fibroblast \
#     --obj2 results/03_subsets/macrophage/processed/macrophage_detailed_annotated.rds --type2 macrophage \
#     --outdir results/04_cellchat/fibroblast_macrophage

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(CellChat)
  library(patchwork)
  library(ggplot2)
})

source("workflow/scripts/00_utils.R")

# ==============================================================================
# 1. SETUP COMMAND LINE ARGUMENTS
# ==============================================================================
option_list <- list(
  make_option(c("--basedir"), type = "character", default = "results/03_subsets",
              help = "Base dir shared with 03/04 (default: results/03_subsets). Used to auto-derive --obj1/--obj2 from --type1/--type2 if those aren't given explicitly."),
  make_option(c("--obj1"), type = "character", default = NULL, help = "Path to first celltype's *_detailed_annotated.rds [default: <basedir>/<type1>/processed/<type1>_detailed_annotated.rds]"),
  make_option(c("--type1"), type = "character", help = "Short name for the first celltype (used for prefixes/labels, e.g. 'fibroblast')"),
  make_option(c("--obj2"), type = "character", default = NULL, help = "Path to second celltype's *_detailed_annotated.rds [default: <basedir>/<type2>/processed/<type2>_detailed_annotated.rds]"),
  make_option(c("--type2"), type = "character", help = "Short name for the second celltype (e.g. 'macrophage')"),
  make_option(c("--outdir"), type = "character", default = NULL, help = "Output directory [default: results/04_cellchat/<type1>_<type2>]"),
  make_option(c("--config"), type = "character", default = "config/config.yaml", help = "Path to config.yaml"),
  make_option(c("--min_cells"), type = "integer", default = NULL, help = "Minimum number of cells required per group [config: cellchat.min_cells]"),
  make_option(c("--pval_thresh"), type = "numeric", default = NULL, help = "P-value threshold for significant interactions [config: cellchat.pval_thresh]"),
  make_option(c("--seed"), type = "integer", default = NULL, help = "Global random seed [config: reproducibility.random_seed]")
)
opt <- parse_args(OptionParser(option_list = option_list))
cfg <- get_config(opt$config)

`%||%` <- function(a, b) if (is.null(a)) b else a
opt$min_cells   <- opt$min_cells   %||% cfg_get(cfg, "cellchat", "min_cells",   default = 10)
opt$pval_thresh <- opt$pval_thresh %||% cfg_get(cfg, "cellchat", "pval_thresh", default = 0.05)
opt$seed        <- opt$seed        %||% cfg_get(cfg, "reproducibility", "random_seed", default = 42)

if (is.null(opt$type1) || is.null(opt$type2)) {
  stop("Missing required arguments: --type1, --type2")
}

# Auto-derive obj1/obj2/outdir from the shared 03_subsets convention if not
# given explicitly, so a caller only needs to name the two celltypes.
opt$obj1   <- opt$obj1   %||% file.path(opt$basedir, opt$type1, "00_data", paste0(opt$type1, "_detailed_annotated.rds"))
opt$obj2   <- opt$obj2   %||% file.path(opt$basedir, opt$type2, "00_data", paste0(opt$type2, "_detailed_annotated.rds"))
opt$outdir <- opt$outdir %||% file.path("results", "04_cellchat", paste0(opt$type1, "_", opt$type2))

if (!file.exists(opt$obj1)) stop("--obj1 not found: ", opt$obj1, " (did 04_detail_annotation.R run for '", opt$type1, "'?)")
if (!file.exists(opt$obj2)) stop("--obj2 not found: ", opt$obj2, " (did 04_detail_annotation.R run for '", opt$type2, "'?)")

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# Short, filesystem/barcode-safe prefixes for disambiguating cell IDs after merge
prefix1 <- toupper(substr(opt$type1, 1, 4))
prefix2 <- toupper(substr(opt$type2, 1, 4))
if (prefix1 == prefix2) {
  # guard against two celltypes sharing the same 4-letter prefix
  prefix1 <- paste0(prefix1, "1")
  prefix2 <- paste0(prefix2, "2")
}

# ==============================================================================
# 2. LOAD OBJECTS
# ==============================================================================
message(sprintf("Loading objects for %s and %s...", opt$type1, opt$type2))
obj1 <- readRDS(opt$obj1)
obj2 <- readRDS(opt$obj2)

for (nm in c("obj1", "obj2")) {
  o <- get(nm)
  if (!"Detailed_Label" %in% colnames(o@meta.data)) {
    stop(sprintf("ERROR: %s object is missing 'Detailed_Label'. Use the file ending in '_detailed_annotated.rds'.", nm))
  }
  if (!"Condition" %in% colnames(o@meta.data)) {
    stop(sprintf("ERROR: %s object is missing 'Condition'. Re-run 03_subset_clusters.R / 04_detail_annotation.R.", nm))
  }
}

# ==============================================================================
# 3. SEURAT V5-NATIVE DATA EXTRACTION & MERGE
# ==============================================================================
message("Extracting raw data to build a clean combined object...")

safe_extract_counts <- function(seu) {
  assay_use <- ifelse("RNA" %in% Assays(seu), "RNA", DefaultAssay(seu))
  if (inherits(seu[[assay_use]], "Assay5")) {
    try({ seu <- JoinLayers(seu) }, silent = TRUE)
  }
  mat <- tryCatch({
    GetAssayData(seu, assay = assay_use, layer = "counts")
  }, error = function(e) {
    GetAssayData(seu, assay = assay_use, slot = "counts")
  })
  return(mat)
}

counts_1 <- safe_extract_counts(obj1)
counts_2 <- safe_extract_counts(obj2)

colnames(counts_1) <- paste0(prefix1, "_", colnames(counts_1))
colnames(counts_2) <- paste0(prefix2, "_", colnames(counts_2))

common_genes <- intersect(rownames(counts_1), rownames(counts_2))
if (length(common_genes) == 0) {
  stop("No common genes found between the two objects - cannot build a combined matrix.")
}
counts_1 <- counts_1[common_genes, ]
counts_2 <- counts_2[common_genes, ]
combined_counts <- cbind(counts_1, counts_2)

meta_1 <- obj1@meta.data
meta_2 <- obj2@meta.data
rownames(meta_1) <- paste0(prefix1, "_", rownames(meta_1))
rownames(meta_2) <- paste0(prefix2, "_", rownames(meta_2))

cols_to_keep <- c("orig.ident1", "orig.ident2", "Condition", "Detailed_Label")
cols_to_keep <- intersect(cols_to_keep, intersect(colnames(meta_1), colnames(meta_2)))
meta_1 <- meta_1[, cols_to_keep, drop = FALSE]
meta_2 <- meta_2[, cols_to_keep, drop = FALSE]

# Track which input object each cell came from - this replaces the old
# label-prefix grepping ("^F"/"^M") as the source of truth for group
# membership, so it works regardless of what the Detailed_Label values
# actually look like (dictionary-based names, "Cluster_N" fallback, etc.)
meta_1$source_celltype <- opt$type1
meta_2$source_celltype <- opt$type2

combined_meta <- rbind(meta_1, meta_2)

message("Building new combined Seurat object...")
combined <- CreateSeuratObject(counts = combined_counts, meta.data = combined_meta)
combined <- NormalizeData(combined, verbose = FALSE)

Idents(combined) <- "Detailed_Label"

rm(obj1, obj2, counts_1, counts_2)
gc()

# ==============================================================================
# 4. CELLCHAT WORKFLOW (BYPASSING INTERNAL BUGS)
# ==============================================================================
message("Initializing CellChat (bypassing Seurat v5 internal conflicts)...")
set.seed(opt$seed)

data.input <- tryCatch({
  GetAssayData(combined, assay = "RNA", layer = "data")
}, error = function(e) {
  GetAssayData(combined, assay = "RNA", slot = "data")
})
meta.data <- combined@meta.data

cellchat <- createCellChat(object = data.input, meta = meta.data, group.by = "Detailed_Label")

cellchat@DB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

message("Computing probabilities (this takes time)...")
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = opt$min_cells)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# ==============================================================================
# 5. VISUALIZATIONS - PMH DISEASE-FOCUSED
# ==============================================================================
message("Saving plots to: ", opt$outdir)

# --- Group membership from object identity, not label-prefix grepping ---
type1_groups <- as.character(unique(combined_meta$Detailed_Label[combined_meta$source_celltype == opt$type1]))
type2_groups <- as.character(unique(combined_meta$Detailed_Label[combined_meta$source_celltype == opt$type2]))

# --- "Disease" subset from the generic Condition column, not label text ---
disease_meta <- combined_meta[combined_meta$Condition == "PMH", ]
type1_disease_groups <- as.character(unique(disease_meta$Detailed_Label[disease_meta$source_celltype == opt$type1]))
type2_disease_groups <- as.character(unique(disease_meta$Detailed_Label[disease_meta$source_celltype == opt$type2]))

label1 <- tools::toTitleCase(opt$type1)
label2 <- tools::toTitleCase(opt$type2)

# --- Original broad plots (keep these) ---
pdf(file.path(opt$outdir, sprintf("%s_to_%s_Crosstalk.pdf", label1, label2)), width = 12, height = 8)
print(netVisual_bubble(cellchat, sources.use = type1_groups, targets.use = type2_groups) +
      ggtitle(sprintf("Signals: %s -> %s", label1, label2)))
dev.off()

pdf(file.path(opt$outdir, sprintf("%s_to_%s_Crosstalk.pdf", label2, label1)), width = 12, height = 8)
print(netVisual_bubble(cellchat, sources.use = type2_groups, targets.use = type1_groups) +
      ggtitle(sprintf("Signals: %s -> %s", label2, label1)))
dev.off()

disease_groups <- c(type1_disease_groups, type2_disease_groups)
disease_idx <- which(rownames(cellchat@net$count) %in% disease_groups)

if (length(disease_idx) >= 2) {
  pdf(file.path(opt$outdir, "PMH_Interaction_Network_Circle_Disease.pdf"), width = 10, height = 10)
  netVisual_circle(
    cellchat@net$count[disease_idx, disease_idx],
    weight.scale = TRUE,
    label.edge   = FALSE,
    title.name   = "Number of interactions (PMH disease groups only)"
  )
  dev.off()

  pdf(file.path(opt$outdir, "PMH_Interaction_Strength_Circle_Disease.pdf"), width = 10, height = 10)
  netVisual_circle(
    cellchat@net$weight[disease_idx, disease_idx],
    weight.scale = TRUE,
    label.edge   = FALSE,
    title.name   = "Interaction strength (PMH disease groups only)"
  )
  dev.off()
} else {
  message("   - SKIPPING disease-only circle plots: fewer than 2 disease-associated groups found.")
}

# --- PMH disease-focused plots (disease-condition clusters only, significant only) ---
if (length(type1_disease_groups) > 0 && length(type2_groups) > 0) {
  pdf(file.path(opt$outdir, sprintf("PMH_%s_Disease_to_%s.pdf", label1, label2)), width = 12, height = 8)
  print(netVisual_bubble(cellchat,
        sources.use  = type1_disease_groups,
        targets.use  = type2_groups,
        remove.isolate = TRUE,
        thresh = opt$pval_thresh) +
        ggtitle(sprintf("PMH: Disease %s -> %s", label1, label2)))
  dev.off()
}

if (length(type2_groups) > 0 && length(type1_disease_groups) > 0) {
  pdf(file.path(opt$outdir, sprintf("PMH_%s_to_%s_Disease.pdf", label2, label1)), width = 12, height = 8)
  print(netVisual_bubble(cellchat,
        sources.use  = type2_groups,
        targets.use  = type1_disease_groups,
        remove.isolate = TRUE,
        thresh = opt$pval_thresh) +
        ggtitle(sprintf("PMH: %s -> Disease %s", label2, label1)))
  dev.off()
}

# --- Top interactions table ---
df_1_to_2 <- subsetCommunication(cellchat,
  sources.use = type1_disease_groups,
  targets.use = type2_groups,
  thresh = opt$pval_thresh)
df_2_to_1 <- subsetCommunication(cellchat,
  sources.use = type2_groups,
  targets.use = type1_disease_groups,
  thresh = opt$pval_thresh)

write.csv(df_1_to_2, file.path(opt$outdir, sprintf("PMH_%sDisease_to_%s_interactions.csv", label1, label2)), row.names = FALSE)
write.csv(df_2_to_1, file.path(opt$outdir, sprintf("PMH_%s_to_%sDisease_interactions.csv", label2, label1)), row.names = FALSE)

saveRDS(cellchat, file.path(opt$outdir, sprintf("%s_%s_cellchat.rds", opt$type1, opt$type2)))
message("Analysis Complete!")