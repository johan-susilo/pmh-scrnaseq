#!/usr/bin/env Rscript
# Usage:
#   Rscript cell_annotation.R -i TN.combined_dim30.rds -s all -o /output/annotations -r 0.3
#
# Steps:
#   read_rds        – load TN.combined_dim30.rds and cache seurat_objects.rds
#   singleR         – SingleR annotation (HPCA + BlueprintEncode)
#   markers         – classical marker DotPlots
#   celliD          – CelliD / PanglaoDB annotation
#   scCATCH         – scCATCH annotation
#   consensus       – voting consensus from all methods → consensus_annotation.tsv
#   apply_labels    – apply consensus to dim30 RDS (updates it in-place) and write
#                     TN.combined_annotated.rds + annotated UMAP plots
#   combined_plots  – cluster-number UMAP / proportion plots (no annotation needed)
#   all             – read_rds → singleR → markers → celliD → scCATCH →
#                     consensus → apply_labels → combined_plots

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(scCATCH)
  library(SingleR)
  library(celldex)
  library(dplyr)
  library(tidyverse)
  library(CelliD)
  library(ggpubr)
  library(RColorBrewer)
  library(parallel)
  library(future)
  library(future.apply)
})

source("workflow/scripts/00_utils.R")

# ==============================================================================
# COMMAND-LINE INTERFACE
# ==============================================================================

option_list <- list(
  make_option(c("-i", "--rds"),        type = "character", default = NULL,
              help = "Path to TN.combined_dim30.rds"),
  make_option(c("-s", "--step"),       type = "character", default = "all",
              help = "Pipeline step: read_rds, singleR, markers, celliD, scCATCH, consensus, apply_labels, combined_plots, all"),
  make_option(c("-o", "--output"),     type = "character", default = "annotations",
              help = "Base output directory [default: annotations]"),
  make_option(c("-p", "--plots"),      type = "character", default = NULL,
              help = "Preprocessing plots directory (source for UMAP copies)"),
  make_option(c("--tissue"),           type = "character", default = "skin",
              help = "Tissue type for scCATCH [default: skin]"),
  make_option(c("-r", "--resolution"), type = "character", default = "0.2",
              help = "Clustering resolution [default: 0.2]")
)

opt    <- parse_args(OptionParser(option_list = option_list))

# Append resolution sub-folder to output base
res_folder  <- paste0("res_", opt$resolution)
output_base <- file.path(opt$output, res_folder)

dir.create(output_base, recursive = TRUE, showWarnings = FALSE)

output_dirs <- list(
  singleR        = file.path(output_base, "singleR"),
  markers        = file.path(output_base, "markers"),
  celliD         = file.path(output_base, "celliD"),
  scCATCH        = file.path(output_base, "scCATCH"),
  consensus      = file.path(output_base, "consensus"),
  annotation_plots = file.path(output_base, "annotation_plots"),
  combined_plots = file.path(output_base, "combined_plots"),
  logs           = file.path(output_base, "logs")
)
lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# LOGGING
# ==============================================================================

log_file <- file.path(output_dirs$logs,
                      paste0("annotation_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_conn <- file(log_file, open = "wt")
sink(log_conn, type = "output", split = TRUE)
sink(log_conn, type = "message")
message("Log: ", log_file, " | Started: ", Sys.time())

# Always sequential to avoid deadlocks with Seurat
plan(sequential)

# ==============================================================================
# SHARED FILE PATHS
# ==============================================================================

annotated_rds <- file.path(output_base, "TN.combined_annotated.rds")

# ==============================================================================
# HELPER UTILITIES
# ==============================================================================

# Normalise cell-type names to a canonical vocabulary.
# Handles plural/singular variants and known aliases (e.g. "M2 Macrophage" -> "Macrophages")
# so that near-duplicate votes get merged before consensus tallying, instead of
# being counted as separate competing cell types.
normalize_cell_type <- function(cell_type) {
  cell_type <- trimws(gsub("_", " ", cell_type))
  key <- tolower(cell_type)

  synonym_map <- c(
    "monocyte" = "Monocytes", "monocytes" = "Monocytes",
    "macrophage" = "Macrophages", "macrophages" = "Macrophages",
    "m1 macrophage" = "Macrophages", "m2 macrophage" = "Macrophages",
    "dc" = "Dendritic cells", "dendritic cell" = "Dendritic cells", "dendritic cells" = "Dendritic cells",
    "plasmacytoid dendritic cell" = "Dendritic cells",
    "killer cell" = "NK cells", "natural killer cell" = "NK cells",
    "nk cell" = "NK cells", "nk cells" = "NK cells",
    "t cell" = "T cells", "t cells" = "T cells",
    "cd4+ t-cells" = "T cells", "cd8+ t-cells" = "T cells",
    "regulatory t cell" = "T cells", "gamma delta t cells" = "T cells",
    "fibroblast" = "Fibroblasts", "fibroblasts" = "Fibroblasts", "myofibroblast" = "Fibroblasts",
    "keratinocyte" = "Keratinocytes", "keratinocytes" = "Keratinocytes",
    "endothelial_cells" = "Endothelial cells", "endothelial cells" = "Endothelial cells",
    "endothelial cell" = "Endothelial cells", "endothelial progenitor cell" = "Endothelial cells",
    "epithelial cell" = "Epithelial cells", "epithelial cells" = "Epithelial cells",
    "mast cell" = "Mast cells", "mast cells" = "Mast cells",
    "smooth muscle cell" = "Smooth muscle cells", "smooth muscle cells" = "Smooth muscle cells",
    "b cell" = "B cells", "b cells" = "B cells", "b cells memory" = "B cells",
    "plasma cell" = "Plasma cells", "plasma cells" = "Plasma cells",
    "adipocyte" = "Adipocytes", "adipocytes" = "Adipocytes"
  )

  if (key %in% names(synonym_map)) return(synonym_map[[key]])
  cell_type
}

tie_break_markers <- list(
  "Mast cells"      = c("TPSAB1", "TPSB2", "CPA3", "MS4A2"),
  "T cells"         = c("CD3D", "CD3E", "CD3G", "TRBC2"),
  "NK cells"        = c("NKG7", "GNLY", "KLRD1"),
  "B cells"         = c("MS4A1", "CD79A", "CD79B"),
  "Macrophages"     = c("CD68", "CD163", "MRC1"),
  "Monocytes"       = c("FCN1", "S100A8", "S100A9"),
  "Dendritic cells" = c("CD1C", "CLEC9A", "LAMP3"),
  "CMP"             = c("MPO", "ELANE", "PRTN3"),   # true bone-marrow progenitor markers, KIT excluded on purpose
  "HSC"             = c("AVP", "CRHBP", "MLLT3")
)

# Broader lineage grouping, used only to resolve near-miss ties that survive
# synonym normalization (e.g. "Monocytes" vs "Macrophages" both being myeloid).
lineage_group_map <- c(
  "Monocytes"           = "Myeloid cells",
  "Macrophages"         = "Myeloid cells",
  "Dendritic cells"     = "Myeloid cells",
  "Neutrophils"         = "Myeloid cells",
  "Mast cells"          = "Myeloid cells",
  "T cells"             = "Lymphoid cells",
  "NK cells"            = "Lymphoid cells",
  "B cells"             = "Lymphoid cells",
  "Plasma cells"        = "Lymphoid cells",
  "Fibroblasts"         = "Stromal cells",
  "Smooth muscle cells" = "Stromal cells",
  "Pericytes"           = "Stromal cells",
  "Keratinocytes"       = "Epithelial cells",
  "Epithelial cells"    = "Epithelial cells"
)

resolve_tie_with_markers <- function(seurat_obj, cluster_id, res_col, candidates,
                                      marker_panel = tie_break_markers,
                                      min_score = 0.1, min_gap = 0.15) {
  candidates <- intersect(candidates, names(marker_panel))
  if (length(candidates) < 2) return(NULL)  # nothing to arbitrate with

  cells <- WhichCells(seurat_obj, idents = cluster_id)
  if (length(cells) == 0) return(NULL)

  expr <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")

  scores <- sapply(candidates, function(ct) {
    genes <- intersect(marker_panel[[ct]], rownames(expr))
    if (length(genes) == 0) return(NA_real_)
    mean(Matrix::colMeans(expr[genes, cells, drop = FALSE]))
  })
  scores <- sort(scores[!is.na(scores)], decreasing = TRUE)
  if (length(scores) < 2) return(NULL)

  # require the winner to actually be expressed, and to clearly beat the runner-up
  if (scores[1] >= min_score && (scores[1] - scores[2]) >= min_gap) {
    return(names(scores)[1])
  }
  NULL  # stays Ambiguous — genuinely inconclusive
}

get_lineage <- function(cell_type) {
  if (cell_type %in% names(lineage_group_map)) lineage_group_map[[cell_type]] else cell_type
}

# Parse "CellType (Count, Percent%); CellType (Count, Percent%); ..." into a data.frame
parse_vote_details <- function(vote_details) {
  entries <- trimws(unlist(strsplit(vote_details, ";")))
  entries <- entries[entries != ""]
  if (length(entries) == 0) return(NULL)

  m <- regmatches(entries, regexec("^(.*)\\s*\\(([0-9]+),\\s*[0-9.]+%\\)$", entries))
  rows <- lapply(m, function(x) {
    if (length(x) < 3) return(NULL)
    data.frame(CellType = trimws(x[2]), Count = as.integer(x[3]), stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

# Re-tally votes from Vote_Details after synonym normalization (so "Mast cell" and
# "Mast cells" count as one vote, etc.), then classify each cluster as:
#   - a confident single cell type (>= min_percent of the re-tallied votes)
#   - "<Lineage> (mixed)" if the still-tied winners share a broader lineage
#   - "Low-confidence: <CellType>" if there's a single winner but below min_percent
#   - "Ambiguous" if the tied winners have nothing in common
resolve_confident_label <- function(vote_details, min_percent = 50) {
  votes <- parse_vote_details(vote_details)
  if (is.null(votes) || nrow(votes) == 0) return("Unknown")

  votes$CellType <- sapply(votes$CellType, normalize_cell_type, USE.NAMES = FALSE)
  agg <- aggregate(Count ~ CellType, data = votes, sum)
  agg <- agg[order(-agg$Count), ]

  total   <- sum(agg$Count)
  top_n   <- agg$Count[1]
  winners <- agg$CellType[agg$Count == top_n]
  percent <- 100 * top_n / total

  if (length(winners) == 1) {
    if (percent >= min_percent) return(winners[1])
    return(paste0("Low-confidence: ", winners[1]))
  }

  # Still tied after synonym merging - check if they at least share a lineage
  lineages <- unique(sapply(winners, get_lineage))
  if (length(lineages) == 1) return(paste0(lineages[1], " (mixed)"))

  "Ambiguous"
}

# Apply resolve_confident_label() across every row of the consensus table
get_confident_labels <- function(consensus_data, min_percent = 50) {
  consensus_data$Confident_Cell_Type <- sapply(
    consensus_data$Vote_Details, resolve_confident_label, min_percent = min_percent
  )
  consensus_data
}

# Resolve resolution column from a Seurat object's metadata
resolve_res_col <- function(obj, resolution) {
  md_cols <- colnames(obj@meta.data)
  res_col_sct <- paste0("SCT_snn_res.", resolution)
  res_col_rna <- paste0("RNA_snn_res.", resolution)
  if      (res_col_sct %in% md_cols) res_col_sct
  else if (res_col_rna %in% md_cols) res_col_rna
  else NULL
}

# ==============================================================================
# OBJECT LOADING — ensure seurat_objects is populated
# ==============================================================================

seurat_objects <- NULL   # global cache for the session

load_seurat_objects <- function() {
  cache_path <- file.path(output_base, "seurat_objects.rds")
  if (!is.null(seurat_objects)) return(invisible(NULL))

  if (file.exists(cache_path)) {
    message("Loading cached seurat_objects from: ", cache_path)
    seurat_objects <<- readRDS(cache_path)
    return(invisible(NULL))
  }

  if (is.null(opt$rds)) stop("--rds path must be specified")
  seurat_objects <<- read_rds(opt$rds)
  saveRDS(seurat_objects, cache_path)
}

# ==============================================================================
# STEP: READ RDS
# ==============================================================================

read_rds <- function(rds_path) {
  message("============================================================")
  message("Reading RDS: ", rds_path)
  message("============================================================")

  TN.combined <- readRDS(rds_path)

  res_col <- resolve_res_col(TN.combined, opt$resolution)
  if (!is.null(res_col)) {
    Idents(TN.combined) <- res_col
    message("Active identity set to: ", res_col)
  } else {
    message("WARNING: resolution column for ", opt$resolution, " not found; using default clusters.")
  }

  DefaultAssay(TN.combined) <- "RNA"
  Joined_TN.combined <- JoinLayers(TN.combined)

  message("Cells: ", ncol(Joined_TN.combined),
          " | Features: ", nrow(Joined_TN.combined),
          " | Clusters: ", length(unique(Idents(Joined_TN.combined))))

  list(TN.combined = TN.combined, Joined_TN.combined = Joined_TN.combined)
}

# ==============================================================================
# STEP: SINGLER
# ==============================================================================

run_singleR <- function(Joined_TN.combined) {
  message("\n============================================================")
  message("Starting SingleR Annotation")
  message("============================================================")

  if (file.exists(file.path(output_dirs$singleR, "SingleR_hpca_summary.tsv")) &&
      file.exists(file.path(output_dirs$singleR, "SingleR_bpe_summary.tsv"))) {
    message("SingleR results already exist. Skipping.")
    return(invisible(NULL))
  }

  counts <- GetAssayData(Joined_TN.combined)

  load_ref <- function(loader, fallback_args) {
    tryCatch(loader(), error = function(e) {
      message("Retrying with workaround: ", conditionMessage(e))
      do.call(loader, fallback_args)
    })
  }

  message("Running SingleR with HumanPrimaryCellAtlas...")
  hpca.se  <- load_ref(HumanPrimaryCellAtlasData,
                        list(ensembl = FALSE, cell.ont = "nonna"))
  pred.hpca <- SingleR(test = counts, ref = hpca.se, assay.type.test = 1,
                       labels = hpca.se$label.main)
  clustering.table_hpca <- table(pred.hpca@listData[["pruned.labels"]],
                                 Joined_TN.combined@active.ident)
  write.csv(clustering.table_hpca,
            file.path(output_dirs$singleR, "SingleR_hpca.csv"), row.names = TRUE)

  message("Running SingleR with BlueprintEncode...")
  bpe.se   <- load_ref(BlueprintEncodeData,
                        list(ensembl = FALSE, cell.ont = "nonna"))
  pred.bpe <- SingleR(test = counts, ref = bpe.se, assay.type.test = 1,
                      labels = bpe.se$label.main)
  clustering.table_bpe <- table(pred.bpe@listData[["pruned.labels"]],
                                Joined_TN.combined@active.ident)
  write.csv(clustering.table_bpe,
            file.path(output_dirs$singleR, "SingleR_bpe.csv"), row.names = TRUE)

  # Summarise: add annotation row (best cell type per cluster)
  summarise_singleR <- function(csv_path, tsv_path) {
    t <- read.csv(csv_path)
    rownames(t) <- t[, 1]; t <- t[, -1]
    t["annotation", ] <- rownames(t)[apply(t, 2, which.max)]
    write.table(t, tsv_path, col.names = TRUE, sep = "\t",
                row.names = TRUE, quote = FALSE)
    t
  }
  summarise_singleR(file.path(output_dirs$singleR, "SingleR_hpca.csv"),
                    file.path(output_dirs$singleR, "SingleR_hpca_summary.tsv"))
  summarise_singleR(file.path(output_dirs$singleR, "SingleR_bpe.csv"),
                    file.path(output_dirs$singleR, "SingleR_bpe_summary.tsv"))

  message("SingleR annotation completed!")
  message("============================================================\n")
}

# ==============================================================================
# STEP: CLASSICAL MARKERS
# ==============================================================================

plot_markers <- function(Joined_TN.combined) {
  message("\n============================================================")
  message("Generating Classical Marker Plots")
  message("============================================================")

  marker_sets <- list(
    "Epithelial"    = c("KRT1","KRT10","KRT5","KRT14","KRT6A","KRT16","KRT17","KRT18","KRT19","KRT7","DSP"),
    "Sweat_gland"   = c("MUCL1","PIP","AQP5"),
    "SMC"           = c("MCAM","ACTA2","MYL9","TAGLN","MYH11"),
    "Pericyte"      = c("NOTCH3","RGS5","PDGFRB","MYL9","TAGLN","MYH11"),
    "Fibroblasts"   = c("PDGFRA","DCN","LUM","POSTN","COL1A1","COL3A1","COL5A1","COL6A3","CD248"),
    "Vascular_EC"   = c("PECAM1","VWF"),
    "Lymphatic_EC"  = c("PROX1","LYVE1"),
    "T_cells"       = c("GZMK","CD3D","CD8A","CD8B","CCR7","GNLY","NKG7"),
    "NK_cells"      = c("GNLY","NKG7"),
    "B_cells"       = c("MS4A1","CD79A","SEC11C","CD79B"),
    "Plasma_cells"  = c("IGJ","MZB1","XBP1","CD79A","CD79B"),
    "Monocytes"     = c("CD14","CD68","CD163","MRC1","CSF1R","IL10RA","FCGR2A","FCGR2B","CD83","LYZ"),
    "Dendritic_cells" = c("IRF7","HLA-DRA","LYZ","S100B","CD1C"),
    "Neutrophils"   = c("ITGAX","ITGAM","FCGR2A","ANPEP"),
    "Mast_cells"    = c("ADCYAP1","CPA3","TPSAB1","VWA5A"),
    "Melanocytes"   = c("DCT","MLANA"),
    "Neuronal_cells"= c("NRXN1","SCN7A","CDH19","S100B","IGFBP5","MIA","EGFL8","NGFR","TYR"),
    "Schwann_cells" = c("NRXN1","CCN3","MPZ","PTN","S100B")
  )

  message("Processing ", length(marker_sets), " marker sets...")
  available_features <- rownames(Joined_TN.combined)

  results <- lapply(names(marker_sets), function(cell_type) {
    markers_filtered <- intersect(marker_sets[[cell_type]], available_features)
    if (length(markers_filtered) == 0) {
      message("Warning: no markers available for ", cell_type, " — skipping")
      return(NULL)
    }
    p <- DotPlot(Joined_TN.combined, features = markers_filtered,
                 cols = c("white", "darkred"), dot.scale = 8) +
      RotatedAxis() +
      labs(title = gsub("_", " ", cell_type)) +
      theme(plot.title = element_text(hjust = 0.5, size = 24))
    save_plot(p, file.path(output_dirs$markers,
                               paste0("Classical_markers_", cell_type, ".pdf")))
    cell_type
  })

  message("Marker plots generated for ",
          length(Filter(Negate(is.null), results)), " cell types")
  message("============================================================\n")
}

# ==============================================================================
# STEP: CELLID
# ==============================================================================

run_celliD <- function(seurat_object) {
  message("\n============================================================")
  message("Starting CelliD Annotation")
  message("============================================================")

  if (file.exists(file.path(output_dirs$celliD, "CelliD_PanglaoDB_summary.tsv"))) {
    message("CelliD results already exist. Skipping.")
    return(invisible(NULL))
  }

  # Downsample if very large
  if (ncol(seurat_object) > 90000) {
    message("Downsampling to 90,000 cells for CelliD...")
    seurat_object <- subset(seurat_object,
                            cells = sample(Cells(seurat_object), 90000))
  }

  message("Joining layers on downsampled object...")
  seurat_joined <- JoinLayers(seurat_object)

  message("Running MCA...")
  DefaultAssay(seurat_joined) <- "RNA"
  Baron <- RunMCA(seurat_joined, features = rownames(seurat_joined))

  message("Downloading PanglaoDB signatures...")
  panglao <- read_tsv("https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz",
                      show_col_types = FALSE)
  all_gs <- panglao %>%
    filter(str_detect(species, "Hs")) %>%
    group_by(`cell type`) %>%
    summarise(geneset = list(`official gene symbol`), .groups = "drop") %>%
    { setNames(.$geneset, .$`cell type`) } %>%
    Filter(function(x) length(x) >= 10, .)

  message("Running RunCellHGT...")
  HGT_all_gs <- RunCellHGT(Baron, pathways = all_gs, dims = 1:50)
  all_gs_prediction <- rownames(HGT_all_gs)[apply(HGT_all_gs, 2, which.max)]
  Baron$all_gs_prediction_signif <- ifelse(apply(HGT_all_gs, 2, max) > 2,
                                           all_gs_prediction, "unassigned")

  save_plot(
    DimPlot(Baron, group.by = "all_gs_prediction_signif", reduction = "umap",
            label = TRUE, label.size = 3, repel = TRUE) +
      theme(legend.text = element_text(size = 7), aspect.ratio = 1),
    file.path(output_dirs$celliD, "Baron_dimplot.pdf")
  )

  message("Summarising CelliD results...")
  clustering.table_CelliD <- table(Baron$all_gs_prediction_signif, Baron@active.ident)
  write.csv(clustering.table_CelliD,
            file.path(output_dirs$celliD, "CelliD_PanglaoDB.csv"))

  # Exclude "unassigned" when picking the best label per cluster
  table_for_summary <- clustering.table_CelliD
  if ("unassigned" %in% rownames(table_for_summary))
    table_for_summary <- table_for_summary[rownames(table_for_summary) != "unassigned", , drop = FALSE]

  annotation_row <- apply(table_for_summary, 2, function(col) {
    if (all(col == 0)) "unassigned" else rownames(table_for_summary)[which.max(col)]
  })

  summary_df <- as.data.frame.matrix(clustering.table_CelliD)
  summary_df["annotation", ] <- annotation_row
  write.table(summary_df, file.path(output_dirs$celliD, "CelliD_PanglaoDB_summary.tsv"),
              col.names = TRUE, sep = "\t", row.names = TRUE, quote = FALSE)

  message("CelliD annotation completed!")
  message("============================================================\n")
}

# ==============================================================================
# STEP: SCCATCH
# ==============================================================================

run_scCATCH <- function(TN.combined, Joined_TN.combined) {
  message("\n============================================================")
  message("Starting scCATCH Annotation")
  message("============================================================")

  if (file.exists(file.path(output_dirs$scCATCH, "scCATCH_summary.tsv"))) {
    message("scCATCH results already exist. Skipping.")
    return(invisible(NULL))
  }

  data.input <- GetAssayData(Joined_TN.combined, assay = "RNA", layer = "data")
  message("Revising gene symbols...")
  data.input <- rev_gene(data = data.input, data_type = "data",
                         species = "Human", geneinfo = geneinfo)

  labels <- Idents(TN.combined)
  obj <- createscCATCH(data = data.input, cluster = as.character(labels))

  tissue_list <- if (opt$tissue == "skin") {
    c('Adipose tissue','Blood','Peripheral blood','Bone','Cartilage',
      'Subcutaneous adipose tissue','Hair follicle','Lung','Muscle','Skin',
      'Dermis','Lymph node','Lymphoid tissue','Pluripotent stem cell',
      'Skeletal muscle','Umbilical cord blood','Plasma','Umbilical cord',
      'Spleen','Serum','Bone marrow','Placenta','Embryonic stem cell',
      'Kidney','Pancreas','Pancreatic islet','Pyloric gland',
      'Pancreatic acinar tissue')
  } else {
    c('Blood','Peripheral blood','Lymph node','Lymphoid tissue',
      'Bone marrow','Spleen')
  }

  message("Finding marker genes (tissue: ", opt$tissue, ")...")
  obj <- findmarkergene(object = obj, species = "Human", marker = cellmatch,
                        tissue = tissue_list, use_method = "1")
  obj <- findcelltype(object = obj)

  write.csv(obj@celltype, file.path(output_dirs$scCATCH, "scCATCH.csv"), row.names = FALSE)

  # Reshape to wide summary for consensus voting
  message("Processing scCATCH results for consensus...")
  data <- read.csv(file.path(output_dirs$scCATCH, "scCATCH.csv"),
                   header = TRUE, stringsAsFactors = FALSE)
  colnames(data) <- tolower(colnames(data))

  data$Cluster   <- if ("cluster"   %in% colnames(data)) data$cluster   else data[[1]]
  data$Cell_Type <- if ("cell_type" %in% colnames(data)) data$cell_type
                    else if ("celltype" %in% colnames(data)) data$celltype
                    else data[[2]]
  data$Cell_Type <- sapply(data$Cell_Type, normalize_cell_type)

  wide_data <- data %>%
    mutate(Cluster = ifelse(is.na(Cluster) | Cluster == "",
                            "", paste0("X", Cluster))) %>%
    select(Cluster, Cell_Type) %>%
    distinct() %>%
    arrange(Cluster) %>%
    pivot_wider(names_from = Cluster, values_from = Cell_Type)

  write.table(wide_data, file.path(output_dirs$scCATCH, "scCATCH_summary.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  message("scCATCH annotation completed!")
  message("============================================================\n")
}

# ==============================================================================
# STEP: CONSENSUS
# ==============================================================================

generate_consensus_annotation <- function() {
  message("\n============================================================")
  message("Generating Consensus Annotations")
  message("============================================================")

  votes_by_cluster <- list()

  add_vote <- function(cluster_id, cell_type_raw, source_name) {
    if (is.null(cell_type_raw) || is.na(cell_type_raw) || cell_type_raw == "") return(NULL)
    clean_cluster <- trimws(as.character(gsub("^X+", "", as.character(cluster_id))))
    types <- trimws(unlist(strsplit(as.character(cell_type_raw), ",")))
    types <- types[types != ""]
    types <- sapply(types, normalize_cell_type, USE.NAMES = FALSE)
    for (t in types) {
      if (is.null(t) || t == "" ||
          tolower(t) == "unassigned" || tolower(t) == "unknown") next
      vote_row <- data.frame(Source = source_name, Vote = t, stringsAsFactors = FALSE)
      if (is.null(votes_by_cluster[[clean_cluster]])) {
        votes_by_cluster[[clean_cluster]] <<- list(vote_row)
      } else {
        votes_by_cluster[[clean_cluster]][[
          length(votes_by_cluster[[clean_cluster]]) + 1]] <<- vote_row
      }
    }
  }

  safe_read_table <- function(path) {
    if (!file.exists(path)) return(NULL)
    df <- tryCatch(
      read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                 check.names = FALSE, row.names = 1),
      error = function(e) NULL
    )
    if (!is.null(df)) return(df)
    tryCatch(
      read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                 check.names = FALSE),
      error = function(e) NULL
    )
  }

  # Process SingleR / CelliD matrix files (annotation as a row or column)
  process_matrix_file <- function(path, tool_name) {
    df <- safe_read_table(path)
    if (is.null(df)) return()

    if ("annotation" %in% rownames(df)) {
      annot_row <- df[rownames(df) == "annotation", , drop = FALSE]
      for (col in colnames(annot_row)) add_vote(col, annot_row[[col]], tool_name)
      message("Added votes from ", tool_name, " (annotation row).")
      return()
    }
    if ("annotation" %in% colnames(df)) {
      for (i in seq_len(nrow(df))) {
        cid <- if (!is.null(rownames(df)) && rownames(df)[i] != "") rownames(df)[i] else df[i, 1]
        add_vote(cid, df[i, "annotation"], tool_name)
      }
      message("Added votes from ", tool_name, " (annotation column).")
      return()
    }
    # Fallback: look for "annotation" in first column
    idx <- which(tolower(df[[1]]) == "annotation")
    if (length(idx) == 1) {
      annot_row <- df[idx, -1, drop = FALSE]
      for (j in seq_along(colnames(df)[-1]))
        add_vote(colnames(df)[-1][j], annot_row[[j]], tool_name)
      message("Added votes from ", tool_name, " (annotation in first column).")
      return()
    }
    message("No annotation row/column found in ", path)
  }

  process_matrix_file(file.path(output_dirs$singleR, "SingleR_hpca_summary.tsv"), "SingleR_HPCA")
  process_matrix_file(file.path(output_dirs$singleR, "SingleR_bpe_summary.tsv"),  "SingleR_BPE")
  process_matrix_file(file.path(output_dirs$celliD,  "CelliD_PanglaoDB_summary.tsv"), "CelliD")

  sccatch_path <- file.path(output_dirs$scCATCH, "scCATCH_summary.tsv")
  if (file.exists(sccatch_path)) {
    scc_df <- safe_read_table(sccatch_path)
    if (!is.null(scc_df)) {
      for (col in colnames(scc_df)) {
        vals <- scc_df[[col]]
        vals <- vals[!is.na(vals) & vals != ""]
        for (v in vals) add_vote(col, v, "scCATCH")
      }
      message("Added votes from scCATCH.")
    }
  }

  # Tally votes
  consensus_results_list <- list()
  for (cluster in sort(names(votes_by_cluster))) {
    vote_df <- votes_by_cluster[[cluster]]
    if (is.list(vote_df) && !is.data.frame(vote_df) && length(vote_df) > 0)
      vote_df <- dplyr::bind_rows(vote_df)

    if (is.null(vote_df) || nrow(vote_df) == 0) {
      consensus_results_list[[cluster]] <- data.frame(
        Cluster = cluster, Top_Cell_Type = "Unknown",
        Count = 0L, Total_Votes = 0L, Percent = 0,
        Sources = "", Vote_Details = "", stringsAsFactors = FALSE
      )
      next
    }

    vc <- as.data.frame(table(vote_df$Vote), stringsAsFactors = FALSE)
    colnames(vc) <- c("CellType", "Votes")
    vc <- vc[order(vc$Votes, decreasing = TRUE), , drop = FALSE]

    total_votes <- sum(vc$Votes)
    max_votes   <- vc$Votes[1]
    winners     <- vc$CellType[vc$Votes == max_votes]
    top_label   <- paste(winners, collapse = "; ")
    percent     <- round(100 * max_votes / total_votes, 1)

    sources_str  <- paste(unique(vote_df$Source[vote_df$Vote %in% winners]), collapse = ", ")
    vc$Percent   <- round(100 * vc$Votes / total_votes, 1)
    vote_details <- paste(paste0(vc$CellType, " (", vc$Votes, ", ", vc$Percent, "%)"),
                          collapse = "; ")

    consensus_results_list[[cluster]] <- data.frame(
      Cluster = cluster, Top_Cell_Type = top_label,
      Count = as.integer(max_votes), Total_Votes = as.integer(total_votes),
      Percent = percent, Sources = sources_str, Vote_Details = vote_details,
      stringsAsFactors = FALSE
    )
  }

  consensus_results <- dplyr::bind_rows(consensus_results_list)
  out_path <- file.path(output_dirs$consensus, "consensus_annotation.tsv")
  write.table(consensus_results, out_path, sep = "\t", row.names = FALSE, quote = FALSE)

  message("Consensus annotation saved to: ", out_path)
  message("Top results:")
  print(head(consensus_results))
  message("============================================================\n")
  consensus_results
}

# ==============================================================================
# STEP: APPLY LABELS
# Reads the consensus TSV, stamps annotation columns onto TN.combined_dim30.rds
# (updates it in-place so downstream tools always see current labels), writes
# TN.combined_annotated.rds in the annotation output dir, generates annotated
# UMAP plots in annotation_plots/, and optionally copies the preprocessing UMAP
# plots from --plots into annotation_plots/ for side-by-side comparison.
# ==============================================================================

apply_labels <- function(TN.combined) {
  message("\n============================================================")
  message("Applying Consensus Annotations")
  message("============================================================")

  consensus_file <- file.path(output_dirs$consensus, "consensus_annotation.tsv")
  if (!file.exists(consensus_file))
    stop("consensus_annotation.tsv not found. Run the 'consensus' step first.")

  consensus_data <- read.delim(consensus_file, sep = "\t", stringsAsFactors = FALSE)
  consensus_data$Cluster <- as.character(consensus_data$Cluster)
  consensus_data <- get_confident_labels(consensus_data, min_percent = 50)

  # Resolve the active resolution
  res_col <- resolve_res_col(TN.combined, opt$resolution)
  if (!is.null(res_col)) {
    Idents(TN.combined) <- res_col
    message("Active identity set to: ", res_col)
  }

  current_ids <- levels(TN.combined)

  # Build label mappings: detailed ("C0_T cells") and clean ("T cells")
  new_names_detailed <- setNames(paste0("C", current_ids, "_Unknown"), current_ids)
  new_names_clean    <- setNames(rep("Unknown", length(current_ids)), current_ids)

  # Use the re-tallied, synonym-merged confidence column (not the raw Top_Cell_Type,
  # which can show false ties like "Mast cell" vs "Mast cells" as separate winners)
  label_col <- "Confident_Cell_Type"

  tiebreak_log <- list()
  final_ctype_by_cluster <- character(0)   # clean_id -> final (post-tiebreak) ctype

for (id in current_ids) {
  clean_id  <- gsub("^X", "", id)
  match_row <- consensus_data[consensus_data$Cluster == clean_id, ]
  if (nrow(match_row) > 0) {
    ctype <- trimws(gsub(";.*", "", match_row[[label_col]][1]))

    if (ctype == "Ambiguous" || grepl("\\(mixed\\)$", ctype)) {
      raw_winners <- trimws(strsplit(match_row$Top_Cell_Type[1], ";")[[1]])
      resolved <- resolve_tie_with_markers(TN.combined, clean_id, res_col, raw_winners)
      if (!is.null(resolved)) {
        message(sprintf("Cluster %s: marker tie-break '%s' -> '%s' (candidates: %s)",
                         clean_id, ctype, resolved, paste(raw_winners, collapse = ", ")))
        tiebreak_log[[clean_id]] <- data.frame(Cluster = clean_id, Was = ctype,
                                                Resolved = resolved,
                                                Candidates = paste(raw_winners, collapse = "; "))
        ctype <- resolved
      }
    }

    new_names_detailed[[id]]         <- paste0("C", id, "_", ctype)
    new_names_clean[[id]]            <- ctype
    final_ctype_by_cluster[clean_id] <- ctype   # remember the final (possibly tie-broken) label
  }
}

if (length(tiebreak_log) > 0) {
  write.table(dplyr::bind_rows(tiebreak_log),
              file.path(output_dirs$consensus, "marker_tiebreaks.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message(length(tiebreak_log), " cluster(s) auto-resolved via marker tie-break — see marker_tiebreaks.tsv")
}
  # --- Annotate object (clean types) ---
  TN.annotated <- RenameIdents(TN.combined, new_names_clean)

  # Add searchable metadata columns
  cluster_nums            <- gsub("^C([0-9]+)_.*", "\\1", unname(new_names_detailed[as.character(Idents(TN.combined))]))
  TN.annotated$cluster_label    <- as.character(Idents(TN.annotated))
  TN.annotated$cell_type_short  <- as.character(Idents(TN.annotated))

  # IMPORTANT: pull from the post-tiebreak lookup, NOT the raw consensus_data column,
  # so cell_type_full stays consistent with cell_type_short after marker tie-breaking.
  # (03_subset_clusters.R prefers cell_type_full when it exists, so this must carry
  # the resolved label or tie-broken clusters like Mast cells become unreachable again.)
  TN.annotated$cell_type_full   <- unname(
    final_ctype_by_cluster[gsub("^X", "", as.character(Idents(TN.combined)))])

  # --- Update TN.combined_dim30.rds in-place with annotation columns ---
  TN.combined$cluster_label   <- TN.annotated$cluster_label
  TN.combined$cell_type_short <- TN.annotated$cell_type_short
  TN.combined$cell_type_full  <- TN.annotated$cell_type_full
  saveRDS(TN.combined, opt$rds)
  message("Updated TN.combined_dim30.rds with annotation columns: ", opt$rds)

  # --- Save final annotated RDS ---
  saveRDS(TN.annotated, annotated_rds)
  message("Annotated object saved: ", annotated_rds)

  # --- Generate annotated UMAP plots ---
  n_types       <- length(unique(Idents(TN.annotated)))
  colors_clean  <- colorRampPalette(brewer.pal(min(n_types, 12), "Set3"))(n_types)

  # Detailed: unique cluster + type label
  TN.detailed <- RenameIdents(TN.combined, new_names_detailed)
  save_plot(
    DimPlot(TN.detailed, reduction = "umap", label = TRUE, repel = TRUE) +
      NoLegend() + ggtitle("UMAP — Cluster + Cell Type"),
    file.path(output_dirs$annotation_plots, "UMAP_annotated_detailed.pdf")
  )

  # Clean labeled
  save_plot(
    DimPlot(TN.annotated, reduction = "umap", label = TRUE, repel = TRUE,
            label.size = 5) +
      NoLegend() + ggtitle("UMAP — Cell Types (labeled)"),
    file.path(output_dirs$annotation_plots, "UMAP_annotated_clean_labelT.pdf")
  )

  # Clean with legend
  save_plot(
    DimPlot(TN.annotated, reduction = "umap", label = FALSE) +
      scale_color_manual(values = colors_clean) +
      ggtitle("UMAP — Cell Types (legend)"),
    file.path(output_dirs$annotation_plots, "UMAP_annotated_clean_labelF.pdf")
  )
  message("Annotated UMAP plots saved to: ", output_dirs$annotation_plots)

  # --- Cell type proportion bar plots (uses the same confident/ambiguous labels
  #     as the UMAP above, so labels stay consistent across every plot type) ---
  n_ct       <- length(unique(TN.annotated$cell_type_short))
  ct_colors  <- colorRampPalette(brewer.pal(min(n_ct, 12), "Set3"))(n_ct)

  if ("orig.ident1" %in% colnames(TN.combined@meta.data)) {
    ct_prop_ident1 <- table(TN.annotated$cell_type_short, TN.combined$orig.ident1)
    ct_prop_ident1 <- round(
      sweep(ct_prop_ident1, MARGIN = 2, STATS = colSums(ct_prop_ident1), FUN = "/") * 100, 2
    )
    write.csv(ct_prop_ident1,
              file.path(output_dirs$annotation_plots, "celltype_proportion_by_group.csv"),
              row.names = TRUE)

    save_plot(
      ggplot(as.data.frame(ct_prop_ident1) %>% setNames(c("CellType", "Group", "Freq")),
             aes(x = Group, y = Freq, fill = CellType)) +
        theme_bw(base_size = 15) +
        geom_col(position = "fill", width = 0.6) +
        scale_fill_manual(values = ct_colors) +
        labs(x = "Sample Group", y = "Proportion", fill = "Cell Type",
             title = "Cell Type Proportion by Sample Group") +
        theme(legend.text = element_text(size = 10)),
      file.path(output_dirs$annotation_plots, "celltype_proportion_by_group.pdf")
    )
  }

  if ("orig.ident2" %in% colnames(TN.combined@meta.data)) {
    ct_prop_ident2 <- table(TN.annotated$cell_type_short, TN.combined$orig.ident2)
    ct_prop_ident2 <- round(
      sweep(ct_prop_ident2, MARGIN = 2, STATS = colSums(ct_prop_ident2), FUN = "/") * 100, 2
    )
    write.csv(ct_prop_ident2,
              file.path(output_dirs$annotation_plots, "celltype_proportion_by_sample.csv"),
              row.names = TRUE)

    save_plot(
      ggplot(as.data.frame(ct_prop_ident2) %>% setNames(c("CellType", "Sample", "Freq")),
             aes(x = Sample, y = Freq, fill = CellType)) +
        theme_bw(base_size = 15) +
        geom_col(position = "fill", width = 0.6) +
        scale_fill_manual(values = ct_colors) +
        labs(x = "Sample", y = "Proportion", fill = "Cell Type",
             title = "Cell Type Proportion by Sample") +
        theme(legend.text = element_text(size = 10),
              axis.text.x = element_text(angle = 45, hjust = 1)),
      file.path(output_dirs$annotation_plots, "celltype_proportion_by_sample.pdf")
    )
  }

  message("Cell type proportion plots saved to: ", output_dirs$annotation_plots)

  # --- Copy preprocessing UMAP plots into annotation_plots/ ---
  if (!is.null(opt$plots) && dir.exists(opt$plots)) {
    message("Copying preprocessing UMAP plots from: ", opt$plots)
    umap_files <- list.files(opt$plots, pattern = "umap.*\\.(pdf|png)$",
                             full.names = TRUE, recursive = FALSE)
    if (length(umap_files) > 0) {
      file.copy(umap_files, output_dirs$annotation_plots, overwrite = TRUE)
      message("Copied ", length(umap_files), " UMAP plot(s) from preprocessing into: ",
              output_dirs$annotation_plots)
    } else {
      message("No UMAP plots found in preprocessing plots directory: ", opt$plots)
    }
  }

  # --- Annotation summary table ---
  summ <- TN.annotated@meta.data %>%
    select(cluster_label, cell_type_short, cell_type_full) %>%
    distinct() %>%
    arrange(cluster_label)
  write.table(summ, file.path(output_dirs$consensus, "annotation_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message("Annotation summary written.")

  message("===== apply_labels complete =====\n")
  TN.annotated
}

# ==============================================================================
# STEP: COMBINED PLOTS (cluster numbers only — no annotation required)
# ==============================================================================

generate_combined_plots <- function(TN.combined) {
  message("\n============================================================")
  message("Generating Combined Plots (Cluster Numbers)")
  message("============================================================")

  n_clusters     <- length(unique(Idents(TN.combined)))
  cluster_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_clusters)

  save_plot(
    DimPlot(TN.combined, reduction = "umap", label = TRUE, label.size = 5,
            repel = TRUE, pt.size = 0.8) +
      NoLegend() + ggtitle("UMAP — Cluster Numbers"),
    file.path(output_dirs$combined_plots, "TNcombined_umap_clusters_labelT.pdf")
  )
  save_plot(
    DimPlot(TN.combined, reduction = "umap", label = FALSE, pt.size = 0.8) +
      ggtitle("UMAP — Cluster Numbers (legend)"),
    file.path(output_dirs$combined_plots, "TNcombined_umap_clusters_labelF.pdf")
  )
  save_plot(
    DimPlot(TN.combined, reduction = "umap", split.by = "orig.ident1",
            label = TRUE, label.size = 3, repel = TRUE, pt.size = 0.5, ncol = 3) +
      NoLegend() + ggtitle("UMAP — Clusters (split by condition)"),
    file.path(output_dirs$combined_plots, "TNcombined_umap_clusters_splitorigident1.pdf")
  )
  save_plot(
    DimPlot(TN.combined, reduction = "umap", split.by = "orig.ident2",
            label = TRUE, label.size = 3, repel = TRUE, pt.size = 0.5, ncol = 3) +
      NoLegend() + ggtitle("UMAP — Clusters (split by sample)"),
    file.path(output_dirs$combined_plots, "TNcombined_umap_clusters_splitorigident2.pdf")
  )

  # Cluster proportion by condition (ident1)
  Cluster_prop_ident1 <- table(Idents(TN.combined), TN.combined$orig.ident1)
  Cluster_prop_ident1 <- round(
    sweep(Cluster_prop_ident1, MARGIN = 2, STATS = colSums(Cluster_prop_ident1), FUN = "/") * 100, 2
  )
  write.csv(Cluster_prop_ident1,
            file.path(output_dirs$combined_plots, "cluster_proportion_by_group.csv"),
            row.names = TRUE)

  save_plot(
    ggplot(as.data.frame(Cluster_prop_ident1) %>% setNames(c("Cluster", "Group", "Freq")),
           aes(x = Group, y = Freq, fill = Cluster)) +
      theme_bw(base_size = 15) +
      geom_col(position = "fill", width = 0.6) +
      scale_fill_manual(values = cluster_colors) +
      labs(x = "Sample Group", y = "Proportion", fill = "Cluster",
           title = "Cluster Proportion by Sample Group") +
      theme(legend.text = element_text(size = 10)),
    file.path(output_dirs$combined_plots, "cluster_proportion_by_group.pdf")
  )

  # Cluster proportion by sample (ident2)
  Cluster_prop_ident2 <- table(Idents(TN.combined), TN.combined$orig.ident2)
  Cluster_prop_ident2 <- round(
    sweep(Cluster_prop_ident2, MARGIN = 2, STATS = colSums(Cluster_prop_ident2), FUN = "/") * 100, 2
  )
  write.csv(Cluster_prop_ident2,
            file.path(output_dirs$combined_plots, "cluster_proportion_by_sample.csv"),
            row.names = TRUE)

  save_plot(
    ggplot(as.data.frame(Cluster_prop_ident2) %>% setNames(c("Cluster", "Sample", "Freq")),
           aes(x = Sample, y = Freq, fill = Cluster)) +
      theme_bw(base_size = 15) +
      geom_col(position = "fill", width = 0.6) +
      scale_fill_manual(values = cluster_colors) +
      labs(x = "Sample", y = "Proportion", fill = "Cluster",
           title = "Cluster Proportion by Sample") +
      theme(legend.text = element_text(size = 10),
            axis.text.x = element_text(angle = 45, hjust = 1)),
    file.path(output_dirs$combined_plots, "cluster_proportion_by_sample.pdf")
  )

  message("Combined plots saved to: ", output_dirs$combined_plots)
  message("============================================================\n")
}

# ==============================================================================
# PIPELINE EXECUTOR
# ==============================================================================

execute_step <- function(step) {
  switch(step,

    read_rds = {
      if (is.null(opt$rds)) stop("--rds path must be specified")
      seurat_objects <<- read_rds(opt$rds)
      saveRDS(seurat_objects, file.path(output_base, "seurat_objects.rds"))
    },

    singleR = {
      load_seurat_objects()
      run_singleR(seurat_objects$Joined_TN.combined)
    },

    markers = {
      load_seurat_objects()
      plot_markers(seurat_objects$Joined_TN.combined)
    },

    celliD = {
      load_seurat_objects()
      run_celliD(seurat_objects$TN.combined)
    },

    scCATCH = {
      load_seurat_objects()
      run_scCATCH(seurat_objects$TN.combined, seurat_objects$Joined_TN.combined)
    },

    consensus = {
      generate_consensus_annotation()
    },

    apply_labels = {
      load_seurat_objects()
      apply_labels(seurat_objects$TN.combined)
    },

    combined_plots = {
      load_seurat_objects()
      generate_combined_plots(seurat_objects$TN.combined)
    },

    all = {
      if (is.null(opt$rds)) stop("--rds path must be specified")
      seurat_objects <<- read_rds(opt$rds)
      saveRDS(seurat_objects, file.path(output_base, "seurat_objects.rds"))
      run_singleR(seurat_objects$Joined_TN.combined)
      plot_markers(seurat_objects$Joined_TN.combined)
      run_celliD(seurat_objects$TN.combined)
      run_scCATCH(seurat_objects$TN.combined, seurat_objects$Joined_TN.combined)
      generate_consensus_annotation()
      apply_labels(seurat_objects$TN.combined)
      generate_combined_plots(seurat_objects$TN.combined)
    },

    stop("Invalid step. Valid options: read_rds, singleR, markers, celliD, scCATCH, ",
         "consensus, apply_labels, combined_plots, all")
  )
}

# ==============================================================================
# MAIN
# ==============================================================================

execute_step(opt$step)
message("\nStep '", opt$step, "' completed at ", Sys.time())

sink(type = "message")
sink(type = "output")
close(log_conn)