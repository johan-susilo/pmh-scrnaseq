#!/usr/bin/env Rscript
# Usage: Rscript preprocessing_clean.R -f input.csv -d /path/to/data -s all -o /path/to/output -c 8
#
# Steps:
#   read_csv   – parse input.csv and cache samples_df.rds
#   process    – QC, doublet removal, per-sample normalisation  (parallel across samples)
#   integrate  – merge + Harmony integration + multi-res clustering → TN.combined_dim30.rds
#   plot       – UMAP / heatmap / proportion plots (uses numeric cluster labels)
#   all        – read_csv → process → integrate → plot
#
# Annotation is handled by cell_annotation.R, which reads TN.combined_dim30.rds,
# applies consensus labels, updates that file in-place, and writes
# TN.combined_annotated.rds plus annotated UMAP plots.

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(DoubletFinder)
  library(dplyr)
  library(ggsci)
  library(Matrix)
  library(ggpubr)
  library(cowplot)
  library(gridExtra)
  library(gplots)
  library(ggplot2)
  library(ggnewscale)
  library(RColorBrewer)
  library(tidyr)
  library(presto)
  library(ggrepel)
  library(stringr)
  library(patchwork)
  library(scales)
  library(parallel)
  library(harmony)
  library(clustree)
})

source("workflow/scripts/00_utils.R")
RNGkind("L'Ecuyer-CMRG")
options(future.globals.maxSize = 100 * 1024^3)

# DoubletFinder inadvertently passes a 1-column data frame to order().
# This helper prevents the "cannot xtfrm data frames" crash.
xtfrm.data.frame <- function(x) {
  if (ncol(x) == 1) return(xtfrm(x[[1]]))
  stop("cannot xtfrm data frames")
}

# Load cell cycle genes
if (!exists("cc.genes")) {
  tryCatch({
    utils::data("cc.genes", package = "Seurat")
    message("Loaded cc.genes from Seurat")
  }, error = function(e) {
    message("Warning: cc.genes unavailable: ", conditionMessage(e))
    cc.genes <<- list(s.genes = character(0), g2m.genes = character(0))
  })
}

# ==============================================================================
# COMMAND-LINE INTERFACE
# ==============================================================================

option_list <- list(
  make_option(c("-f", "--file"),       type = "character", help = "Input CSV (columns: sample_names, ident1, ident2)"),
  make_option(c("-r", "--resolution"), type = "numeric",   default = 0.2,          help = "Clustering resolution [default: 0.2, matches config.yaml preprocessing.resolution]"),
  make_option(c("-d", "--datadir"),    type = "character", default = ".",           help = "Base directory with sample folders [default: .]"),
  make_option(c("-s", "--step"),       type = "character", help = "Pipeline step: read_csv, process, integrate, plot, all"),
  make_option(c("-o", "--output"),     type = "character", default = NULL,          help = "Base output directory"),
  make_option(c("-c", "--cores"),      type = "integer",   default = NULL,          help = "Cores for parallel sample processing [default: all available - 1]"),
  make_option("--doublet_rate",        type = "numeric",   default = 0.08,          help = "Expected doublet rate [default: 0.08]"),
  make_option("--min_features",        type = "integer",   default = 200,           help = "Min features per cell [default: 200]"),
  make_option("--max_features",        type = "integer",   default = 5000,          help = "Max features per cell [default: 5000]"),
  make_option("--max_mt",              type = "numeric",   default = 30,            help = "Max mitochondrial % [default: 30]"),
  make_option("--use_sct",             type = "logical",   default = TRUE,          help = "Use SCTransform normalization [default: TRUE]"),
  make_option("--batch_var", type = "character", default = "orig.ident2", help = "Batch variable for Harmony [default: orig.ident2]"),
  make_option(c("--pcs_local"), type = "integer", default = 30, help = "Number of PCs for per-sample normalization"),
  make_option(c("--pcs_global"), type = "integer", default = 50, help = "Number of PCs for global integration"),
  make_option(c("--seed"), type = "integer", default = 42, help = "Global random seed for reproducibility")
)

opt        <- parse_args(OptionParser(option_list = option_list))
res_folder <- paste0("res_", opt$resolution)

# ==============================================================================
# OUTPUT DIRECTORIES
# ==============================================================================

setup_output_dirs <- function(base) {
  dirs <- list(
    processed = file.path(base, "processed"),
    plots     = file.path(base, "plots"),
    tables    = file.path(base, "tables"),
    logs      = file.path(base, "logs")
  )
  lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
  dirs
}

output_base <- if (!is.null(opt$output) && nzchar(opt$output)) opt$output else "output"
dir.create(output_base, recursive = TRUE, showWarnings = FALSE)
output_dirs <- setup_output_dirs(output_base)

# ==============================================================================
# LOGGING
# ==============================================================================

log_file <- file.path(output_dirs$logs,
                      paste0("pipeline_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_conn <- file(log_file, open = "wt")
sink(log_conn, type = "output", split = TRUE)
sink(log_conn, type = "message")
message("Log: ", log_file, " | Started: ", Sys.time())

# ==============================================================================
# PARALLEL SETUP
# Sample processing is parallelised with parallel::mclapply.
# Integration and plotting stay sequential (single shared object; not fork-safe).
# ==============================================================================

n_cores <- if (!is.null(opt$cores)) {
  opt$cores
} else {
  max(1L, parallel::detectCores(logical = FALSE) - 1L)
}
message("Parallel sample processing: ", n_cores, " core(s)")

# ==============================================================================
# SHARED FILE PATHS
# ==============================================================================

integrated_rds <- file.path(output_base, "TN.combined_dim30.rds")

# ==============================================================================
# HELPER UTILITIES
# ==============================================================================

pick_colors <- function(n) {
  if (n <= length(mycolor)) mycolor[1:n] else colorRampPalette(mycolor)(n)
}


# Save a QC violin plot to PDF
save_qc_violin <- function(seur_obj, sample_id, file_prefix, plot_title) {
  qc_feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
  if (!all(qc_feats %in% colnames(seur_obj@meta.data)) || ncol(seur_obj) == 0) {
    message("Warning: skipping QC violin for ", sample_id, " (missing columns / empty)")
    return(invisible(NULL))
  }
  tryCatch({
    p <- create_qc_violin_plot(seur_obj, qc_feats, plot_title)
    if (!is.null(p)) {
      pdf(file_prefix, width = 15, height = 5)
      print(p)
      dev.off()
      message("QC violin saved: ", file_prefix)
    }
  }, error = function(e) {
    message("Warning: QC violin failed for ", sample_id, ": ", conditionMessage(e))
    if (length(dev.list()) > 0) tryCatch(dev.off(), error = function(e2) {})
  })
}

# Build a QC violin ggplot, grouped by cluster when available
create_qc_violin_plot <- function(seurat_obj, features, title) {
  if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    cluster_id <- "seurat_clusters"
  } else if (length(levels(Idents(seurat_obj))) > 1) {
    seurat_obj@meta.data$temp_ident <- Idents(seurat_obj)
    cluster_id <- "temp_ident"
  } else {
    cluster_id <- NULL
  }

  if (!is.null(cluster_id)) {
    qc_data <- seurat_obj@meta.data %>%
      select(all_of(c(features, cluster_id))) %>%
      mutate(cell = rownames(seurat_obj@meta.data)) %>%
      pivot_longer(cols = all_of(features), names_to = "metric", values_to = "value") %>%
      rename(Identity = all_of(cluster_id)) %>%
      mutate(Identity = as.factor(Identity))

    n      <- length(unique(qc_data$Identity))
    pal    <- if (n <= 12) brewer.pal(min(n, 12), "Paired") else brewer.pal(12, "Set3")
    colors <- colorRampPalette(pal)(n)

    ggplot(qc_data, aes(x = Identity, y = value, fill = Identity)) +
      geom_violin(trim = FALSE, scale = "width") +
      geom_jitter(size = 0.1, alpha = 0.1, width = 0.2) +
      facet_wrap(~metric, scales = "free", ncol = length(features)) +
      scale_fill_manual(values = colors) +
      theme_bw(base_size = 12) +
      theme(legend.position  = "none",
            plot.title       = element_text(hjust = 0.5, size = 16, face = "bold"),
            axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
            strip.text       = element_text(size = 12, face = "bold"),
            strip.background = element_rect(fill = "lightgray")) +
      labs(title = title, x = "Identity", y = "Value")
  } else {
    qc_data <- seurat_obj@meta.data %>%
      select(all_of(features)) %>%
      mutate(cell = rownames(seurat_obj@meta.data), Identity = "All") %>%
      pivot_longer(cols = all_of(features), names_to = "metric", values_to = "value")

    ggplot(qc_data, aes(x = Identity, y = value, fill = metric)) +
      geom_violin(trim = FALSE, scale = "width") +
      geom_jitter(size = 0.1, alpha = 0.2, width = 0.2) +
      facet_wrap(~metric, scales = "free", ncol = length(features)) +
      theme_bw(base_size = 12) +
      theme(legend.position  = "none",
            plot.title       = element_text(hjust = 0.5, size = 16, face = "bold"),
            axis.text.x      = element_text(angle = 45, hjust = 1),
            strip.text       = element_text(size = 12, face = "bold"),
            strip.background = element_rect(fill = "lightgray")) +
      labs(title = title, x = "", y = "Value")
  }
}

# Normalise + run PCA (SCTransform v2 or LogNormalize)
normalize_and_pca <- function(seur_obj, use_sct, n_pcs = 30) {
  if (isTRUE(use_sct)) {
    message("Running SCTransform v2...")
    seur_obj <- SCTransform(seur_obj, vars.to.regress = "percent.mt",
                            method = "glmGamPoi", vst.flavor = "v2", verbose = FALSE)
    seur_obj <- RunPCA(seur_obj, assay = "SCT", npcs = opt$pcs_global, verbose = FALSE)
  } else {
    message("Running LogNormalize pipeline...")
    seur_obj <- seur_obj %>%
      NormalizeData() %>%
      FindVariableFeatures() %>%
      ScaleData(vars.to.regress = "percent.mt") %>%
      RunPCA(npcs = opt$pcs_global, verbose = FALSE)
  }
  seur_obj
}

# Call JoinLayers when available (Seurat v5), silently skip otherwise
try_join_layers <- function(obj) {
  if (exists("JoinLayers", where = asNamespace("Seurat"), mode = "function")) {
    tryCatch(JoinLayers(obj), error = function(e) {
      message("Warning: JoinLayers failed: ", conditionMessage(e)); obj
    })
  } else {
    message("JoinLayers() not available; skipping")
    obj
  }
}

# ==============================================================================
# STEP 1 — READ CSV
# ==============================================================================

read_samples_csv <- function(csv_file) {
  message("Reading sample CSV: ", csv_file)
  if (!file.exists(csv_file))         stop("CSV not found: ",  csv_file, call. = FALSE)
  fi <- file.info(csv_file)
  if (is.na(fi$size) || fi$size == 0) stop("CSV is empty: ",   csv_file, call. = FALSE)

  samples_df <- tryCatch(
    read.csv(csv_file, header = TRUE, stringsAsFactors = FALSE),
    error = function(e) stop("Cannot read CSV: ", e$message, call. = FALSE)
  )
  if (nrow(samples_df) == 0) stop("CSV has no data rows: ", csv_file, call. = FALSE)

  required_cols <- c("sample_names", "ident1", "ident2")
  if (!all(required_cols %in% colnames(samples_df)))
    stop("CSV must contain: ", paste(required_cols, collapse = ", "), call. = FALSE)

  saveRDS(samples_df, file.path(output_base, "samples_df.rds"))
  message("Found ", nrow(samples_df), " samples")
  samples_df
}

# ==============================================================================
# STEP 2 — PROCESS ONE SAMPLE  (called inside mclapply)
# ==============================================================================

process_sample <- function(sample_name, sample_ident1, sample_ident2,
                           base_data_dir, out_dirs, use_sct, params) {
  tryCatch({
    output_rds      <- file.path(out_dirs$processed, paste0(sample_name, "_processed.rds"))
    sample_plot_dir <- file.path(out_dirs$plots, sample_name)
    plot_path       <- function(suffix) file.path(sample_plot_dir, paste0(sample_name, suffix))

    # -- Already done? Ensure SCT flag is present then skip --
    if (file.exists(output_rds)) {
      message("Already processed: ", sample_name)
      tryCatch({
        obj <- readRDS(output_rds)
        if (is.null(obj@misc$processed_with_sct)) {
          obj@misc$processed_with_sct <- isTRUE(use_sct)
          saveRDS(obj, output_rds)
          message("  Backfilled processed_with_sct flag for ", sample_name)
        }
      }, error = function(e)
        message("Warning: could not update RDS for ", sample_name, ": ", conditionMessage(e))
      )
      return(output_rds)
    }

    message("\n===== Processing: ", sample_name, " =====")
    dir.create(sample_plot_dir, recursive = TRUE, showWarnings = FALSE)

    # -- Create Seurat object --
    seur_obj <- CreateSeuratObject(
      counts       = Read10X(data.dir = file.path(base_data_dir, sample_name)),
      project      = sample_name,
      min.cells    = 3,
      min.features = 10
    )
    message("Dimensions: ", dim(seur_obj)[1], " features × ", dim(seur_obj)[2], " cells")
    seur_obj[["percent.mt"]] <- PercentageFeatureSet(seur_obj, pattern = "^MT-")

    save_qc_violin(seur_obj, sample_name, plot_path("_01_qc_unfiltered"),
                   paste0("QC (Unfiltered) — ", sample_name))

    # -- Pre-filter (remove empty droplets before DoubletFinder) --
    n_before <- ncol(seur_obj)
    seur_obj <- subset(seur_obj,
                       subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 50)
    message("Pre-filter: removed ", n_before - ncol(seur_obj),
            " empty droplets; keeping ", ncol(seur_obj))

    save_qc_violin(seur_obj, sample_name, plot_path("_02_qc_prefiltered"),
                   paste0("QC (Pre-filtered) — ", sample_name))

    # -- Normalisation + PCA --
    seur_obj <- normalize_and_pca(seur_obj, use_sct, n_pcs = opt$pcs_local)
    save_plot(ElbowPlot(seur_obj, ndims = opt$pcs_local), plot_path("_03_elbow"))

    # -- Clustering + UMAP (needed for DoubletFinder) --
    seur_obj <- seur_obj %>%
      FindNeighbors(dims = 1:opt$pcs_local) %>%
      FindClusters(random.seed = opt$seed) %>%
      RunUMAP(dims = 1:opt$pcs_local, seed.use = opt$seed, umap.method = "uwot", metric = "cosine")
    save_plot(DimPlot(seur_obj, reduction = "umap", label = TRUE),
                   plot_path("_04_umap_initial"))

    # -- DoubletFinder --
    message("Running DoubletFinder pK sweep...")
    sweep.res <- paramSweep(seur_obj, PCs = 1:20, sct = isTRUE(use_sct))
    bcmvn    <- find.pK(summarizeSweep(sweep.res, GT = FALSE))

    save_plot(
      ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
        geom_point() + geom_line() +
        ggtitle(paste0("pK — ", sample_name)) + theme_bw(),
      plot_path("_05_pk")
    )

    pK             <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
    homotypic.prop <- modelHomotypic(seur_obj@meta.data$seurat_clusters)
    nExp_adj       <- round(params$doublet_rate * nrow(seur_obj@meta.data) *
                            (1 - homotypic.prop))
    message("Optimal pK: ", pK, "  |  Adjusted expected doublets: ", nExp_adj)

    # New code
    seur_obj <- doubletFinder(seur_obj, PCs = 1:20, pN = 0.25, pK = pK,
                              nExp = nExp_adj, sct = isTRUE(use_sct))

    DF_col <- grep("DF.classifications", colnames(seur_obj@meta.data), value = TRUE)[1]
    if (is.na(DF_col))
      stop("DoubletFinder produced no classification column for ", sample_name)

    save_plot(
      DimPlot(seur_obj, reduction = "umap", group.by = DF_col, cols = doublet_color) +
        ggtitle(paste0("Doublets (before removal) — ", sample_name)),
      plot_path("_06_doublets_before")
    )
    message("Doublet counts:"); print(table(seur_obj@meta.data[[DF_col]]))

    seur_obj <- seur_obj[, seur_obj@meta.data[[DF_col]] == "Singlet"]
    message("Cells after doublet removal: ", ncol(seur_obj))

    save_plot(
      DimPlot(seur_obj, reduction = "umap", group.by = DF_col, cols = doublet_color) +
        ggtitle(paste0("Doublets (after removal) — ", sample_name)),
      plot_path("_07_doublets_after")
    )
    save_qc_violin(seur_obj, sample_name, plot_path("_08_qc_post_doublet"),
                   paste0("QC (Post-Doublet) — ", sample_name))

    # -- Final QC filter --
    seur_obj$orig.ident1 <- sample_ident1
    seur_obj$orig.ident2 <- sample_ident2
    n_before <- ncol(seur_obj)
    seur_obj <- subset(seur_obj,
                       subset = nFeature_RNA > params$min_features &
                                nFeature_RNA < params$max_features &
                                percent.mt   < params$max_mt)
    message("Final QC: removed ", n_before - ncol(seur_obj),
            " cells; keeping ", ncol(seur_obj))

    save_qc_violin(seur_obj, sample_name, plot_path("_09_qc_final"),
                   paste0("QC (Final) — ", sample_name))

    # -- Save --
    seur_obj@misc$processed_with_sct <- isTRUE(use_sct)
    saveRDS(seur_obj, output_rds)
    message("Saved: ", output_rds, "\n===== Done: ", sample_name, " =====\n")
    output_rds

  }, error = function(e) {
    message("\nERROR processing ", sample_name, ": ", e$message)
    NULL
  })
}

# ==============================================================================
# STEP 3 — INTEGRATE SAMPLES
# ==============================================================================

integrate_samples <- function(sample_list, chosen_res = 0.4) {
  sample_list <- Filter(Negate(is.null), sample_list)
  if (length(sample_list) == 0) stop("No valid samples for integration.")

  # Load file paths into objects
  sample_objs <- lapply(sample_list, function(x) {
    if (is.character(x) && file.exists(x)) readRDS(x) else x
  })

  # De-duplicate barcodes by prefixing with sample tag
  if (any(duplicated(unlist(lapply(sample_objs, colnames))))) {
    message("Duplicated barcodes detected — prefixing with sample tags")
    sample_objs <- lapply(seq_along(sample_objs), function(i) {
      obj <- sample_objs[[i]]
      tag <- NA
      for (src in list(
        function(o) as.character(unique(o$orig.ident2)[1]),
        function(o) as.character(unique(o$orig.ident1)[1]),
        function(o) as.character(o@project.name)
      )) {
        try({ tag <- src(obj) }, silent = TRUE)
        if (!is.null(tag) && !is.na(tag) && nchar(tag) > 0) break
      }
      if (is.na(tag)) tag <- paste0("sample", i)
      colnames(obj) <- paste0(make.names(tag), "_", colnames(obj))
      obj
    })
  }

  # Verify consistent SCT usage across samples
  used_sct <- sapply(sample_objs, function(obj) {
    flag <- try(obj@misc$processed_with_sct, silent = TRUE)
    if (!is.null(flag) && !inherits(flag, "try-error")) as.logical(flag)
    else "SCT" %in% names(obj)
  })
  if (length(unique(used_sct)) > 1)
    stop("Inconsistent SCTransform usage across samples. Reprocess with a consistent --use_sct.")

  inferred_sct <- isTRUE(unique(used_sct))
  if (inferred_sct != isTRUE(opt$use_sct)) {
    message("Overriding opt$use_sct → ", inferred_sct, " to match processed samples")
    opt$use_sct <<- inferred_sct
  }

  message("\n===== Starting Harmony Integration =====")

  # Strip per-sample SCT assays — global SCTransform will be re-run on the merge
  sample_objs <- lapply(sample_objs, function(obj) {
    DefaultAssay(obj) <- "RNA"
    if ("SCT" %in% names(obj)) obj[["SCT"]] <- NULL
    obj
  })

  # Merge
  TN.combined <- if (length(sample_objs) > 1)
    merge(sample_objs[[1]], y = sample_objs[-1])
  else
    sample_objs[[1]]

  if (!inherits(TN.combined, "Seurat"))
    stop("Merged result is not a Seurat object. Check your input samples.")
  TN.combined <- try_join_layers(TN.combined)

  # Cell cycle scoring on the RNA assay before global normalisation
  message("Cell cycle scoring on merged object...")
  DefaultAssay(TN.combined) <- "RNA"
  tryCatch(TN.combined <- NormalizeData(TN.combined, verbose = FALSE),
           error = function(e) message("Warning: NormalizeData failed: ", conditionMessage(e)))

  if (exists("cc.genes") && length(cc.genes$s.genes) > 0 && length(cc.genes$g2m.genes) > 0) {
    s_g   <- intersect(cc.genes$s.genes,   rownames(TN.combined))
    g2m_g <- intersect(cc.genes$g2m.genes, rownames(TN.combined))
    if (length(s_g) > 0 && length(g2m_g) > 0) {
      tryCatch({
        TN.combined <- CellCycleScoring(TN.combined, s.features = s_g,
                                        g2m.features = g2m_g, set.ident = FALSE)
        message("CellCycleScoring complete")
      }, error = function(e) message("Warning: CellCycleScoring failed: ", conditionMessage(e)))
    } else {
      message("Skipping CellCycleScoring: cc.genes not present in merged object")
    }
  } else {
    message("Skipping CellCycleScoring: cc.genes unavailable")
  }

  # Batch column for Harmony
  if (!(opt$batch_var %in% colnames(TN.combined@meta.data))) {
    fb <- if ("orig.ident2" %in% colnames(TN.combined@meta.data)) "orig.ident2" else "orig.ident1"
    message("batch_var '", opt$batch_var, "' not found; using '", fb, "'")
    TN.combined$batch <- factor(as.character(TN.combined@meta.data[[fb]]))
  } else {
    TN.combined$batch <- TN.combined[[opt$batch_var]]
  }

  # Global normalisation + PCA
  if (isTRUE(opt$use_sct)) {
    vars_reg <- if (all(c("S.Score", "G2M.Score") %in% colnames(TN.combined@meta.data)))
      c("percent.mt", "S.Score", "G2M.Score") else "percent.mt"
    message("Global SCTransform v2 (regressing: ", paste(vars_reg, collapse = ", "), ")")
    TN.combined <- SCTransform(TN.combined, assay = "RNA", vars.to.regress = vars_reg,
                               method = "glmGamPoi", vst.flavor = "v2", verbose = FALSE)
    TN.combined <- PrepSCTFindMarkers(TN.combined, assay = "SCT", verbose = FALSE)
    TN.combined <- RunPCA(TN.combined, assay = "SCT", npcs = opt$pcs_global, seed.use = opt$seed, verbose = FALSE)
    harmony_assay <- "SCT"
  } else {
    message("Global LogNormalize pipeline...")
    TN.combined <- NormalizeData(TN.combined, verbose = FALSE) %>%
      FindVariableFeatures(nfeatures = 2000, verbose = FALSE) %>%
      ScaleData(vars.to.regress = "percent.mt", features = rownames(TN.combined), verbose = FALSE) %>%
      RunPCA(npcs = opt$pcs_global, verbose = FALSE)
    harmony_assay <- "RNA"
  }

  save_plot(ElbowPlot(TN.combined, ndims = 50),
                 file.path(output_dirs$plots, "TNcombined_elbow"))

  # Harmony → UMAP → clustering
  # Harmony → UMAP → clustering
  message("Running Harmony (assay: ", harmony_assay, ")...")
  TN.combined <- RunHarmony(TN.combined, group.by.vars = "batch",
                            assay.use = harmony_assay, verbose = FALSE)
  
  TN.combined <- RunUMAP(TN.combined, reduction = "harmony", dims = 1:opt$pcs_global, seed.use = opt$seed,
                         umap.method = "uwot", metric = "cosine", verbose = FALSE)
                         
  TN.combined <- FindNeighbors(TN.combined, reduction = "harmony", dims = 1:opt$pcs_global, verbose = FALSE)

  resolutions <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8)
  message("Clustering at resolutions: ", paste(resolutions, collapse = ", "))
  TN.combined <- FindClusters(TN.combined, resolution = resolutions, random.seed = opt$seed, verbose = FALSE)
  
  cluster_prefix <- paste0(harmony_assay, "_snn_res.")

  # Clustree
  p_tree <- clustree(TN.combined, prefix = cluster_prefix,
                     node_text_size = 3, edge_arrow = FALSE) +
    ggtitle("Clustree Resolution Tracker") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  save_plot(p_tree, file.path(output_dirs$plots, "TNcombined_clustree"), w = 15, h = 10)

  # Set default resolution
  default_col <- paste0(cluster_prefix, chosen_res)
  if (default_col %in% colnames(TN.combined@meta.data)) {
    Idents(TN.combined) <- default_col
    TN.combined$seurat_clusters <- TN.combined[[default_col]]
  } else {
    message("Warning: resolution ", chosen_res, " not found; idents unchanged")
  }

  write.csv(table(Idents(TN.combined), TN.combined$orig.ident1),
            file.path(output_dirs$tables, "CellNumber_bygroup.csv"))

  # --- 1. Save the Full Object First ---
  full_rds_path <- file.path(output_base, "TN.combined_dim30_full.rds")
  message("Saving full integrated object to: ", full_rds_path)
  saveRDS(TN.combined, full_rds_path)

  # --- 2. Create and Save the Diet Object for Downstream Analysis ---
  message("Trimming a lightweight copy with DietSeurat...")
  TN.diet <- DietSeurat(
    TN.combined,
    counts = TRUE,       # Keep raw counts (required for DGE / pseudobulk)
    data = TRUE,         # Keep normalized data (required for CellChat / SingleR)
    scale.data = FALSE,  # Drop heavy scaled matrices to save RAM and disk space
    assays = c("RNA", "SCT"),
    dimreducs = c("pca", "harmony", "umap")
  )

  message("Saving diet integrated object to: ", integrated_rds)
  saveRDS(TN.diet, integrated_rds) # integrated_rds points to "TN.combined_dim30.rds"
  
  message("===== Harmony integration and dual-saving complete =====\n")
  TN.combined # Return the full object to memory so subsequent steps inside the R session don't lose anything
}

# ==============================================================================
# STEP 4 — GENERATE PLOTS
# Uses numeric cluster labels only.  Run cell_annotation.R afterwards to
# produce annotated UMAP plots with cell-type labels.
# ==============================================================================

generate_plots <- function(chosen_res = 0.4) {
  if (!file.exists(integrated_rds))
    stop("Integrated object not found. Run the 'integrate' step first.")

  message("Loading integrated object for plotting: ", integrated_rds)
  TN.combined <- readRDS(integrated_rds)

  message("\n===== Generating Plots =====")

  # Resolve clustering column
  md_cols        <- colnames(TN.combined@meta.data)
  cluster_prefix <- if      (any(startsWith(md_cols, "SCT_snn_res."))) "SCT_snn_res."
                    else if (any(startsWith(md_cols, "RNA_snn_res."))) "RNA_snn_res."
                    else {
                      fb <- if (isTRUE(opt$use_sct)) "SCT_snn_res." else "RNA_snn_res."
                      message("Warning: no clustering columns found; falling back to ", fb)
                      fb
                    }

  res_col <- paste0(cluster_prefix, chosen_res)
  if (!res_col %in% md_cols)
    stop("Resolution ", chosen_res, " not found (prefix: ", cluster_prefix, ")")

  Idents(TN.combined) <- res_col
  message("Using numeric cluster labels (run cell_annotation.R for cell-type labels)")

  plot_colors <- pick_colors(length(levels(TN.combined)))
  res_title   <- paste0("Global Integration (res: ", chosen_res, ")")
  umap_theme  <- theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  # Prepare RNA assay for DE / visualisation
  DefaultAssay(TN.combined) <- "RNA"
  TN.combined <- try_join_layers(TN.combined)
  TN.combined <- NormalizeData(TN.combined, assay = "RNA", verbose = FALSE)

  # ---- UMAP plots ----
  save_plot(
    DimPlot(TN.combined, reduction = "umap", label = FALSE,
            pt.size = 0.8, cols = plot_colors) +
      ggtitle(paste(res_title, "— Unlabeled")) + umap_theme,
    file.path(output_dirs$plots, "TNcombined_umap_labelF")
  )
  save_plot(
    DimPlot(TN.combined, reduction = "umap", label = TRUE, label.size = 3,
            repel = TRUE, pt.size = 0.8, cols = plot_colors) +
      ggtitle(paste(res_title, "— Labeled")) + umap_theme,
    file.path(output_dirs$plots, "TNcombined_umap_labelT")
  )
  save_plot(
    DimPlot(TN.combined, group.by = "orig.ident2", pt.size = 0.8, cols = plot_colors) +
      ggtitle("UMAP by Sample") + umap_theme,
    file.path(output_dirs$plots, "TNcombined_umap_by_sample")
  )
  save_plot(
    DimPlot(TN.combined, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE,
            split.by = "orig.ident1", pt.size = 0.8, ncol = 2, cols = plot_colors) +
      ggtitle(paste(res_title, "— Split by Condition")) + umap_theme,
    file.path(output_dirs$plots, "TNcombined_umap_split_condition")
  )
  save_plot(
    DimPlot(TN.combined, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE,
            split.by = "orig.ident2", pt.size = 0.8, ncol = 2, cols = plot_colors) +
      ggtitle(paste(res_title, "— Split by Sample")) + umap_theme,
    file.path(output_dirs$plots, "TNcombined_umap_split_sample")
  )

  # ---- Marker detection + heatmap ----
  message("Finding markers...")
  markers <- tryCatch(
    suppressWarnings(FindAllMarkers(TN.combined, assay = "RNA", only.pos = TRUE,
                                    min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)),
    error = function(e) { message("FindAllMarkers error: ", conditionMessage(e)); NULL }
  )

  if (!is.null(markers) && nrow(markers) > 0) {
    write.csv(markers, file.path(output_dirs$tables, "Findallmarkers.csv"), row.names = FALSE)

    cluster_col <- intersect(c("cluster", "group"), colnames(markers))[1]
    if (!is.na(cluster_col)) {
      top_genes <- markers %>%
        dplyr::group_by_at(vars(all_of(cluster_col))) %>%
        dplyr::slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
        dplyr::pull(gene) %>% unique()

      if (length(top_genes) > 0) {
        obj_heat <- ScaleData(TN.combined, features = top_genes, assay = "RNA", verbose = FALSE)
        p_heat   <- DoHeatmap(obj_heat, features = top_genes,
                              group.colors = plot_colors, assay = "RNA") +
          scale_fill_gradient2(low = "magenta", mid = "black", high = "yellow",
                               midpoint = 0, name = "Z-Score") +
          ggtitle(paste0("Top 10 Markers per Cluster (res: ", chosen_res, ")")) +
          theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
        save_plot(p_heat, file.path(output_dirs$plots, "heatmap_top10"), w = 25, h = 25)
      }
    }
  } else {
    message("No DE markers found; skipping heatmap")
  }

  # ---- Stacked bar proportion plots ----
  make_stacked_bar <- function(counts_table, x_label, save_name, w = 10, h = 8) {
    df <- as.data.frame(prop.table(counts_table, margin = 2)) %>%
      setNames(c("Cluster", x_label, "Proportion")) %>%
      filter(Proportion > 0)
    p <- ggplot(df, aes_string(x = x_label, y = "Proportion", fill = "Cluster")) +
      geom_bar(stat = "identity", color = "white", linewidth = 0.2) +
      geom_text(aes(label = as.character(Cluster), size = Proportion),
                position = position_stack(vjust = 0.5), color = "black") +
      scale_size_continuous(range = c(0.5, 4), guide = "none") +
      scale_y_continuous(labels = scales::percent) +
      scale_fill_manual(values = plot_colors) +
      theme_minimal(base_size = 15) +
      labs(x = x_label, y = "% of Total Cells") +
      theme(legend.title = element_blank(), legend.text = element_text(size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1))
    write.csv(prop.table(counts_table, margin = 2) * 100,
              file.path(output_dirs$tables, paste0(save_name, ".csv")), row.names = TRUE)
    save_plot(p, file.path(output_dirs$plots, save_name), w = w, h = h)
  }

  make_stacked_bar(table(Idents(TN.combined), TN.combined$orig.ident1),
                   "Condition", "proportion_by_condition")
  make_stacked_bar(table(Idents(TN.combined), TN.combined$orig.ident2),
                   "Sample", "proportion_by_sample", w = 12)

  # ---- Faceted proportion boxplots (fixed + free, condition + sample) ----
  message("Generating faceted proportion plots...")
  prop_data <- as.data.frame(
    table(Idents(TN.combined), TN.combined$orig.ident2, TN.combined$orig.ident1)
  )
  colnames(prop_data) <- c("Cluster", "Sample", "Condition", "Count")
  prop_data <- prop_data %>%
    filter(Count > 0 |
           Condition == TN.combined$orig.ident1[match(Sample, TN.combined$orig.ident2)]) %>%
    group_by(Sample) %>%
    mutate(Percentage = Count / sum(Count) * 100) %>%
    ungroup() %>%
    mutate(Cluster = trimws(gsub(";.*", "", as.character(Cluster))))

  u_clust           <- unique(prop_data$Cluster)
  prop_data$Cluster <- factor(prop_data$Cluster,
                               levels = u_clust[order(as.numeric(sub(":.*", "", u_clust)))])

  facet_theme <- theme_bw(base_size = 14) +
    theme(legend.position  = "none",
          axis.text.x      = element_text(angle = 45, hjust = 1, face = "bold"),
          strip.text        = element_text(face = "bold", size = 10),
          strip.background  = element_rect(fill = "lightgray"))

  for (scale_arg in c("fixed", "free_y")) {
    sfx <- toupper(sub("_y$", "", scale_arg))
    p_box <- ggplot(prop_data, aes(x = Condition, y = Percentage, fill = Condition)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6) +
      geom_jitter(width = 0.2, size = 2, color = "black", alpha = 0.8) +
      facet_wrap(~Cluster, scales = scale_arg) +
      labs(title = paste0("Cluster Proportions by Condition (", sfx, " axis)"),
           x = "Condition", y = "% of Cells") + facet_theme
    save_plot(p_box,
                   file.path(output_dirs$plots, paste0("prop_boxplot_condition_", sfx)),
                   w = 16, h = 12)

    p_bar <- ggplot(prop_data, aes(x = Sample, y = Percentage, fill = Sample)) +
      geom_col(color = "black", alpha = 0.8) +
      facet_wrap(~Cluster, scales = scale_arg) +
      labs(title = paste0("Cluster Proportions by Sample (", sfx, " axis)"),
           x = "Sample", y = "% of Cells") + facet_theme
    save_plot(p_bar,
                   file.path(output_dirs$plots, paste0("prop_barplot_sample_", sfx)),
                   w = 16, h = 12)
  }

  # ---- Pseudo-bulk PCA ----
  message("Generating pseudo-bulk PCA...")
  avg_expr <- AggregateExpression(TN.combined, group.by = "orig.ident2", assays = "RNA",
                                  normalization.method = "LogNormalize",
                                  return.seurat = FALSE)
  mat      <- t(avg_expr$RNA)
  gene_var <- apply(mat, 2, var)
  mat      <- mat[, gene_var > 0]
  message("Removed ", sum(gene_var == 0), " zero-variance genes for PCA")

  pca_res  <- prcomp(mat, scale. = TRUE)
  pca_data <- as.data.frame(pca_res$x); pca_data$sample <- rownames(pca_data)
  pct_var  <- pca_res$sdev^2 / sum(pca_res$sdev^2)

  p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = sample), size = 6, alpha = 0.9) +
    geom_text_repel(aes(label = sample), size = 4, box.padding = 0.5) +
    guides(color = "none") + theme_bw() +
    ggtitle("PCA of Sample Similarity (Pseudo-bulk)") +
    labs(x = paste0("PC1 (", round(pct_var[1] * 100, 2), "%)"),
         y = paste0("PC2 (", round(pct_var[2] * 100, 2), "%)")) +
    theme(plot.title = element_text(hjust = 0.5, size = 16))
  save_plot(p_pca, file.path(output_dirs$plots, "pca_sample_similarity"), w = 10, h = 8)

  # ---- Summary DotPlot ----
  markers_file <- file.path(output_dirs$tables, "Findallmarkers.csv")
  dot_markers  <- NULL
  if (file.exists(markers_file))
    try({ dot_markers <- read.csv(markers_file, stringsAsFactors = FALSE) }, silent = TRUE)

  if (!is.null(dot_markers) && nrow(dot_markers) > 0) {
    top_genes <- dot_markers %>%
      group_by(cluster) %>%
      slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE) %>%
      pull(gene) %>% unique()

    if (length(top_genes) > 0) {
      p_dot <- tryCatch(
        DotPlot(TN.combined, features = top_genes, assay = "RNA") +
          ggtitle("Top 5 markers per cluster") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)),
        error = function(e) { message("DotPlot failed: ", conditionMessage(e)); NULL }
      )
      if (!is.null(p_dot))
        save_plot(p_dot, file.path(output_dirs$plots, "summary_dotplot"), w = 14, h = 8)
    }
  } else {
    message("No markers available for DotPlot; skipping")
  }

  message("===== Plots complete =====\n")
}

# ==============================================================================
# PIPELINE EXECUTOR
# ==============================================================================

execute_step <- function(step) {
  switch(step,

    # ---- read_csv ----
    read_csv = {
      read_samples_csv(opt$file)
    },

    # ---- process ----
    # Samples are processed in PARALLEL using parallel::mclapply.
    # Each call to process_sample is fully self-contained and writes its own
    # RDS file, so workers never share state.  On Windows mclapply falls back
    # to lapply (no forking), which is still correct.
    process = {
      samples_df <- readRDS(file.path(output_base, "samples_df.rds"))

      pending <- which(!file.exists(
        file.path(output_dirs$processed,
                  paste0(samples_df$sample_names, "_processed.rds"))
      ))

      if (length(pending) == 0) {
        message("All samples already processed; skipping")
        return(invisible(NULL))
      }
      message("Processing ", length(pending), " sample(s) across ", n_cores, " core(s)")

      params <- list(
        doublet_rate = opt$doublet_rate,
        min_features = opt$min_features,
        max_features = opt$max_features,
        max_mt       = opt$max_mt
      )
      use_sct <- isTRUE(opt$use_sct)

      results <- parallel::mclapply(
        seq_len(nrow(samples_df)),
        function(i) {
          process_sample(
            sample_name   = samples_df$sample_names[i],
            sample_ident1 = samples_df$ident1[i],
            sample_ident2 = samples_df$ident2[i],
            base_data_dir = opt$datadir,
            out_dirs      = output_dirs,
            use_sct       = use_sct,
            params        = params
          )
        },
        mc.cores       = n_cores,
        mc.preschedule = FALSE,
        mc.set.seed    = TRUE
      )

      failed <- which(sapply(results, is.null))
      if (length(failed) > 0)
        message("WARNING: failed samples — ",
                paste(samples_df$sample_names[failed], collapse = ", "))
      else
        message("All samples processed successfully")
      invisible(results)
    },

    # ---- integrate ----
    integrate = {
      if (file.exists(integrated_rds)) {
        message("Loading existing integrated object: ", integrated_rds)
        return(readRDS(integrated_rds))
      }
      sample_files <- list.files(output_dirs$processed,
                                 pattern = "_processed.rds$", full.names = TRUE)
      if (length(sample_files) == 0)
        stop("No processed samples in ", output_dirs$processed)
      message("Integrating ", length(sample_files), " sample(s)")
      integrate_samples(lapply(sample_files, readRDS), chosen_res = opt$resolution)
    },

    # ---- plot ----
    plot = {
      output_dirs$plots  <<- file.path(output_base, "plots",  res_folder)
      output_dirs$tables <<- file.path(output_base, "tables", res_folder)
      dir.create(output_dirs$plots,  recursive = TRUE, showWarnings = FALSE)
      dir.create(output_dirs$tables, recursive = TRUE, showWarnings = FALSE)
      generate_plots(chosen_res = opt$resolution)
    },

    # ---- all ----
    all = {
      execute_step("read_csv")
      execute_step("process")
      execute_step("integrate")
      execute_step("plot")
    },

    stop("Invalid step. Choose: read_csv, process, integrate, plot, all")
  )
}

# ==============================================================================
# MAIN
# ==============================================================================

if (is.null(opt$file)) stop("Specify input file with -f")
if (is.null(opt$step)) stop("Specify pipeline step with -s")

execute_step(opt$step)
message("Step '", opt$step, "' completed at ", Sys.time())

sink(type = "message")
sink(type = "output")
close(log_conn)