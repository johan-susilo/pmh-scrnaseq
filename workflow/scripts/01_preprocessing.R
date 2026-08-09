#!/usr/bin/env Rscript
# Usage: Rscript 01_preprocessing.R -f input.csv -d /path/to/data -s all -o results --config config/config.yaml
#
# Steps:
#   read_csv   – parse input.csv and cache samples_df.rds
#   process    – QC, doublet removal, per-sample normalisation  (parallel across samples)
#   integrate  - merge + Harmony integration + multi-res clustering -> TN.combined_dim30.rds
#   plot       – UMAP / heatmap / proportion plots (uses numeric cluster labels)
#   all        – read_csv -> process -> integrate -> plot
#
# Annotation is handled by 02_global_annotation.R, which reads
# TN.combined_dim30.rds, applies consensus labels, and writes TN.combined_annotated.rds.
#
# CONTRACT WITH 02: this script must always produce
#   <outdir>/TN.combined_dim30.rds  (<outdir> = whatever -o this rule is given, used as-is)
# with clustering columns named "<assay>_snn_res.<r>" for every resolution in
# preprocessing.resolution / the --resolution flag, resolvable via
# resolve_res_col() in 00_utils.R.

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
# All numeric/logical defaults are NULL here and resolved AFTER config.yaml is
# loaded (see "RESOLVE DEFAULTS" below), so config.yaml is the single source
# of truth and a CLI flag is only needed to override it.
# ==============================================================================

option_list <- list(
  make_option(c("-f", "--file"),       type = "character", help = "Input CSV (columns: sample_names, ident1, ident2)"),
  make_option(c("-d", "--datadir"),    type = "character", default = "data", help = "Base directory with sample folders"),
  make_option(c("-s", "--step"),       type = "character", help = "Pipeline step: read_csv, process, integrate, plot, all"),
  make_option(c("-o", "--output"),     type = "character", default = "results", help = "Base output directory"),
  make_option(c("--config"),           type = "character", default = "config/config.yaml", help = "Path to config.yaml"),
  make_option(c("-c", "--cores"),      type = "integer",   default = NULL, help = "Cores for parallel sample processing [default: all available - 1]"),
  make_option(c("-r", "--resolution"), type = "numeric",   default = NULL, help = "Default clustering resolution to use for plotting [config: preprocessing.resolution]"),
  make_option("--doublet_rate",        type = "numeric",   default = NULL, help = "Expected doublet rate [config: preprocessing.doublet_rate]"),
  make_option("--min_features",        type = "integer",   default = NULL, help = "Min features per cell [config: preprocessing.min_features]"),
  make_option("--max_features",        type = "integer",   default = NULL, help = "Max features per cell [config: preprocessing.max_features]"),
  make_option("--max_mt",              type = "numeric",   default = NULL, help = "Max mitochondrial % [config: preprocessing.max_mt]"),
  make_option("--use_sct",             type = "logical",   default = NULL, help = "Use SCTransform normalization [config: preprocessing.use_sct]"),
  make_option("--batch_var",           type = "character", default = NULL, help = "Batch variable for Harmony [config: preprocessing.batch_var]"),
  make_option(c("--pcs_local"),        type = "integer",   default = NULL, help = "PCs for per-sample normalization [config: preprocessing.pcs_local]"),
  make_option(c("--pcs_global"),       type = "integer",   default = NULL, help = "PCs for global integration [config: preprocessing.pcs_global]"),
  make_option(c("--seed"),             type = "integer",   default = NULL, help = "Global random seed [config: reproducibility.random_seed]")
)

opt <- parse_args(OptionParser(option_list = option_list))
cfg <- get_config(opt$config)

# ---- RESOLVE DEFAULTS: config.yaml first, CLI flag overrides if supplied ----
`%||%` <- function(a, b) if (is.null(a)) b else a
opt$resolution    <- opt$resolution    %||% cfg_get(cfg, "preprocessing", "resolution",    default = 0.2)
opt$doublet_rate  <- opt$doublet_rate  %||% cfg_get(cfg, "preprocessing", "doublet_rate",  default = 0.08)
opt$min_features  <- opt$min_features  %||% cfg_get(cfg, "preprocessing", "min_features",  default = 200)
opt$max_features  <- opt$max_features  %||% cfg_get(cfg, "preprocessing", "max_features",  default = 5000)
opt$max_mt        <- opt$max_mt        %||% cfg_get(cfg, "preprocessing", "max_mt",        default = 30)
opt$use_sct       <- opt$use_sct       %||% cfg_get(cfg, "preprocessing", "use_sct",       default = TRUE)
opt$batch_var     <- opt$batch_var     %||% cfg_get(cfg, "preprocessing", "batch_var",     default = "orig.ident2")
opt$pcs_local     <- opt$pcs_local     %||% cfg_get(cfg, "preprocessing", "pcs_local",     default = 30)
opt$pcs_global    <- opt$pcs_global    %||% cfg_get(cfg, "preprocessing", "pcs_global",    default = 50)
opt$seed          <- opt$seed          %||% cfg_get(cfg, "reproducibility", "random_seed", default = 42)
if (is.null(opt$file)) opt$file <- cfg_get(cfg, "input_csv", default = NULL)

res_folder <- paste0("res_", opt$resolution)
set.seed(opt$seed)

# ==============================================================================
# OUTPUT DIRECTORIES
#
# `-o` is used AS-IS -- no extra subfolder is appended here. The Snakemake
# rule for this script already points `-o` at a stage-specific directory
# (e.g. .../01_preprocessing); nesting a second "01_integration" folder
# underneath that silently breaks the fixed output path the rule's `output:`
# block expects (this caused a MissingOutputException previously -- see the
# SHARED FILE PATHS note below for the matching filename fix).
# ==============================================================================

output_base <- if (!is.null(opt$output) && nzchar(opt$output)) opt$output else "results"
dir.create(output_base, recursive = TRUE, showWarnings = FALSE)
output_dirs <- make_stage_dirs(output_base, c("01_per_sample_processed", "02_plots", "03_tables", "logs"))
# This maps your new neat folders back to the original variables the script expects
names(output_dirs) <- c("processed", "plots", "tables", "logs")

# ==============================================================================
# PARALLEL SETUP
# ==============================================================================

n_cores <- if (!is.null(opt$cores)) opt$cores else max(1L, parallel::detectCores(logical = FALSE) - 1L)
message("Parallel sample processing: ", n_cores, " core(s)")

# ==============================================================================
# SHARED FILE PATHS
#
# Kept as TN.combined_dim30.rds / _full.rds -- these are the exact names the
# Snakefile's `output:` blocks and 02_global_annotation.R's `-i` default
# expect. "dim30" no longer literally reflects opt$pcs_global (which is now
# config-driven), but renaming is a coordinated Snakefile + downstream-script
# change, not something to do silently inside this one file.
# ==============================================================================

integrated_full_rds <- file.path(output_base, "TN.combined_dim30_full.rds")
integrated_rds       <- file.path(output_base, "TN.combined_dim30.rds")

# ==============================================================================
# HELPER UTILITIES
# ==============================================================================

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
      theme_pipeline(12) +
      theme(legend.position = "none") +
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
      theme_pipeline(12) +
      theme(legend.position = "none") +
      labs(title = title, x = "", y = "Value")
  }
}

normalize_and_pca <- function(seur_obj, use_sct) {
  if (isTRUE(use_sct)) {
    message("Running SCTransform v2...")
    seur_obj <- SCTransform(seur_obj, vars.to.regress = "percent.mt",
                            method = "glmGamPoi", vst.flavor = "v2", verbose = FALSE)
    seur_obj <- RunPCA(seur_obj, assay = "SCT", npcs = opt$pcs_global, seed.use = opt$seed, verbose = FALSE)
  } else {
    message("Running LogNormalize pipeline...")
    seur_obj <- seur_obj %>%
      NormalizeData() %>%
      FindVariableFeatures() %>%
      ScaleData(vars.to.regress = "percent.mt") %>%
      RunPCA(npcs = opt$pcs_global, seed.use = opt$seed, verbose = FALSE)
  }
  seur_obj
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
                           base_data_dir, out_dirs, use_sct, params, seed) {
  tryCatch({
    output_rds      <- file.path(out_dirs$processed, paste0(sample_name, "_processed.rds"))
    sample_plot_dir <- file.path(out_dirs$plots, sample_name)
    plot_path       <- function(suffix) file.path(sample_plot_dir, paste0(sample_name, suffix))

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

    seur_obj <- CreateSeuratObject(
      counts       = Read10X(data.dir = file.path(base_data_dir, sample_name)),
      project      = sample_name,
      min.cells    = 3,
      min.features = 10
    )
    message("Dimensions: ", dim(seur_obj)[1], " features x ", dim(seur_obj)[2], " cells")
    seur_obj[["percent.mt"]] <- PercentageFeatureSet(seur_obj, pattern = "^MT-")

    save_qc_violin(seur_obj, sample_name, plot_path("_01_qc_unfiltered"),
                   paste0("QC (Unfiltered) - ", sample_name))

    n_before <- ncol(seur_obj)
    seur_obj <- subset(seur_obj,
                       subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 50)
    message("Pre-filter: removed ", n_before - ncol(seur_obj),
            " empty droplets; keeping ", ncol(seur_obj))

    save_qc_violin(seur_obj, sample_name, plot_path("_02_qc_prefiltered"),
                   paste0("QC (Pre-filtered) - ", sample_name))

    seur_obj <- normalize_and_pca(seur_obj, use_sct)
    save_plot(ElbowPlot(seur_obj, ndims = opt$pcs_local), plot_path("_03_elbow"))

    seur_obj <- seur_obj %>%
      FindNeighbors(dims = 1:opt$pcs_local) %>%
      FindClusters(random.seed = seed) %>%
      RunUMAP(dims = 1:opt$pcs_local, seed.use = seed, umap.method = "uwot", metric = "cosine")
    save_plot(DimPlot(seur_obj, reduction = "umap", label = TRUE),
                   plot_path("_04_umap_initial"))

    # -- DoubletFinder --
    # set.seed() here is the fix: paramSweep()/find.pK() previously had no
    # seed control, so the chosen pK (and therefore which cells get called
    # doublets, and everything computed downstream of that) could vary
    # between runs even with opt$seed fixed everywhere else.
    message("Running DoubletFinder pK sweep (seed = ", seed, ")...")
    set.seed(seed)
    sweep.res <- paramSweep(seur_obj, PCs = 1:20, sct = isTRUE(use_sct))
    bcmvn    <- find.pK(summarizeSweep(sweep.res, GT = FALSE))

    save_plot(
      ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
        geom_point() + geom_line() +
        ggtitle(paste0("pK - ", sample_name)) + theme_bw(),
      plot_path("_05_pk")
    )

    pK             <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
    homotypic.prop <- modelHomotypic(seur_obj@meta.data$seurat_clusters)
    nExp_adj       <- round(params$doublet_rate * nrow(seur_obj@meta.data) *
                            (1 - homotypic.prop))
    message("Optimal pK: ", pK, "  |  Adjusted expected doublets: ", nExp_adj)

    set.seed(seed)
    seur_obj <- doubletFinder(seur_obj, PCs = 1:20, pN = 0.25, pK = pK,
                              nExp = nExp_adj, sct = isTRUE(use_sct))

    DF_col <- grep("DF.classifications", colnames(seur_obj@meta.data), value = TRUE)[1]
    if (is.na(DF_col))
      stop("DoubletFinder produced no classification column for ", sample_name)

    save_plot(
      DimPlot(seur_obj, reduction = "umap", group.by = DF_col, cols = doublet_color) +
        ggtitle(paste0("Doublets (before removal) - ", sample_name)),
      plot_path("_06_doublets_before")
    )
    message("Doublet counts:"); print(table(seur_obj@meta.data[[DF_col]]))

    seur_obj <- seur_obj[, seur_obj@meta.data[[DF_col]] == "Singlet"]
    message("Cells after doublet removal: ", ncol(seur_obj))

    save_plot(
      DimPlot(seur_obj, reduction = "umap", group.by = DF_col, cols = doublet_color) +
        ggtitle(paste0("Doublets (after removal) - ", sample_name)),
      plot_path("_07_doublets_after")
    )
    save_qc_violin(seur_obj, sample_name, plot_path("_08_qc_post_doublet"),
                   paste0("QC (Post-Doublet) - ", sample_name))

    seur_obj$orig.ident1 <- sample_ident1
    seur_obj$orig.ident2 <- sample_ident2
    
    seur_obj$Detailed_Condition <- case_when(
      grepl("HTY|UA", sample_ident2, ignore.case = TRUE) ~ "Healthy",
      grepl("AC", sample_ident2, ignore.case = TRUE) ~ "Acute",
      grepl("CH", sample_ident2, ignore.case = TRUE) ~ "Chronic",
      TRUE ~ "Unknown"
    )
    seur_obj$Condition <- ifelse(seur_obj$Detailed_Condition == "Healthy", "Healthy", "PMH")

    n_before <- ncol(seur_obj)
    seur_obj <- subset(seur_obj,
                       subset = nFeature_RNA > params$min_features &
                                nFeature_RNA < params$max_features &
                                percent.mt   < params$max_mt)
    message("Final QC: removed ", n_before - ncol(seur_obj),
            " cells; keeping ", ncol(seur_obj))

    save_qc_violin(seur_obj, sample_name, plot_path("_09_qc_final"),
                   paste0("QC (Final) - ", sample_name))

    seur_obj@misc$processed_with_sct <- isTRUE(use_sct)
    seur_obj@misc$pipeline_seed      <- seed
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

integrate_samples <- function(sample_list, chosen_res) {
  sample_list <- Filter(Negate(is.null), sample_list)
  if (length(sample_list) == 0) stop("No valid samples for integration.")

  sample_objs <- lapply(sample_list, function(x) {
    if (is.character(x) && file.exists(x)) readRDS(x) else x
  })

  if (any(duplicated(unlist(lapply(sample_objs, colnames))))) {
    message("Duplicated barcodes detected - prefixing with sample tags")
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

  used_sct <- sapply(sample_objs, function(obj) {
    flag <- try(obj@misc$processed_with_sct, silent = TRUE)
    if (!is.null(flag) && !inherits(flag, "try-error")) as.logical(flag)
    else "SCT" %in% names(obj)
  })
  if (length(unique(used_sct)) > 1)
    stop("Inconsistent SCTransform usage across samples. Reprocess with a consistent --use_sct.")

  inferred_sct <- isTRUE(unique(used_sct))
  if (inferred_sct != isTRUE(opt$use_sct)) {
    message("Overriding opt$use_sct -> ", inferred_sct, " to match processed samples")
    opt$use_sct <<- inferred_sct
  }

  message("\n===== Starting Harmony Integration =====")

  sample_objs <- lapply(sample_objs, function(obj) {
    DefaultAssay(obj) <- "RNA"
    if ("SCT" %in% names(obj)) obj[["SCT"]] <- NULL
    obj
  })

  TN.combined <- if (length(sample_objs) > 1)
    merge(sample_objs[[1]], y = sample_objs[-1])
  else
    sample_objs[[1]]

  if (!inherits(TN.combined, "Seurat"))
    stop("Merged result is not a Seurat object. Check your input samples.")
  TN.combined <- try_join_layers(TN.combined)

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

  if (!(opt$batch_var %in% colnames(TN.combined@meta.data))) {
    fb <- if ("orig.ident2" %in% colnames(TN.combined@meta.data)) "orig.ident2" else "orig.ident1"
    message("batch_var '", opt$batch_var, "' not found; using '", fb, "'")
    TN.combined$batch <- factor(as.character(TN.combined@meta.data[[fb]]))
  } else {
    TN.combined$batch <- TN.combined[[opt$batch_var]]
  }

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
      RunPCA(npcs = opt$pcs_global, seed.use = opt$seed, verbose = FALSE)
    harmony_assay <- "RNA"
  }

  save_plot(ElbowPlot(TN.combined, ndims = opt$pcs_global),
                 file.path(output_dirs$plots, "TNcombined_elbow"))

  message("Running Harmony (assay: ", harmony_assay, ")...")
  TN.combined <- RunHarmony(TN.combined, group.by.vars = "batch",
                            assay.use = harmony_assay, verbose = FALSE)

  TN.combined <- RunUMAP(TN.combined, reduction = "harmony", dims = 1:opt$pcs_global, seed.use = opt$seed,
                         umap.method = "uwot", metric = "cosine", verbose = FALSE)

  TN.combined <- FindNeighbors(TN.combined, reduction = "harmony", dims = 1:opt$pcs_global, verbose = FALSE)

  # Resolutions to sweep for clustree come from config; chosen_res (the
  # single "active" resolution) must always be included in the sweep.
  resolutions <- unique(sort(c(
    cfg_get(cfg, "subsetting", "resolutions", default = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8)),
    chosen_res
  )))
  message("Clustering at resolutions: ", paste(resolutions, collapse = ", "))
  TN.combined <- FindClusters(TN.combined, resolution = resolutions, random.seed = opt$seed, verbose = FALSE)

  cluster_prefix <- paste0(harmony_assay, "_snn_res.")

  p_tree <- clustree(TN.combined, prefix = cluster_prefix,
                     node_text_size = 3, edge_arrow = FALSE) +
    ggtitle("Clustree Resolution Tracker") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  save_plot(p_tree, file.path(output_dirs$plots, "TNcombined_clustree"), w = 15, h = 10)

  default_col <- paste0(cluster_prefix, chosen_res)
  if (default_col %in% colnames(TN.combined@meta.data)) {
    Idents(TN.combined) <- default_col
    TN.combined$seurat_clusters <- TN.combined[[default_col]]
  } else {
    message("Warning: resolution ", chosen_res, " not found; idents unchanged")
  }

  write.csv(table(Idents(TN.combined), TN.combined$orig.ident1),
            file.path(output_dirs$tables, "CellNumber_bygroup.csv"))

  TN.combined@misc$pipeline_seed <- opt$seed
  TN.combined@misc$pipeline_config_path <- opt$config

  message("Saving full integrated object to: ", integrated_full_rds)
  saveRDS(TN.combined, integrated_full_rds)

  message("Trimming a lightweight copy with DietSeurat...")
  TN.diet <- DietSeurat(
    TN.combined,
    counts = TRUE,
    data = TRUE,
    scale.data = FALSE,
    assays = c("RNA", "SCT"),
    dimreducs = c("pca", "harmony", "umap")
  )

  message("Saving diet integrated object to: ", integrated_rds)
  saveRDS(TN.diet, integrated_rds)

  message("===== Harmony integration and dual-saving complete =====\n")
  TN.combined
}

# ==============================================================================
# STEP 4 — GENERATE PLOTS
# ==============================================================================

generate_plots <- function(chosen_res) {
  if (!file.exists(integrated_rds))
    stop("Integrated object not found. Run the 'integrate' step first.")

  message("Loading integrated object for plotting: ", integrated_rds)
  TN.combined <- readRDS(integrated_rds)

  message("\n===== Generating Plots =====")

  res_col <- resolve_res_col(TN.combined, chosen_res)
  if (is.null(res_col)) {
    fb <- if (isTRUE(opt$use_sct)) "SCT_snn_res." else "RNA_snn_res."
    stop("Resolution ", chosen_res, " not found (expected prefix: ", fb, ")")
  }

  Idents(TN.combined) <- res_col
  message("Using numeric cluster labels (run 02_global_annotation.R for cell-type labels)")

  plot_colors <- pick_colors(length(levels(TN.combined)))
  res_title   <- paste0("Global Integration (res: ", chosen_res, ")")
  umap_theme  <- theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  DefaultAssay(TN.combined) <- "RNA"
  TN.combined <- try_join_layers(TN.combined)
  TN.combined <- NormalizeData(TN.combined, assay = "RNA", verbose = FALSE)

  save_plot(
    DimPlot(TN.combined, reduction = "umap", label = FALSE,
            pt.size = 0.8, cols = plot_colors) +
      ggtitle(paste(res_title, "- Unlabeled")) + umap_theme,
    file.path(output_dirs$plots, "TNcombined_umap_labelF")
  )
  save_plot(
    DimPlot(TN.combined, reduction = "umap", label = TRUE, label.size = 3,
            repel = TRUE, pt.size = 0.8, cols = plot_colors) +
      ggtitle(paste(res_title, "- Labeled")) + umap_theme,
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
      ggtitle(paste(res_title, "- Split by Condition")) + umap_theme,
    file.path(output_dirs$plots, "TNcombined_umap_split_condition")
  )
  save_plot(
    DimPlot(TN.combined, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE,
            split.by = "orig.ident2", pt.size = 0.8, ncol = 2, cols = plot_colors) +
      ggtitle(paste(res_title, "- Split by Sample")) + umap_theme,
    file.path(output_dirs$plots, "TNcombined_umap_split_sample")
  )

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

    read_csv = { read_samples_csv(opt$file) },

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
            params        = params,
            seed          = opt$seed
          )
        },
        mc.cores       = n_cores,
        mc.preschedule = FALSE,
        mc.set.seed    = TRUE
      )

      failed <- which(sapply(results, is.null))
      if (length(failed) > 0)
        message("WARNING: failed samples - ",
                paste(samples_df$sample_names[failed], collapse = ", "))
      else
        message("All samples processed successfully")
      invisible(results)
    },

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

    plot = {
      output_dirs$plots  <<- file.path(output_base, "plots",  res_folder)
      output_dirs$tables <<- file.path(output_base, "tables", res_folder)
      dir.create(output_dirs$plots,  recursive = TRUE, showWarnings = FALSE)
      dir.create(output_dirs$tables, recursive = TRUE, showWarnings = FALSE)
      generate_plots(chosen_res = opt$resolution)
    },

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

if (is.null(opt$file)) stop("Specify input file with -f, or set input_csv in config.yaml")
if (is.null(opt$step)) stop("Specify pipeline step with -s")

execute_step(opt$step)
message("Step '", opt$step, "' completed at ", Sys.time())

