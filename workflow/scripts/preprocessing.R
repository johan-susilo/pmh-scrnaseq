#!/usr/bin/env Rscript
# Usage: Rscript /home/johan/pipeline/scRNA/skin/preprocessing.R -f /home/johan/pipeline/scRNA/skin/input.csv -d /home/johan/data/PMH_scRNA-seq -s all -o /home/johan/output/skin_pmh/

# Example: Rscript /home/johan/pipeline/scRNA/skin/preprocessing_harmony.R -f /home/johan/pipeline/scRNA/skin/input.csv -d /home/johan/data/PMH_scRNA-seq -s all -o /home/johan/output/skin_pmh_harmony_0.4/ -s plot -r 0.1
# Example: Rscript /home/johan/pipeline/scRNA/skin/preprocessing.R -f input.csv -d /path/to/data -s read_csv
# Example: Rscript /home/johan/pipeline/scRNA/skin/preprocessing.R -f input.csv -d /path/to/data -s process -c 4
# Example: Rscript /home/johan/pipeline/scRNA/skin/preprocessing.R -f input.csv -d . -s integrate
# Example: Rscript /home/johan/pipeline/scRNA/skin/preprocessing.R -f input.csv -d /home/johan/data/PMH_scRNA-seq -s plot
# Example: Rscript /home/johan/pipeline/scRNA/skin/preprocessing.R -f /home/johan/pipeline/scRNA/skin/input.csv -d /home/johan/data/PMH_scRNA-seq -s all -o /home/johan/output/skin_pmh -c 8
# Example: Rscript /home/johan/pipeline/scRNA/skin/preprocessing.R -f /home/johan/pipeline/scRNA/skin/input_side_task.csv -d /mnt/80T/pipeline/single_cell/side_task/HC+DF+K -s all -o /home/johan/output/skin_pmh/side_task/HC+DF+K
# Rscript /home/johan/pipeline/scRNA/skin/preprocessing.R -f /home/johan/pipeline/scRNA/skin/input_side_task2.csv -d /mnt/80T/pipeline/single_cell/side_task/HC+DF+SS -s all -o /home/johan/output/skin_pmh/side_task/HC+DF+SS
# Rscript /home/johan/pipeline/scRNA/skin/preprocessing.R -f /home/johan/pipeline/scRNA/skin/input_side_task3.csv -d /mnt/80T/pipeline/single_cell/side_task/HC+DF+PMH -s all -o /home/johan/output/skin_pmh/side_task/HC+DF+PMH
# Rscript /home/johan/pipeline/scRNA/skin/preprocessing.R -f /home/johan/pipeline/scRNA/skin/input_side_task4.csv -d /mnt/80T/pipeline/single_cell/side_task/HC+DF+HS -s all -o /home/johan/output/skin_pmh/side_task/HC+DF+HS
#Rscript /home/johan/pipeline/scRNA/skin/preprocessing.R -f /home/johan/pipeline/scRNA/skin/input_side_task5.csv -d /mnt/80T/pipeline/single_cell/side_task/HC+DF+DFonly1 -s all -o /home/johan/output/skin_pmh/side_task/HC+DF+DFonly1
# Rscript /home/johan/pipeline/scRNA/skin/preprocessing.R -f /home/johan/pipeline/scRNA/skin/input.csv -d /home/johan/data/PMH_scRNA-seq -o /home/johan/output/skin_pmh_harmony_sctransform2/ -s all -r 0.2


Sys.time()



if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
  message("Downloading DoubletFinder directly from GitHub...")
  remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", upgrade = "never", quiet = TRUE)
}

if (!requireNamespace("scCATCH", quietly = TRUE)) {
  message("Downloading scCATCH directly from CRAN...")
  install.packages("scCATCH", repos = "https://cloud.r-project.org", quiet = TRUE)
}


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
  # library(clusterProfiler)  # Not used in this pipeline, commented out to avoid dependency issues
  library(gplots)
  library(ggplot2)
  library(ggnewscale)
  library(RColorBrewer)
  library(tidyr)
  library(ggrepel)
  library(stringr)
  library(patchwork)
  library(scales)
  library(parallel)
  library(future)
  library(future.apply)
  library(scCATCH)
  library(harmony)  
  library(clustree)
})


set.seed(42)

# Ensure cc.genes is available (used for CellCycleScoring). Try to load from Seurat data; provide safe fallback.
if (!exists("cc.genes")) {
  tryCatch({
    utils::data("cc.genes", package = "Seurat")
    if (!exists("cc.genes")) stop("cc.genes not found after data()")
    message("Loaded cc.genes from Seurat package")
  }, error = function(e) {
    message("Warning: failed to load 'cc.genes' from Seurat: ", conditionMessage(e))
    message("Proceeding with empty cc.genes lists; CellCycleScoring will be skipped if lists are empty.")
    cc.genes <- list(s.genes = character(0), g2m.genes = character(0))
  })
}

# --- Increase Global Memory Limit for Parallel Processing ---
options(future.globals.maxSize = 100 * 1024^3)

# --- Helper Function for QC Violin Plots (workaround for Seurat v5.3.1 VlnPlot bug) ---
create_qc_violin_plot <- function(seurat_obj, features, title) {
  # Get cluster identity if available
  if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    cluster_id <- "seurat_clusters"
  } else if (length(levels(Idents(seurat_obj))) > 1) {
    # Add current idents as a column
    seurat_obj@meta.data$temp_ident <- Idents(seurat_obj)
    cluster_id <- "temp_ident"
  } else {
    cluster_id <- NULL
  }

  # Extract metadata for specified features
  if (!is.null(cluster_id)) {
    qc_data <- seurat_obj@meta.data %>%
      select(all_of(c(features, cluster_id))) %>%
      mutate(cell = rownames(seurat_obj@meta.data)) %>%
      pivot_longer(cols = all_of(features), names_to = 'metric', values_to = 'value')

    # Rename cluster column for consistency
    colnames(qc_data)[colnames(qc_data) == cluster_id] <- "Identity"
    qc_data$Identity <- as.factor(qc_data$Identity)

    # Generate colors based on number of clusters
    n_clusters <- length(unique(qc_data$Identity))
    if (n_clusters <= 12) {
      colors <- colorRampPalette(brewer.pal(min(n_clusters, 12), "Paired"))(n_clusters)
    } else {
      colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_clusters)
    }

    # Create violin plot grouped by cluster
    p <- ggplot(qc_data, aes(x = Identity, y = value, fill = Identity)) +
      geom_violin(trim = FALSE, scale = "width") +
      geom_jitter(size = 0.1, alpha = 0.1, width = 0.2) +
      facet_wrap(~metric, scales = 'free', ncol = length(features)) +
      scale_fill_manual(values = colors) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "lightgray")
      ) +
      labs(title = title, x = 'Identity', y = 'Value')
  } else {
    # Fallback: no clustering information available
    qc_data <- seurat_obj@meta.data %>%
      select(all_of(features)) %>%
      mutate(cell = rownames(seurat_obj@meta.data),
             Identity = "All") %>%
      pivot_longer(cols = all_of(features), names_to = 'metric', values_to = 'value')

    # Create simple violin plot
    p <- ggplot(qc_data, aes(x = Identity, y = value, fill = metric)) +
      geom_violin(trim = FALSE, scale = "width") +
      geom_jitter(size = 0.1, alpha = 0.2, width = 0.2) +
      facet_wrap(~metric, scales = 'free', ncol = length(features)) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "lightgray")
      ) +
      labs(title = title, x = '', y = 'Value')
  }

  return(p)
}

# Command-line Interface ----------------------------------------------------
option_list <- list(
  make_option(c("-f", "--file"), type = "character", help = "Input CSV file (columns: sample_names, ident1, ident2)"),
  make_option(c("-r", "--resolution"), type = "numeric", default = 0.4,
              help = "Default clustering resolution for global plots [default: 0.4]"),
  make_option(c("-d", "--datadir"), type = "character", default = ".",
              help = "Base directory containing sample folders [default: current directory]"),
  make_option(c("-s", "--step"), type = "character",
              help = "Pipeline step: read_csv, process, integrate, plot, all"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Base output directory [overrides default output_base]"),
  make_option(c("-c", "--cores"), type = "integer", default = NULL,
              help = "Number of cores for parallel processing [default: detect available cores]"),
  make_option(c("--doublet_rate"), type = "numeric", default = 0.08,
              help = "Expected doublet formation rate [default: 0.08 for ~10k cells]"),
  make_option(c("--min_features"), type = "integer", default = 200,
              help = "Minimum number of features per cell [default: 200]"),
  make_option(c("--max_features"), type = "integer", default = 5000,
              help = "Maximum number of features per cell [default: 5000]"),
  make_option(c("--max_mt"), type = "numeric", default = 30,
              help = "Maximum mitochondrial percentage [default: 30]"),
  make_option(c("--use_sct"), type = "logical", default = TRUE,
              help = "Use per-sample SCTransform before merging/integration [default: TRUE]"),
  make_option(c("--batch_var"), type = "character", default = "orig.ident2",
              help = "Metadata column used as technical batch for Harmony [default: orig.ident2]")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)
res_folder <- opt$resolution

# Provide a sensible default output base and directories
output_base <- "output"
output_dirs <- list(
  processed = file.path(output_base, "processed"),
  plots = file.path(output_base, "plots"),
  tables = file.path(output_base, "tables"),
  logs = file.path(output_base, "logs")
)
lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

# If -o provided, override output_base and recreate directories
if (!is.null(opt$output) && nzchar(opt$output)) {
  output_base <- opt$output
  dir.create(output_base, recursive = TRUE, showWarnings = FALSE)
  output_dirs <- list(
    processed = file.path(output_base, "processed"),
    plots = file.path(output_base, "plots"),
    tables = file.path(output_base, "tables"),
    logs = file.path(output_base, "logs")
  )
  lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
}

# Set up processing mode - ALWAYS use sequential to avoid deadlocks
message("Using SEQUENTIAL processing mode to avoid deadlocks")
plan(sequential)

# Define color schemes
mycolor <- c("#CCCCCC", # gray for other
  "#A6CEE3","#FF7F00","#09519C","#FFFB33","#556b2f",
  "#FF00FF","#377EB8","#bb8fce","#666666","#90EE90",
  "#ff4500","#A6761D","#E67E22","#323695","#E81E32",
  "#006837","#CBC3E3","#F1C40F","#3498DB","#34495E",
  "#FA9399","#48C9B0","#7D3C98","#ff4500","#8b4513",
  "#8a2be2","#f0e68c","#00ffff","#32CD32","#b03060")

doublet_color <- c("Doublet" = "#e35473", "Singlet" = "#54cdeb")

read_samples_csv <- function(csv_file) {
  Sys.time()
  # Read CSV file with sample information
  message("Reading sample information from ", csv_file)

  if (!file.exists(csv_file)) {
    stop("Input CSV file not found: ", csv_file, call. = FALSE)
  }
  fi <- file.info(csv_file)
  if (is.na(fi$size) || fi$size == 0) {
    stop("Input CSV exists but is empty: ", csv_file,
         "\nPlease provide a CSV with a header and at least one sample row.", call. = FALSE)
  }

  samples_df <- tryCatch({
    read.csv(csv_file, header = TRUE, stringsAsFactors = FALSE)
  }, error = function(e) {
    stop("Failed to read CSV file: ", csv_file, "\nError: ", e$message, call. = FALSE)
  })

  if (nrow(samples_df) == 0) {
    stop("Input CSV has a header but no data rows. Add at least one sample row to ", csv_file, call. = FALSE)
  }

  # Verify expected columns are present
  required_cols <- c("sample_names", "ident1", "ident2")
  if (!all(required_cols %in% colnames(samples_df))) {
    stop("CSV file must contain columns: ", paste(required_cols, collapse=", "), call. = FALSE)
  }

  # Save the dataframe for future steps
  saveRDS(samples_df, file.path(output_base, "samples_df.rds"))

  message("Found ", nrow(samples_df), " samples in the CSV file")
  return(samples_df)
}

# Helper function to safely save PDF & png plots
# --- Helper function to safely save BOTH PDF and PNG ---
safe_save_plot <- function(plot_obj, base_filepath, w = 15, h = 15) {
  # 1. Save PDF
  pdf_path <- paste0(base_filepath, ".pdf")
  tryCatch({
    pdf(pdf_path, width = w, height = h)
    print(plot_obj)
    dev.off()
    message("Saved PDF: ", pdf_path)
  }, error = function(e) {
    message("Warning: Failed to save PDF ", pdf_path, ": ", conditionMessage(e))
    if (length(dev.list()) > 0) dev.off()
  })
  
  # 2. Save High-Res PNG
  png_path <- paste0(base_filepath, ".png")
  tryCatch({
    png(png_path, width = w, height = h, units = "in", res = 300)
    print(plot_obj)
    dev.off()
    message("Saved PNG: ", png_path)
  }, error = function(e) {
    message("Warning: Failed to save PNG ", png_path, ": ", conditionMessage(e))
    if (length(dev.list()) > 0) dev.off()
  })
}

process_sample <- function(sample_name, sample_ident1, sample_ident2, base_data_dir, out_dirs) {
  # Wrap the ENTIRE function in tryCatch so one bad sample doesn't kill the parallel run
  Sys.time()
  tryCatch({
    output_rds <- file.path(out_dirs$processed, paste0(sample_name, "_processed.rds"))
    sample_id <- sample_name

    if(file.exists(output_rds)) {
      message("Loading preprocessed: ", sample_name)
      # Load existing object to ensure it has the processed_with_sct flag (for older files)
      tryCatch({
        existing_obj <- readRDS(output_rds)
        flag_present <- FALSE
        try({
          flag_present <- !is.null(existing_obj@misc$processed_with_sct)
        }, silent = TRUE)
        if (!flag_present) {
          try({
            existing_obj@misc$processed_with_sct <- isTRUE(opt$use_sct)
            saveRDS(existing_obj, output_rds)
            message("Updated existing processed RDS with 'processed_with_sct' flag for ", sample_name)
          }, silent = TRUE)
        }
      }, error = function(e) {
        message("Warning: failed to inspect/update existing RDS for ", sample_name, ": ", conditionMessage(e))
      })
      return(output_rds)  # Return path, not object
    }

    message("\n==================== Processing sample: ", sample_name, " ====================")
    data_dir <- file.path(base_data_dir, sample_name)

    # Create sample-specific plot directories
    sample_plot_dir <- file.path(out_dirs$plots, sample_name)
    dir.create(sample_plot_dir, recursive = TRUE, showWarnings = FALSE)

    # Create a Seurat object
    seur_obj <- CreateSeuratObject(
      counts = Read10X(data.dir = data_dir),
      project = sample_name,
      min.cells = 3,
      min.features = 10
    )
    message("Created Seurat object for ", sample_name, " with dimensions: ", dim(seur_obj)[1], " features and ", dim(seur_obj)[2], " cells")

    # ============ QC METRICS (BEFORE ANY FILTERING) ============
    message("Calculating QC metrics...")
    seur_obj[["percent.mt"]] <- PercentageFeatureSet(seur_obj, pattern = "^MT-")

    # QC violin plots BEFORE filtering
    pdf_file <- file.path(sample_plot_dir, paste0(sample_id, "_01_qc_violins_unfiltered"))
    tryCatch({
      qc_feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

      if (all(qc_feats %in% colnames(seur_obj@meta.data)) && ncol(seur_obj) > 0) {
        pdf(pdf_file, width = 15, height = 5)

        vln_plot <- tryCatch({
          create_qc_violin_plot(seur_obj, qc_feats, paste0("QC Metrics (Unfiltered) - ", sample_name))
        }, error = function(e) {
          message("Warning: QC plot creation failed: ", conditionMessage(e))
          return(NULL)
        })

        if (!is.null(vln_plot)) {
          print(vln_plot)
          message("QC violin plots (unfiltered) saved: ", pdf_file)
        }

        dev.off()
      } else {
        message("Warning: Cannot create VlnPlot - missing columns or empty object for ", sample_id)
      }
    }, error = function(e) {
      message("Warning: Error plotting QC violins for ", sample_id, ": ", conditionMessage(e))
      if (length(dev.list()) > 0) {
        tryCatch(dev.off(), error = function(e2) {})
      }
    })

    # ============ PRE-FILTERING (BEFORE DOUBLET FINDER) ============
    # This significantly speeds up DoubletFinder by removing obvious empty droplets
    message("Pre-filtering to remove empty droplets...")

    # Apply basic pre-filtering to remove obvious junk
    cells_before_prefilter <- ncol(seur_obj)
    seur_obj <- subset(seur_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 50)
    cells_after_prefilter <- ncol(seur_obj)

    message("Pre-filtering removed ", cells_before_prefilter - cells_after_prefilter,
            " empty droplets (", round((cells_before_prefilter - cells_after_prefilter) / cells_before_prefilter * 100, 2),
            "%), keeping ", cells_after_prefilter, " cells for DoubletFinder")

    # QC violin plots AFTER pre-filtering
    pdf_file <- file.path(sample_plot_dir, paste0(sample_id, "_02_qc_violins_prefiltered"))
    tryCatch({
      qc_feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

      if (all(qc_feats %in% colnames(seur_obj@meta.data)) && ncol(seur_obj) > 0) {
        pdf(pdf_file, width = 15, height = 5)

        vln_plot <- tryCatch({
          create_qc_violin_plot(seur_obj, qc_feats, paste0("QC Metrics (Pre-filtered) - ", sample_name))
        }, error = function(e) {
          message("Warning: QC plot creation failed: ", conditionMessage(e))
          return(NULL)
        })

        if (!is.null(vln_plot)) {
          print(vln_plot)
          message("QC violin plots (pre-filtered) saved: ", pdf_file)
        }

        dev.off()
      } else {
        message("Warning: Cannot create VlnPlot - missing columns or empty object for ", sample_id)
      }
    }, error = function(e) {
      message("Warning: Error plotting QC violins for ", sample_id, ": ", conditionMessage(e))
      if (length(dev.list()) > 0) {
        tryCatch(dev.off(), error = function(e2) {})
      }
    })

    # ============ NORMALIZATION & DIM REDUCTION ============
    message("Starting DoubletFinder workflow...")

    if (isTRUE(opt$use_sct)) {
      message("Running SCTransform v2 + glmGamPoi (regressing percent.mt only)...")
      seur_obj <- SCTransform(
        seur_obj, 
        vars.to.regress = "percent.mt", 
        method = "glmGamPoi", 
        vst.flavor = "v2", 
        verbose = FALSE
      )
      seur_obj <- RunPCA(seur_obj, assay = "SCT", npcs = 30, verbose = FALSE)
    } else {
      message("Running LogNormalize -> FindVariableFeatures -> ScaleData pipeline...")
      seur_obj <- seur_obj %>%
        NormalizeData() %>%
        FindVariableFeatures() %>%
        ScaleData(vars.to.regress = "percent.mt") %>%
        RunPCA(npcs = 30, verbose = FALSE)
    }

    # Plot elbow plot
    safe_save_plot(ElbowPlot(seur_obj, ndims = 30),
                  file.path(sample_plot_dir, paste0(sample_id, "_03_elbow")))

    # Find neighbors and clusters (needed for DoubletFinder)
    seur_obj <- seur_obj %>%
      FindNeighbors(dims = 1:30) %>%
      FindClusters() %>%
      RunUMAP(dims = 1:30, umap.method = "uwot", metric = "cosine")

    # Plot initial UMAP
    safe_save_plot(DimPlot(seur_obj, reduction = "umap", label = TRUE),
                  file.path(sample_plot_dir, paste0(sample_id, "_04_umap_initial")))

    # ============ DOUBLET FINDER ============
    # pK Identification (no ground-truth)
    message("Identifying optimal pK parameter...")
    sweep.res.list <- paramSweep(seur_obj, PCs = 1:20, sct = isTRUE(opt$use_sct))
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)

    safe_save_plot(ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
                    geom_point() + geom_line() +
                    ggtitle(paste0("pK Identification - ", sample_name)) +
                    theme_bw(),
                  file.path(sample_plot_dir, paste0(sample_id, "_05_doubletfinder_pk")))

    # Select optimal pK
    pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
    message("Optimal pK: ", pK)

    # Homotypic Doublet Proportion Estimate
    annotations <- seur_obj@meta.data$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations)
    nExp_poi <- round(opt$doublet_rate * nrow(seur_obj@meta.data)) 
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
    message("Expected doublets: ", nExp_poi, " (adjusted: ", nExp_poi.adj, ")")

    # Identify doublet cells
    message("Running DoubletFinder...")
    seur_obj <- doubletFinder(seur_obj, PCs = 1:20, pN = 0.25, pK = pK,
                              nExp = nExp_poi.adj, sct = isTRUE(opt$use_sct))

    DF_cols <- grep("DF.classifications", colnames(seur_obj@meta.data), value = TRUE)

    if (length(DF_cols) == 0) {
      stop("DoubletFinder failed to add classification column for sample: ", sample_name)
    }
    DF_classification <- DF_cols[1]
    message("DoubletFinder classification column: ", DF_classification)

    # Plot unfiltered doublet results
    safe_save_plot(DimPlot(seur_obj, reduction = 'umap', group.by = DF_classification, cols = doublet_color) +
                    ggtitle(paste0("Doublets (Before Removal) - ", sample_name)),
                  file.path(sample_plot_dir, paste0(sample_id, "_06_doublets_unfiltered")))

    # Check doublet statistics
    doublet_table <- table(seur_obj@meta.data[[DF_classification]])
    message("Doublet statistics:")
    print(doublet_table)

    # Remove doublet cells (keep only Singlets)
    singlet_indices <- which(seur_obj@meta.data[[DF_classification]] == "Singlet")
    seur_obj <- seur_obj[, singlet_indices]
    message("Cells after doublet removal: ", ncol(seur_obj))

    # Plot filtered doublet results
    safe_save_plot(DimPlot(seur_obj, reduction = "umap", group.by = DF_classification, cols = doublet_color) +
                    ggtitle(paste0("Doublets (After Removal) - ", sample_name)),
                  file.path(sample_plot_dir, paste0(sample_id, "_07_doublets_filtered")))

    # QC violin plots AFTER doublet removal
    pdf_file <- file.path(sample_plot_dir, paste0(sample_id, "_08_qc_violins_post_doublet"))
    tryCatch({
      qc_feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

      if (all(qc_feats %in% colnames(seur_obj@meta.data)) && ncol(seur_obj) > 0) {
        pdf(pdf_file, width = 15, height = 5)

        vln_plot <- tryCatch({
          create_qc_violin_plot(seur_obj, qc_feats, paste0("QC Metrics (Post-Doublet Removal) - ", sample_name))
        }, error = function(e) {
          message("Warning: QC plot creation failed: ", conditionMessage(e))
          return(NULL)
        })

        if (!is.null(vln_plot)) {
          print(vln_plot)
          message("QC violin plots (post-doublet) saved: ", pdf_file)
        }

        dev.off()
      } else {
        message("Warning: Cannot create VlnPlot - missing columns or empty object for ", sample_id)
      }
    }, error = function(e) {
      message("Warning: Error plotting QC violins for ", sample_id, ": ", conditionMessage(e))
      if (length(dev.list()) > 0) {
        tryCatch(dev.off(), error = function(e2) {})
      }
    })

    # ============ QUALITY CONTROL  ============
    message("Performing final quality control filtering...")

    # Add metadata identifiers
    seur_obj$orig.ident1 <- sample_ident1
    seur_obj$orig.ident2 <- sample_ident2

    # Quality filtering using configurable parameters
    message("Filtering cells with criteria: nFeature_RNA > ", opt$min_features,
            " & < ", opt$max_features, ", percent.mt < ", opt$max_mt)
    cells_before_final_qc <- ncol(seur_obj)
    seur_obj <- subset(seur_obj, subset = nFeature_RNA > opt$min_features &
                                          nFeature_RNA < opt$max_features &
                                          percent.mt < opt$max_mt)
    cells_after_final_qc <- ncol(seur_obj)
    message("Final QC filtering removed ", cells_before_final_qc - cells_after_final_qc,
            " cells, keeping ", cells_after_final_qc, " high-quality cells")

    # QC violin plots AFTER final filtering
    pdf_file <- file.path(sample_plot_dir, paste0(sample_id, "_09_qc_violins_final"))
    tryCatch({
      qc_feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

      if (all(qc_feats %in% colnames(seur_obj@meta.data)) && ncol(seur_obj) > 0) {
        pdf(pdf_file, width = 15, height = 5)

        vln_plot <- tryCatch({
          create_qc_violin_plot(seur_obj, qc_feats, paste0("QC Metrics (Final Filtered) - ", sample_name))
        }, error = function(e) {
          message("Warning: QC plot creation failed: ", conditionMessage(e))
          return(NULL)
        })

        if (!is.null(vln_plot)) {
          print(vln_plot)
          message("QC violin plots (final) saved: ", pdf_file)
        }

        dev.off()
      } else {
        message("Warning: Cannot create VlnPlot - missing columns or empty object for ", sample_id)
      }
    }, error = function(e) {
      message("Warning: Error plotting QC violins for ", sample_id, ": ", conditionMessage(e))
      if (length(dev.list()) > 0) {
        tryCatch(dev.off(), error = function(e2) {})
      }
    })

    # Save processed object
    tryCatch({
      seur_obj@misc$processed_with_sct <- isTRUE(opt$use_sct)
    }, error = function(e) {
      message("Warning: Unable to set 'processed_with_sct' flag in @misc: ", conditionMessage(e))
    })
    saveRDS(seur_obj, output_rds)
    message("Processed object saved: ", output_rds)
    message("==================== Completed: ", sample_name, " ====================\n")

    return(output_rds)

  }, error = function(e) {
    message("\nCRITICAL ERROR processing sample ", sample_name, ": ", e$message)
    return(NULL) 
  })
}

integrate_samples <- function(sample_list, chosen_res = 0.4) {
  Sys.time()
  # Remove NULLs from failed samples
  sample_list <- Filter(Negate(is.null), sample_list)

  if (length(sample_list) == 0) {
    stop("No valid samples available for integration.")
  }

  # Ensure we are working with Seurat objects (read file paths if present)
  sample_objs <- lapply(sample_list, function(elem) {
    if (is.character(elem) && file.exists(elem)) {
      readRDS(elem)
    } else {
      elem
    }
  })

  # If duplicated barcodes exist across samples, prefix barcodes with a sample-specific tag
  all_barcodes <- unlist(lapply(sample_objs, colnames))
  if (any(duplicated(all_barcodes))) {
    message("Detected duplicated cell barcodes across input samples. Prefixing barcodes with sample-specific tags.")
    sample_objs <- lapply(seq_along(sample_objs), function(i) {
      obj <- sample_objs[[i]]
      # Prefer orig.ident2, then orig.ident1, then project name, fallback to index
      tag <- NA
      try({ tag <- as.character(unique(obj$orig.ident2)[1]) }, silent = TRUE)
      if (is.null(tag) || is.na(tag) || nchar(tag) == 0) {
        try({ tag <- as.character(unique(obj$orig.ident1)[1]) }, silent = TRUE)
      }
      if (is.null(tag) || is.na(tag) || nchar(tag) == 0) {
        try({ tag <- as.character(obj@project.name) }, silent = TRUE)
      }
      if (is.null(tag) || is.na(tag) || nchar(tag) == 0) tag <- paste0("sample", i)
      tag_clean <- make.names(tag)
      colnames(obj) <- paste0(tag_clean, "_", colnames(obj))
      return(obj)
    })
  }

  # Replace sample_list with the prepared objects
  sample_list <- sample_objs

  # --- CONSISTENCY CHECK: ensure per-sample --use_sct matches integration setting ---
  # sample_list may contain Seurat objects (or, in other contexts, file paths). Infer processed_with_sct per sample.
  processed_sct_vec <- sapply(sample_list, function(elem) {
    obj <- elem
    if (is.character(obj) && file.exists(obj)) {
      obj <- readRDS(obj)
    }
    # Prefer explicit flag stored in @misc
    flag <- NULL
    try({
      flag <- obj@misc$processed_with_sct
    }, silent = TRUE)
    if (!is.null(flag)) {
      return(as.logical(flag))
    }
    # Fallback: infer from presence of "SCT" assay
    return("SCT" %in% names(obj))
  })

  unique_flags <- unique(processed_sct_vec)
  if (length(unique_flags) > 1) {
    stop("Inconsistent processing detected: some samples were processed with SCTransform and others were not. Reprocess samples to use a consistent --use_sct setting.")
  }

  inferred_use_sct <- as.logical(unique_flags)
  if (is.na(inferred_use_sct)) inferred_use_sct <- FALSE

  if (isTRUE(inferred_use_sct) && !isTRUE(opt$use_sct)) {
    message("Detected processed samples were created with SCTransform but --use_sct is FALSE. Overriding opt$use_sct -> TRUE for integration.")
    opt$use_sct <- TRUE
  } else if (!isTRUE(inferred_use_sct) && isTRUE(opt$use_sct)) {
    message("Detected processed samples were NOT created with SCTransform but --use_sct is TRUE. Overriding opt$use_sct -> FALSE for integration.")
    opt$use_sct <- FALSE
  }

  message("\n==================== Starting Harmony Integration ====================")

 # --- 1. Merge Samples (Required for Harmony) ---
  message("Merging ", length(sample_list), " samples into a single object...")

 # We will re-run SCTransform globally on the merged object (Best practice for Harmony).
  sample_list <- lapply(sample_list, function(obj) {
    DefaultAssay(obj) <- "RNA"
    if ("SCT" %in% names(obj)) {
      obj[["SCT"]] <- NULL
    }
    return(obj)
  })

  if (length(sample_list) > 1) {
    TN.combined <- merge(x = sample_list[[1]], 
                         y = sample_list[2:length(sample_list)])
  } else {
    TN.combined <- sample_list[[1]]
  }
  
  # Try to ensure we have a Seurat object and (when possible) join layers.
  # Some Seurat versions / unexpected input objects can lead to JoinLayers being called on
  # an Assay (or other) object which causes the pipeline to fail. Guard against that.
  if (!inherits(TN.combined, "Seurat")) {
    stop("Merged object is not a Seurat object. Class(es) found: ",
         paste(class(TN.combined), collapse = ", "),
         "\nCheck that each entry in sample_list is a Seurat object (readRDS of processed samples).")
  }

  # Attempt JoinLayers only if the function is available; if it fails, log and continue.
  join_layers_available <- exists("JoinLayers", where = asNamespace("Seurat"), mode = "function")
  if (isTRUE(join_layers_available)) {
    tryCatch({
      TN.combined <- JoinLayers(TN.combined)
    }, error = function(e) {
      message("Warning: JoinLayers failed with message: ", conditionMessage(e))
      message("Proceeding without JoinLayers. This may be fine depending on your Seurat version;",
              " if you encounter downstream issues, consider updating Seurat or inspecting sample objects.")
    })
  } else {
    message("JoinLayers() not found in the loaded Seurat namespace; skipping JoinLayers step.")
  }

  # Final sanity check
  if (!inherits(TN.combined, "Seurat")) {
    stop("After attempted JoinLayers fallback, combined object is not a Seurat object. Aborting.")
  }

  # =========================================================================
  # --- CRITICAL FIX: Add Cell Cycle Scoring back to the merged object safely ---
  # =========================================================================
  message("Preparing merged object for Cell Cycle Scoring...")
  
  # Ensure basic RNA assay is active for normalization
  DefaultAssay(TN.combined) <- "RNA"
  
  # Normalize RNA (Required before CellCycleScoring can compute S/G2M scores)
  tryCatch({
    TN.combined <- NormalizeData(TN.combined, assay = "RNA", verbose = FALSE)
  }, error = function(e) {
    message("Warning: NormalizeData failed (skipping): ", conditionMessage(e))
  })

  # Cell cycle scoring: run safely only if genes are present
  if (exists("cc.genes") && length(cc.genes$s.genes) > 0 && length(cc.genes$g2m.genes) > 0) {
    present_s <- intersect(cc.genes$s.genes, rownames(TN.combined))
    present_g2m <- intersect(cc.genes$g2m.genes, rownames(TN.combined))
    
    if (length(present_s) > 0 && length(present_g2m) > 0) {
      message("Running CellCycleScoring with ", length(present_s), " S genes and ", length(present_g2m), " G2M genes...")
      tryCatch({
        TN.combined <- CellCycleScoring(TN.combined, s.features = present_s, g2m.features = present_g2m, set.ident = FALSE)
        message("Successfully added S.Score and G2M.Score to merged object.")
      }, error = function(e) {
         message("Warning: CellCycleScoring failed on merged object (skipping): ", conditionMessage(e))
      })
    } else {
      message("Skipping CellCycleScoring: canonical cc.genes were not found in the merged object.")
    }
  } else {
    message("Skipping CellCycleScoring: 'cc.genes' list is not available in the environment.")
  }
  # =========================================================================

  # --- 2. Global Normalization & Cell Cycle Scoring ---
  message("Preparing merged object for PCA/Harmony...")
  
  # Create an explicit batch column from the requested metadata (safer for Harmony)
  if (!(opt$batch_var %in% colnames(TN.combined@meta.data))) {
    message("Warning: requested batch_var '", opt$batch_var, "' not found in meta.data; falling back to orig.ident2 if present.")
    if ("orig.ident2" %in% colnames(TN.combined@meta.data)) {
      TN.combined$batch <- TN.combined$orig.ident2
    } else {
      TN.combined$batch <- TN.combined$orig.ident1
    }
  } else {
    TN.combined$batch <- TN.combined[[opt$batch_var]]
  }
  message("Using metadata column 'batch' for Harmony grouping (derived from ", opt$batch_var, ")")

  if (isTRUE(opt$use_sct)) {
    message("Running Global SCTransform v2 on the merged object (Best practice for Harmony)...")
    
    # Check if cell cycle scores merged successfully from the per-sample step
    vars_regress <- c("percent.mt")
    if ("S.Score" %in% colnames(TN.combined@meta.data) && "G2M.Score" %in% colnames(TN.combined@meta.data)) {
      vars_regress <- c("percent.mt", "S.Score", "G2M.Score")
      message("Regressing out percent.mt, S.Score, and G2M.Score...")
    } else {
      message("Regressing out percent.mt only...")
    }

    TN.combined <- SCTransform(
      TN.combined, 
      assay = "RNA",
      vars.to.regress = vars_regress, 
      method = "glmGamPoi", 
      vst.flavor = "v2", 
      verbose = FALSE
    )
    
    message("Preparing global SCT model for downstream differential expression...")
    TN.combined <- PrepSCTFindMarkers(TN.combined, assay = "SCT", verbose = FALSE)

    message("Running PCA on unified SCT assay...")
    TN.combined <- RunPCA(TN.combined, assay = "SCT", npcs = 50, verbose = FALSE)
    harmony_assay <- "SCT"
  } else {
    message("Running global LogNormalize/ScaleData on merged RNA assay (no per-sample SCT requested).")
    TN.combined <- NormalizeData(TN.combined, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
    TN.combined <- FindVariableFeatures(TN.combined, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    TN.combined <- ScaleData(TN.combined, vars.to.regress = c("percent.mt"), features = rownames(TN.combined), verbose = TRUE)
    TN.combined <- RunPCA(TN.combined, npcs = 50, verbose = FALSE)
    harmony_assay <- "RNA"
  }

  safe_save_plot(ElbowPlot(TN.combined, ndims = 50),
                file.path(output_dirs$plots, "TNcombined_elbow"))

  # --- 4. HARMONY INTEGRATION ---
  message("Running Harmony Integration on PCA embedding (assay: ", harmony_assay, ") ...")
  # Use the explicit 'batch' column created above
  TN.combined <- RunHarmony(TN.combined, group.by.vars = "batch", verbose = FALSE)
  # --- 5. UMAP & Multi-Resolution Clustering ---
  message("Running UMAP and FindNeighbors on Harmony reduction...")
  # Use R-native UWOT implementation and cosine metric to avoid reticulate warnings
  TN.combined <- RunUMAP(TN.combined, reduction = "harmony", dims = 1:30, verbose = FALSE, umap.method = "uwot", metric = "cosine")
  TN.combined <- FindNeighbors(TN.combined, reduction = "harmony", dims = 1:30, verbose = FALSE)

  message("Running multiple clustering resolutions...")
  resolutions <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8)
  TN.combined <- FindClusters(TN.combined, resolution = resolutions, verbose = FALSE)

  # Clustering columns depend on the assay used for clustering (SCT -> SCT_snn_res., RNA -> RNA_snn_res.)
  cluster_prefix <- paste0(harmony_assay, "_snn_res.")

  # --- 6. Clustree Generation ---
  message("Generating Clustree visualization...")
  p_tree <- clustree(TN.combined, prefix = cluster_prefix, node_text_size = 3, edge_arrow = FALSE) +
    ggtitle("Global Dataset - Clustree Resolution Tracker") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  safe_save_plot(p_tree, file.path(output_dirs$plots, "TNcombined_clustree"), w = 15, h = 10)

  # --- 7. Finalize Object ---
  # Set a default resolution to ensure downstream plotting scripts don't break
  default_res <- chosen_res
  default_res_col <- paste0(cluster_prefix, default_res)
  
  if(default_res_col %in% colnames(TN.combined@meta.data)) {
    Idents(TN.combined) <- default_res_col
    TN.combined$seurat_clusters <- TN.combined[[default_res_col]]
  } else {
    message("Warning: Resolution ", default_res, " not found. Using highest available.")
  }

  message("Saving cell count tables...")
  CellNumber <- table(Idents(TN.combined), TN.combined$orig.ident1)
  write.csv(CellNumber, file = file.path(output_dirs$tables, "CellNumber_bygroup.csv"))

  message("Saving globally integrated object...")
  saveRDS(TN.combined, file = file.path(output_base, "TN.combined_dim30.rds"))
  message("==================== Harmony Integration complete! ====================\n")

  return(TN.combined)
}

generate_plots <- function(TN.combined, chosen_res = 0.4) {
  Sys.time()
  message("\n==================== Generating Plots ====================")

  # --- Determine clustering prefix dynamically (SCT vs RNA) ---
  md_cols <- colnames(TN.combined@meta.data)
  if (any(startsWith(md_cols, "SCT_snn_res."))) {
    cluster_prefix <- "SCT_snn_res."
  } else if (any(startsWith(md_cols, "RNA_snn_res."))) {
    cluster_prefix <- "RNA_snn_res."
  } else {
    # Fallback: if neither prefix present, prefer SCT when use_sct TRUE, otherwise RNA
    cluster_prefix <- if (isTRUE(opt$use_sct)) "SCT_snn_res." else "RNA_snn_res."
    message("Warning: Could not find existing *_snn_res.* columns. Falling back to prefix: ", cluster_prefix)
  }

  res_col <- paste0(cluster_prefix, chosen_res)
  if (res_col %in% md_cols) {
    Idents(TN.combined) <- res_col
    message("Successfully switched active clusters to: ", res_col)
  } else {
    stop("Error: Resolution ", chosen_res, " was not calculated during the integrate step! Checked prefix: ", cluster_prefix)
  }

  # --- 2. Attempt to Apply Annotations ---
  # Look for the consensus file in the standard location relative to output_base
  consensus_file <- file.path(
  output_base,
  "annotations",
  paste0("res_", res_folder),
  "consensus",
  "consensus_annotation.tsv"
  )

  if (!file.exists(consensus_file) && !is.null(opt$output)) {
     consensus_file <- file.path(opt$output, "annotations", "consensus", "consensus_annotation.tsv")
  }

  # NOW create the plotting copy (it will inherit the new resolution!)
  TN.plotting <- TN.combined
  is_annotated <- FALSE

  if (file.exists(consensus_file)) {
    message("Found consensus annotation file. Applying labels...")
    consensus_data <- read.delim(consensus_file, sep = "\t", stringsAsFactors = FALSE)
    
    consensus_data$Cluster <- as.character(gsub("^X", "", consensus_data$Cluster))
    
    current_ids <- levels(TN.combined)
    new_names <- character(length(current_ids))
    names(new_names) <- current_ids
    
    # --- FLEXIBLE COLUMN DETECTION ---
    possible_cols <- c("Top_Cell_Type", "Cell_Type", "cell_type", "celltype", "Annotation", "annotation")
    actual_col <- intersect(possible_cols, colnames(consensus_data))[1]
    
    if (is.na(actual_col) && ncol(consensus_data) >= 2) {
      # Fallback to the second column if standard names aren't found
      actual_col <- colnames(consensus_data)[2]
      message(paste("Warning: Annotation column name not recognized. Defaulting to column:", actual_col))
    }
    
    # --- RENAME LOOP (With Global Shortening) ---
    for (id in current_ids) {
      clean_id <- as.character(gsub("^X", "", id))
      match_row <- consensus_data[consensus_data$Cluster == clean_id, ]
      
      if (nrow(match_row) > 0 && !is.na(actual_col) && !is.na(match_row[[actual_col]][1])) {
        
        # Extract the raw string, then delete the semicolon and everything after it
        raw_name <- as.character(match_row[[actual_col]][1])
        short_name <- trimws(gsub(";.*", "", raw_name))
        
        new_names[id] <- paste0(clean_id, ": ", short_name)
      } else {
        new_names[id] <- paste0(clean_id, ": Unknown")
      }
    }
    
    TN.plotting <- RenameIdents(TN.combined, new_names)
    is_annotated <- TRUE
    
    # --- FIX: Mathematically sort the new factor levels ---
    # Prevents R from alphabetically sorting "10: Cell" before "2: Cell"
    current_levels <- levels(TN.plotting)
    level_nums <- as.numeric(sub(":.*", "", current_levels))
    sorted_levels <- current_levels[order(level_nums)]
    Idents(TN.plotting) <- factor(Idents(TN.plotting), levels = sorted_levels)
    # ------------------------------------------------------
    
  } else {
    message("Consensus file not found. Using numeric cluster labels.")
  }

  n_groups <- length(levels(TN.plotting))
  if (n_groups <= length(mycolor)) {
    plot_colors <- mycolor[1:n_groups]
  } else {
    plot_colors <- colorRampPalette(mycolor)(n_groups)
  }

  # --- FIX: Force 'RNA' assay (LogNormalize) for DE and visualization ---
  message("Preparing RNA assay for marker/DE calls (NormalizeData on RNA)...")
  DefaultAssay(TN.plotting) <- "RNA"
  
  # Attempt JoinLayers for Seurat v5 compatibility
  try({ TN.plotting <- JoinLayers(TN.plotting) }, silent = TRUE)
  
  # Force standard LogNormalize for visualization and marker detection
  TN.plotting <- NormalizeData(TN.plotting, assay = "RNA", normalization.method = "LogNormalize", verbose = FALSE)
  de_assay <- "RNA"
  # ----------------------------------------------------------------------

  # --- 2. UMAP Plots ---
  # Note: label.size set to 3 to accommodate longer text labels
  
  # Create a dynamic title string
  res_title <- paste("Global Integration (Resolution:", chosen_res, ")")
  p_umap_f <- DimPlot(TN.plotting, reduction = "umap", label = FALSE, pt.size = 0.8, cols = plot_colors) +
    ggtitle(paste(res_title, "- Unlabeled")) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  safe_save_plot(p_umap_f, file.path(output_dirs$plots, "TNcombined_umap_labelF"))

  p_umap_t <- DimPlot(TN.plotting, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE, pt.size = 0.8, cols = plot_colors) +
    ggtitle(paste(res_title, "- Labeled")) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  safe_save_plot(p_umap_t, file.path(output_dirs$plots, "TNcombined_umap_labelT"))

  p_umap_orig <- DimPlot(TN.plotting, group.by = "orig.ident2", pt.size = 0.8, cols = plot_colors) +
    ggtitle("UMAP by Sample") + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  safe_save_plot(p_umap_orig, file.path(output_dirs$plots, "TNcombined_umap_groupbyorigident"))

  p_umap_split1 <- DimPlot(TN.plotting, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE, split.by = "orig.ident1", pt.size = 0.8, ncol = 2, cols = plot_colors) +
    ggtitle(paste(res_title, "- Split by Condition")) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  safe_save_plot(p_umap_split1, file.path(output_dirs$plots, "TNcombined_umap_labelT_splitorigident1"))

  p_umap_split2 <- DimPlot(TN.plotting, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE, split.by = "orig.ident2", pt.size = 0.8, ncol = 2, cols = plot_colors) +
    ggtitle(paste(res_title, "- Split by Sample")) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  safe_save_plot(p_umap_split2, file.path(output_dirs$plots, "TNcombined_umap_labelT_splitsample"))

  message("Generating Before/After Harmony Integration Plots...")

  # 1. PCA Plot (Before Integration)
  p_pca_before <- DimPlot(TN.plotting, reduction = "pca", group.by = "orig.ident2", pt.size = 0.5, cols = plot_colors) +
    ggtitle("PCA (Before Harmony)") + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  safe_save_plot(p_pca_before, file.path(output_dirs$plots, "TNcombined_01_PCA_Before_Harmony"))

  # 2. Harmony Plot (After Integration)
  p_harmony_after <- DimPlot(TN.plotting, reduction = "harmony", group.by = "orig.ident2", pt.size = 0.5, cols = plot_colors) +
    ggtitle("Harmony Embeddings (After Integration)") + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  safe_save_plot(p_harmony_after, file.path(output_dirs$plots, "TNcombined_02_Harmony_After"))

  # 3. Side-by-Side Comparison Plot
  # Uses patchwork library (already loaded in your script) to stitch them together
  p_integration_comparison <- p_pca_before + p_harmony_after
  safe_save_plot(p_integration_comparison, file.path(output_dirs$plots, "TNcombined_03_Integration_Comparison"), w = 16, h = 8)

  # --- 3. Heatmap ---
  CellNumber <- table(Idents(TN.plotting), TN.plotting$orig.ident1)
  cluster_count <- nrow(CellNumber)

  message("Finding markers for heatmap...")
  Heatmapall <- subset(TN.plotting, idents = levels(TN.plotting))

  # Safely run FindAllMarkers: capture errors/warnings and handle empty results gracefully
  Heatmapall.markers <- tryCatch({
    suppressWarnings(
      FindAllMarkers(
        Heatmapall,
        assay = de_assay,
        only.pos = TRUE,
        min.pct = 0.1,
        logfc.threshold = 0.25,
        verbose = FALSE
      )
    )
  }, error = function(e) {
    message("Error running FindAllMarkers: ", conditionMessage(e))
    return(NULL)
  })

  # If no markers found or function failed, skip heatmap/marker-based plotting
  if (is.null(Heatmapall.markers) || nrow(Heatmapall.markers) == 0) {
    message("No DE markers detected by FindAllMarkers. Skipping heatmap and marker-dependent plots.")
  } else {
    # Save full markers table
    write.csv(Heatmapall.markers, file = file.path(output_dirs$tables, "Findallmarkers.csv"), row.names = FALSE)

    # Determine which column names contain cluster labels (robust across Seurat versions)
    cluster_col <- if ("cluster" %in% colnames(Heatmapall.markers)) {
      "cluster"
    } else if ("group" %in% colnames(Heatmapall.markers)) {
      "group"
    } else {
      NA
    }

    if (is.na(cluster_col)) {
      message("FindAllMarkers returned results but no 'cluster' or 'group' column was found. Skipping heatmap generation.")
    } else {
      # Compute top genes per cluster (use slice_max which is the modern replacement for top_n)
      top_genes <- Heatmapall.markers %>%
        dplyr::group_by_at(vars(all_of(cluster_col))) %>%
        dplyr::slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
        dplyr::pull(gene) %>%
        unique()

      if (length(top_genes) == 0) {
        message("No top marker genes identified after grouping. Skipping heatmap.")
      } else {
        message("Scaling data for heatmap features...")
        # Make sure the object is scaled on the correct assay before plotting
        DefaultAssay(Heatmapall) <- de_assay 
        Heatmapall <- ScaleData(Heatmapall, features = top_genes, assay = de_assay, verbose = FALSE)
        
        p_heat <- DoHeatmap(Heatmapall, features = top_genes, group.colors = plot_colors, assay = de_assay) +          scale_fill_gradient2(low = "magenta", mid = "black", high = "yellow", midpoint = 0, name = "Z-Score") + 
          ggtitle(paste("Top 10 Markers per Cluster (Resolution:", chosen_res, ")")) +
          theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
          
        safe_save_plot(p_heat, file.path(output_dirs$plots, "heatmap_top10"), w = 25, h = 25)
      }
    }
  }

 # --- 4. Proportion Plots ---
  message("Generating stacked barplots with contained micro-labels...")
  
  # By Condition (ident1)
  Cellproportion <- prop.table(table(Idents(TN.plotting), TN.plotting$orig.ident1), margin = 2)
  write.csv(Cellproportion * 100, file = file.path(output_dirs$tables, "Cellproportion.csv"), row.names = TRUE)

  Cellproportion_df <- as.data.frame(Cellproportion)
  colnames(Cellproportion_df) <- c("Cluster", "Condition", "Proportion")
  Cellprop_plot <- Cellproportion_df %>% filter(Proportion > 0)
  
  proportion_plot <- ggplot(Cellprop_plot, aes(x = Condition, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity", color = "white", linewidth = 0.2) + 
    # Print all labels, regardless of how small the slice is
    geom_text(aes(label = as.character(Cluster), size = Proportion), 
              position = position_stack(vjust = 0.5), color = "black") +
    # TRICK: Stretch the range down to 0.5. 1% slices will get 0.5 size text (tiny specks).
    scale_size_continuous(range = c(0.5, 4), guide = "none") + 
    theme_minimal(base_size = 15) +
    xlab("Condition") +
    ylab("Percentage of Total Cells") +
    scale_y_continuous(labels = scales::percent) +
    theme(legend.title = element_blank(), 
          legend.text = element_text(size = 8)) +
    scale_fill_manual(values = plot_colors)

  safe_save_plot(proportion_plot, file.path(output_dirs$plots, "proportion_plot"), w = 10, h = 8)

  # By Sample (ident2)
  Sampleproportion <- prop.table(table(Idents(TN.plotting), TN.plotting$orig.ident2), margin = 2)
  write.csv(Sampleproportion * 100, file = file.path(output_dirs$tables, "Sampleproportion.csv"), row.names = TRUE)

  Sampleproportion_df <- as.data.frame(Sampleproportion)
  colnames(Sampleproportion_df) <- c("Cluster", "Sample", "Proportion")
  Sampleprop_plot <- Sampleproportion_df %>% filter(Proportion > 0)

  sample_proportion_plot <- ggplot(Sampleprop_plot, aes(x = Sample, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity", color = "white", linewidth = 0.2) + 
    # Print all labels, regardless of how small the slice is
    geom_text(aes(label = as.character(Cluster), size = Proportion), 
              position = position_stack(vjust = 0.5), color = "black") +
    # TRICK: Stretch the range down to 0.5. 1% slices will get 0.5 size text (tiny specks).
    scale_size_continuous(range = c(0.5, 4), guide = "none") + 
    theme_minimal(base_size = 15) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    xlab("Sample") +
    ylab("Percentage of Total Cells") +
    scale_y_continuous(labels = scales::percent) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 8)) +
    scale_fill_manual(values = plot_colors)

  safe_save_plot(sample_proportion_plot, file.path(output_dirs$plots, "sample_proportion_plot"), w = 12, h = 8)


  # --- 5. Proportion Plots (Variance & Distribution) ---
  message("Generating all 4 faceted proportion plots (Fixed and Free axes)...")
  
  # Calculate raw counts per cluster, per sample, per condition
  prop_data <- as.data.frame(table(Idents(TN.plotting), TN.plotting$orig.ident2, TN.plotting$orig.ident1))
  colnames(prop_data) <- c("Cluster", "Sample", "Condition", "Count")
  
  # Remove phantom combinations (where a sample doesn't belong to a condition)
  prop_data <- prop_data %>% filter(Count > 0 | Condition == TN.plotting$orig.ident1[match(Sample, TN.plotting$orig.ident2)])
  
  # Calculate percentages per sample
  prop_data <- prop_data %>%
    group_by(Sample) %>%
    mutate(Percentage = (Count / sum(Count)) * 100) %>%
    ungroup()
    
  # Shorten the cluster names for cleaner facet headers
  prop_data$Cluster <- as.character(prop_data$Cluster)
  prop_data$Cluster <- trimws(gsub(";.*", "", prop_data$Cluster))
  
  # --- NUMERICAL SORTING TRICK ---
  # Force R to sort mathematically (0,1,2,3...) instead of alphabetically (0,1,10,11...)
  unique_clusters <- unique(prop_data$Cluster)
  cluster_nums <- as.numeric(sub(":.*", "", unique_clusters)) # Grab the number before the colon
  sorted_clusters <- unique_clusters[order(cluster_nums)]     # Sort based on that number
  prop_data$Cluster <- factor(prop_data$Cluster, levels = sorted_clusters) # Lock in the correct order
    
  # --- PLOT 1: Boxplot by Condition (FIXED Axis) ---
  box_plot_cond_fixed <- ggplot(prop_data, aes(x = Condition, y = Percentage, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.2, size = 2, color = "black", alpha = 0.8) +
    facet_wrap(~ Cluster, scales = "fixed") +
    theme_bw(base_size = 14) +
    labs(title = "Cluster Proportions Across Conditions (Fixed Axis)", x = "Condition", y = "Percentage of Cells (%)") +
    theme(legend.position = "none", 
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          strip.text = element_text(face = "bold", size = 10), 
          strip.background = element_rect(fill = "lightgray"))
          
  safe_save_plot(box_plot_cond_fixed, file.path(output_dirs$plots, "proportion_boxplot_condition_FIXED"), w = 16, h = 12)

  # --- PLOT 2: Boxplot by Condition (FREE Axis) ---
  box_plot_cond_free <- ggplot(prop_data, aes(x = Condition, y = Percentage, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.2, size = 2, color = "black", alpha = 0.8) +
    facet_wrap(~ Cluster, scales = "free_y") +
    theme_bw(base_size = 14) +
    labs(title = "Cluster Proportions Across Conditions (Free Axis)", x = "Condition", y = "Percentage of Cells (%)") +
    theme(legend.position = "none", 
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          strip.text = element_text(face = "bold", size = 10), 
          strip.background = element_rect(fill = "lightgray"))
          
  safe_save_plot(box_plot_cond_free, file.path(output_dirs$plots, "proportion_boxplot_condition_FREE"), w = 16, h = 12)

  # --- PLOT 3: Column Plot by Sample (FIXED Axis) ---
  bar_plot_samp_fixed <- ggplot(prop_data, aes(x = Sample, y = Percentage, fill = Sample)) +
    geom_col(color = "black", alpha = 0.8) +
    facet_wrap(~ Cluster, scales = "fixed") +
    theme_bw(base_size = 14) +
    labs(title = "Cluster Proportions per Sample (Fixed Axis)", x = "Sample", y = "Percentage of Cells (%)") +
    theme(legend.position = "none", 
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          strip.text = element_text(face = "bold", size = 10), 
          strip.background = element_rect(fill = "lightgray"))
          
  safe_save_plot(bar_plot_samp_fixed, file.path(output_dirs$plots, "proportion_barplot_sample_FIXED"), w = 16, h = 12)

  # --- PLOT 4: Column Plot by Sample (FREE Axis) ---
  bar_plot_samp_free <- ggplot(prop_data, aes(x = Sample, y = Percentage, fill = Sample)) +
    geom_col(color = "black", alpha = 0.8) +
    facet_wrap(~ Cluster, scales = "free_y") +
    theme_bw(base_size = 14) +
    labs(title = "Cluster Proportions per Sample (Free Axis)", x = "Sample", y = "Percentage of Cells (%)") +
    theme(legend.position = "none", 
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          strip.text = element_text(face = "bold", size = 10), 
          strip.background = element_rect(fill = "lightgray"))
          
  safe_save_plot(bar_plot_samp_free, file.path(output_dirs$plots, "proportion_barplot_sample_FREE"), w = 16, h = 12)
  
  message("==================== Plots completed! ====================\n")
}

run_and_plot_pca <- function(seurat_object, output_plots_dir) {
  Sys.time()
  message("Generating sample similarity plots (Expression Pseudo-bulk & Embedding Centroids)...")

  # ========================================================================
  # 1. ORIGINAL: Pseudo-bulk RNA PCA (Raw Biological/Technical Variance)
  # ========================================================================
  avg_expr <- AggregateExpression(
    seurat_object,
    group.by = "orig.ident2",
    assays = "RNA",
    normalization.method = "LogNormalize",
    return.seurat = FALSE
  )

  avg_expr_matrix <- t(avg_expr$RNA)
  gene_vars <- apply(avg_expr_matrix, 2, var)
  avg_expr_matrix_filtered <- avg_expr_matrix[, gene_vars > 0]
  
  pca_results <- prcomp(avg_expr_matrix_filtered, scale. = TRUE)

  pca_data <- as.data.frame(pca_results$x)
  pca_data$sample <- rownames(pca_data)
  pca_variance <- pca_results$sdev^2 / sum(pca_results$sdev^2)

  p_rna_pca <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = sample), size = 6, alpha = 0.9) +
    geom_text_repel(aes(label = sample), size = 4, box.padding = 0.5) +
    theme_bw() +
    guides(color = "none") +
    ggtitle("PCA of Sample Similarity\n(RNA Pseudo-bulk)") +
    labs(
      x = paste0("PC1 (", round(pca_variance[1] * 100, 2), "%)"),
      y = paste0("PC2 (", round(pca_variance[2] * 100, 2), "%)")
    ) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

  safe_save_plot(p_rna_pca, file.path(output_plots_dir, "pca_sample_similarity_01_RNA"), w = 8, h = 7)

  # ========================================================================
  # 2. NEW: Embedding Centroids (Before vs After Harmony Comparison)
  # ========================================================================
  meta <- seurat_object@meta.data
  
  if ("pca" %in% names(seurat_object@reductions) && "harmony" %in% names(seurat_object@reductions)) {
    
    # Extract the first 2 dimensions of each cell's spatial embeddings
    pca_embed <- Embeddings(seurat_object, "pca")[, 1:2]
    harmony_embed <- Embeddings(seurat_object, "harmony")[, 1:2]
    
    # Calculate the centroid (average spatial position) for each sample inside the PCA space
    pca_centroids <- aggregate(pca_embed, by = list(sample = meta$orig.ident2), FUN = mean)
    colnames(pca_centroids)[2:3] <- c("Dim1", "Dim2")
    
    # Calculate the centroid for each sample inside the Harmony space
    harmony_centroids <- aggregate(harmony_embed, by = list(sample = meta$orig.ident2), FUN = mean)
    colnames(harmony_centroids)[2:3] <- c("Dim1", "Dim2")
    
    # Plot 1: Uncorrected PCA Centroids
    p_embed_pca <- ggplot(pca_centroids, aes(x = Dim1, y = Dim2)) +
      geom_point(aes(color = sample), size = 6, alpha = 0.9) +
      geom_text_repel(aes(label = sample), size = 4, box.padding = 0.5) +
      theme_bw() +
      guides(color = "none") +
      ggtitle("Sample Centroids\nBefore Harmony (PCA Space)") +
      labs(x = "PC 1", y = "PC 2") +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

    # Plot 2: Corrected Harmony Centroids
    p_embed_harmony <- ggplot(harmony_centroids, aes(x = Dim1, y = Dim2)) +
      geom_point(aes(color = sample), size = 6, alpha = 0.9) +
      geom_text_repel(aes(label = sample), size = 4, box.padding = 0.5) +
      theme_bw() +
      guides(color = "none") +
      ggtitle("Sample Centroids\nAfter Harmony (Harmony Space)") +
      labs(x = "Harmony 1", y = "Harmony 2") +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

    # Stitch them together using patchwork
    p_centroid_comparison <- p_embed_pca + p_embed_harmony
    
    safe_save_plot(p_centroid_comparison, file.path(output_plots_dir, "pca_sample_similarity_02_Before_vs_After_Harmony"), w = 14, h = 7)
    message("Before vs After Harmony similarity comparison saved.")
    
  } else {
    message("Warning: 'pca' or 'harmony' reductions not found. Skipping centroid comparison.")
  }
}

create_summary_dot_plot <- function(TN.combined, output_plots_dir, top_n = 5) {
  # Ensure plots dir exists
  dir.create(output_plots_dir, recursive = TRUE, showWarnings = FALSE)

  # Try to reuse previously computed markers if available
  markers_file <- file.path(output_dirs$tables, "Findallmarkers.csv")
  markers <- NULL
  if (file.exists(markers_file)) {
    try({
      markers <- read.csv(markers_file, stringsAsFactors = FALSE)
      message("Loaded markers from ", markers_file)
    }, silent = TRUE)
  }

  # If markers not loaded, compute them (may be slow)
  # If markers not loaded, compute them (may be slow)
  if (is.null(markers)) {
    message("Computing FindAllMarkers (may take a while)...")
    tryCatch({
      # --- FIX: Force 'RNA' assay for fallback marker computation ---
      DefaultAssay(TN.combined) <- "RNA"
      try({ TN.combined <- JoinLayers(TN.combined) }, silent = TRUE)
      TN.combined <- NormalizeData(TN.combined, assay = "RNA", normalization.method = "LogNormalize", verbose = FALSE)
      de_assay <- "RNA"
      
      markers <- FindAllMarkers(TN.combined, assay = de_assay, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)

      write.csv(markers, file = markers_file, row.names = FALSE)
      message("Saved FindAllMarkers to ", markers_file)
    }, error = function(e) {
      message("Failed to compute markers: ", conditionMessage(e))
      return(NULL)
    })
  }

  if (is.null(markers) || nrow(markers) == 0) {
    message("No marker genes available to plot.")
    return(invisible(NULL))
  }

  # Select top N genes per cluster
  library(dplyr)
  top_genes <- markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = top_n, with_ties = FALSE) %>%
    pull(gene) %>%
    unique()

  if (length(top_genes) == 0) {
    message("No top genes found for DotPlot.")
    return(invisible(NULL))
  }

  plot_assay <- "RNA"
  DefaultAssay(TN.combined) <- plot_assay

  p_dot <- tryCatch({
    DotPlot(TN.combined, features = top_genes, assay = plot_assay) +
      ggtitle(paste0("Top ", top_n, " markers per cluster (DotPlot)")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }, error = function(e) {
    message("Failed to create DotPlot: ", conditionMessage(e))
    return(NULL)
  })

  if (!is.null(p_dot)) {
    safe_save_plot(p_dot, file.path(output_plots_dir, "summary_dotplot"), w = 14, h = 8)
    message("Saved summary dot plot to: ", output_plots_dir)
  }
  invisible(p_dot)
}

execute_step <- function(step) {
  switch(step,
    read_csv = {
      read_samples_csv(opt$file)
    },
    process = {
      samples_df <- readRDS(file.path(output_base, "samples_df.rds"))
      message("\n==================== Checking Sample Processing Status ====================")

      # Check which samples already have processed RDS files
      existing_samples <- character()
      new_samples <- character()

      for (i in 1:nrow(samples_df)) {
        sample_name <- samples_df$sample_names[i]
        output_rds <- file.path(output_dirs$processed, paste0(sample_name, "_processed.rds"))

        if (file.exists(output_rds)) {
          existing_samples <- c(existing_samples, sample_name)
        } else {
          new_samples <- c(new_samples, sample_name)
        }
      }

      if (length(existing_samples) > 0) {
        message("Found ", length(existing_samples), " already processed samples:")
        message("  ", paste(existing_samples, collapse=", "))
      }

      if (length(new_samples) > 0) {
        message("Will process ", length(new_samples), " new samples:")
        message("  ", paste(new_samples, collapse=", "))
      } else {
        message("All samples already processed! Skipping processing step.")
        return(invisible(NULL))
      }

      message("\n==================== Starting Sample Processing ====================")

      # Process samples sequentially using lapply
      results <- lapply(1:nrow(samples_df), function(i) {
        process_sample(samples_df$sample_names[i],
                      samples_df$ident1[i],
                      samples_df$ident2[i],
                      opt$datadir,
                      output_dirs)
      })

      # Check for failures
      failed_indices <- which(sapply(results, is.null))
      if (length(failed_indices) > 0) {
        message("\nWARNING: The following samples FAILED to process:")
        message(paste(samples_df$sample_names[failed_indices], collapse=", "))
        message("Check the logs for details.")
      } else {
        message("All samples processed successfully")
      }
      invisible(results)
    },
    integrate = {
      # Check if integrated object already exists
      integrated_rds <- file.path(output_base, "TN.combined_dim30.rds")

      if (file.exists(integrated_rds)) {
        message("Loading existing integrated object from: ", integrated_rds)
        TN_combined <- readRDS(integrated_rds)
        message("Integrated object loaded successfully with ", ncol(TN_combined), " cells")
        return(TN_combined)
      }

      # If not, proceed with integration
      sample_files <- list.files(output_dirs$processed,
                               pattern = "_processed.rds$",
                               full.names = TRUE)

      if (length(sample_files) == 0) stop("No processed sample files found in ", output_dirs$processed)

      message("Found ", length(sample_files), " processed samples for integration")
      sample_list <- lapply(sample_files, readRDS)
      integrate_samples(sample_list, chosen_res = opt$resolution)
    },
    plot = {
      TN_combined <- readRDS(file.path(output_base, "TN.combined_dim30.rds"))
      
      # --- NEW: Dynamically update output directories for specific resolutions ---
      res_folder <- paste0("res_", opt$resolution)
      output_dirs$plots <<- file.path(output_base, "plots", res_folder)
      output_dirs$tables <<- file.path(output_base, "tables", res_folder)
      
      # Create the new resolution-specific directories
      dir.create(output_dirs$plots, recursive = TRUE, showWarnings = FALSE)
      dir.create(output_dirs$tables, recursive = TRUE, showWarnings = FALSE)
      # -------------------------------------------------------------------------

      generate_plots(TN_combined, chosen_res = opt$resolution)
      run_and_plot_pca(TN_combined, output_dirs$plots)
      create_summary_dot_plot(TN_combined, output_dirs$plots)
    },
    all = {
      execute_step("read_csv")
      execute_step("process")
      execute_step("integrate")
      execute_step("plot")
    },
    stop("Invalid step. Valid options: read_csv, process, integrate, plot, all")
  )
}

# Main Execution
if(is.null(opt$file)) stop("Must specify input file with -f")
if(is.null(opt$step)) stop("Must specify execution step with -s")

execute_step(opt$step)
message("Pipeline step '", opt$step, "' completed successfully")

# Close log file
message("Completed at: ", Sys.time())
# Main Execution
if(is.null(opt$file)) stop("Must specify input file with -f")
if(is.null(opt$step)) stop("Must specify execution step with -s")

execute_step(opt$step)
message("Pipeline step '", opt$step, "' completed successfully")
message("Completed at: ", Sys.time())