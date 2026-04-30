#!/usr/bin/env Rscript
# Usage: Rscript ./cell_annotation.R -r TN.combined_dim30.rds -s all
# Example: Rscript ./cell_annotation.R -r /home/johan/output/skin_pmh/TN.combined_dim30.rds -s scCATCH 
# Example: Rscript ./cell_annotation.R -r /home/johan/output/skin_pmh/TN.combined_dim30.rds -s all -o /home/johan/output/skin_pmh/annotations --consensus
# Example: Rscript ./cell_annotation.R -r TN.combined_dim30.rds -s singleR,scCATCH -c 2
# Example: Rscript /home/johan/pipeline/scRNA/skin/cell_annotation.R  -i /home/johan/output/skin_pmh_harmony_test/TN.combined_dim30.rds -s all -o /home/johan/output/skin_pmh_harmony_test/ --consensus -r 0.3
# Example: Rscript /home/johan/pipeline/scRNA/skin/cell_annotation.R -i /home/johan/output/skin_pmh_harmony_test/TN.combined_dim30.rds -s all -o /home/johan/output/skin_pmh_harmony_sctransform/annotations --consensus -r 0.3

Sys.time()

if (!requireNamespace("scCATCH", quietly = TRUE)) {
  message("Downloading scCATCH directly from CRAN...")
  install.packages("scCATCH", repos = "https://cloud.r-project.org", quiet = TRUE)
}

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
  library(parallel)
  library(future)
  library(future.apply)
})

# Command-line Interface ----------------------------------------------------
option_list <- list(
  make_option(c("-i", "--rds"), type = "character", default = NULL,
              help = "Path to RDS file (TN.combined_dim30.rds)"),
  make_option(c("-s", "--step"), type = "character", default = "all",
              help = "Pipeline step: read_rds, singleR, markers, celliD, scCATCH, consensus, annotated_plots, combined_plots, all"),
  make_option(c("-o", "--output"), type = "character", default = "annotations",
              help = "Base output directory [default: annotations]"),
  make_option(c("-p", "--plots"), type = "character", default = NULL,
              help = "Path to plots directory (for annotated_plots step)"),
  make_option(c("--consensus"), action = "store_true", default = FALSE,
              help = "Generate consensus annotations from all methods"),
  make_option(c("--tissue"), type = "character", default = "skin",
              help = "Tissue type for scCATCH [default: skin]"),
  make_option(c("-r", "--resolution"), type = "character", default = "0.2",
              help = "Clustering resolution to use for annotation [default: 0.2]")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

# --- NEW: Append the resolution to the base output path ---
res_folder <- paste0("res_", opt$res)
output_base <- file.path(opt$output, res_folder)
# ----------------------------------------------------------

dir.create(output_base, recursive = TRUE, showWarnings = FALSE)

output_dirs <- list(
  singleR = file.path(output_base, "singleR"),
  markers = file.path(output_base, "markers"),
  celliD = file.path(output_base, "celliD"),
  scCATCH = file.path(output_base, "scCATCH"),
  consensus = file.path(output_base, "consensus"),
  logs = file.path(output_base, "logs"),
  annotated_plots = if (!is.null(opt$plots)) file.path(opt$plots, "annotated") else file.path(output_base, "annotated_plots"),
  combined_plots = if (!is.null(opt$plots)) file.path(opt$plots, "combined_plots") else file.path(output_base, "combined_plots")
)
lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

# Set up logging
log_file <- file.path(output_dirs$logs, paste0("annotation_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_conn <- file(log_file, open = "wt")
sink(log_conn, type = "output", split = TRUE)  # split=TRUE means output to both console and file
sink(log_conn, type = "message")

message("Log file: ", log_file)
message("Started at: ", Sys.time())

message("Using MULTISESSION processing mode to optimize HPC resources")
plan(multisession, workers = parallel::detectCores() - 1)

# Helper function to safely save PDF plots
safe_save_pdf <- function(plot_obj, filepath, w = 15, h = 15) {
  tryCatch({
    pdf(filepath, width = w, height = h)
    on.exit(dev.off(), add = TRUE)  # Ensure device closes even if print fails
    print(plot_obj)
    message("Saved plot: ", filepath)
  }, error = function(e) {
    message("Warning: Failed to save plot ", filepath, ": ", conditionMessage(e))
    # Make sure device is closed
    if (length(dev.list()) > 0) dev.off()
  })
}

# Normalization function for cell type names (from Cell_anno_consensus.pl)
normalize_cell_type <- function(cell_type) {
  # Normalization mappings
  normalization_map <- c(
    "Monocyte" = "Monocytes",
    "Endothelial_cells" = "Endothelial cells",
    "Endothelial cells" = "Endothelial cells",
    "Macrophage" = "Macrophages",
    "DC" = "Dendritic cells",
    "Killer Cell" = "NK cells",
    "Natural Killer Cell" = "NK cells",
    "NK cell" = "NK cells",
    "T cells" = "T cells",
    "Fibroblasts" = "Fibroblasts",
    "Keratinocytes" = "Keratinocytes",
    "Epithelial cells" = "Epithelial cells"
  )

  # Replace underscores with spaces
  cell_type <- gsub("_", " ", cell_type)

  # Apply normalization mapping if exists
  if (cell_type %in% names(normalization_map)) {
    return(normalization_map[cell_type])
  }

  return(cell_type)
}

# Pipeline Functions --------------------------------------------------------
# Pipeline Functions --------------------------------------------------------
read_rds <- function(rds_path) {
  message("============================================================")
  message("Reading RDS file: ", rds_path)
  message("============================================================")
  TN.combined <- readRDS(file = rds_path)
  
  # --- NEW: Switch active identity to the requested resolution ---
  md_cols <- colnames(TN.combined@meta.data)
  res_col_sct <- paste0("SCT_snn_res.", opt$res)
  res_col_rna <- paste0("RNA_snn_res.", opt$res)
  
  if (res_col_sct %in% md_cols) {
    Idents(TN.combined) <- res_col_sct
    message("Successfully switched active identity to: ", res_col_sct)
  } else if (res_col_rna %in% md_cols) {
    Idents(TN.combined) <- res_col_rna
    message("Successfully switched active identity to: ", res_col_rna)
  } else {
    message("WARNING: Resolution column for ", opt$res, " not found in object! Using default clusters.")
  }
  # -------------------------------------------------------------
  
  DefaultAssay(TN.combined) <- "RNA"
  Joined_TN.combined <- JoinLayers(TN.combined)
  message("Done reading RDS")
  message("Cells: ", ncol(Joined_TN.combined), ", Features: ", nrow(Joined_TN.combined))
  message("Clusters: ", length(unique(Idents(Joined_TN.combined))))
  return(list(TN.combined = TN.combined, Joined_TN.combined = Joined_TN.combined))
}

#-----------------------------------------------SingleR------------------------------------------------------------------
run_singleR <- function(Joined_TN.combined) {
  message("\n============================================================")
  message("Starting SingleR Annotation")
  message("============================================================")
  # --- ADDED: Skip logic ---
  if (file.exists(file.path(output_dirs$singleR, "SingleR_hpca_summary.tsv")) &&
      file.exists(file.path(output_dirs$singleR, "SingleR_bpe_summary.tsv"))) {
    message("SingleR results already exist. Skipping calculation...")
    return(NULL)
  }
  
  counts <- GetAssayData(Joined_TN.combined)

  # hpca database
  message("Running SingleR with HumanPrimaryCellAtlas database...")
  hpca.se <- tryCatch({
    HumanPrimaryCellAtlasData()
  }, error = function(e) {
    message("Error loading HumanPrimaryCellAtlas: ", e$message)
    message("Attempting workaround...")
    # Clear cache and retry
    eh <- ExperimentHub::ExperimentHub()
    celldex::HumanPrimaryCellAtlasData(ensembl = FALSE, cell.ont = "nonna")
  })

  pred.hpca <- SingleR(test = counts, ref = hpca.se, assay.type.test=1, labels = hpca.se$label.main)
  clustering.table_hpca <- table(pred.hpca@listData[["pruned.labels"]], Joined_TN.combined@active.ident)
  write.csv(clustering.table_hpca, file = file.path(output_dirs$singleR, "SingleR_hpca.csv"), row.names = TRUE)

  # bpe database
  message("Running SingleR with BlueprintEncode database...")
  bpe.se <- tryCatch({
    BlueprintEncodeData()
  }, error = function(e) {
    message("Error loading BlueprintEncode: ", e$message)
    message("Attempting workaround...")
    # Clear cache and retry
    eh <- ExperimentHub::ExperimentHub()
    celldex::BlueprintEncodeData(ensembl = FALSE, cell.ont = "nonna")
  })

  pred.bpe <- SingleR(test = counts, ref = bpe.se, assay.type.test=1, labels = bpe.se$label.main)
  clustering.table_bpe <- table(pred.bpe@listData[["pruned.labels"]], Joined_TN.combined@active.ident)
  write.csv(clustering.table_bpe, file = file.path(output_dirs$singleR, "SingleR_bpe.csv"), row.names = TRUE)

  # Process and summarize results
  message("Processing SingleR results...")
  clustering.table_hpca <- read.csv(file.path(output_dirs$singleR, "SingleR_hpca.csv"))
  rownames(clustering.table_hpca) <- clustering.table_hpca[,1]
  clustering.table_hpca <- clustering.table_hpca[,-1]
  clustering.table_hpca["annotation",] <- rownames(clustering.table_hpca)[apply(clustering.table_hpca,2,which.max)]
  write.table(clustering.table_hpca, file = file.path(output_dirs$singleR, "SingleR_hpca_summary.tsv"),
              col.names = TRUE, sep= "\t", row.names = TRUE, quote = FALSE)

  clustering.table_bpe <- read.csv(file.path(output_dirs$singleR, "SingleR_bpe.csv"))
  rownames(clustering.table_bpe) <- clustering.table_bpe[,1]
  clustering.table_bpe <- clustering.table_bpe[,-1]
  clustering.table_bpe["annotation",] <- rownames(clustering.table_bpe)[apply(clustering.table_bpe,2,which.max)]
  write.table(clustering.table_bpe, file = file.path(output_dirs$singleR, "SingleR_bpe_summary.tsv"),
              sep= "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)

  message("SingleR annotation completed!")
  message("============================================================\n")
  return(list(hpca = clustering.table_hpca, bpe = clustering.table_bpe))
}
#-------------------------------------------SingleR end-----------------------------------------------------------

#-------------------------------------------Classical Markers------------------------------------------------------
plot_markers <- function(Joined_TN.combined) {
  message("\n============================================================")
  message("Generating Classical Marker Plots")
  message("============================================================")

  # Define marker sets for different cell types 
  marker_sets <- list(
    "Epithelial" = c("KRT1","KRT10","KRT5","KRT14","KRT6A","KRT16","KRT17","KRT18","KRT19","KRT7","DSP"),
    "Sweat_gland" = c("MUCL1","PIP","AQP5"),
    "SMC" = c("MCAM","ACTA2","MYL9","TAGLN","MYH11"),
    "Pericyte" = c("NOTCH3","RGS5","PDGFRB","MYL9","TAGLN","MYH11"),
    "Fibroblasts" = c("PDGFRA","DCN","LUM","POSTN","COL1A1","COL3A1","COL5A1","COL6A3","CD248"),
    "Vascular_EC" = c("PECAM1","VWF"),
    "Lymphatic_EC" = c("PROX1","LYVE1"),
    "T_cells" = c("GZMK","CD3D","CD8A","CD8B","CCR7","GNLY","NKG7"),
    "NK_cells" = c("GNLY","NKG7"),
    "B_cells" = c("MS4A1","CD79A","SEC11C","CD79B"),
    "Plasma_cells" = c("IGJ","MZB1","XBP1","CD79A","CD79B"),
    "Monocytes" = c("CD14","CD68","CD163","MRC1","CSF1R","IL10RA","FCGR2A","FCGR2B","CD83","LYZ"),
    "Dendritic_cells" = c("IRF7","HLA-DRA","LYZ","S100B","CD1C"),
    "Neutrophils" = c("ITGAX","ITGAM","FCGR2A","ANPEP"),
    "Mast_cells" = c("ADCYAP1","CPA3","TPSAB1","VWA5A"),
    "Melanocytes" = c("DCT","MLANA"),
    "Neuronal_cells" = c("NRXN1","SCN7A","CDH19", "S100B", "IGFBP5", "MIA", "EGFL8", "NGFR", "TYR"),
    "Schwann_cells" = c("NRXN1","CCN3","MPZ","PTN","S100B")
  )

  message("Processing ", length(marker_sets), " marker sets...")

  # Process each marker set sequentially
  results <- lapply(names(marker_sets), function(cell_type) {


    markers <- marker_sets[[cell_type]]
    plot_title <- gsub("_", " ", cell_type)

    # Filter to available features to avoid "requested variables were not found" warnings
    available_features <- rownames(Joined_TN.combined)
    markers_filtered <- intersect(markers, available_features)

    if (length(markers_filtered) == 0) {
      message("Warning: No markers available in object for ", cell_type, " - skipping plot.")
      return(NULL)
    }

    p <- DotPlot(Joined_TN.combined, features = markers_filtered, cols = c("white", "darkred"), dot.scale = 8) +
      RotatedAxis() + labs(title=plot_title) +
      theme(plot.title = element_text(hjust = 0.5, size=24))

    plot_file <- file.path(output_dirs$markers, paste0("Classical_markers_", cell_type, ".pdf"))
    # Use safe_save_pdf helper to ensure device cleanup and consistent messages
    safe_save_pdf(p, plot_file, w = 15, h = 15)

    return(cell_type)
  })

  message("Marker plots generated for ", length(Filter(Negate(is.null), results)), " cell types")
  message("============================================================\n")
}


#-------------------------------------------CelliD----------------------------------------------------------------
run_celliD <- function(seurat_object) {
  message("\n============================================================")
  message("Starting CelliD Annotation")
  message("============================================================")

  if (file.exists(file.path(output_dirs$celliD, "CelliD_PanglaoDB_summary.tsv"))) {
    message("CelliD results already exist. Skipping calculation...")
    return(NULL)
  }

  # Step 1: Downsample the object to a manageable size to prevent memory errors.
  if (ncol(seurat_object) > 90000) {
    message("Object has > 90,000 cells. Downsampling to 90,000 for CelliD analysis.")
    seurat_object_subset <- subset(seurat_object, cells = sample(Cells(seurat_object), 90000))
  } else {
    seurat_object_subset <- seurat_object
  }

  # Step 2: Join layers ON THE SMALL, SUBSETTED OBJECT.
  message("Joining layers on the downsampled object...")
  seurat_subset_joined <- JoinLayers(seurat_object_subset)

  # Step 3: CelliD dimensionality reduction through MCA
  message("Running MCA (Multiple Correspondence Analysis)...")
  Baron <- RunMCA(seurat_subset_joined, slot = "RNA")

  # Download all cell-type gene signatures from PanglaoDB
  message("Downloading PanglaoDB gene signatures...")
  panglao <- read_tsv("https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz", show_col_types = FALSE)

  # Filter to get human specific genes
  panglao_all <- panglao %>% filter(str_detect(species,"Hs"))

  # Convert dataframes to a list of named vectors which is the format for CelliD input
  panglao_all <- panglao_all %>%
    group_by(`cell type`) %>%
    summarise(geneset = list(`official gene symbol`))

  all_gs <- setNames(panglao_all$geneset, panglao_all$`cell type`)

  # Remove very short signatures
  all_gs <- all_gs[sapply(all_gs, length) >= 10]

  message("Running RunCellHGT for cell type prediction...")
  HGT_all_gs <- RunCellHGT(Baron, pathways = all_gs, dims = 1:50)
  all_gs_prediction <- rownames(HGT_all_gs)[apply(HGT_all_gs, 2, which.max)]
  Baron$all_gs_prediction_signif <- ifelse(apply(HGT_all_gs, 2, max)>2, yes = all_gs_prediction, "unassigned")

  # Safe plotting with device cleanup
  safe_save_pdf(
    DimPlot(Baron, group.by = "all_gs_prediction_signif", reduction = "umap",
            label = TRUE, label.size = 3, repel = TRUE) +
      theme(legend.text = element_text(size = 7), aspect.ratio = 1),
    file.path(output_dirs$celliD, "Baron_dimplot.pdf")
  )

  # Output data
  message("Summarizing CelliD results...")
  clustering.table_CelliD <- table(Baron@meta.data[["all_gs_prediction_signif"]], Baron@active.ident)
  write.csv(clustering.table_CelliD, file = file.path(output_dirs$celliD, "CelliD_PanglaoDB.csv"))

  table_for_summary <- clustering.table_CelliD
  if ("unassigned" %in% rownames(table_for_summary)) {
    table_for_summary <- table_for_summary[rownames(table_for_summary) != "unassigned", , drop = FALSE]
  }

  get_annotation <- function(col) {
    if (all(col == 0)) {
      return("unassigned")
    } else {
      return(rownames(table_for_summary)[which.max(col)])
    }
  }

  annotation_row <- apply(table_for_summary, 2, get_annotation)
  summary_df <- as.data.frame.matrix(clustering.table_CelliD)
  summary_df["annotation",] <- annotation_row

  write.table(summary_df, file = file.path(output_dirs$celliD, "CelliD_PanglaoDB_summary.tsv"),
              col.names = TRUE, sep= "\t", row.names =TRUE, quote = FALSE)

  message("CelliD annotation completed!")
  message("============================================================\n")
  return(Baron)
}

#-------------------------------------------scCATCH----------------------------------------------------------------
run_scCATCH <- function(TN.combined, Joined_TN.combined) {
  message("\n============================================================")
  message("Starting scCATCH Annotation")
  message("============================================================")

  # --- ADDED: Skip logic ---
  if (file.exists(file.path(output_dirs$scCATCH, "scCATCH_summary.tsv"))) {
    message("scCATCH results already exist. Skipping calculation...")
    return(NULL)
  }

  # Get normalized data matrix
  data.input <- GetAssayData(Joined_TN.combined, assay = "RNA", layer = "data")

  # Revise gene symbols
  message("Revising gene symbols...")
  data.input <- rev_gene(data = data.input, data_type = "data", species = "Human", geneinfo = geneinfo)

  # Create scCATCH object
  labels <- Idents(TN.combined)
  meta <- data.frame(group = labels, row.names = names(labels))
  obj <- createscCATCH(data = data.input, cluster = as.character(meta$group))

  # Define tissue types based on input
  tissue_list <- if (opt$tissue == "skin") {
    c('Adipose tissue','Blood','Peripheral blood','Bone','Cartilage','Subcutaneous adipose tissue',
      'Hair follicle','Lung','Muscle','Skin','Dermis','Lymph node','Lymphoid tissue',
      'Pluripotent stem cell','Skeletal muscle','Umbilical cord blood','Plasma',
      'Umbilical cord','Spleen','Serum','Bone marrow','Placenta','Embryonic stem cell','Kidney',
      'Pancreas','Pancreatic islet','Pyloric gland','Pancreatic acinar tissue')
  } else {
    c('Blood','Peripheral blood','Lymph node','Lymphoid tissue','Bone marrow','Spleen')
  }

  message("Finding marker genes with tissue filter: ", opt$tissue)
  obj <- findmarkergene(object = obj,
                        species = "Human",
                        marker = cellmatch,
                        tissue = tissue_list,
                        use_method = "1")

  # Find cell types
  message("Finding cell types...")
  obj <- findcelltype(object = obj)

  # Save results
  write.csv(obj@celltype, file = file.path(output_dirs$scCATCH, "scCATCH.csv"), row.names = FALSE)

  # Process scCATCH output for consensus
  message("Processing scCATCH results for consensus...")
  data <- read.csv(file.path(output_dirs$scCATCH, "scCATCH.csv"), header = TRUE, stringsAsFactors = FALSE)

  # --- FIXED: Dynamic column mapping to avoid length mismatch errors ---
  colnames(data) <- tolower(colnames(data))
  
  if ("cluster" %in% colnames(data)) data$Cluster <- data$cluster else data$Cluster <- data[[1]]
  
  if ("cell_type" %in% colnames(data)) {
    data$Cell_Type <- data$cell_type
  } else if ("celltype" %in% colnames(data)) {
    data$Cell_Type <- data$celltype
  } else {
    data$Cell_Type <- data[[2]] # Fallback
  }

  # Normalize cell types
  data$Cell_Type <- sapply(data$Cell_Type, normalize_cell_type)

  # Create summary with cluster prefix
  unique_mapping <- data %>%
    mutate(Cluster = ifelse(is.na(Cluster) | Cluster == "", "", paste0("X", Cluster))) %>%
    select(Cluster, Cell_Type) %>%
    distinct() %>%
    arrange(Cluster)

  wide_data <- unique_mapping %>%
    pivot_wider(names_from = Cluster, values_from = Cell_Type)

  write.table(wide_data, file = file.path(output_dirs$scCATCH, "scCATCH_summary.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  message("scCATCH annotation completed!")
  message("============================================================\n")
  return(obj)
}


#-------------------------------------------scCATCH end-----------------------------------------------------------

#-------------------------------------------Consensus Annotation---------------------------------------------------
generate_consensus_annotation <- function() {
  message("\n============================================================")
  message("Generating Consensus Annotations (Improved Voting)")
  message("============================================================")

  votes_by_cluster <- list()

  # add_vote: register one or more types (comma-separated) as separate votes
  add_vote <- function(cluster_id, cell_type_raw, source_name) {
    if (is.null(cell_type_raw) || is.na(cell_type_raw) || cell_type_raw == "") return(NULL)

    # Normalize cluster id: remove leading X or non-digit chars, keep as character
    clean_cluster <- as.character(gsub("^X+", "", as.character(cluster_id)))
    clean_cluster <- trimws(clean_cluster)

    # Split multi-type annotations and normalize each
    types <- unlist(strsplit(as.character(cell_type_raw), ","))
    types <- trimws(types)
    types <- types[types != ""]
    types <- sapply(types, normalize_cell_type, USE.NAMES = FALSE)

    for (t in types) {
      if (is.null(t) || t == "" || tolower(t) == "unassigned" || tolower(t) == "unknown") next
      vote_row <- data.frame(Source = source_name, Vote = t, stringsAsFactors = FALSE)
      if (is.null(votes_by_cluster[[clean_cluster]])) {
        votes_by_cluster[[clean_cluster]] <<- vote_row
      } else {
        votes_by_cluster[[clean_cluster]] <<- rbind(votes_by_cluster[[clean_cluster]], vote_row)
      }
    }
  }

  # robust reader for tabular summary files
  safe_read_table <- function(path) {
    if (!file.exists(path)) return(NULL)
    # try reading with row.names = 1 first (common when rownames were written)
    df <- tryCatch(read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1),
                   error = function(e) NULL)
    if (!is.null(df)) return(df)
    # fallback: read without row.names
    df2 <- tryCatch(read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE),
                    error = function(e) NULL)
    return(df2)
  }

  # process matrix-like summary (SingleR/CelliD) — supports annotation as a row or column
  process_matrix_file <- function(path, tool_name) {
    df <- safe_read_table(path)
    if (is.null(df)) return()

    # case 1: annotation as a row (rownames include 'annotation')
    if ("annotation" %in% rownames(df)) {
      annot_row <- df[rownames(df) == "annotation", , drop = FALSE]
      for (col in colnames(annot_row)) {
        add_vote(col, annot_row[[col]], tool_name)
      }
      message("Added votes from ", tool_name, " (annotation row).")
      return()
    }

    # case 2: annotation as a column
    if ("annotation" %in% colnames(df)) {
      # rownames may be cluster identifiers
      for (i in seq_len(nrow(df))) {
        cluster_id <- if (!is.null(rownames(df)) && rownames(df)[i] != "") rownames(df)[i] else df[i, 1]
        add_vote(cluster_id, df[i, "annotation"], tool_name)
      }
      message("Added votes from ", tool_name, " (annotation column).")
      return()
    }

    # case 3: sometimes annotation appears as the last row but without proper rownames
    # try to detect a row named like "annotation" in the first column
    first_col <- df[[1]]
    idx <- which(tolower(first_col) == "annotation")
    if (length(idx) == 1) {
      annot_row <- df[idx, -1, drop = FALSE]
      cols <- colnames(df)[-1]
      for (j in seq_along(cols)) {
        add_vote(cols[j], annot_row[[j]], tool_name)
      }
      message("Added votes from ", tool_name, " (annotation detected in first column).")
      return()
    }

    message("No annotation row/column found in ", path)
  }

  # Collect votes
  process_matrix_file(file.path(output_dirs$singleR, "SingleR_hpca_summary.tsv"), "SingleR_HPCA")
  process_matrix_file(file.path(output_dirs$singleR, "SingleR_bpe_summary.tsv"), "SingleR_BPE")
  process_matrix_file(file.path(output_dirs$celliD, "CelliD_PanglaoDB_summary.tsv"), "CelliD")

  # scCATCH wide summary: columns are clusters (names like X0); values are cell types (possibly multiple rows)
  sccatch_path <- file.path(output_dirs$scCATCH, "scCATCH_summary.tsv")
  if (file.exists(sccatch_path)) {
    scc_df <- safe_read_table(sccatch_path)
    if (!is.null(scc_df)) {
      # iterate columns (each column corresponds to a cluster)
      for (col in colnames(scc_df)) {
        vals <- scc_df[[col]]
        vals <- vals[!is.na(vals) & vals != ""]
        for (v in vals) {
          # scCATCH may list comma-separated types in cells; add_vote will split them
          add_vote(col, v, "scCATCH")
        }
      }
      message("Added votes from scCATCH.")
    }
  }

  # Tally votes and build consensus table
  consensus_results <- data.frame(
    Cluster = character(),
    Top_Cell_Type = character(),
    Count = integer(),
    Total_Votes = integer(),
    Percent = numeric(),
    Sources = character(),
    Vote_Details = character(),
    stringsAsFactors = FALSE
  )

  clusters <- sort(names(votes_by_cluster), decreasing = FALSE)
  for (cluster in clusters) {
    vote_df <- votes_by_cluster[[cluster]]
    if (is.null(vote_df) || nrow(vote_df) == 0) {
      consensus_results <- rbind(consensus_results, data.frame(
        Cluster = cluster,
        Top_Cell_Type = "Unknown",
        Count = 0,
        Total_Votes = 0,
        Percent = 0,
        Sources = "",
        Vote_Details = "",
        stringsAsFactors = FALSE
      ))
      next
    }

    # vote counts per cell type
    vc <- as.data.frame(table(vote_df$Vote), stringsAsFactors = FALSE)
    colnames(vc) <- c("CellType", "Votes")
    vc <- vc[order(vc$Votes, decreasing = TRUE, vc$CellType), , drop = FALSE]

    total_votes <- sum(vc$Votes)
    # Winner(s): highest vote count; include ties
    max_votes <- vc$Votes[1]
    winners <- vc$CellType[vc$Votes == max_votes]
    top_label <- paste(winners, collapse = "; ")

    percent <- round(100 * max_votes / total_votes, 1)

    # sources that voted for winner(s)
    voters_for_winner <- unique(vote_df$Source[vote_df$Vote %in% winners])
    sources_str <- paste(voters_for_winner, collapse = ", ")

    # vote details: type(count,%) ; sorted
    vc$Percent <- round(100 * vc$Votes / total_votes, 1)
    vote_details <- paste(paste0(vc$CellType, " (", vc$Votes, ", ", vc$Percent, "%)"), collapse = "; ")

    consensus_results <- rbind(consensus_results, data.frame(
      Cluster = cluster,
      Top_Cell_Type = top_label,
      Count = as.integer(max_votes),
      Total_Votes = as.integer(total_votes),
      Percent = percent,
      Sources = sources_str,
      Vote_Details = vote_details,
      stringsAsFactors = FALSE
    ))
  }

  # write neat table (Cluster as numeric if possible)
  consensus_results$Cluster <- as.character(consensus_results$Cluster)
  out_path <- file.path(output_dirs$consensus, "consensus_annotation.tsv")
  write.table(consensus_results, out_path, sep = "\t", row.names = FALSE, quote = FALSE)

  message("Consensus annotation completed! Saved to: ", out_path)
  message("Top results:")
  print(head(consensus_results))

  return(consensus_results)
}
#-------------------------------------------Consensus Annotation end-----------------------------------------------

#-------------------------------------------Generate Annotated Plots-----------------------------------------------
generate_annotated_plots <- function(TN.combined) {
  message("\n============================================================")
  message("Generating Annotated Plots (Robust)")
  message("============================================================")

  consensus_file <- file.path(output_dirs$consensus, "consensus_annotation.tsv")
  if (!file.exists(consensus_file)) stop("Consensus file missing.")
  
  consensus_data <- read.delim(consensus_file, sep = "\t", stringsAsFactors = FALSE)
  
  # Ensure Cluster column is character for matching
  consensus_data$Cluster <- as.character(consensus_data$Cluster)

  # Get current Cluster IDs from Seurat object
  current_ids <- levels(TN.combined)
  
  # --- Create Mapping Vectors ---
  # 1. Detailed: "C1_T cell"
  # 2. Clean: "T cell"
  
  new_names_detailed <- character(length(current_ids))
  new_names_clean <- character(length(current_ids))
  names(new_names_detailed) <- current_ids
  names(new_names_clean) <- current_ids
  
  for (id in current_ids) {
    # Match ID (Strip X if it exists in Seurat object just in case, though unlikely)
    clean_id <- gsub("^X", "", id) 
    
    # Look up in consensus file
    match_row <- consensus_data[consensus_data$Cluster == clean_id, ]
    
    if (nrow(match_row) > 0) {
      ctype <- match_row$Cell_Type[1]
      new_names_detailed[id] <- paste0("C", id, "_", ctype)
      new_names_clean[id] <- ctype
    } else {
      new_names_detailed[id] <- paste0("C", id, "_Unknown")
      new_names_clean[id] <- "Unknown"
    }
  }
  
  # --- Apply Renaming & Plotting ---
  
  # 1. Detailed Plot (Unique Clusters)
  TN.detailed <- RenameIdents(TN.combined, new_names_detailed)
  
  p1 <- DimPlot(TN.detailed, reduction = "umap", label = TRUE, repel = TRUE) + 
        NoLegend() + ggtitle("Annotation (Cluster + Type)")
        
  safe_save_pdf(p1, file.path(output_dirs$annotated_plots, "UMAP_Annotated_Detailed.pdf"))

  # 2. Clean Plot (Merged Cell Types)
  TN.clean <- RenameIdents(TN.combined, new_names_clean)
  
  # Generate colors based on number of CELL TYPES (not clusters)
  n_types <- length(unique(Idents(TN.clean)))
  colors_clean <- colorRampPalette(RColorBrewer::brewer.pal(min(n_types, 12), "Set3"))(n_types)

  p2 <- DimPlot(TN.clean, reduction = "umap", label = TRUE, repel = TRUE, label.size = 5) + 
        NoLegend() + ggtitle("Annotation (Cell Types Only)")
  
  p3 <- DimPlot(TN.clean, reduction = "umap", label = FALSE) + 
        scale_color_manual(values = colors_clean) + ggtitle("Annotation (Legend)")

  safe_save_pdf(p2, file.path(output_dirs$annotated_plots, "UMAP_Annotated_Clean_LabelT.pdf"))
  safe_save_pdf(p3, file.path(output_dirs$annotated_plots, "UMAP_Annotated_Clean_LabelF.pdf"))

  # Save the CLEAN object as the final annotated RDS
  saveRDS(TN.clean, file.path(output_dirs$annotated_plots, "TN.combined_annotated.rds"))
  message("Saved annotated object (Clean types) to RDS.")
  
  return(TN.clean)
}
#-------------------------------------------Generate Annotated Plots end-------------------------------------------

#-------------------------------------------Generate Combined Plots (Cluster Numbers Only)--------------------
generate_combined_plots <- function(TN.combined) {
  message("\n============================================================")
  message("Generating Combined Plots (Cluster Numbers Only)")
  message("============================================================")

  # Define colors for 18 clusters (0-17)
  cluster_colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(18)

  # UMAP plots with cluster numbers
  message("Generating UMAP plots with cluster numbers...")

  # UMAP with cluster labels
  safe_save_pdf(
    DimPlot(TN.combined, reduction = "umap", label = TRUE, label.size = 5,
            repel = TRUE, pt.size = 0.8) +
      NoLegend() +
      ggtitle("UMAP - Cluster Numbers"),
    file.path(output_dirs$combined_plots, "TNcombined_umap_clusters_labelT.pdf")
  )

  # UMAP with legend
  safe_save_pdf(
    DimPlot(TN.combined, reduction = "umap", label = FALSE, pt.size = 0.8) +
      ggtitle("UMAP - Cluster Numbers (with legend)"),
    file.path(output_dirs$combined_plots, "TNcombined_umap_clusters_labelF.pdf")
  )

  # UMAP split by condition (orig.ident1)
  safe_save_pdf(
    DimPlot(TN.combined, reduction = "umap", split.by = "orig.ident1",
            label = TRUE, label.size = 3, repel = TRUE, pt.size = 0.5, ncol = 3) +
      NoLegend() +
      ggtitle("UMAP - Cluster Numbers (split by condition)"),
    file.path(output_dirs$combined_plots, "TNcombined_umap_clusters_splitorigident1.pdf")
  )

  # UMAP split by sample (orig.ident2)
  safe_save_pdf(
    DimPlot(TN.combined, reduction = "umap", split.by = "orig.ident2",
            label = TRUE, label.size = 3, repel = TRUE, pt.size = 0.5, ncol = 3) +
      NoLegend() +
      ggtitle("UMAP - Cluster Numbers (split by sample)"),
    file.path(output_dirs$combined_plots, "TNcombined_umap_clusters_splitorigident2.pdf")
  )

  # Cluster proportion plots
  message("Generating cluster proportion plots...")

  # Cluster proportion by sample group (ident1)
  Cluster_proportion_ident1 <- table(Idents(TN.combined), TN.combined$orig.ident1)
  Cluster_proportion_ident1 <- round(sweep(Cluster_proportion_ident1, MARGIN = 2,
                                            STATS = colSums(Cluster_proportion_ident1), FUN = "/") * 100, 2)
  write.csv(Cluster_proportion_ident1,
            file = file.path(output_dirs$combined_plots, "cluster_proportion_by_group.csv"),
            row.names = TRUE)

  Cluster_proportion_ident1_df <- as.data.frame(Cluster_proportion_ident1)
  colnames(Cluster_proportion_ident1_df) <- c("Cluster", "Group", "Freq")

  cluster_proportion_ident1_plot <- ggplot(Cluster_proportion_ident1_df,
                                           aes(x = Group, y = Freq, fill = Cluster)) +
    theme_bw(base_size = 15) +
    geom_col(position = "fill", width = 0.6) +
    xlab("Sample Group") +
    ylab("Proportion") +
    labs(fill = "Cluster") +
    theme(legend.text = element_text(size = 10)) +
    scale_fill_manual(values = cluster_colors) +
    ggtitle("Cluster Proportion by Sample Group")

  safe_save_pdf(cluster_proportion_ident1_plot,
                file.path(output_dirs$combined_plots, "cluster_proportion_by_group.pdf"))

  # Cluster proportion by individual sample (ident2)
  Cluster_proportion_ident2 <- table(Idents(TN.combined), TN.combined$orig.ident2)
  Cluster_proportion_ident2 <- round(sweep(Cluster_proportion_ident2, MARGIN = 2,
                                            STATS = colSums(Cluster_proportion_ident2), FUN = "/") * 100, 2)
  write.csv(Cluster_proportion_ident2,
            file = file.path(output_dirs$combined_plots, "cluster_proportion_by_sample.csv"),
            row.names = TRUE)

  Cluster_proportion_ident2_df <- as.data.frame(Cluster_proportion_ident2)
  colnames(Cluster_proportion_ident2_df) <- c("Cluster", "Sample", "Freq")

  cluster_proportion_ident2_plot <- ggplot(Cluster_proportion_ident2_df,
                                           aes(x = Sample, y = Freq, fill = Cluster)) +
    theme_bw(base_size = 15) +
    geom_col(position = "fill", width = 0.6) +
    xlab("Sample") +
    ylab("Proportion") +
    labs(fill = "Cluster") +
    theme(legend.text = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = cluster_colors) +
    ggtitle("Cluster Proportion by Sample")

  safe_save_pdf(cluster_proportion_ident2_plot,
                file.path(output_dirs$combined_plots, "cluster_proportion_by_sample.pdf"))

  message("Combined plots (cluster numbers only) generation completed!")
  message("Plots saved to: ", output_dirs$combined_plots)
  message("============================================================\n")

  return(TN.combined)
}
#-------------------------------------------Generate Combined Plots end-------------------------------------------


# Global object for sharing between steps
seurat_objects <- NULL

# Step execution function
execute_step <- function(step) {
  switch(step,
         read_rds = {
           if (is.null(opt$rds)) stop("RDS file path must be specified with --rds")
           seurat_objects <<- read_rds(opt$rds)
           saveRDS(seurat_objects, file.path(output_base, "seurat_objects.rds"))
         },
         singleR = {
           if (is.null(seurat_objects)) {
             if (file.exists(file.path(output_base, "seurat_objects.rds"))) {
               seurat_objects <<- readRDS(file.path(output_base, "seurat_objects.rds"))
             } else {
               if (is.null(opt$rds)) stop("RDS file path must be specified with --rds")
               seurat_objects <<- read_rds(opt$rds)
               saveRDS(seurat_objects, file.path(output_base, "seurat_objects.rds"))
             }
           }
           run_singleR(seurat_objects$Joined_TN.combined)
         },
         markers = {
           if (is.null(seurat_objects)) {
             if (file.exists(file.path(output_base, "seurat_objects.rds"))) {
               seurat_objects <<- readRDS(file.path(output_base, "seurat_objects.rds"))
             } else {
               if (is.null(opt$rds)) stop("RDS file path must be specified with --rds")
               seurat_objects <<- read_rds(opt$rds)
               saveRDS(seurat_objects, file.path(output_base, "seurat_objects.rds"))
             }
           }
           plot_markers(seurat_objects$Joined_TN.combined)
         },
         celliD = {
           if (is.null(seurat_objects)) {
             if (file.exists(file.path(output_base, "seurat_objects.rds"))) {
               seurat_objects <<- readRDS(file.path(output_base, "seurat_objects.rds"))
             } else {
               if (is.null(opt$rds)) stop("RDS file path must be specified with --rds")
               seurat_objects <<- read_rds(opt$rds)
               saveRDS(seurat_objects, file.path(output_base, "seurat_objects.rds"))
             }
           }
           run_celliD(seurat_objects$TN.combined)
         },
         scCATCH = {
           if (is.null(seurat_objects)) {
             if (file.exists(file.path(output_base, "seurat_objects.rds"))) {
               seurat_objects <<- readRDS(file.path(output_base, "seurat_objects.rds"))
             } else {
               if (is.null(opt$rds)) stop("RDS file path must be specified with --rds")
               seurat_objects <<- read_rds(opt$rds)
               saveRDS(seurat_objects, file.path(output_base, "seurat_objects.rds"))
             }
           }
           run_scCATCH(seurat_objects$TN.combined, seurat_objects$Joined_TN.combined)
         },
         consensus = {
           generate_consensus_annotation()
         },
         annotated_plots = {
           if (is.null(seurat_objects)) {
             if (file.exists(file.path(output_base, "seurat_objects.rds"))) {
               seurat_objects <<- readRDS(file.path(output_base, "seurat_objects.rds"))
             } else {
               if (is.null(opt$rds)) stop("RDS file path must be specified with --rds")
               seurat_objects <<- read_rds(opt$rds)
               saveRDS(seurat_objects, file.path(output_base, "seurat_objects.rds"))
             }
           }
           generate_annotated_plots(seurat_objects$TN.combined)
         },
         combined_plots = {
           if (is.null(seurat_objects)) {
             if (file.exists(file.path(output_base, "seurat_objects.rds"))) {
               seurat_objects <<- readRDS(file.path(output_base, "seurat_objects.rds"))
             } else {
               if (is.null(opt$rds)) stop("RDS file path must be specified with --rds")
               seurat_objects <<- read_rds(opt$rds)
               saveRDS(seurat_objects, file.path(output_base, "seurat_objects.rds"))
             }
           }
           generate_combined_plots(seurat_objects$TN.combined)
         },
         all = {
           if (is.null(opt$rds)) stop("RDS file path must be specified with --rds")
           seurat_objects <<- read_rds(opt$rds)
           saveRDS(seurat_objects, file.path(output_base, "seurat_objects.rds"))
           run_singleR(seurat_objects$Joined_TN.combined)
           plot_markers(seurat_objects$Joined_TN.combined)
           run_celliD(seurat_objects$TN.combined)
           run_scCATCH(seurat_objects$TN.combined, seurat_objects$Joined_TN.combined)
           # Always generate consensus after running all annotation steps to avoid accidental omission
           message("Generating consensus annotations (forced for 'all' step)...")
           generate_consensus_annotation()
         },
         stop("Invalid step. Valid options: read_rds, singleR, markers, celliD, scCATCH, consensus, annotated_plots, combined_plots, all")
  )
}

# Main Execution ------------------------------------------------------------
execute_step(opt$step)
message("\n============================================================")
message("Cell annotation pipeline step '", opt$step, "' completed successfully")
message("============================================================")

# Close log file
message("Completed at: ", Sys.time())
sink(type = "message")
sink(type = "output")
close(log_conn)
