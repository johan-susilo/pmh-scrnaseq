#!/usr/bin/env Rscript
# GO:BP + KEGG pathway enrichment on the DESeq2 output of 05_dge.R
#
# FIX vs original go.R: cell_type is now passed down from the outer folder-loop
# instead of being re-derived from the CSV filename with
# str_extract(clean_name, "^[^_]+"). That regex broke on the "ALL_COMBINED_*"
# macro-view files (e.g. "deseq2_ALL_COMBINED_fibroblast.csv" -> extracted "ALL"
# instead of "fibroblast"), which silently merged every celltype's macro-view
# GO/KEGG results into a single stray "ALL/pathways/" folder.
#
# CONTRACT WITH 05: --indir must be the SAME base directory 05_dge.R used,
# since this script scans <indir>/<celltype>/dge_pseudobulk/*.csv for its
# input and writes to <indir>/<celltype>/pathways/.
#
# Usage:
#   Rscript 06_go.R --indir results/03_subsets

suppressPackageStartupMessages({
  library(optparse)
  library(clusterProfiler)
  library(org.Hs.eg.db) # Human database for gene translation
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

source("workflow/scripts/00_utils.R")

# ==============================================================================
# 0. COMMAND-LINE INTERFACE
# ==============================================================================
option_list <- list(
  make_option(c("-i", "--indir"), type = "character", default = "results/03_subsets",
              help = "Base directory containing one folder per celltype (SAME dir as 05_dge.R --indir) [default: results/03_subsets]"),
  make_option(c("--config"), type = "character", default = "config/config.yaml", help = "Path to config.yaml"),
  make_option(c("--logfc_cutoff"), type = "numeric", default = NULL,
              help = "Absolute log2FoldChange cutoff for calling a DE gene significant [config: pathway.logfc_cutoff]"),
  make_option(c("--pval_cutoff"), type = "numeric", default = NULL,
              help = "Adjusted p-value cutoff for enrichGO/enrichKEGG and gene selection [config: pathway.pval_cutoff]"),
  make_option(c("--seed"), type = "integer", default = NULL, help = "Global random seed [config: reproducibility.random_seed]")
)
opt <- parse_args(OptionParser(option_list = option_list))
cfg <- get_config(opt$config)

`%||%` <- function(a, b) if (is.null(a)) b else a
opt$logfc_cutoff <- opt$logfc_cutoff %||% cfg_get(cfg, "pathway", "logfc_cutoff", default = 0.5)
opt$pval_cutoff  <- opt$pval_cutoff  %||% cfg_get(cfg, "pathway", "pval_cutoff",  default = 0.05)
opt$seed         <- opt$seed         %||% cfg_get(cfg, "reproducibility", "random_seed", default = 42)
set.seed(opt$seed)

if (is.null(opt$indir)) {
  stop("Missing required argument: --indir (base 03_subsets directory)")
}
if (!dir.exists(opt$indir)) {
  stop(sprintf("--indir does not exist: %s", opt$indir))
}

base_subset_dir <- opt$indir

# ==============================================================================
# 1. Define the clusterProfiler Analysis Function
# ==============================================================================
run_pathway_analysis <- function(csv_path, out_base_dir, cell_type,
                                  logfc_cutoff = 0.5, pval_cutoff = 0.05) {

  # Extract just the comparison name for titling/saving.
  # cell_type is now supplied by the caller (the outer folder loop) rather
  # than re-derived from the filename - see header note above.
  file_name <- basename(csv_path)
  clean_name <- str_remove(file_name, "^deseq2_") %>% str_remove("\\.csv$")

  message(paste("\n-> Running clusterProfiler Analysis for:", clean_name))

  out_dir <- file.path(out_base_dir, cell_type, "04_pathways")
  if (!dir.exists(out_dir)) { dir.create(out_dir, recursive = TRUE) }

  # Load DESeq2 DGE data
  dge_data <- read.csv(csv_path)

  # --- STEP 1: Gene ID Translation (Symbol to Entrez) ---
  # clusterProfiler requires Entrez IDs. We must map your gene symbols.
  dge_data <- dge_data %>% filter(!is.na(padj))

  if (nrow(dge_data) == 0) {
    message("   - SKIPPING: No genes with a valid padj in this file.")
    return(NULL)
  }

  mapped_genes <- bitr(dge_data$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

  # Merge the Entrez IDs back into your DGE results
  dge_data <- left_join(dge_data, mapped_genes, by = c("gene" = "SYMBOL")) %>%
    filter(!is.na(ENTREZID))

  # Define the custom background (Universe)
  universe_entrez <- dge_data$ENTREZID

  # --- STEP 2: Extract Significant Upregulated and Downregulated Genes ---
  genes_up <- dge_data %>%
    filter(log2FoldChange > logfc_cutoff & padj < pval_cutoff) %>%
    pull(ENTREZID)

  genes_down <- dge_data %>%
    filter(log2FoldChange < -logfc_cutoff & padj < pval_cutoff) %>%
    pull(ENTREZID)

  # Create a named list for clusterProfiler's compareCluster
  query_list <- list()
  if (length(genes_up) > 0) query_list[["Upregulated_in_PMH"]] <- genes_up
  if (length(genes_down) > 0) query_list[["Downregulated_in_PMH"]] <- genes_down

  if (length(query_list) == 0) {
    message("   - SKIPPING: No significant DE genes found to analyze.")
    return(NULL)
  }

  # ============================================================================
  # 3A. Run GO: Biological Process (GO:BP)
  # ============================================================================
  message("   - Running GO:BP...")
  go_res <- tryCatch({
    compareCluster(
      geneCluster = query_list,
      fun = "enrichGO",
      universe = universe_entrez,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = pval_cutoff,
      qvalueCutoff = pval_cutoff,
      readable = TRUE # Automatically translates Entrez back to human-readable Symbols for the plot!
    )
  }, error = function(e) {
    message("   - GO:BP failed: ", conditionMessage(e))
    NULL
  })

  if (!is.null(go_res) && nrow(as.data.frame(go_res)) > 0) {

    # --- Replace slashes with commas so gene lists survive CSV round-tripping ---
    go_df <- as.data.frame(go_res)
    go_df$geneID <- str_replace_all(go_df$geneID, "/", ", ")
    write.csv(go_df, file.path(out_dir, paste0(clean_name, "_GOBP_results.csv")), row.names = FALSE)

    # Generate native clusterProfiler DotPlot
    p_go <- dotplot(go_res, showCategory = 10) +
      ggtitle(paste("Biological Processes (GO:BP) -", tools::toTitleCase(str_replace_all(clean_name, "_", " ")))) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))

    ggsave(filename = file.path(out_dir, paste0(clean_name, "_GOBP_bubbleplot.png")), plot = p_go, width = 10, height = 7, dpi = 300)
  } else {
    message("   - No significant GO:BP terms found.")
  }

  # ============================================================================
  # 3B. Run KEGG Pathways
  # ============================================================================
  message("   - Running KEGG Pathways...")
  kegg_res <- tryCatch({
    compareCluster(
      geneCluster = query_list,
      fun = "enrichKEGG",
      universe = universe_entrez,
      organism = "hsa",
      pAdjustMethod = "BH",
      pvalueCutoff = pval_cutoff,
      qvalueCutoff = pval_cutoff
    )
  }, error = function(e) {
    message("   - KEGG failed: ", conditionMessage(e))
    NULL
  })

  if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
    # Translate KEGG Entrez IDs back to Symbols for readable CSVs and plots
    kegg_res <- setReadable(kegg_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

    # --- Replace slashes with commas so gene lists survive CSV round-tripping ---
    kegg_df <- as.data.frame(kegg_res)
    kegg_df$geneID <- str_replace_all(kegg_df$geneID, "/", ", ")
    write.csv(kegg_df, file.path(out_dir, paste0(clean_name, "_KEGG_results.csv")), row.names = FALSE)

    # Generate native clusterProfiler DotPlot
    p_kegg <- dotplot(kegg_res, showCategory = 10) +
      ggtitle(paste("KEGG Pathways -", tools::toTitleCase(str_replace_all(clean_name, "_", " ")))) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))

    ggsave(filename = file.path(out_dir, paste0(clean_name, "_KEGG_bubbleplot.png")), plot = p_kegg, width = 10, height = 7, dpi = 300)
  } else {
    message("   - No significant KEGG pathways found.")
  }

  message(paste("   - Success! Saved separate GO and KEGG plots to:", out_dir))
}

# ==============================================================================
# 2. Automated Execution Loop with Automated Logging
# ==============================================================================

cell_type_folders <- list.dirs(base_subset_dir, recursive = FALSE, full.names = FALSE)
# Exclude non-celltype scratch/marker folders (logs, dge/go .done flags, etc.)
cell_type_folders <- cell_type_folders[cell_type_folders != "logs" & !startsWith(cell_type_folders, "_")]

for (cell_type in cell_type_folders) {

  dge_dir <- file.path(base_subset_dir, cell_type, "03_dge_pseudobulk")
  if (!dir.exists(dge_dir)) { next }

  dge_files <- list.files(dge_dir, pattern = "^deseq2_.*\\.csv$", full.names = TRUE, recursive = TRUE)

  if (length(dge_files) == 0) {
    message(paste("No DGE CSV files found for", cell_type))
    next
  }

  for (csv_file in dge_files) {
    # If the analysis crashes on one file, tryCatch ensures the loop continues to the next!
    tryCatch({
      run_pathway_analysis(csv_file, base_subset_dir, cell_type,
                            logfc_cutoff = opt$logfc_cutoff,
                            pval_cutoff  = opt$pval_cutoff)
    }, error = function(e) {
      message(paste("   [ERROR] Failed to process", basename(csv_file), ":", conditionMessage(e)))
    })
  }
}

message("\n==================================================================")
message("=== All Pathway Analyses Completed Successfully! ===")
message(paste("Run finished at:", Sys.time()))
message("==================================================================")
