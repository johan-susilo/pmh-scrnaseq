#!/usr/bin/env Rscript
# Usage: Rscript 02c_global_go.R -i results/02_global_annotation/res_0.1/dge_pseudobulk -o results/02_global_annotation/res_0.1/pathways

suppressPackageStartupMessages({
  library(optparse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
  library(stringr)
})
source("workflow/scripts/00_utils.R")

option_list <- list(
  make_option(c("-i", "--indir"), type = "character", help = "Path to global DGE CSVs (e.g. dge_pseudobulk)"),
  make_option(c("-o", "--outdir"), type = "character", help = "Output directory for pathways"),
  make_option(c("--config"), type = "character", default = "config/config.yaml", help = "Path to config.yaml"),
  make_option(c("--logfc_cutoff"), type = "numeric", default = 0.5),
  make_option(c("--pval_cutoff"), type = "numeric", default = 0.05)
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
dge_files <- list.files(opt$indir, pattern = "^deseq2_.*\\.csv$", full.names = TRUE)

if (length(dge_files) == 0) stop("No DGE CSV files found in ", opt$indir)

for (csv_path in dge_files) {
  file_name <- basename(csv_path)
  clean_name <- str_remove(file_name, "^deseq2_") %>% str_remove("\\.csv$")
  message("\n-> Running clusterProfiler for: ", clean_name)
  
  dge_data <- read.csv(csv_path) %>% filter(!is.na(padj))
  if (nrow(dge_data) == 0) next
  
  mapped_genes <- bitr(dge_data$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  dge_data <- left_join(dge_data, mapped_genes, by = c("gene" = "SYMBOL")) %>% filter(!is.na(ENTREZID))
  universe_entrez <- dge_data$ENTREZID
  
  genes_up <- dge_data %>% filter(log2FoldChange > opt$logfc_cutoff & padj < opt$pval_cutoff) %>% pull(ENTREZID)
  genes_down <- dge_data %>% filter(log2FoldChange < -opt$logfc_cutoff & padj < opt$pval_cutoff) %>% pull(ENTREZID)
  
  query_list <- list()
  if (length(genes_up) > 0) query_list[["Upregulated"]] <- genes_up
  if (length(genes_down) > 0) query_list[["Downregulated"]] <- genes_down
  if (length(query_list) == 0) next
  
  # GO:BP
  message("   - Running GO:BP...")
  go_res <- tryCatch(compareCluster(geneCluster = query_list, fun = "enrichGO", universe = universe_entrez, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = opt$pval_cutoff, qvalueCutoff = opt$pval_cutoff, readable = TRUE), error = function(e) NULL)
  if (!is.null(go_res) && nrow(as.data.frame(go_res)) > 0) {
    go_df <- as.data.frame(go_res)
    go_df$geneID <- str_replace_all(go_df$geneID, "/", ", ")
    write.csv(go_df, file.path(opt$outdir, paste0(clean_name, "_GOBP_results.csv")), row.names = FALSE)
    p_go <- dotplot(go_res, showCategory = 10) + ggtitle(paste("GO:BP -", clean_name))
    ggsave(file.path(opt$outdir, paste0(clean_name, "_GOBP_bubbleplot.png")), plot = p_go, width = 10, height = 7)
  }
  
  # KEGG
  message("   - Running KEGG...")
  kegg_res <- tryCatch(compareCluster(geneCluster = query_list, fun = "enrichKEGG", universe = universe_entrez, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = opt$pval_cutoff, qvalueCutoff = opt$pval_cutoff), error = function(e) NULL)
  if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
    kegg_res <- setReadable(kegg_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    kegg_df <- as.data.frame(kegg_res)
    kegg_df$geneID <- str_replace_all(kegg_df$geneID, "/", ", ")
    write.csv(kegg_df, file.path(opt$outdir, paste0(clean_name, "_KEGG_results.csv")), row.names = FALSE)
    p_kegg <- dotplot(kegg_res, showCategory = 10) + ggtitle(paste("KEGG -", clean_name))
    ggsave(file.path(opt$outdir, paste0(clean_name, "_KEGG_bubbleplot.png")), plot = p_kegg, width = 10, height = 7)
  }
}
message("=== Global Pathway Analysis Complete! ===")