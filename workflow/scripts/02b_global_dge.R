#!/usr/bin/env Rscript
# Usage: Rscript 02b_global_dge.R -i results/02_global_annotation/res_0.2/TN.combined_annotated.rds -o results/02_global_annotation/res_0.2

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(DESeq2)
  library(ggrepel)
  library(stringr)
})
source("workflow/scripts/00_utils.R")

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Path to TN.combined_annotated.rds"),
  make_option(c("-o", "--outdir"), type = "character", help = "Output directory for global DGE"),
  make_option(c("--config"), type = "character", default = "config/config.yaml", help = "Path to config.yaml"),
  make_option(c("--ident_1"), type = "character", default = "PMH"),
  make_option(c("--ident_2"), type = "character", default = "Healthy")
)
opt <- parse_args(OptionParser(option_list = option_list))

message("Loading annotated global object...")
TN.annotated <- readRDS(opt$input)

out_dir <- file.path(opt$outdir, "05_dge_pseudobulk")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# DESeq2 Helper Function (Identical to 05_dge.R)
run_and_save_deseq2 <- function(counts_matrix, meta, comparison_name, title_prefix) {
  ident1_count <- sum(meta$condition == opt$ident_1)
  ident2_count <- sum(meta$condition == opt$ident_2)
  
  if (ident1_count < 2 | ident2_count < 2) {
    message(sprintf(" -> SKIPPING %s: Insufficient replicates (%s: %d, %s: %d)", comparison_name, opt$ident_1, ident1_count, opt$ident_2, ident2_count))
    return(NULL)
  }
  
  message(sprintf(" -> Running %s...", comparison_name))
  dds <- DESeqDataSetFromMatrix(countData = counts_matrix, colData = meta, design = ~ condition)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", opt$ident_1, opt$ident_2))
  
  res_df <- as.data.frame(res) %>%
    filter(!is.na(padj)) %>%
    mutate(gene = rownames(.)) %>%
    mutate(significance = case_when(
        padj < 0.05 & log2FoldChange > 1 ~ paste0("Up in ", opt$ident_1),
        padj < 0.05 & log2FoldChange < -1 ~ paste0("Down in ", opt$ident_1),
        TRUE ~ "Not Significant"
    )) %>%
    dplyr::select(gene, everything()) %>%
    arrange(padj)
    
  write.csv(res_df, file.path(out_dir, paste0("deseq2_", comparison_name, ".csv")), row.names = FALSE)
  
  top_genes <- res_df %>% filter(significance != "Not Significant") %>% group_by(significance) %>% slice_head(n = 20) %>% ungroup()
  pal <- setNames(c("red", "blue", "grey80"), c(paste0("Up in ", opt$ident_1), paste0("Down in ", opt$ident_1), "Not Significant"))
  
  p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_text_repel(data = top_genes, aes(label = gene), color = "black", box.padding = 0.5, max.overlaps = Inf) +
    scale_color_manual(values = pal) + theme_minimal() +
    labs(title = paste(title_prefix, "(", opt$ident_1, "vs", opt$ident_2, ")"))
    
  ggsave(filename = file.path(out_dir, paste0("volcano_", comparison_name, ".png")), plot = p_volcano, width = 7, height = 6)
}

# 1. MACRO VIEW: All cells combined
message("\n--- Aggregating ALL cells (Macro View) ---")
pb_all <- AggregateExpression(TN.annotated, assays = "RNA", slot = "counts", group.by = c("orig.ident2", "Condition"), return.seurat = FALSE)$RNA
ident_regex <- paste0("(", opt$ident_1, "|", opt$ident_2, ")$")
meta_all <- data.frame(pseudobulk_id = colnames(pb_all)) %>%
  mutate(condition = str_extract(pseudobulk_id, ident_regex))
rownames(meta_all) <- meta_all$pseudobulk_id
meta_all$condition <- factor(meta_all$condition, levels = c(opt$ident_2, opt$ident_1))

run_and_save_deseq2(pb_all, meta_all, "GLOBAL_ALL_COMBINED", "All Skin Cells Combined")

# 2. MICRO VIEW: DGE per Global Cell Type
message("\n--- Aggregating by Global Cell Type (Micro View) ---")
pb_sub <- AggregateExpression(TN.annotated, assays = "RNA", slot = "counts", group.by = c("cell_type_full", "orig.ident2", "Condition"), return.seurat = FALSE)$RNA
clusters <- levels(as.factor(TN.annotated$cell_type_full))
all_sub_cols <- colnames(pb_sub)

for (cluster_id in clusters) {
  if (cluster_id %in% c("Unknown", "Ambiguous")) next # Skip unresolved clusters
  
  seurat_safe_id <- str_replace_all(cluster_id, "_", "-")
  cluster_cols <- all_sub_cols[grepl(paste0("^g?", seurat_safe_id, "_"), all_sub_cols)]
  if (length(cluster_cols) == 0) next
  
  counts_sub <- pb_sub[, cluster_cols, drop = FALSE]
  meta_sub <- data.frame(pseudobulk_id = cluster_cols) %>%
    mutate(condition = str_extract(pseudobulk_id, ident_regex))
  rownames(meta_sub) <- meta_sub$pseudobulk_id
  meta_sub$condition <- factor(meta_sub$condition, levels = c(opt$ident_2, opt$ident_1))
  
  run_and_save_deseq2(counts_sub, meta_sub, paste0("GLOBAL_", make.names(cluster_id)), paste("Global:", cluster_id))
}

message("=== Global DGE Complete ===")