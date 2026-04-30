suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
})

# ==============================================================================
# 1. COMMAND LINE ARGUMENTS
# ==============================================================================
option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Path to the input RDS object (e.g., annotated fibroblast subset)"),
  make_option(c("-c", "--compare_by"), type="character", help="Metadata column to group by (e.g., 'Condition', 'SingleR_label')"),
  make_option(c("-1", "--ident_1"), type="character", help="The primary group to test (e.g., 'PMH', 'Myofibroblast_Disease_Signature')"),
  make_option(c("-2", "--ident_2"), type="character", default=NULL, help="The control group to test against (e.g., 'Healthy'. If NULL, compares against all other cells)"),
  make_option(c("-p", "--prefix"), type="character", help="Prefix for naming outputs (e.g., 'fibroblast', 'macrophage')"),
  make_option(c("-o", "--outdir"), type="character", help="Directory to save the DGE tables and Volcano plots")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (!is.null(opt$ident_2)) {
  if (opt$ident_2 == "" || toupper(opt$ident_2) == "NULL" || toupper(opt$ident_2) == "NONE") {
    opt$ident_2 <- NULL
  }
}

# Ensure output directory exists
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 2. LOAD & PREPARE OBJECT
# ==============================================================================
message(paste("\n=== Running Differential Expression for:", toupper(opt$prefix), "==="))
seurat_obj <- readRDS(opt$input)

# Check if the requested metadata column actually exists
if (!opt$compare_by %in% colnames(seurat_obj@meta.data)) {
  stop(paste("ERROR: Column", opt$compare_by, "not found in seurat metadata!"))
}

# Set the active identity dynamically based on your parameter
Idents(seurat_obj) <- opt$compare_by

message(paste("Grouping cells by:", opt$compare_by))
message(paste("Comparing:", opt$ident_1, "vs", ifelse(is.null(opt$ident_2), "All Other Cells", opt$ident_2)))

# ==============================================================================
# 3. RUN DIFFERENTIAL EXPRESSION
# ==============================================================================
# Ensure we use RNA assay for DGE
DefaultAssay(seurat_obj) <- "RNA"
try({ seurat_obj <- JoinLayers(seurat_obj) }, silent = TRUE)

# Run FindMarkers
markers <- FindMarkers(
  seurat_obj,
  ident.1 = opt$ident_1,
  ident.2 = opt$ident_2,
  assay = "RNA",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

# Move gene names to a proper column
markers$gene <- rownames(markers)

# ==============================================================================
# 4. FORMAT VOLCANO PLOT DATA
# ==============================================================================
message("Formatting data for Volcano Plot...")

# Determine significance and add a labeling column
p_val_cutoff <- 0.05
logfc_cutoff <- 0.5

# Dynamic labels based on what you are comparing
up_label <- paste("Upregulated in", opt$ident_1)
down_label <- paste("Upregulated in", ifelse(is.null(opt$ident_2), "Others", opt$ident_2))

markers <- markers %>%
  mutate(
    Significance = case_when(
      p_val_adj < p_val_cutoff & avg_log2FC > logfc_cutoff ~ up_label,
      p_val_adj < p_val_cutoff & avg_log2FC < -logfc_cutoff ~ down_label,
      TRUE ~ "Not Significant"
    )
  )

# Select the top genes to label on the plot (Top 15 up, Top 15 down by p-value)
top_up <- markers %>% filter(Significance == up_label) %>% top_n(15, wt = -p_val_adj)
top_down <- markers %>% filter(Significance == down_label) %>% top_n(15, wt = -p_val_adj)
genes_to_label <- bind_rows(top_up, top_down)

# ==============================================================================
# 5. GENERATE ELITE VOLCANO PLOT
# ==============================================================================
# Set custom colors dynamically
volcano_colors <- setNames(c("#E31A1C", "#1F78B4", "#D9D9D9"), c(up_label, down_label, "Not Significant"))

# Create dynamic title
plot_title <- paste0(toupper(opt$prefix), " - ", opt$ident_1, " vs ", ifelse(is.null(opt$ident_2), "Others", opt$ident_2))

p_volcano <- ggplot(markers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = volcano_colors) +
  theme_bw(base_size = 14) +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(p_val_cutoff), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_text_repel(data = genes_to_label, aes(label = gene), size = 4, color = "black", max.overlaps = 20, box.padding = 0.5) +
  labs(
    title = plot_title,
    subtitle = paste("Grouped by:", opt$compare_by),
    x = paste("Log2 Fold Change (", opt$ident_1, "/", ifelse(is.null(opt$ident_2), "Others", opt$ident_2), ")"),
    y = "-Log10 Adjusted P-Value"
  ) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom", legend.title = element_blank())

# ==============================================================================
# 6. SAVE OUTPUTS
# ==============================================================================
# Create a safe filename without spaces or weird characters
safe_filename <- paste0(opt$prefix, "_", opt$ident_1, "_vs_", ifelse(is.null(opt$ident_2), "Others", opt$ident_2))
safe_filename <- gsub(" ", "_", safe_filename)

message(paste("Saving outputs to:", opt$outdir))

# Save CSV
write.csv(markers, file.path(opt$outdir, paste0(safe_filename, "_DGE.csv")), row.names = FALSE)

# Save Plot
ggsave(file.path(opt$outdir, paste0("Volcano_", safe_filename, ".png")), plot = p_volcano, width = 10, height = 8, dpi = 300)
ggsave(file.path(opt$outdir, paste0("Volcano_", safe_filename, ".pdf")), plot = p_volcano, width = 10, height = 8)

message("=== DGE Complete ===")