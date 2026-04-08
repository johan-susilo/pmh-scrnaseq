#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
  library(ggrepel)
  library(patchwork)
})

option_list <- list(
  make_option(c("--project_dir"), type="character", default=NULL),
  make_option(c("--rds"), type="character", default=NULL),
  make_option(c("--consensus"), type="character", default=NULL),
  make_option(c("--output_dir"), type="character", default="./plots"),
  make_option(c("--resolution"), type="numeric", default=0.4)
)

opt <- parse_args(OptionParser(option_list=option_list))
source(file.path(opt$project_dir, "utils", "plot.R"))

out_umap <- file.path(opt$output_dir, "1_UMAPs")
out_markers <- file.path(opt$output_dir, "2_Markers")
out_props <- file.path(opt$output_dir, "3_Proportions")
lapply(c(out_umap, out_markers, out_props), dir.create, recursive=TRUE, showWarnings=FALSE)

# --- 1. Load & Apply Annotations ---
message("Loading Integrated RDS...")
TN.combined <- readRDS(opt$rds)

res_col <- paste0("SCT_snn_res.", opt$resolution)
if(!res_col %in% colnames(TN.combined@meta.data)) res_col <- paste0("RNA_snn_res.", opt$resolution)
Idents(TN.combined) <- res_col

message("Applying Consensus Annotations...")
consensus <- read.delim(opt$consensus, sep="\t")
consensus$Cluster <- as.character(consensus$Cluster)

new_ids_clean <- character(length(levels(TN.combined)))
names(new_ids_clean) <- levels(TN.combined)

for(id in levels(TN.combined)) {
  match_row <- consensus[consensus$Cluster == id, ]
  new_ids_clean[id] <- if(nrow(match_row) > 0) match_row$Top_Cell_Type[1] else "Unknown"
}

TN.annotated <- RenameIdents(TN.combined, new_ids_clean)
n_types <- length(unique(Idents(TN.annotated)))
plot_colors <- colorRampPalette(brewer.pal(min(n_types, 12), "Set3"))(n_types)

# --- 2. UMAPs ---
message("Plotting UMAPs...")
p_umap_clean <- DimPlot(TN.annotated, reduction="umap", label=TRUE, repel=TRUE, pt.size=0.8) + 
  NoLegend() + ggtitle("Annotated Cell Types")
safe_save_plot(p_umap_clean, file.path(out_umap, "UMAP_Annotated_Clean"))

p_umap_split <- DimPlot(TN.annotated, reduction="umap", split.by="orig.ident1", label=TRUE, repel=TRUE, pt.size=0.5) + 
  NoLegend() + ggtitle("Annotated Types by Condition")
safe_save_plot(p_umap_split, file.path(out_umap, "UMAP_Split_Condition"), w=20, h=8)

# --- 3. Classical Markers ---
message("Plotting Classical Markers...")
marker_sets <- list(
  "T_cells" = c("CD3D","CD8A","CD8B","CCR7","NKG7"),
  "Monocytes" = c("CD14","CD68","CD163","LYZ"),
  "Fibroblasts" = c("PDGFRA","DCN","COL1A1","LUM"),
  "Epithelial" = c("KRT14","KRT5","KRT10")
)
DefaultAssay(TN.annotated) <- "RNA"
try({ TN.annotated <- JoinLayers(TN.annotated) }, silent=TRUE)

for(ctype in names(marker_sets)) {
  feats <- intersect(marker_sets[[ctype]], rownames(TN.annotated))
  if(length(feats) > 0) {
    p <- DotPlot(TN.annotated, features=feats, cols=c("white", "darkred")) + RotatedAxis() + ggtitle(ctype)
    safe_save_plot(p, file.path(out_markers, paste0("Classical_", ctype)))
  }
}

# --- 4. Proportions ---
message("Plotting Proportions...")
prop_table <- as.data.frame(prop.table(table(Idents(TN.annotated), TN.annotated$orig.ident1), margin=2) * 100)
colnames(prop_table) <- c("CellType", "Condition", "Percentage")

p_prop <- ggplot(prop_table, aes(x=Condition, y=Percentage, fill=CellType)) +
  geom_col(color="black") + theme_bw() + scale_fill_manual(values=plot_colors) +
  ggtitle("Cell Type Proportions by Condition")
safe_save_plot(p_prop, file.path(out_props, "Proportion_Barplot"))

# --- 5. Save Final Object ---
saveRDS(TN.annotated, file.path(opt$output_dir, "TN.combined_ANNOTATED.rds"))
message("Successfully saved final Annotated RDS object!")