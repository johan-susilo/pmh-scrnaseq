#!/usr/bin/env Rscript

Sys.time()

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(scCATCH)
  library(SingleR)
  library(celldex)
  library(dplyr)
  library(tidyverse)
  library(CelliD)
  library(parallel)
  library(future)
})

option_list <- list(
  make_option(c("--project_dir"), type = "character", default = NULL, help = "Path to Nextflow project"),
  make_option(c("--rds"), type = "character", default = NULL, help = "Path to integrated RDS file"),
  make_option(c("--output_dir"), type = "character", default = "./annotations", help = "Base output directory"),
  make_option(c("--tissue"), type = "character", default = "skin", help = "Tissue type for scCATCH"),
  make_option(c("--resolution"), type = "numeric", default = 0.4, help = "Clustering resolution to use")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

res_folder <- paste0("res_", opt$resolution)
output_base <- file.path(opt$output_dir, res_folder)

output_dirs <- list(
  singleR = file.path(output_base, "singleR"),
  celliD = file.path(output_base, "celliD"),
  scCATCH = file.path(output_base, "scCATCH"),
  consensus = file.path(output_base, "consensus")
)
lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

plan(sequential)

# Normalization function
normalize_cell_type <- function(cell_type) {
  normalization_map <- c(
    "Monocyte" = "Monocytes", "Endothelial_cells" = "Endothelial cells",
    "Macrophage" = "Macrophages", "DC" = "Dendritic cells",
    "Killer Cell" = "NK cells", "Natural Killer Cell" = "NK cells", "NK cell" = "NK cells",
    "T cells" = "T cells", "Fibroblasts" = "Fibroblasts", "Keratinocytes" = "Keratinocytes"
  )
  cell_type <- gsub("_", " ", cell_type)
  if (cell_type %in% names(normalization_map)) return(normalization_map[cell_type])
  return(cell_type)
}

# Load Data
message("Reading RDS file: ", opt$rds)
TN.combined <- readRDS(opt$rds)

res_col <- paste0("SCT_snn_res.", opt$resolution)
if (!res_col %in% colnames(TN.combined@meta.data)) res_col <- paste0("RNA_snn_res.", opt$resolution)
Idents(TN.combined) <- res_col

DefaultAssay(TN.combined) <- "RNA"
try({ Joined_TN.combined <- JoinLayers(TN.combined) }, silent = TRUE)

# --- 1. SingleR ---
message("\n--- Running SingleR ---")
counts <- GetAssayData(Joined_TN.combined)

hpca.se <- tryCatch(HumanPrimaryCellAtlasData(), error = function(e) celldex::HumanPrimaryCellAtlasData(ensembl = FALSE, cell.ont = "nonna"))
pred.hpca <- SingleR(test = counts, ref = hpca.se, assay.type.test=1, labels = hpca.se$label.main)
table_hpca <- table(pred.hpca@listData[["pruned.labels"]], Joined_TN.combined@active.ident)
write.csv(table_hpca, file.path(output_dirs$singleR, "SingleR_hpca.csv"))
table_hpca_df <- as.data.frame.matrix(table_hpca)
table_hpca_df["annotation",] <- rownames(table_hpca_df)[apply(table_hpca_df, 2, which.max)]
write.table(table_hpca_df, file.path(output_dirs$singleR, "SingleR_hpca_summary.tsv"), sep="\t", quote=FALSE)

bpe.se <- tryCatch(BlueprintEncodeData(), error = function(e) celldex::BlueprintEncodeData(ensembl = FALSE, cell.ont = "nonna"))
pred.bpe <- SingleR(test = counts, ref = bpe.se, assay.type.test=1, labels = bpe.se$label.main)
table_bpe <- table(pred.bpe@listData[["pruned.labels"]], Joined_TN.combined@active.ident)
write.csv(table_bpe, file.path(output_dirs$singleR, "SingleR_bpe.csv"))
table_bpe_df <- as.data.frame.matrix(table_bpe)
table_bpe_df["annotation",] <- rownames(table_bpe_df)[apply(table_bpe_df, 2, which.max)]
write.table(table_bpe_df, file.path(output_dirs$singleR, "SingleR_bpe_summary.tsv"), sep="\t", quote=FALSE)

# --- 2. CelliD ---
message("\n--- Running CelliD ---")
seurat_subset <- if(ncol(TN.combined) > 90000) subset(TN.combined, cells = sample(Cells(TN.combined), 90000)) else TN.combined
try({ seurat_subset <- JoinLayers(seurat_subset) }, silent=TRUE)


message("Extracting dense matrix to bypass Seurat v5 defunct API...")

# 1. Extract the sparse matrix
sparse_mat <- GetAssayData(seurat_subset, assay = "RNA", layer = "data")

# 2. Convert to dense matrix (CelliD natively requires this)
dense_mat <- as.matrix(sparse_mat)

# 3. Run MCA pure math function
message("Running MCA on dense matrix...")
mca_res <- RunMCA(dense_mat)

# 4. Inject results cleanly back into the Seurat object
message("Injecting MCA reduction back into Seurat object...")
seurat_subset[["mca"]] <- CreateDimReducObject(
  embeddings = mca_res$obscoord, 
  loadings = mca_res$varcoord, 
  assay = "RNA", 
  key = "MCA_"
)
Baron <- seurat_subset

# 5. Free up memory!
rm(sparse_mat, dense_mat, mca_res)
gc()


panglao <- read_tsv("https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz", show_col_types = FALSE) %>% 
  filter(str_detect(species,"Hs")) %>% group_by(`cell type`) %>% summarise(geneset = list(`official gene symbol`))
all_gs <- setNames(panglao$geneset, panglao$`cell type`)
all_gs <- all_gs[sapply(all_gs, length) >= 10]

HGT_all_gs <- RunCellHGT(Baron, pathways = all_gs, dims = 1:50)
Baron$cellid_pred <- ifelse(apply(HGT_all_gs, 2, max) > 2, rownames(HGT_all_gs)[apply(HGT_all_gs, 2, which.max)], "unassigned")
table_cellid <- table(Baron$cellid_pred, Baron@active.ident)
table_cellid_df <- as.data.frame.matrix(table_cellid)
table_cellid_df["annotation",] <- apply(table_cellid_df, 2, function(x) if(all(x==0)) "unassigned" else rownames(table_cellid_df)[which.max(x)])
write.table(table_cellid_df, file.path(output_dirs$celliD, "CelliD_summary.tsv"), sep="\t", quote=FALSE)

# --- 3. scCATCH ---
message("\n--- Running scCATCH ---")
data.input <- GetAssayData(Joined_TN.combined, assay = "RNA", layer = "data")
data.input <- rev_gene(data = data.input, data_type = "data", species = "Human", geneinfo = geneinfo)
obj <- createscCATCH(data = data.input, cluster = as.character(Idents(TN.combined)))
tissue_list <- if (opt$tissue == "skin") c('Skin','Dermis','Blood','Peripheral blood','Fibroblast','Endothelial') else c('Blood')
obj <- findmarkergene(object = obj, species = "Human", marker = cellmatch, tissue = tissue_list, use_method = "1")
obj <- findcelltype(object = obj)

catch_res <- obj@celltype
colnames(catch_res) <- tolower(colnames(catch_res))
catch_res$Cluster <- if("cluster" %in% colnames(catch_res)) catch_res$cluster else catch_res[[1]]
catch_res$Cell_Type <- if("cell_type" %in% colnames(catch_res)) catch_res$cell_type else catch_res[[2]]
catch_res$Cell_Type <- sapply(catch_res$Cell_Type, normalize_cell_type)

catch_wide <- catch_res %>% mutate(Cluster = paste0("X", Cluster)) %>% select(Cluster, Cell_Type) %>% distinct() %>% pivot_wider(names_from = Cluster, values_from = Cell_Type)
write.table(catch_wide, file.path(output_dirs$scCATCH, "scCATCH_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# --- 4. Consensus ---
message("\n--- Generating Consensus ---")
votes <- list()
add_vote <- function(cluster, type, source) {
  if(is.null(type) || is.na(type) || type == "") return()
  cl <- trimws(gsub("^X+", "", as.character(cluster)))
  t_list <- sapply(unlist(strsplit(as.character(type), ",")), function(x) normalize_cell_type(trimws(x)))
  for(t in t_list) {
    if(t == "" || tolower(t) %in% c("unassigned", "unknown")) next
    votes[[cl]] <<- rbind(votes[[cl]], data.frame(Source=source, Vote=t, stringsAsFactors=FALSE))
  }
}

for(col in colnames(table_hpca_df)) add_vote(col, table_hpca_df["annotation", col], "SingleR_HPCA")
for(col in colnames(table_bpe_df)) add_vote(col, table_bpe_df["annotation", col], "SingleR_BPE")
for(col in colnames(table_cellid_df)) add_vote(col, table_cellid_df["annotation", col], "CelliD")
for(col in colnames(catch_wide)) add_vote(col, catch_wide[[1, col]], "scCATCH")

consensus_df <- data.frame(Cluster=character(), Top_Cell_Type=character(), Count=integer(), Percent=numeric(), stringsAsFactors=FALSE)
for(cl in sort(names(votes))) {
  vc <- as.data.frame(table(votes[[cl]]$Vote))
  vc <- vc[order(-vc$Freq), ]
  consensus_df <- rbind(consensus_df, data.frame(
    Cluster = cl, Top_Cell_Type = as.character(vc$Var1[1]), 
    Count = vc$Freq[1], Percent = round(100 * vc$Freq[1] / sum(vc$Freq), 1)
  ))
}

out_path <- file.path(output_dirs$consensus, "consensus_annotation.tsv")
write.table(consensus_df, out_path, sep="\t", row.names=FALSE, quote=FALSE)
message("Saved consensus to: ", out_path)