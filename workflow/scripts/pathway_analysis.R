suppressPackageStartupMessages({
  library(optparse)
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
})

# ==============================================================================
# 1. COMMAND LINE ARGUMENTS
# ==============================================================================
option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Path to the DGE CSV file from Step 4"),
  make_option(c("-p", "--prefix"), type="character", help="Prefix for the plot titles/files (e.g., 'Macrophage_Chronic_vs_Acute')"),
  make_option(c("-o", "--outdir"), type="character", help="Directory to save the pathway plots"),
  make_option(c("-l", "--logfc"), type="numeric", default=0.5, help="Log2 Fold Change cutoff [default: 0.5]"),
  make_option(c("-v", "--pval"), type="numeric", default=0.05, help="Adjusted P-Value cutoff [default: 0.05]")
)
opt <- parse_args(OptionParser(option_list=option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
message(paste("\n=== Running Pathway Analysis for:", toupper(opt$prefix), "==="))

# ==============================================================================
# 2. LOAD & FILTER DGE DATA
# ==============================================================================
dge_data <- read.csv(opt$input)

# Check if required columns exist
if (!all(c("gene", "avg_log2FC", "p_val_adj") %in% colnames(dge_data))) {
  stop("ERROR: Input CSV must contain 'gene', 'avg_log2FC', and 'p_val_adj' columns.")
}

# Filter Upregulated and Downregulated genes
up_genes <- dge_data %>% filter(p_val_adj < opt$pval & avg_log2FC > opt$logfc) %>% pull(gene)
down_genes <- dge_data %>% filter(p_val_adj < opt$pval & avg_log2FC < -opt$logfc) %>% pull(gene)

message(paste("Found", length(up_genes), "upregulated genes."))
message(paste("Found", length(down_genes), "downregulated genes."))

if (length(up_genes) == 0 && length(down_genes) == 0) {
  stop("No significant genes found based on cutoffs. Exiting pathway analysis.")
}

# ==============================================================================
# 3. GENE ID CONVERSION (Symbol -> Entrez)
# ==============================================================================
message("Converting Gene Symbols to Entrez IDs...")

# Helper function to convert gene lists safely
convert_genes <- function(gene_vector) {
  if (length(gene_vector) == 0) return(character(0))
  bitr_res <- bitr(gene_vector, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  return(bitr_res$ENTREZID)
}

up_entrez <- convert_genes(up_genes)
down_entrez <- convert_genes(down_genes)

# Create a named list for compareCluster
gene_list <- list()
if (length(up_entrez) > 0) gene_list[[paste0("Upregulated\n(", length(up_entrez), ")")]] <- up_entrez
if (length(down_entrez) > 0) gene_list[[paste0("Downregulated\n(", length(down_entrez), ")")]] <- down_entrez

# ==============================================================================
# 4. RUN COMPARATIVE ENRICHMENT (GO & REACTOME)
# ==============================================================================
message("Running GO Biological Process Enrichment...")
comp_go <- tryCatch({
  compareCluster(geneCluster = gene_list, fun = "enrichGO", OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05)
}, error = function(e) { message("GO Enrichment failed or found no terms."); return(NULL) })

message("Running Reactome Pathway Enrichment...")
comp_reactome <- tryCatch({
  compareCluster(geneCluster = gene_list, fun = "enrichPathway", pvalueCutoff = 0.05)
}, error = function(e) { message("Reactome Enrichment failed or found no terms."); return(NULL) })

# ==============================================================================
# 5. GENERATE AND SAVE PLOTS
# ==============================================================================
message("Generating Plots...")

# Clean filename string
safe_prefix <- gsub(" ", "_", opt$prefix)

if (!is.null(comp_go)) {
  p_go <- dotplot(comp_go, showCategory = 15) + 
    ggtitle(paste("GO Biological Processes:", opt$prefix)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))
  
  ggsave(file.path(opt$outdir, paste0(safe_prefix, "_GO_Comparative.pdf")), plot = p_go, width = 10, height = 12)
  ggsave(file.path(opt$outdir, paste0(safe_prefix, "_GO_Comparative.png")), plot = p_go, width = 10, height = 12, dpi = 300)
}

if (!is.null(comp_reactome)) {
  p_reactome <- dotplot(comp_reactome, showCategory = 15) + 
    ggtitle(paste("Reactome Pathways:", opt$prefix)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))
  
  ggsave(file.path(opt$outdir, paste0(safe_prefix, "_Reactome_Comparative.pdf")), plot = p_reactome, width = 10, height = 12)
  ggsave(file.path(opt$outdir, paste0(safe_prefix, "_Reactome_Comparative.png")), plot = p_reactome, width = 10, height = 12, dpi = 300)
}

message("=== Pathway Analysis Complete ===")