#!/usr/bin/env Rscript

# 1. Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
})

option_list <- list(
  make_option(c("--project_dir"), type="character", default=NULL, help="Path to the Nextflow project directory"),
  make_option(c("--sample_name"), type="character", default=NULL, help="Name of the sample"),
  make_option(c("--ident1"), type="character", default=NULL, help="Condition metadata"),
  make_option(c("--ident2"), type="character", default=NULL, help="Patient/Sample ID"),
  make_option(c("--data_path"), type="character", default=NULL, help="Path to 10x data"),
  make_option(c("--output_dir"), type="character", default="./output", help="Output directory")
)

# parse argument
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser) #reads what the user actually passed in

io_path <- file.path(opt$project_dir, "utils", "io.R")
source(io_path)

# check if essential arguments are missing
if (is.null(opt$sample_name) || is.null(opt$data_path)) {
  print_help(opt_parser)
  stop("Execution Halted: --sample_name and --data_path must be supplied.", call.=FALSE)
}


###########################################################

message(paste("Starting pipeline for:", opt$sample_name))

# 1. Load data AND metadata at the same time
obj <- load_10x_data(
    data_path   = opt$data_path, 
    sample_name = opt$sample_name,
    ident1      = opt$ident1, 
    ident2      = opt$ident2
)


save_seurat_object(obj, opt$output_dir, paste0("raw_load", opt$sample_name))
message(paste("Finished preprocessing for:", opt$sample_name))