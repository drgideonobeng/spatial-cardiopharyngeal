#!/usr/bin/env Rscript
# =============================================================================
# Step 01: SEURAT_SPATIAL
# -----------------------------------------------------------------------------
# Converts .h5ad to Seurat, identifies Spatially Variable Features, and maps 
# cardiopharyngeal markers over spatial coordinates.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(optparse)
  library(tidyverse)
  library(fs)
  library(glue)
  library(patchwork)
})

option_list <- list(
  make_option("--in_h5ad",    type = "character"),
  make_option("--sample_id",  type = "character"),
  make_option("--timepoint",  type = "character"),
  make_option("--replicate",  type = "character"),
  make_option("--puck_id",    type = "character"),
  make_option("--genotype",   type = "character"),
  make_option("--markers",    type = "character"),
  make_option("--out_rds",    type = "character"),
  make_option("--out_pdf",    type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- 1. Convert AnnData to Seurat -------------------------------------------
message(glue("[01] Converting {opt$in_h5ad} to SeuratDisk format..."))
temp_h5seurat <- glue("{opt$sample_id}_temp.h5seurat")
Convert(opt$in_h5ad, dest = temp_h5seurat, overwrite = TRUE)

message(glue("[01] Loading Seurat object for {opt$sample_id}..."))
obj <- LoadH5Seurat(temp_h5seurat)
file_delete(temp_h5seurat)

# ---- Inject Global Metadata into the Cell Observation Table -----------------
obj$sample_id  <- opt$sample_id
obj$timepoint  <- opt$timepoint
obj$replicate  <- opt$replicate
obj$puck_id    <- opt$puck_id
obj$genotype   <- opt$genotype
obj$pipeline_version <- "v1.0-spatial-cardiopharyngeal"

# Confirm injection in logs
print(head(obj@meta.data))

# ---- 2. Spatial Visualization -----------------------------------------------
genes_to_plot <- str_split(opt$markers, ",")[[1]]
available_genes <- intersect(genes_to_plot, rownames(obj))

message(glue("[01] Plotting spatial expression for: {str_c(available_genes, collapse=', ')}"))
pdf(opt$out_pdf, width = 12, height = 5)
if(length(available_genes) > 0) {
    p <- SpatialFeaturePlot(obj, features = available_genes, pt.size.factor = 1.5, ncol = length(available_genes))
    print(p)
} else {
    message("Warning: None of the requested markers were found in the dataset.")
}
dev.off()

# ---- 3. Save Object ---------------------------------------------------------
saveRDS(obj, opt$out_rds)
message(glue(
  "\n[01] Built spatial object for '{opt$sample_id}':\n",
  "       cells              : {ncol(obj)}\n",
  "       genes              : {nrow(obj)}\n",
  "       Wrote {opt$out_rds}\n"
))