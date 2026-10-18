#!/usr/bin/env Rscript
suppressMessages({library(Seurat); library(Matrix)})

obj <- readRDS("/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/vangalen.seurat.rds")
cat("dim:", dim(obj), "\n")

cnt <- tryCatch(GetAssayData(obj, assay = "RNA", layer = "counts"),
                 error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts"))
cat("class before conversion:", class(cnt)[1], "\n")
if (!inherits(cnt, "sparseMatrix")) cnt <- as(cnt, "CsparseMatrix") else cnt <- as(cnt, "CsparseMatrix")
vals <- cnt@x
vals <- head(vals[vals != 0], 2e6)
cat("fraction non-integer (nonzero, sampled):", round(mean(vals %% 1 != 0), 6), "\n")
cat("min/max:", min(vals), max(vals), "\n")

out_dir <- "/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder"
writeMM(cnt, file.path(out_dir, "vangalen_counts.mtx"))
writeLines(rownames(cnt), file.path(out_dir, "vangalen_genes.tsv"))
writeLines(colnames(cnt), file.path(out_dir, "vangalen_barcodes.tsv"))

md <- obj@meta.data[, "PredictionRefined", drop = FALSE]
md$barcode <- rownames(md)
write.csv(md, file.path(out_dir, "vangalen_prediction_refined.csv"), row.names = FALSE)

cat("done\n")
