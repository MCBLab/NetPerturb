#!/usr/bin/env Rscript
suppressMessages({library(Seurat); library(Matrix)})

B <- "/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/aml_bonemarrowmap"
OUT <- "/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder"

export_one <- function(study, rds_name) {
  obj <- readRDS(file.path(B, rds_name))
  cat(study, "dim:", dim(obj), "\n")

  cnt <- tryCatch(GetAssayData(obj, assay = "RNA", layer = "counts"),
                   error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts"))
  if (!inherits(cnt, "sparseMatrix")) cnt <- as(cnt, "CsparseMatrix") else cnt <- as(cnt, "CsparseMatrix")

  # Every cell, labeled or not -- this is the whole point: fill in the
  # scMalignantFinder prediction for cells whose original NormalvsMalignant
  # value is NA/ND/NoCall/unsure, alongside the ones that already had a
  # curated label (kept untouched, just carried through for reference).
  writeMM(cnt, file.path(OUT, paste0("annotate_", study, "_counts.mtx")))
  writeLines(rownames(cnt), file.path(OUT, paste0("annotate_", study, "_genes.tsv")))
  writeLines(colnames(cnt), file.path(OUT, paste0("annotate_", study, "_barcodes.tsv")))
  write.csv(data.frame(barcode = colnames(cnt), NormalvsMalignant = obj$NormalvsMalignant),
            file.path(OUT, paste0("annotate_", study, "_original_label.csv")), row.names = FALSE)
  cat(study, "exported", ncol(cnt), "cells,", nrow(cnt), "genes\n")
}

export_one("Naldini_2023", "filtered_Naldini_2023.seurat.rds")
export_one("BeneytoCalabuig_2023", "filtered_BeneytoCalabuig_2023.seurat.rds")
export_one("Lasry_2022", "filtered_Lasry_2022.seurat.rds")
cat("done\n")
