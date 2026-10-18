#!/usr/bin/env Rscript
suppressMessages({library(Seurat); library(Matrix)})

B <- "/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/aml_bonemarrowmap"
OUT <- "/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder"

meta <- read.csv(file.path(B, "scAML_21cohorts_Annotations_byCell.csv"),
                  check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("NA", ""))
rownames(meta) <- meta$Cell

# Per-study label mapping, decided from inspecting Diagnosis x NormalvsMalignant:
#  - BeneytoCalabuig_2023: "leukemic" -> Malignant; "healthy" (non-malignant
#    cells inside an AML sample) and "Normal" (healthy-donor sample) -> Normal;
#    "unsure"/NA dropped.
#  - Naldini_2023: genotype calls, not a Malignant/Normal vocabulary at all:
#    "MUT" (mutation detected) -> Malignant; "WT" (no mutation) -> Normal;
#    "ND"/"NoCall"/NA dropped (indeterminate genotyping).
map_label <- function(study, nvm) {
  out <- rep(NA_character_, length(nvm))
  if (study == "BeneytoCalabuig_2023") {
    out[nvm == "leukemic"] <- "Malignant"
    out[nvm %in% c("healthy", "Normal")] <- "Normal"
  } else if (study == "Naldini_2023") {
    out[nvm == "MUT"] <- "Malignant"
    out[nvm == "WT"] <- "Normal"
  }
  out
}

export_one <- function(study, rds_name) {
  obj <- readRDS(file.path(B, rds_name))
  cat(study, "input dim:", dim(obj), "\n")

  cnt <- tryCatch(GetAssayData(obj, assay = "RNA", layer = "counts"),
                   error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts"))
  if (!inherits(cnt, "sparseMatrix")) cnt <- as(cnt, "CsparseMatrix") else cnt <- as(cnt, "CsparseMatrix")

  m <- meta[colnames(cnt), , drop = FALSE]
  lab <- map_label(study, m$NormalvsMalignant)
  keep <- !is.na(lab)
  cat(study, "labeled cells kept:", sum(keep), "of", ncol(cnt), "\n")
  print(table(lab[keep]))

  cnt <- cnt[, keep, drop = FALSE]
  lab <- lab[keep]

  writeMM(cnt, file.path(OUT, paste0("train_", study, "_counts.mtx")))
  writeLines(rownames(cnt), file.path(OUT, paste0("train_", study, "_genes.tsv")))
  writeLines(colnames(cnt), file.path(OUT, paste0("train_", study, "_barcodes.tsv")))
  write.csv(data.frame(barcode = colnames(cnt), Raw_annotation = lab, study = study),
            file.path(OUT, paste0("train_", study, "_labels.csv")), row.names = FALSE)
}

export_one("Naldini_2023", "scAML_Naldini_2023_projected_cleaned.rds")
export_one("BeneytoCalabuig_2023", "scAML_BeneytoCalabuig_2023_projected_cleaned.rds")
cat("done\n")
