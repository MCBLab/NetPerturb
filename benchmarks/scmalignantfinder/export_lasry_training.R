#!/usr/bin/env Rscript
suppressMessages({library(Matrix)})

ext_dir <- "/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/aml_bonemarrowmap/geo_lasry/extracted"
meta_path <- "/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/aml_bonemarrowmap/scAML_21cohorts_Annotations_byCell.csv"
out_dir <- "/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder"

meta <- read.csv(meta_path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("NA", ""))
meta <- meta[meta$Study == "Lasry_2022", ]
rownames(meta) <- meta$Cell
cat("total Lasry cells:", nrow(meta), "\n")

# Lasry's NormalvsMalignant vocabulary is its own (neither Naldini's genotype
# calls nor Beneyto's healthy/leukemic/Normal): "malignant" (AML or AMKL
# blasts), "Control" (cells from a Healthy-diagnosis sample), and
# "microenvironment" (non-malignant cells inside an AML/AMKL sample).
map_label <- function(nvm) {
  out <- rep(NA_character_, length(nvm))
  out[nvm == "malignant"] <- "Malignant"
  out[nvm %in% c("Control", "microenvironment")] <- "Normal"
  out
}
meta$Raw_annotation <- map_label(meta$NormalvsMalignant)
cat("label counts:\n"); print(table(meta$Raw_annotation, useNA = "always"))

prefixes <- unique(sub(":.*", "", meta$Cell))
cat("samples:", length(prefixes), "\n")

read_sample <- function(p) {
  mtx <- Sys.glob(file.path(ext_dir, paste0("*_", p, "_matrix.mtx.gz")))
  bcs <- Sys.glob(file.path(ext_dir, paste0("*_", p, "_barcodes.tsv.gz")))
  fts <- Sys.glob(file.path(ext_dir, paste0("*_", p, "_features.tsv.gz")))
  stopifnot(length(mtx) == 1, length(bcs) == 1, length(fts) == 1)

  m <- as(readMM(gzfile(mtx)), "CsparseMatrix")
  bc <- sub("-1$", "", readLines(gzfile(bcs)))
  ft <- read.delim(gzfile(fts), header = FALSE, stringsAsFactors = FALSE)

  stopifnot(nrow(m) == nrow(ft), ncol(m) == length(bc))
  rownames(m) <- make.unique(ft$V2)
  colnames(m) <- paste0(p, ":", bc)
  m
}

mats <- vector("list", length(prefixes))
names(mats) <- prefixes
for (p in prefixes) {
  m <- read_sample(p)
  wanted <- intersect(colnames(m), rownames(meta)[!is.na(meta$Raw_annotation)])
  mats[[p]] <- m[, wanted, drop = FALSE]
  cat(p, ":", ncol(mats[[p]]), "labeled cells kept\n")
}

common_genes <- Reduce(intersect, lapply(mats, rownames))
cat("common genes:", length(common_genes), "\n")
counts <- do.call(cbind, lapply(mats, function(m) m[common_genes, , drop = FALSE]))
cat("final counts dim:", dim(counts), "\n")

writeMM(counts, file.path(out_dir, "train_Lasry_2022_counts.mtx"))
writeLines(rownames(counts), file.path(out_dir, "train_Lasry_2022_genes.tsv"))
writeLines(colnames(counts), file.path(out_dir, "train_Lasry_2022_barcodes.tsv"))
write.csv(
  data.frame(barcode = colnames(counts), Raw_annotation = meta[colnames(counts), "Raw_annotation"], study = "Lasry_2022"),
  file.path(out_dir, "train_Lasry_2022_labels.csv"), row.names = FALSE
)
cat("done\n")
