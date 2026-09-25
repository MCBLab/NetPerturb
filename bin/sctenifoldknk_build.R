#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(scTenifoldNet)
})

# First half of scTenifoldKnk, run once per cell type: the wild-type network.
# scTenifoldKnk() builds this network and then knocks out one gene on it, so
# a run with several targets used to rebuild the same network once per target
# -- the network does not depend on the target, and it is seeded, so every
# rebuild came out identical. It is built here once and saved, and
# SCTENIFOLDKNK_KO then runs only the knockout, alignment and differential
# regulation on it, once per (cell type, target), in parallel.
#
# The steps and their arguments are scTenifoldKnk()'s own defaults, in the
# same order and seeded the same way -- exactly so with --seed 1, the seed
# scTenifoldKnk() hardcodes and the pipeline default. Of its quality control,
# only the gene filter is applied: DOWNSAMPLE has already chosen the cells.
#
# The network is built on scTenifoldKnk's own gene set, not on gene4use, the
# set DOWNSAMPLE picks for the rank-score methods. scTenifoldKnk keeps every
# gene detected in more than qc_minPCT (5%) of cells; most of gene4use falls
# below that in a given cell type, and most of what clears it is not in
# gene4use, so the knockout would otherwise be modelled on sparse genes and
# miss the well-measured ones. The knockout track never feeds RANK_SCORE, so
# nothing requires it to share that universe. One departure: the targets are
# kept below the threshold, where scTenifoldKnk() would refuse to knock them
# out, so every target has a knockout wherever it is expressed at all.

args <- commandArgs(trailingOnly = TRUE)

seuratObj <- args[1]
n_cores   <- as.integer(args[2])
# --sctknk_min_pct: the fraction of the cell type's cells a gene has to be
# detected in. Missing or out of [0, 1) falls back to 0.05, scTenifoldKnk's
# own qc_minPCT.
min_pct   <- suppressWarnings(as.numeric(if (length(args) >= 3) args[3] else NA))
if (is.na(min_pct) || min_pct < 0 || min_pct >= 1) min_pct <- 0.05
# --seed; missing or unparseable falls back to 1, the pipeline default
seed      <- suppressWarnings(as.integer(if (length(args) >= 4) args[4] else NA))
if (is.na(seed)) seed <- 1L

cell_type <- sub("\\.RDS$", "", basename(seuratObj))
out_file  <- paste0(cell_type, "_sctknk_wt.rds")

# Number of principal components pcNet regresses on, named here because the
# degenerate-shape guard below has to compare the gene count against it.
n_comp <- 3

sc_obj <- readRDS(seuratObj)

# scTenifoldKnk's gene filter, scQC's X[rowMeans(X != 0) > minPCT, ], over
# every gene the object carries, plus the targets that passed DOWNSAMPLE's QC.
# The rowSums drop then takes out a target with no count in this cell type,
# which has no edges to knock out. DOWNSAMPLE counts genes_sctknk for the
# report by this same rule.
counts   <- sc_obj@assays$RNA$counts
detected <- rownames(counts)[Matrix::rowMeans(counts > 0) > min_pct]
targets  <- intersect(sc_obj@misc$targets, rownames(counts))
genes    <- union(detected, targets)

mat <- as.matrix(counts[genes, , drop = FALSE])
mat <- mat[rowSums(mat) > 0, , drop = FALSE]

below <- intersect(setdiff(targets, detected), rownames(mat))
message(nrow(mat), " genes for ", cell_type, ": detected in more than ",
        100 * min_pct, "% of its ", ncol(mat), " cells",
        if (length(below) > 0) {
          paste0(", plus target(s) kept below that: ", paste(below, collapse = ", "))
        } else "", ".")

# scTenifoldKnk CPM-normalises by dividing every cell by its own total count,
# and scTenifoldNet::cpmNormalization is a bare t(t(X)/colSums(X)) with no
# guard for a total of zero. A cell with no counts left becomes a column of
# NaN, and makeNetworks' per-bootstrap `Z[apply(Z, 1, sum) > 0, ]` then
# subsets with NA (NaN > 0 is NA) and dies with "missing value where
# TRUE/FALSE needed", before the first network is built. The genes kept above
# are a subset of the transcriptome, so a cell that is healthy overall can
# still carry zero counts across all of them, rare as that is with thousands of
# well-detected genes. This cannot re-zero a gene: every gene kept above has a
# count in some cell, and that cell is kept here.
empty_cells <- sum(colSums(mat) == 0)
if (empty_cells > 0) {
  message("Dropping ", empty_cells, " cell(s) with no counts across the network's genes in ",
          cell_type, " (CPM cannot normalise them).")
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
}

# A cell type with no network still writes its file, holding NULL, so each of
# its knockouts leaves an empty table and says why, rather than the cell type
# silently dropping out of the knockout track.
write_wt <- function(wt) {
  saveRDS(wt, out_file)
}

# pcNet requires nComp < nGenes, and scale() needs more than one cell to get a
# standard deviation from. A cell type left this thin after the filters above
# has nothing to build a network from.
if (nrow(mat) <= n_comp || ncol(mat) < 3) {
  message("Too little data left in ", cell_type, " (", nrow(mat), " genes x ",
          ncol(mat), " cells) to build a network; skipping.")
  write_wt(NULL)
  quit(save = "no", status = 0)
}

# scTenifoldKnk's strictDirection(), which it does not export: of each pair of
# opposite edges, keep the stronger. Copied rather than reached for with :::,
# so a change to the package's internals cannot break this step silently.
strict_direction <- function(X, lambda = 1) {
  S <- as.matrix(X)
  S[abs(S) < abs(t(S))] <- 0
  Matrix::Matrix(((1 - lambda) * X) + (lambda * S))
}

# A failure here is a fact about the data, not the target, so the whole cell
# type is skipped rather than the task failed -- the same call the knockout
# step used to make per pair. nc_nCells is capped at the cell type's own size,
# since makeNetworks subsamples that many cells per bootstrap network.
wt <- tryCatch({
  X <- cpmNormalization(mat)

  set.seed(seed)
  nets <- makeNetworks(X = X, q = 0.9, priorNetwork = NULL, nNet = 10,
                       nCells = min(500, ncol(X)), scaleScores = TRUE,
                       symmetric = FALSE, nComp = n_comp, nCores = n_cores)

  set.seed(seed)
  td <- tensorDecomposition(xList = nets, K = 3, maxError = 1e-05,
                            maxIter = 1000, nDecimal = 3)

  W <- as.matrix(strict_direction(td$X, lambda = 0))
  diag(W) <- 0
  t(W)
}, error = function(e) {
  message("Building the wild-type network failed for ", cell_type, ": ",
          conditionMessage(e))
  NULL
})

# Stored sparse: most entries of a thresholded network are zero, and every
# knockout task stages this file.
write_wt(if (is.null(wt)) NULL else Matrix::Matrix(wt, sparse = TRUE))
