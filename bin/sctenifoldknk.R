#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(scTenifoldKnk)
})

# Unlike GENIE3/SCRANK/HDWGCNA, this method does not produce a per-cell-type
# edge weight for RANK_SCORE to read a target's score off of. scTenifoldKnk
# runs its own virtual-KO-and-compare internally (build the WT network, zero
# the target's outgoing edges, align the two networks, chi-squared test) and
# returns a genome-wide differentially-regulated (DR) gene table directly --
# that table *is* the result, so this script's output goes straight to MERGE.

args <- commandArgs(trailingOnly = TRUE)

seuratObj <- args[1]
target    <- args[2]
n_cores   <- as.integer(args[3])

cell_type <- sub("\\.RDS$", "", basename(seuratObj))

# Number of principal components pcNet regresses on, named here because the
# degenerate-shape guard below has to compare the gene count against it.
n_comp <- 3

sc_obj <- readRDS(seuratObj)

# gene4use is the same gene universe DOWNSAMPLE hands to every other network
# method, so the KO target is guaranteed to be in it whenever it is expressed
# at all in this cell type.
genes_4_use <- intersect(sc_obj@misc$gene4use, rownames(sc_obj))

mat <- as.matrix(sc_obj[genes_4_use, ]@assays$RNA$counts)
mat <- mat[rowSums(mat) > 0, , drop = FALSE]

# scTenifoldKnk CPM-normalises by dividing every cell by its own total count,
# and scTenifoldNet::cpmNormalization is a bare t(t(X)/colSums(X)) with no
# guard for a total of zero. A cell with no counts left becomes a column of
# NaN, and makeNetworks' per-bootstrap `Z[apply(Z, 1, sum) > 0, ]` then
# subsets with NA (NaN > 0 is NA) and dies with "missing value where
# TRUE/FALSE needed", before the first network is built. gene4use is only a
# few hundred genes, so a cell that is perfectly healthy transcriptome-wide
# can easily carry zero counts across just those. This cannot re-zero a gene:
# every gene kept above has a count in some cell, and that cell is kept here.
empty_cells <- sum(colSums(mat) == 0)
if (empty_cells > 0) {
  message("Dropping ", empty_cells, " cell(s) with no counts across gene4use in ",
          cell_type, " (CPM cannot normalise them).")
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
}

target_id <- gsub("[^A-Za-z0-9_.-]+", "_", target)
out_file  <- paste0(cell_type, "_sctenifoldknk_", target_id, ".txt")

# Same 6 columns dRegulation() returns, plus cell_type/target up front.
# Explicit, empty-but-named columns so a target that gets skipped still
# leaves a properly headered (if empty) table -- see rank_score.R for why a
# bare empty data.frame() breaks MERGE's plain `head -n 1`.
empty_dr <- data.frame(
  cell_type = character(), target = character(), gene = character(),
  distance = numeric(), Z = numeric(), FC = numeric(),
  p.value = numeric(), p.adj = numeric(), stringsAsFactors = FALSE
)

write_dr <- function(dr) {
  write.table(dr, out_file, quote = FALSE, row.names = FALSE,
              col.names = TRUE, sep = "\t")
}

# scTenifoldKnk's gKO is a SINGLE gene symbol in its default (non
# transcriptome-wide) mode -- confirmed against the package source
# (R/scTenifoldKnk.R: `length(gKO) != 1` is a hard error). It has no mode
# that knocks out several genes at once the way a ';'-joined scRank target
# (e.g. "BCL2;EIF4A1") asks for; transcriptomeWide loops each gene through
# its own separate single-gene KO rather than combining them. So a combined
# target is skipped here for now rather than silently doing something else.
if (grepl(";", target, fixed = TRUE)) {
  message("scTenifoldKnk has no simultaneous multi-gene knockout mode; ",
          "skipping combined target '", target, "' in ", cell_type, ".")
  write_dr(empty_dr)
  quit(save = "no", status = 0)
}

if (!target %in% rownames(mat)) {
  message("Target '", target, "' is not expressed in ", cell_type,
          "; skipping.")
  write_dr(empty_dr)
  quit(save = "no", status = 0)
}

# pcNet requires nc_nComp < nGenes, and scale() needs more than one cell to
# get a standard deviation from. A cell type left this thin after the filters
# above has nothing to build a network from, so skip it rather than let
# scTenifoldKnk fail partway through.
if (nrow(mat) <= n_comp || ncol(mat) < 3) {
  message("Too little data left in ", cell_type, " (", nrow(mat), " genes x ",
          ncol(mat), " cells) to build a network; skipping.")
  write_dr(empty_dr)
  quit(save = "no", status = 0)
}

# qc = FALSE: DOWNSAMPLE already curated cells and genes upstream (the same
# input every other network method trusts as-is); scQC()'s defaults
# (minLibSize = 1000 UMI, 5th-percentile gene filtering, outlier-cell
# removal) are tuned for raw, unfiltered data and would re-filter an
# already-small per-cell-type split on top of that, inconsistently with how
# every other method here treats it. nc_nCells is capped at the cell type's
# own size the same way bin/sctenifoldnet.R capped it, since nc_nCells asks
# scTenifoldNet::makeNetworks() to subsample that many cells per bootstrap
# network.
result <- tryCatch({
  scTenifoldKnk(
    countMatrix = mat,
    gKO         = target,
    qc          = FALSE,
    nc_nNet     = 10,
    nc_nCells   = min(500, ncol(mat)),
    nc_nComp    = n_comp,
    nCores      = n_cores
  )
}, error = function(e) {
  message("scTenifoldKnk failed for ", target, " in ", cell_type, ": ",
          conditionMessage(e))
  NULL
})

dr <- if (is.null(result)) NULL else result$diffRegulation

# A pair that still fails is skipped rather than aborted on, the same way the
# guards above skip one they can rule out in advance. This runs once per
# (cell type, target) pair, so failing the task would take a whole sweep down
# over a single degenerate combination; the reason is already on the log
# above, and the pair simply contributes no rows to the merged table.
if (is.null(dr) || nrow(dr) == 0) {
  message("No differentially-regulated genes returned for ", target, " in ",
          cell_type, "; writing an empty table for this pair.")
  write_dr(empty_dr)
  quit(save = "no", status = 0)
}

# dRegulation() already returns dr sorted by p.value.
dr$cell_type <- cell_type
dr$target    <- target
dr <- dr[, c("cell_type", "target", "gene", "distance", "Z", "FC",
             "p.value", "p.adj")]

write_dr(dr)
