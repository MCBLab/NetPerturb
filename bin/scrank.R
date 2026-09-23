#!/usr/bin/env Rscript
library(Seurat)
library(dplyr)
library(scRank)

args <- commandArgs(trailingOnly = TRUE)

seuratObj <- args[1]
species <- args[2]
targets <- args[3]
column <- args[4]
n_cores <- args[5]

targets <- readLines(targets)
target <- targets[1]

n_cores <- as.integer(n_cores)

cell_type <- sub(".RDS", "", seuratObj)

sc_obj <- readRDS(seuratObj)

obj <- CreateScRank(input = sc_obj,
                    species = species,
                    cell_type = column,
                    target = target)

# CreateScRank only forces the one gene it is handed into gene4use; every other
# target reaches the network only if it happens to be an HVG, TF or drug gene,
# and RANK_SCORE then fails with "Drug target gene is not in the network". Use
# the set DOWNSAMPLE built instead -- the same one GENIE3 and hdWGCNA use, which
# already carries every requested target. Constr_net keeps a row for each of
# these genes, so a target is in the network even with no edges in this type.
gene4use <- sc_obj@misc$gene4use
if (is.null(gene4use)) {
  gene4use <- unique(c(obj@para$gene4use,
                       unlist(strsplit(targets, split = ";"))))
}
obj@para$gene4use <- gene4use[gene4use %in% rownames(sc_obj)]

obj <- Constr_net(obj, n.core = n_cores)

# obj@net is keyed by the raw, unsanitized value of the identity column
# (e.g. "8, endothelial cells"), not by `cell_type` (derived from the
# already-sanitized input filename, e.g. "8__endothelial_cells") — those
# never match whenever the raw label has a space/comma/etc., which
# silently fell through to the zero-matrix fallback below for every
# affected cell type. Each invocation processes exactly one cell type's
# split object, so obj@net normally has exactly one element regardless
# of its name; index by position instead of by name. Constr_net() can
# also return a genuinely EMPTY list (population too small to build any
# network — e.g. a 25-cell cluster), so guard the length before indexing
# rather than indexing first and checking for NULL: `[[1]]` on an empty
# list errors instead of returning NULL.
weight <- if (length(obj@net) >= 1) obj@net[[1]] else NULL

# If NULL, create gene x gene zero matrix
if (is.null(weight)) {

  genes <- obj@para$gene4use
  n <- length(genes)

  weight <- matrix(0, nrow = n, ncol = n)
  rownames(weight) <- genes
  colnames(weight) <- genes

  message(paste0("Weight was NULL for ", cell_type, " — replaced with zero matrix"))
}

n_cells <- dim(sc_obj)[2]

# Save the object
saveRDS(weight, file = paste0(cell_type, "_weight_SCRANK_", n_cells, ".rds"))
