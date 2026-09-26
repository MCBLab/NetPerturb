#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Matrix)
  library(scTenifoldNet)
})

# The knockout track's null model, once per cell type: the same knockout,
# alignment and distance SCTENIFOLDKNK_KO runs for a target, run for random
# genes of the wild-type network instead. Every knockout of this network moves
# the same genes -- its hubs -- by an amount set by how many edges the
# knocked-out gene had, so a target's distances say little about the target
# until they are compared with what any gene's knockout does to the same
# network. sctenifoldknk_ko.R reads this file and tests each gene of a
# target's knockout against that gene's own distances here.
#
# The random genes are drawn among those with at least one outgoing edge,
# since a gene without one removes nothing and its knockout is a no-op (the
# KO step skips those), and never among the run's own targets. They are
# spread evenly over the range of outgoing edge strength -- equal numbers
# from each fifth of it -- because how far a knockout moves each gene depends
# on how much it removed, and the KO step fits that dependence across these
# knockouts before testing a target against it. A target's strength then sits
# inside the fitted range, whether it is a hub or a gene kept below
# --sctknk_min_pct with a few weak edges.
#
# Usage: sctenifoldknk_null.R <wt.rds> <n_null> <n_cores> <seed> <ndim> [targets_qc.txt]
# Writes <cell type>_sctknk_null.rds: a list of the genes, the genes knocked
# out, their out-strength, each knockout's median distance, and a genes x
# knockouts matrix of log10 distances.

args <- commandArgs(trailingOnly = TRUE)

wt_file <- args[1]
n_null  <- suppressWarnings(as.integer(args[2]))
n_cores <- as.integer(args[3])
seed    <- suppressWarnings(as.integer(if (length(args) >= 4) args[4] else NA))
if (is.na(seed)) seed <- 1L
ndim    <- suppressWarnings(as.integer(if (length(args) >= 5) args[5] else NA))
if (is.na(ndim) || ndim < 1) ndim <- 2L
targets <- if (length(args) >= 6 && file.exists(args[6])) {
  unique(trimws(unlist(strsplit(readLines(args[6], warn = FALSE), ";"))))
} else character()

cell_type <- sub("_sctknk_wt\\.rds$", "", basename(wt_file))
out_file  <- paste0(cell_type, "_sctknk_null.rds")

empty_null <- function(genes = character()) {
  list(cell_type = cell_type, ndim = ndim, genes = genes, null_genes = character(),
       out_strength = numeric(), median_distance = numeric(),
       log_distance = matrix(numeric(), nrow = length(genes), ncol = 0,
                             dimnames = list(genes, NULL)))
}

wt <- readRDS(wt_file)
if (is.null(wt) || is.na(n_null) || n_null < 1) {
  message("No null knockouts for ", cell_type,
          if (is.null(wt)) ": no wild-type network." else ": --sctknk_null is 0.")
  saveRDS(empty_null(), out_file)
  quit(save = "no", status = 0)
}

wt <- as.matrix(wt)
genes <- rownames(wt)
n <- length(genes)
out_strength <- rowSums(abs(wt))

candidates <- setdiff(genes[out_strength > 0], targets)
if (length(candidates) == 0) {
  message("No gene of ", cell_type, "'s network has an outgoing edge; no null knockouts.")
  saveRDS(empty_null(genes), out_file)
  quit(save = "no", status = 0)
}

set.seed(seed)
strata <- cut(rank(log10(out_strength[candidates]), ties.method = "first"), breaks = 5, labels = FALSE)
per_stratum <- ceiling(n_null / 5)
null_genes <- unlist(lapply(split(candidates, strata), function(g) sample(g, min(per_stratum, length(g)))))
null_genes <- sample(null_genes, min(n_null, length(null_genes)))
# The two ends of the range are always in, so no target's strength falls
# outside what the line was fitted on: the weakest genes with an edge are
# exactly the targets kept below --sctknk_min_pct.
ends <- candidates[c(which.min(out_strength[candidates]), which.max(out_strength[candidates]))]
null_genes <- unique(c(ends, null_genes))[seq_len(min(n_null, length(candidates)))]
message("Null knockouts for ", cell_type, ": ", length(null_genes), " random gene(s) with outgoing edges, ",
        "spread over out-strength ", signif(min(out_strength[null_genes]), 2), "-",
        signif(max(out_strength[null_genes]), 2), ", d = ", ndim, ".")

# The same distance sctenifoldknk_ko.R computes: each gene's Euclidean
# distance between its wild-type (X_) and knocked-out (Y_) positions in the
# aligned manifold.
knockout_distance <- function(ko_genes) {
  ko <- wt
  ko[ko_genes, ] <- 0
  set.seed(seed)
  ma <- manifoldAlignment(wt, ko, d = ndim, nCores = n_cores)
  sqrt(rowSums((ma[seq_len(n), , drop = FALSE] - ma[n + seq_len(n), , drop = FALSE])^2))
}

log_distance <- matrix(NA_real_, nrow = n, ncol = length(null_genes),
                       dimnames = list(genes, null_genes))
ok <- logical(length(null_genes))
for (j in seq_along(null_genes)) {
  d <- tryCatch(knockout_distance(null_genes[j]), error = function(e) {
    message("Null knockout of ", null_genes[j], " failed: ", conditionMessage(e))
    NULL
  })
  if (!is.null(d)) {
    log_distance[, j] <- log10(pmax(d, 1e-300))
    ok[j] <- TRUE
  }
}

null_genes   <- null_genes[ok]
log_distance <- log_distance[, ok, drop = FALSE]
message(length(null_genes), " null knockout(s) completed.")

saveRDS(list(
  cell_type       = cell_type,
  ndim            = ndim,
  genes           = genes,
  null_genes      = null_genes,
  out_strength    = unname(out_strength[null_genes]),
  median_distance = apply(10^log_distance, 2, median),
  log_distance    = log_distance
), out_file)
