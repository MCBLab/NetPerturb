#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Matrix)
  library(scTenifoldNet)
  library(scTenifoldKnk)
})

# Second half of scTenifoldKnk, run once per (cell type, target): knock the
# target out of the wild-type network SCTENIFOLDKNK_BUILD saved for that cell
# type, align the two networks and score every gene on how far it moved. The
# result is scTenifoldKnk's differentially-regulated (DR) gene table.
#
# Unlike GENIE3/SCRANK/HDWGCNA, this method does not produce a per-cell-type
# edge weight for RANK_SCORE to read a target's score off of; the DR table
# *is* the result, so this script's output goes straight to MERGE.
#
# The knockout and alignment are scTenifoldKnk()'s own with its defaults,
# seeded the same way (exactly so with --seed 1, the seed scTenifoldKnk()
# hardcodes and the pipeline default). The differential regulation is
# scTenifoldKnk 1.0.3's, computed here rather than by the installed package;
# see d_regulation() below for why. On top of it, each gene is tested against
# the null model SCTENIFOLDKNK_NULL built for this cell type -- the same
# knockout run for random genes -- since scTenifoldKnk's own statistic only
# ranks a knockout's genes against each other; see null_calibration() below.

args <- commandArgs(trailingOnly = TRUE)

wt_file <- args[1]
target  <- args[2]
n_cores <- as.integer(args[3])
# --sctknk_plot: also draw scTenifoldKnk's plotKO() network for this pair.
# Anything but "true" leaves it off, so the script run by hand without the
# argument behaves as it always did.
make_plot <- length(args) >= 4 && tolower(args[4]) == "true"
# --seed; missing or unparseable falls back to 1, the pipeline default
seed <- suppressWarnings(as.integer(if (length(args) >= 5) args[5] else NA))
if (is.na(seed)) seed <- 1L
# The cell type's null knockouts from SCTENIFOLDKNK_NULL. Missing, "null" or
# an empty file means none, and the table then carries only scTenifoldKnk's
# own statistic.
null_file <- if (length(args) >= 6 && nzchar(args[6]) && args[6] != "null" &&
                 file.exists(args[6]) && file.size(args[6]) > 0) args[6] else NA_character_
# --sctknk_ndim: dimensions of the aligned manifold; scTenifoldKnk()'s own
# default is 2.
ndim <- suppressWarnings(as.integer(if (length(args) >= 7) args[7] else NA))
if (is.na(ndim) || ndim < 1) ndim <- 2L

cell_type <- sub("_sctknk_wt\\.rds$", "", basename(wt_file))

target_id <- gsub("[^A-Za-z0-9_.-]+", "_", target)
out_file  <- paste0(cell_type, "_sctenifoldknk_", target_id, ".txt")

# Same 6 columns d_regulation() returns, plus cell_type/target up front.
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

# One line about this pair for the merged status table: whether it ran, what
# the knockout removed and how far it moved the network next to the null
# knockouts. Written for skipped pairs too, so the report can say why a pair
# has no rows rather than leaving it out silently. (".status.tsv" so the
# module's "*.txt" glob for the DR table does not pick it up.)
status_file <- paste0(cell_type, "_sctenifoldknk_", target_id, ".status.tsv")
write_status <- function(status, genes_knocked = NA, out_degree = NA, out_strength = NA,
                         out_strength_pct = NA, median_distance = NA, effect_pct = NA,
                         n_null = NA, n_tested = NA, hits_raw = NA, hits_null = NA) {
  write.table(data.frame(cell_type = cell_type, target = target, status = status,
                         genes_knocked = genes_knocked, out_degree = out_degree,
                         out_strength = out_strength, out_strength_pct = out_strength_pct,
                         median_distance = median_distance, effect_pct = effect_pct,
                         n_null = n_null, n_tested = n_tested, hits_raw = hits_raw,
                         hits_null = hits_null, stringsAsFactors = FALSE),
              status_file, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
}

# Every skip below is a fact about this pair that a retry would not change, so
# it leaves an empty table and exits clean rather than failing the task: this
# runs once per (cell type, target), and one degenerate pair should not take a
# whole sweep down. The reason is on the log, and the pair simply contributes
# no rows to the merged table.
skip <- function(status, ...) {
  message(...)
  write_dr(empty_dr)
  write_status(status)
  quit(save = "no", status = 0)
}

wt <- readRDS(wt_file)

if (is.null(wt)) {
  skip("no wild-type network", "No wild-type network was built for ", cell_type,
       " (see SCTENIFOLDKNK_BUILD); skipping ", target, ".")
}

wt <- as.matrix(wt)

# A ';'-joined target is a joint knockout: every one of its genes loses its
# outgoing edges in the same knocked-out network, the way scRank perturbs a
# combined target. scTenifoldKnk() itself only takes one gene, but the
# knockout is nothing more than zeroing rows, so several are zeroed at once.
# Targets are kept in the network whatever their detection rate (see
# sctenifoldknk_build.R), so a gene of the target is missing only when it has
# no counts at all in this cell type. It has no edges to remove, so the rest of
# the target is knocked
# out without it -- "A;B" with B absent is the same knockout as "A" -- and the
# label stays as given.
genes   <- unique(trimws(strsplit(target, split = ";", fixed = TRUE)[[1]]))
genes   <- genes[nzchar(genes)]
present <- intersect(genes, rownames(wt))

if (length(present) == 0) {
  skip("not expressed", "No gene of '", target, "' is expressed in ", cell_type, "; skipping.")
}

if (length(present) < length(genes)) {
  message("Not expressed in ", cell_type, ", left out of the knockout: ",
          paste(setdiff(genes, present), collapse = ", "), ".")
}

# What the knockout removes is the target's outgoing edges, so a target with
# none -- a gene the network kept below --sctknk_min_pct that got no edge, or
# one whose edges were all cut -- is a no-op: the knocked-out network is the
# wild-type one, the alignment returns floating-point noise (distances of
# 1e-16), and scTenifoldKnk's statistic, which only compares genes with each
# other, still reports "significant" genes from that noise. Skipped here,
# with the reason on the status line.
out_strength <- sum(abs(wt[present, , drop = FALSE]))
out_degree   <- sum(colSums(abs(wt[present, , drop = FALSE])) > 0)
if (out_strength == 0) {
  skip("no outgoing edges", "'", target, "' has no outgoing edge in the ", cell_type,
       " network, so knocking it out changes nothing; skipping.")
}

ko <- wt
ko[present, ] <- 0

# Differential regulation as scTenifoldKnk 1.0.3 computes it. Each gene's
# distance between its wild-type and knocked-out positions in the aligned
# manifold becomes FC = distance^2 / mean(distance^2), tested against a
# chi-square with 1 df. In 1.0.3 that mean leaves the knocked-out genes out.
# 1.1, the current CRAN release, dropped its gKO argument and takes the mean
# over every gene -- and the knocked-out gene, whose distance is enormous by
# construction, then carries almost all of it: on a 5,600-gene network one
# knocked-out gene held ~98% of the summed FC, which shrinks every other
# gene's FC ~50-fold and leaves nothing below p.adj 1. So the statistic is
# taken from 1.0.3, with the knocked-out genes excluded by name rather than
# by rank, and they are also left out of the FDR correction, since they are
# not among the genes being tested. The Box-Cox Z column is 1.0.3's as well.
d_regulation <- function(manifoldOutput, gKO) {
  genes   <- gsub("^X_", "", grep("^X_", rownames(manifoldOutput), value = TRUE))
  y_genes <- gsub("^Y_", "", grep("^Y_", rownames(manifoldOutput), value = TRUE))
  n <- length(genes)
  if (n != nrow(manifoldOutput) / 2 || !all(y_genes == genes)) {
    stop("manifold output does not pair X_ and Y_ genes in the same order")
  }

  distance <- vapply(seq_len(n), function(i) {
    as.numeric(stats::dist(rbind(manifoldOutput[i, ], manifoldOutput[i + n, ])))
  }, numeric(1))

  lambdas <- seq(-2, 2, length.out = 1000)
  lambdas <- lambdas[lambdas != 0]
  bc <- try(MASS::boxcox(distance[distance > 0] ~ 1, plot = FALSE, lambda = lambdas),
            silent = TRUE)
  nD <- if (inherits(bc, "try-error")) {
    distance
  } else {
    lambda <- bc$x[which.max(bc$y)]
    if (lambda < 0) 1 / (distance^lambda) else distance^lambda
  }

  tested <- !genes %in% gKO
  FC <- distance^2 / mean(distance[tested]^2)
  p_value <- stats::pchisq(FC, df = 1, lower.tail = FALSE)

  out <- data.frame(gene = genes, distance = distance, Z = as.numeric(scale(nD)),
                    FC = FC, p.value = p_value, p.adj = NA_real_)[tested, , drop = FALSE]
  out$p.adj <- stats::p.adjust(out$p.value, method = "fdr")
  out[order(out$p.value), , drop = FALSE]
}

dr <- tryCatch({
  set.seed(seed)
  ma <- manifoldAlignment(wt, ko, d = ndim, nCores = n_cores)
  d_regulation(ma, gKO = present)
}, error = function(e) {
  message("scTenifoldKnk failed for ", target, " in ", cell_type, ": ",
          conditionMessage(e))
  NULL
})

# d_regulation() has already left the knocked-out gene(s) out: zeroing a
# gene's edges is what the distance measures, so it would sit at the top by
# construction and say nothing about the knockout's effect on the rest of the
# network. A network holding nothing but the knocked-out gene(s) leaves no
# rows at all.
if (is.null(dr) || nrow(dr) == 0) {
  skip("alignment failed", "No differentially-regulated genes returned for ", target, " in ",
       cell_type, "; writing an empty table for this pair.")
}

# Calibration against the cell type's null knockouts. Whatever gene is
# knocked out of this network, the genes that move most are its hubs, by an
# amount set by the knocked-out gene's outgoing edges: the pattern is the
# network's, not the target's. So each gene's log10 distance, centred on the
# knockout's median so that the overall scale drops out, is compared with
# what the same gene does under the random knockouts. That is not one value
# per gene: the centred pattern itself shifts with how much a knockout
# removed (a hub's knockout and a weak gene's leave different residues), so
# per gene a line is fitted to its centred distance against log10
# out-strength across the null knockouts, and the target is compared with the
# line at its own out-strength. z_null is that residual in MAD units of the
# null residuals, p_null the two-sided normal tail and p_null_adj its BH
# adjustment over the tested genes. A gene that moves further for this target
# than the random knockouts of comparable strength move it is what the target
# did; a hub that moves for every knockout is not. The MAD rather than the SD
# because a few null knockouts can move a gene a long way, and the SD would
# let them hide a real effect. The scores are then standardised once more by
# their own median and MAD over the knockout's genes, an empirical null per
# knockout in the sense of Efron (and of scTenifoldKnk 1.1's empiricalNull):
# at the ends of the strength range the fitted lines rest on few null
# knockouts, and a weak target's scores came out spread twice as wide as
# they should be, which turned into hundreds of hits. The MAD ignores a
# minority of genes, so a knockout that really moved a few genes keeps them.
# Two knockout-level numbers go on the status line: effect_pct, the share of
# null knockouts that moved the network less than this one (median
# distance), and out_strength_pct, the same for the edges removed, so a
# target whose effect is small next to random genes' can be read as such.
null_calibration <- function(dr, null_model) {
  genes_ok <- identical(sort(null_model$genes), sort(rownames(wt)))
  if (!genes_ok || ncol(null_model$log_distance) < 8) {
    message("Null model for ", cell_type, if (!genes_ok) " is for a different gene set" else
              paste0(" has only ", ncol(null_model$log_distance), " knockout(s), fewer than 8"),
            "; no null calibration.")
    return(NULL)
  }
  L <- null_model$log_distance[dr$gene, , drop = FALSE]
  L <- sweep(L, 2, apply(L, 2, median, na.rm = TRUE))
  r <- log10(pmax(dr$distance, 1e-300))
  r <- r - median(r)
  x <- log10(pmax(null_model$out_strength, 1e-300))
  x0 <- log10(out_strength)
  # One least-squares line per gene, all genes at once: L is genes x
  # knockouts, so t(L) ~ X gives every gene's intercept and slope. The line
  # needs the null to vary in strength; if it does not, its mean is used.
  if (stats::sd(x) > 0.1) {
    X    <- cbind(1, x)
    beta <- solve(crossprod(X), crossprod(X, t(L)))
    pred <- as.numeric(cbind(1, x0) %*% beta)
    res  <- t(L) - X %*% beta
  } else {
    pred <- rowMeans(L)
    res  <- t(L) - rep(pred, each = nrow(L))
  }
  s <- pmax(apply(res, 2, function(v) stats::mad(v, na.rm = TRUE)), 1e-6)
  z <- (r - pred) / s
  z <- (z - stats::median(z)) / max(stats::mad(z), 1e-6)
  p <- 2 * stats::pnorm(-abs(z))
  list(z_null = z, p_null = p, p_null_adj = stats::p.adjust(p, method = "fdr"),
       n_null = ncol(L),
       effect_pct = 100 * mean(null_model$median_distance < median(dr$distance)),
       out_strength_pct = 100 * mean(null_model$out_strength < out_strength))
}

null_model <- if (!is.na(null_file)) readRDS(null_file) else NULL
cal <- if (!is.null(null_model)) null_calibration(dr, null_model) else NULL
if (!is.null(cal)) {
  dr$z_null     <- cal$z_null
  dr$p_null     <- cal$p_null
  dr$p_null_adj <- cal$p_null_adj
}

# d_regulation() already returns dr sorted by p.value.
dr$cell_type <- cell_type
dr$target    <- target
dr <- dr[, c("cell_type", "target", "gene", "distance", "Z", "FC", "p.value", "p.adj",
             if (!is.null(cal)) c("z_null", "p_null", "p_null_adj"))]

write_dr(dr)
write_status("ok", genes_knocked = paste(present, collapse = ";"), out_degree = out_degree,
             out_strength = signif(out_strength, 4),
             out_strength_pct = if (!is.null(cal)) round(cal$out_strength_pct, 1) else NA,
             median_distance = signif(median(dr$distance), 3),
             effect_pct = if (!is.null(cal)) round(cal$effect_pct, 1) else NA,
             n_null = if (!is.null(cal)) cal$n_null else 0L, n_tested = nrow(dr),
             hits_raw = sum(dr$p.adj < 0.05, na.rm = TRUE),
             hits_null = if (!is.null(cal)) sum(dr$p_null_adj < 0.05 & dr$z_null > 0, na.rm = TRUE) else NA)

# plotKO() draws the knocked-out gene(s) with the genes the knockout moved
# (p.adj < 0.05) and the wild-type edges among them. It only reads
# tensorNetworks$WT and diffRegulation off scTenifoldKnk()'s result, so that
# much of the result is rebuilt here. annotate = FALSE: annotation queries the
# Enrichr web service against human-only libraries, which a compute node with
# no internet cannot reach and which would mislabel a mouse run. igraph is
# attached because plotKO() calls its E() without a namespace. A plot that
# fails is logged and dropped; it never costs the table above.
if (make_plot) {
  plot_file <- paste0(cell_type, "_sctknk_plot_", target_id, ".pdf")
  suppressPackageStartupMessages(library(igraph))
  grDevices::pdf(plot_file, width = 8, height = 8)
  ok <- tryCatch({
    plotKO(list(tensorNetworks = list(WT = Matrix::Matrix(wt)),
                diffRegulation = dr),
           gKO = present, annotate = FALSE)
    TRUE
  }, error = function(e) {
    message("plotKO failed for ", target, " in ", cell_type, ": ",
            conditionMessage(e))
    FALSE
  })
  grDevices::dev.off()
  if (!ok) {
    unlink(plot_file)
  }
}
