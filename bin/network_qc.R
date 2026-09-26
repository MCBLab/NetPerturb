#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Matrix)
})

# Structure of every network the run built, for the report's Data quality
# section. Perturbation scores are only comparable across identities when the
# networks they come from are comparable: a denser network, or one whose
# strength sits in a few hubs, can shift every target's score in that identity
# at once. So this summarises each network on its own terms, and places each
# target gene inside it. It reads the networks REPORT never sees (it only
# receives tables), once for the whole run rather than once per target, which
# is what computing this in RANK_SCORE would have meant.
#
# Usage: network_qc.R <network method> <targets_qc.txt> <network .rds ...>
#   Rank-score networks are <identity>_weight_<METHOD>_<n>.rds (GENIE3, scRank,
#   hdWGCNA); knockout networks are <identity>_sctknk_wt.rds. hdWGCNA's
#   <identity>_hdwgcna_qc.tsv files, when given among them, are joined in.
#
# Writes network_qc.tsv (one row per network) and target_network_qc.tsv (one
# row per network x single target gene).

args    <- commandArgs(trailingOnly = TRUE)
method  <- args[1]
targets <- unique(unlist(strsplit(readLines(args[2], warn = FALSE), ";")))
targets <- trimws(targets[nzchar(trimws(targets))])
files   <- args[-(1:2)]

net_files <- files[grepl("\\.rds$", files, ignore.case = TRUE)]
hd_files  <- files[grepl("_hdwgcna_qc\\.tsv$", files)]

# Every metric is scale-free or compared within one network, since the three
# rank-score methods and scTenifoldKnk put weights on unrelated scales: GENIE3
# importances, regression betas, rank-normalised TOM. Strength is a gene's
# summed absolute weight, in and out, so it is defined for a dense network
# (GENIE3, where every gene touches every other and degree says nothing) as
# well as a sparse one.
gini <- function(x) {
  x <- sort(x[is.finite(x)])
  n <- length(x)
  if (n == 0 || sum(x) == 0) return(NA_real_)
  sum((2 * seq_len(n) - n - 1) * x) / (n * sum(x))
}

describe <- function(file) {
  base  <- basename(file)
  track <- if (grepl("_sctknk_wt\\.rds$", base)) "knockout" else "rank_score"
  ident <- sub("(_weight_.*|_sctknk_wt)\\.rds$", "", base)
  net   <- tryCatch(readRDS(file), error = function(e) NULL)

  empty <- is.null(net) || length(dim(net)) != 2 || nrow(net) == 0
  if (!empty) {
    net <- abs(Matrix(as.matrix(net), sparse = TRUE))
    diag(net) <- 0
    net <- drop0(net)
  }
  n     <- if (empty) 0L else nrow(net)
  edges <- if (empty) 0L else nnzero(net)

  row <- data.frame(
    track = track, method = if (track == "knockout") "sctenifoldknk" else method,
    identity = ident, n_genes = n, n_edges = edges,
    density = if (n > 1) edges / (n * (n - 1)) else NA_real_,
    mean_abs_weight = if (edges > 0) sum(net@x) / edges else NA_real_,
    isolated_frac = NA_real_, hub_share_top1pct = NA_real_, strength_gini = NA_real_,
    stringsAsFactors = FALSE
  )
  tgt <- data.frame(track = character(), method = character(), identity = character(),
                    gene = character(), in_network = logical(), degree = integer(),
                    strength = numeric(), strength_percentile = numeric(),
                    stringsAsFactors = FALSE)

  if (edges > 0) {
    strength <- rowSums(net) + colSums(net)
    degree   <- rowSums(net > 0) + colSums(net > 0)
    names(strength) <- names(degree) <- rownames(net)
    top <- max(1, ceiling(0.01 * n))
    row$isolated_frac     <- mean(strength == 0)
    row$hub_share_top1pct <- sum(sort(strength, decreasing = TRUE)[seq_len(top)]) / sum(strength)
    row$strength_gini     <- gini(strength)

    # Where each target sits among the genes of this network: the share of
    # genes weaker than it, so an isolated target is at 0. This percentile is
    # comparable across identities and methods where its raw
    # strength is not. A target absent from the network is listed as such.
    tgt <- do.call(rbind, lapply(targets, function(g) {
      inn <- g %in% rownames(net)
      data.frame(track = track, method = row$method, identity = ident, gene = g,
                 in_network = inn,
                 degree = if (inn) as.integer(degree[[g]]) else NA_integer_,
                 strength = if (inn) strength[[g]] else NA_real_,
                 strength_percentile = if (inn) 100 * mean(strength < strength[[g]]) else NA_real_,
                 stringsAsFactors = FALSE)
    }))
  } else {
    message("Network for ", ident, " (", track, ") has no edges.")
  }
  list(row = row, tgt = tgt)
}

res <- lapply(net_files, describe)
network_qc <- do.call(rbind, lapply(res, `[[`, "row"))
target_qc  <- do.call(rbind, lapply(res, `[[`, "tgt"))

# hdWGCNA's own numbers -- how many metacells each network was correlated
# over and how scale-free the fit was -- join its rank-score rows.
hd_cols <- c("n_metacells", "soft_power", "power_estimate", "sft_r2")
for (col in hd_cols) network_qc[[col]] <- NA_real_
if (length(hd_files) > 0) {
  hd <- do.call(rbind, lapply(hd_files, read.delim, stringsAsFactors = FALSE))
  m  <- match(network_qc$identity, hd$identity)
  on <- network_qc$track == "rank_score" & !is.na(m)
  for (col in hd_cols) network_qc[[col]][on] <- hd[[col]][m[on]]
}
network_qc$metacells_per_gene <- network_qc$n_metacells / network_qc$n_genes

write.table(network_qc, "network_qc.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(target_qc, "target_network_qc.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
message("Summarised ", nrow(network_qc), " network(s).")
