#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Matrix)
})

# Which --extra_target targets the networks this run already built can take.
# --extra_target adds targets to a finished run without rebuilding anything:
# the extras never reach DOWNSAMPLE, so its gene set (gene4use) and the
# networks built on it are the ones the run already has. A gene the networks
# were not built on has nothing to score or knock out, so it is checked here,
# against the networks themselves, and only what they carry goes on to
# RANK_SCORE or SCTENIFOLDKNK_KO. The rest is listed for the report, which is
# where the run says an extra target was not run.
#
# Run once per track, since the two build different networks over different
# genes. In the rank-score track a gene counts as in a network when it has an
# edge there, either way; in the knockout track it needs an outgoing edge,
# since that is all a knockout removes (SCTENIFOLDKNK_KO skips a target with
# none). A gene in at least one identity's network is enough: the identities
# lacking it are handled the way they are for any target. A combined target
# loses only the genes no network carries, as DOWNSAMPLE's target QC does, and
# is not run when none is left.
#
# Usage: extra_target_qc.R <track> <method> <extra_target.txt> <target.txt> <network .rds ...>
#   track is rank_score or knockout; target.txt is the run's --target, whose
#   targets are already analysed and are dropped from the extras.
#
# Writes extra_targets_<track>.txt (the extras to run, one per line, possibly
# none) and extra_target_qc_<track>.tsv (one row per extra target x gene).

args      <- commandArgs(trailingOnly = TRUE)
track     <- args[1]
method    <- args[2]
extra     <- args[3]
main      <- args[4]
net_files <- args[-(1:4)]

if (!track %in% c("rank_score", "knockout")) {
  stop("track must be rank_score or knockout, got: ", track)
}

read_targets <- function(path) {
  x <- trimws(readLines(path, warn = FALSE))
  unique(x[nzchar(x)])
}
split_genes <- function(line) {
  g <- trimws(strsplit(line, ";", fixed = TRUE)[[1]])
  unique(g[nzchar(g)])
}
# A target is the set of genes knocked out together, so "A;B" and "B;A" are
# the same one.
gene_set_key <- function(line) paste(sort(split_genes(line)), collapse = ";")

extra_lines <- read_targets(extra)
main_keys   <- vapply(read_targets(main), gene_set_key, character(1))
already     <- vapply(extra_lines, gene_set_key, character(1)) %in% main_keys
if (any(already)) {
  message("Already a --target, so already analysed: ",
          paste(extra_lines[already], collapse = ", "), ".")
}
extra_lines <- extra_lines[!already]

ident_of <- function(file) sub("(_weight_.*|_sctknk_wt)\\.rds$", "", basename(file))

# Per identity, the genes with an edge (any edge in the rank-score track, an
# outgoing one in the knockout track). Only that set is kept, so one network
# is held in memory at a time. A network that failed to build -- a NULL
# SCTENIFOLDKNK_BUILD saves, or an empty matrix -- carries no gene.
with_edges <- lapply(net_files, function(file) {
  net <- tryCatch(readRDS(file), error = function(e) NULL)
  if (is.null(net) || length(dim(net)) != 2 || nrow(net) == 0 || is.null(rownames(net))) {
    message("No usable network in ", basename(file), ".")
    return(character())
  }
  net <- abs(Matrix(as.matrix(net), sparse = TRUE))
  diag(net) <- 0
  out <- rowSums(net) > 0
  if (track == "rank_score") {
    out <- out | colSums(net)[match(rownames(net), colnames(net))] > 0
  }
  rownames(net)[!is.na(out) & out]
})
names(with_edges) <- vapply(net_files, ident_of, character(1))
n_nets <- length(with_edges)

network_word <- if (track == "knockout") "knockout" else method

rows <- list()
runs <- character()
for (line in extra_lines) {
  genes  <- split_genes(line)
  n_in   <- vapply(genes, function(g) sum(vapply(with_edges, function(s) g %in% s, logical(1))),
                   integer(1))
  kept   <- genes[n_in > 0]
  as_run <- paste(kept, collapse = ";")
  # "A;B" with B in no network is the knockout of A, which may already be a
  # --target: it is then analysed there, and not run a second time.
  dup    <- length(kept) > 0 && gene_set_key(as_run) %in% main_keys
  action <- if (length(kept) == 0) {
    "not run"
  } else if (dup) {
    paste0("already analysed as --target ", as_run)
  } else if (length(kept) < length(genes)) {
    paste0("run as ", as_run)
  } else {
    "run"
  }
  if (length(kept) > 0 && !dup) runs <- c(runs, as_run)
  rows[[length(rows) + 1]] <- data.frame(
    track = track, target = line, gene = genes,
    networks = paste0(n_in, "/", n_nets),
    reason = ifelse(n_in > 0, "",
                    paste0("not in any ", network_word, " network this run built",
                           if (track == "knockout") " with an outgoing edge" else "")),
    action = action, stringsAsFactors = FALSE
  )
}

qc <- if (length(rows) > 0) do.call(rbind, rows) else
  data.frame(track = character(), target = character(), gene = character(),
             networks = character(), reason = character(), action = character())

not_run <- unique(qc$target[qc$action == "not run"])
if (length(not_run) > 0) {
  message("Extra target(s) not run in the ", track, " track, no gene in its networks: ",
          paste(not_run, collapse = ", "), ". Add them to --target and run from scratch ",
          "to build networks that carry them.")
}
message(length(unique(runs)), " extra target(s) to run in the ", track, " track.")

write.table(qc, paste0("extra_target_qc_", track, ".tsv"), quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)
# Written even when empty, so the channel downstream has a file to split into
# no targets rather than a missing output.
writeLines(unique(runs), paste0("extra_targets_", track, ".txt"))
