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
# The knockout, alignment and dRegulation() call are scTenifoldKnk()'s own
# with its defaults, under the same seed, so a single-gene target gives the
# same table scTenifoldKnk() would have.

args <- commandArgs(trailingOnly = TRUE)

wt_file <- args[1]
target  <- args[2]
n_cores <- as.integer(args[3])
# --sctknk_plot: also draw scTenifoldKnk's plotKO() network for this pair.
# Anything but "true" leaves it off, so the script run by hand without the
# argument behaves as it always did.
make_plot <- length(args) >= 4 && tolower(args[4]) == "true"

cell_type <- sub("_sctknk_wt\\.rds$", "", basename(wt_file))

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

# Every skip below is a fact about this pair that a retry would not change, so
# it leaves an empty table and exits clean rather than failing the task: this
# runs once per (cell type, target), and one degenerate pair should not take a
# whole sweep down. The reason is on the log, and the pair simply contributes
# no rows to the merged table.
skip <- function(...) {
  message(...)
  write_dr(empty_dr)
  quit(save = "no", status = 0)
}

wt <- readRDS(wt_file)

if (is.null(wt)) {
  skip("No wild-type network was built for ", cell_type,
       " (see SCTENIFOLDKNK_BUILD); skipping ", target, ".")
}

wt <- as.matrix(wt)

# A ';'-joined target is a joint knockout: every one of its genes loses its
# outgoing edges in the same knocked-out network, the way scRank perturbs a
# combined target. scTenifoldKnk() itself only takes one gene, but the
# knockout is nothing more than zeroing rows, so several are zeroed at once.
# A gene with no counts in this cell type was filtered out before the network
# was built and has no edges to remove, so the rest of the target is knocked
# out without it -- "A;B" with B absent is the same knockout as "A" -- and the
# label stays as given.
genes   <- unique(trimws(strsplit(target, split = ";", fixed = TRUE)[[1]]))
genes   <- genes[nzchar(genes)]
present <- intersect(genes, rownames(wt))

if (length(present) == 0) {
  skip("No gene of '", target, "' is expressed in ", cell_type, "; skipping.")
}

if (length(present) < length(genes)) {
  message("Not expressed in ", cell_type, ", left out of the knockout: ",
          paste(setdiff(genes, present), collapse = ", "), ".")
}

ko <- wt
ko[present, ] <- 0

dr <- tryCatch({
  set.seed(1)
  ma <- manifoldAlignment(wt, ko, d = 2, nCores = n_cores)
  dRegulation(ma, empiricalNull = FALSE)
}, error = function(e) {
  message("scTenifoldKnk failed for ", target, " in ", cell_type, ": ",
          conditionMessage(e))
  NULL
})

if (is.null(dr) || nrow(dr) == 0) {
  skip("No differentially-regulated genes returned for ", target, " in ",
       cell_type, "; writing an empty table for this pair.")
}

# The knocked-out gene(s) are dropped from the table. Zeroing a gene's edges is
# what the distance measures, so it always sits at the top by construction and
# says nothing about the knockout's effect on the rest of the network. The
# other rows are unaffected: dRegulation() has already computed their p.adj,
# over every gene, before this.
dr <- dr[!dr$gene %in% present, , drop = FALSE]

if (nrow(dr) == 0) {
  skip("Only the knocked-out gene(s) came back for ", target, " in ",
       cell_type, "; writing an empty table for this pair.")
}

# dRegulation() already returns dr sorted by p.value.
dr$cell_type <- cell_type
dr$target    <- target
dr <- dr[, c("cell_type", "target", "gene", "distance", "Z", "FC",
             "p.value", "p.adj")]

write_dr(dr)

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
