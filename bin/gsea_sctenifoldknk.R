#!/usr/bin/env Rscript
library(fgsea)

args <- commandArgs(trailingOnly = TRUE)

dr_file  <- args[1]
gmt_file <- args[2]
min_size <- args[3]
max_size <- args[4]

out_file <- "gsea_all_targets.txt"

# fgsea() dispatches to fgseaMultilevel, an adaptive sampler: the same ranked
# list gives slightly different p-values on every call. Nextflow's cache hides
# that on -resume, which makes pinning it more important rather than less --
# without a seed, nobody re-running this by hand can reproduce the numbers in
# the published table.
set.seed(20260921)

# Same defaulting rank_score.R applies to its own numeric args: a missing or
# unparseable value falls back rather than propagating an NA into fgsea.
as_int <- function(x, fallback) {
  value <- suppressWarnings(as.integer(x))
  if (length(value) != 1 || is.na(value) || value < 1) fallback else value
}

min_size <- as_int(min_size, 10)
max_size <- as_int(max_size, 500)

# The 8 columns fgsea::fgsea() returns, with cell_type/target up front, the
# same way sctenifoldknk_ko.R fronts dRegulation()'s columns. Explicit and
# empty-but-named for the reason empty_dr gives there: a bare data.frame()
# writes no header, and every reader downstream -- MERGE-style concatenation,
# report.qmd's filters -- needs the columns to exist even when there are no
# rows.
empty_gsea <- data.frame(
  cell_type = character(), target = character(), pathway = character(),
  pval = numeric(), padj = numeric(), log2err = numeric(),
  ES = numeric(), NES = numeric(), size = integer(),
  leadingEdge = character(), stringsAsFactors = FALSE
)

write_gsea <- function(x) {
  write.table(x, out_file, quote = FALSE, row.names = FALSE,
              col.names = TRUE, sep = "\t")
}

# Nothing here is worth failing a run over: this step is downstream of every
# knockout the run did, so aborting would throw away the DR table too. Each
# guard below writes the header-only table and exits clean, which report.qmd
# renders as "nothing to show" -- identical to a real run that found nothing.
bail <- function(...) {
  message(...)
  write_gsea(empty_gsea)
  quit(save = "no", status = 0)
}

# A zero-byte file is what MERGE_SCTENIFOLDKNK's stub touches, and read.delim
# refuses it with "no lines available in input"; the header-only sentinel
# assets/NO_SCTKNK_TABLE is what a run without --sctknk hands over.
if (!file.exists(dr_file) || file.size(dr_file) == 0) {
  bail("No differentially-regulated gene table at '", dr_file, "'; nothing to enrich.")
}

dr <- read.delim(dr_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

needed <- c("cell_type", "target", "gene", "FC")
missing <- setdiff(needed, names(dr))
if (length(missing) > 0) {
  bail("Table is missing column(s): ", paste(missing, collapse = ", "), ".")
}

if (nrow(dr) == 0) {
  bail("The differentially-regulated gene table is empty; nothing to enrich.")
}

# The ranking statistic. dRegulation() reports FC as a gene's squared
# manifold-alignment distance over the mean squared distance of the run, so
# log2(FC) is signed around that mean: positive for a gene the knockout moved
# further than the average gene, negative for one it moved less. Ranking on it
# rather than on the distance is what gives fgsea a two-tailed list and lets
# NES come back negative -- a set can now be concentrated at either end.
#
# Read the sign for what it is: this is distance from the wild-type network,
# not expression, so a negative NES means "this set sat still while the rest of
# the network moved", not "this set went down". There is no direction of
# regulation anywhere in the scTenifoldKnk output to recover.
#
# FC is a ratio of squares, so it is never negative; zero distance gives
# log2(0) = -Inf, which is dropped with the other non-finite values below.
dr$log2FC <- ifelse(is.finite(dr$FC) & dr$FC > 0, log2(dr$FC), NA_real_)

pathways <- tryCatch(
  fgsea::gmtPathways(gmt_file),
  error = function(e) {
    message("Could not read the GMT at '", gmt_file, "': ", conditionMessage(e))
    NULL
  }
)

if (is.null(pathways) || length(pathways) == 0) {
  bail("No gene sets read from '", gmt_file, "'.")
}

# The commonest way this step comes back empty is a species mismatch: MSigDB
# mouse collections carry MGI symbols (Title case) and human ones uppercase
# HGNC, so a human GMT against mouse data overlaps almost nothing. Counting the
# overlap up front turns that into a line in the log instead of a table of
# zero rows with no explanation.
universe <- unique(dr$gene)
overlap  <- intersect(universe, unique(unlist(pathways, use.names = FALSE)))

message(length(pathways), " gene set(s) read; ", length(overlap), " of ",
        length(universe), " gene(s) in the table are found in them.")

if (length(overlap) < min_size) {
  message("  table symbols look like: ",
          paste(utils::head(universe, 5), collapse = ", "))
  message("  gene set symbols look like: ",
          paste(utils::head(unique(unlist(pathways, use.names = FALSE)), 5),
                collapse = ", "))

  # Almost always a species mismatch. Saying how many would match case-folded
  # names the cause precisely -- but the script does not fold: Myc -> MYC is
  # right and Trp53 -> TP53 is not, so folding half-works, which is worse than
  # not working, since the half that lands looks like a real result.
  folded <- length(intersect(toupper(universe),
                             toupper(unique(unlist(pathways, use.names = FALSE)))))
  if (folded > 0) {
    message("  ", folded, " symbol(s) would match if case were ignored. That is ",
            "not done here: mouse and human symbols differ by more than case ",
            "(Trp53 vs TP53), so folding would produce plausible but wrong hits.")
  }

  bail("Only ", length(overlap), " gene(s) overlap the gene sets, fewer than ",
       "the minimum set size of ", min_size, ". Supply a GMT for --species: ",
       "MSigDB publishes .Mm. and .Hs. symbol files.")
}

# Serial on purpose, on two counts. fgsea's default (nproc = 0) hands the work
# to BiocParallel::bpparam(), which is MulticoreParam over every core on the
# machine -- in an allocation of 2 cpus that is exactly the oversubscription
# conf/modules.config documents for SCTENIFOLDKNK_KO. And forked workers would
# undo the set.seed() above. A few hundred ranked genes against a few hundred
# sets is sub-second work, so there is nothing to win here anyway.
#
# progressbar = FALSE because fgsea draws one bar per call and this loop runs
# once per cell type x gene pair, which buries every message above in redraws.
bp_param <- BiocParallel::SerialParam(progressbar = FALSE)

combos <- unique(dr[, c("cell_type", "target")])
results <- list()

for (i in seq_len(nrow(combos))) {
  this_cell   <- combos$cell_type[i]
  this_target <- combos$target[i]

  sub <- dr[dr$cell_type == this_cell & dr$target == this_target, ]

  # fgsea stops with "Not all stats values are finite numbers", so non-finite
  # values are dropped here rather than taken to it.
  sub <- sub[!is.na(sub$gene) & nzchar(sub$gene) & is.finite(sub$log2FC), ]

  # Genes sorted on log2FC, furthest-moved first. The order is the paper's --
  # log2FC is a monotone transform of the manifold-alignment distance, so
  # sorting on either puts the same genes in the same sequence -- but the
  # values handed to fgsea straddle zero rather than piling up above it, which
  # is what makes the score two-tailed. A gene appearing twice keeps its
  # largest value, since a duplicated name in the stats vector would make the
  # ranking ambiguous.
  sub <- sub[order(sub$log2FC, decreasing = TRUE), ]
  sub <- sub[!duplicated(sub$gene), ]

  # The knocked-out gene is dropped from its own ranking. It is first by
  # construction -- zeroing its edges is what the distance measures -- so
  # leaving it in hands a guaranteed top-of-list hit to every gene set that
  # annotates it, which is exactly the set a reader would most want to believe.
  # SCTENIFOLDKNK_KO already drops it from the DR table for the same reason, so
  # this only matters for a table from a run before that. A ';'-joined target
  # is a joint knockout, so every one of its genes goes.
  sub <- sub[!sub$gene %in% strsplit(this_target, ";", fixed = TRUE)[[1]], ]

  if (nrow(sub) < min_size) {
    message("Only ", nrow(sub), " ranked gene(s) for ", this_target, " in ",
            this_cell, "; skipping.")
    next
  }

  stats <- stats::setNames(sub$log2FC, sub$gene)

  # scoreType = "std", fgsea's default, because log2FC has a real negative tail
  # -- genes the knockout moved less than the run's average gene. Both ends of
  # the ranking are therefore meaningful and both get scored, which is where a
  # negative NES comes from. ("pos", which this used while the ranking was the
  # raw distance, would now throw that tail away.)
  #
  # BPPARAM rather than nproc: fgsea's own progress bar is drawn once per call
  # and this loop runs once per cell type x gene pair, which buries the
  # messages above in bar redraws. Passing both is a warning ("Both nproc and
  # BPPARAM arguments were set"), so only BPPARAM goes in.
  #
  # p-values floor at fgsea's default eps of 1e-50, where log2err comes back
  # NA. That is a strong hit reported conservatively, not a failure; the report
  # only ever asks whether p.adj clears 0.25, GSEA's own significance
  # convention.
  res <- tryCatch(
    fgsea::fgsea(
      pathways  = pathways,
      stats     = stats,
      minSize   = min_size,
      maxSize   = max_size,
      scoreType = "std",
      BPPARAM   = bp_param
    ),
    # error only, never warning: ties in the distances are routine (a
    # dRegulation table repeats values) and fgsea says so with a warning. A
    # warning handler here would abort the call and silently skip a perfectly
    # good combination.
    error = function(e) {
      message("fgsea failed for ", this_target, " in ", this_cell, ": ",
              conditionMessage(e))
      NULL
    }
  )

  if (is.null(res) || nrow(res) == 0) {
    message("No gene set passed the size filters for ", this_target, " in ",
            this_cell, ".")
    next
  }

  res <- as.data.frame(res, stringsAsFactors = FALSE)

  # leadingEdge comes back as a list column, which write.table cannot put in a
  # cell. ';' joins it the same way a combined --target line names its genes.
  res$leadingEdge <- vapply(
    res$leadingEdge,
    function(genes) paste(genes, collapse = ";"),
    character(1)
  )

  results[[length(results) + 1]] <- data.frame(
    cell_type = this_cell,
    target    = this_target,
    res,
    stringsAsFactors = FALSE
  )
}

if (length(results) == 0) {
  bail("No combination produced an enrichment result.")
}

gsea <- do.call(rbind, results)
gsea <- gsea[order(gsea$target, gsea$cell_type, gsea$pval), names(empty_gsea)]

message("Wrote ", nrow(gsea), " row(s) for ", length(results),
        " target x cell type combination(s); ",
        sum(gsea$padj < 0.25, na.rm = TRUE), " at FDR < 0.25.")

write_gsea(gsea)
