#!/usr/bin/env Rscript
library(Seurat)
library(dplyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

# Inputs
seuratObj <- args[1]
targets <- args[2]
column <- args[3]
species <- args[4]
n_cells <- as.integer(args[5])

# 0 disables the filter and keeps every identity; anything unparseable falls
# back to the default rather than propagating an NA into the comparison below.
min_cells <- suppressWarnings(as.integer(args[6]))
if (is.na(min_cells) || min_cells < 0) {
  min_cells <- 150L
}

# Assay holding the counts; see select_assay() below.
assay <- if (length(args) >= 7 && nzchar(args[7])) args[7] else "RNA"

# Number of highly variable genes carried into gene4use; see the HVG block
# below. Unparseable or below 1 falls back to the default.
n_hvg <- suppressWarnings(as.integer(if (length(args) >= 8) args[8] else NA))
if (is.na(n_hvg) || n_hvg < 1) {
  n_hvg <- 2000L
}

# --seed: which cells slice_sample() keeps below, and the UMAP computed when
# the object has none. Missing or unparseable falls back to 1, the pipeline
# default.
seed <- suppressWarnings(as.integer(if (length(args) >= 9) args[9] else NA))
if (is.na(seed)) {
  seed <- 1L
}

# --sctknk_min_pct: the detection rate a gene needs in an identity's cells to
# enter its scTenifoldKnk network (see sctenifoldknk_build.R). Only counted
# here, for the report. Missing or out of [0, 1) falls back to 0.05,
# scTenifoldKnk's own qc_minPCT.
sctknk_min_pct <- suppressWarnings(as.numeric(if (length(args) >= 10) args[10] else NA))
if (is.na(sctknk_min_pct) || sctknk_min_pct < 0 || sctknk_min_pct >= 1) {
  sctknk_min_pct <- 0.05
}
# --batch: metadata column naming the donor, sample or batch of each cell,
# used only for the Data quality section. Missing, empty or "null" (what
# Nextflow passes for an unset parameter) means none.
batch <- if (length(args) >= 11 && nzchar(args[11]) && args[11] != "null") args[11] else NA_character_
# Mitochondrial and ribosomal protein genes, by symbol. They are removed from
# the object as soon as it is loaded (below), so no step of the pipeline ever
# sees them. Symbols are uppercased first, so mouse (mt-, Rps, Rpl) matches the
# same patterns as human (MT-, RPS, RPL).
#   mitochondrial: MT- plus the "." or "_" some converters write in its place
#     (make.names() turns MT-CO1 into MT.CO1). A bare ^MT would also take
#     nuclear genes such as MT2A, MTCH2 and MTIF3.
#   ribosomal: the cytosolic ribosomal proteins only -- RPS/RPL followed by a
#     number and an optional A/X/Y variant (RPL13A, RPS4X, RPS4Y1), plus RPLP0-2
#     and RPSA. A bare ^RPS/^RPL also took the RPS6K kinases (RPS6KA1/3/6
#     are drug targets in scRank's table), the RPL*L paralogues, RPS27L,
#     RPS19BP1 and pseudogenes like RPSAP58.
is_mt_rb <- function(genes) {
  g <- toupper(genes)
  grepl("^MT[-._]", g) |
    grepl("^RP[SL][0-9]+[AXY]?[0-9]*$|^RPLP[0-2]$|^RPSA$", g)
}

if (is.na(n_cells) || n_cells < 1) {
  stop("--n_cells must be a positive whole number, got: ", args[5])
}

# One target per line, ';' joining genes perturbed together. Blank lines and
# stray whitespace are dropped, as main.nf used to do before reading the
# targets DOWNSAMPLE writes instead. warn = FALSE: a missing final newline is
# not worth a warning.
target_lines <- trimws(readLines(targets, warn = FALSE))
target_lines <- target_lines[nzchar(target_lines)]
target_genes <- lapply(strsplit(target_lines, split = ";"),
                       function(g) unique(trimws(g[nzchar(trimws(g))])))

# The container pairs Seurat v4 with SeuratObject v5, so an object saved with a
# v5 assay is invisible to every Seurat v4 entry point: they resolve assays with
# FilterObjects(classes.keep = "Assay"), Assay5 does not inherit from Assay, and
# the search comes back empty ("RNA is not an assay present in the given object.
# Available assays are:"). Converting the assay once here keeps the rest of this
# script, and the per-cell-type objects it writes, on ground Seurat v4 handles.
#
# The assay cannot be replaced under its own name while it is still an Assay5,
# so the converted copy goes in beside it and is renamed once the original is
# gone. as() joins split layers on the way, which reading a counts layer
# directly would not: LayerData() returns only the first layer of an object
# whose layers are split by sample.
as_v3_assay <- function(obj) {
  assay <- DefaultAssay(obj)
  if (!inherits(obj[[assay]], "Assay5")) {
    return(obj)
  }
  message("Converting v5 assay '", assay, "' to a v3 assay for Seurat v4.")
  key <- Key(obj[[assay]])
  tmp <- paste0(assay, ".v3")
  obj[[tmp]] <- as(obj[[assay]], "Assay")
  DefaultAssay(obj) <- tmp
  obj[[assay]] <- NULL
  obj <- do.call(RenameAssays, c(list(object = obj), setNames(assay, tmp)))
  # The original key is still taken while both assays coexist, so the copy is
  # given a derived one; put the original back now that it is free again.
  Key(obj[[assay]]) <- key
  obj
}

# scRank, GENIE3, hdWGCNA and scTenifoldKnk all read the assay named "RNA"
# (scRank hardcodes GetAssayData(assay = "RNA")), so an object whose counts sit
# under another name -- e.g. "originalexp" from a SingleCellExperiment
# conversion -- fails in all of them. The --assay chosen is made the default,
# converted like any other, and renamed to "RNA", replacing whatever assay
# already went by that name, so everything downstream reads the right counts.
select_assay <- function(obj, assay) {
  if (!assay %in% Assays(obj)) {
    stop("--assay '", assay, "' is not an assay of this object. ",
         "Available: ", paste(Assays(obj), collapse = ", "))
  }
  DefaultAssay(obj) <- assay
  obj <- as_v3_assay(obj)
  if (assay != "RNA") {
    message("Using assay '", assay, "' as 'RNA'.")
    if ("RNA" %in% Assays(obj)) {
      obj[["RNA"]] <- NULL
    }
    obj <- do.call(RenameAssays, c(list(object = obj), setNames(list("RNA"), assay)))
  }
  # An object converted from AnnData often carries X alone, which lands in
  # `data` and leaves `counts` empty -- and scRank and scTenifoldKnk read
  # counts. Integer values there are raw counts filed in the wrong slot, so
  # they are moved over and normalised; anything else is already transformed
  # and no counts can be recovered from it.
  if (length(obj[["RNA"]]@counts) == 0) {
    x <- obj[["RNA"]]@data
    if (length(x) == 0 || !all(x@x == round(x@x))) {
      stop("Assay '", assay, "' has no counts, and its data layer is not raw ",
           "counts either; pass --assay pointing at an assay with raw counts.")
    }
    message("Assay '", assay, "' has no counts; using its integer data layer ",
            "as counts and log-normalising it.")
    obj[["RNA"]] <- CreateAssayObject(counts = x)
    obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
  }
  obj
}

if (seuratObj == 'AML_object.rda') {
    load(seuratObj)
    seuratObj <- seuratObj[unique(c(VariableFeatures(seuratObj)[1:200],
                                  intersect(unlist(target_genes), rownames(seuratObj)))),]
} else {
    seuratObj <- readRDS(seuratObj)
}

seuratObj <- select_assay(seuratObj, assay)

# Mitochondrial and ribosomal protein genes are removed here, before anything
# else reads the object, so no network, gene set, QC check or figure downstream
# is built on them. They are highly expressed and tightly co-expressed, and
# they pull edges and alignment distance towards themselves whatever the
# question: on Kang 2018 PBMCs ribosomal genes were 3% of the knockout network
# but 37% of its differentially-regulated genes, the same ones for every target.
#
# Two things are kept from before the removal. A gene requested as a target is
# never removed for its name, so it can still be scored and knocked out. And
# each cell's total count over every gene is kept in the metadata, because
# ribosomal genes alone can be a fifth or more of a cell's counts: normalising
# over what is left would shift every gene by a different amount in every
# cell. Subsetting genes leaves the data layer's values as they were, so an
# object normalised before it reached the pipeline keeps that normalisation.
seuratObj$netperturb_lib_size <- Matrix::colSums(seuratObj[["RNA"]]@counts)
mt_rb_genes <- setdiff(rownames(seuratObj)[is_mt_rb(rownames(seuratObj))],
                       unlist(target_genes))
n_genes_input <- nrow(seuratObj)
if (length(mt_rb_genes) > 0) {
  seuratObj <- seuratObj[setdiff(rownames(seuratObj), mt_rb_genes), ]
}
message("Removed ", length(mt_rb_genes), " mitochondrial/ribosomal protein ",
        "gene(s) of ", n_genes_input, "; ", nrow(seuratObj), " left.")

if (!column %in% colnames(seuratObj@meta.data)) {
  stop("--column '", column, "' is not a metadata column of this object. ",
       "Available: ", paste(colnames(seuratObj@meta.data), collapse = ", "))
}
if (!is.na(batch) && !batch %in% colnames(seuratObj@meta.data)) {
  stop("--batch '", batch, "' is not a metadata column of this object. ",
       "Available: ", paste(colnames(seuratObj@meta.data), collapse = ", "))
}

# Identities too small to build a network from are dropped here, before the
# split, so nothing downstream ever receives one: no network is inferred from
# it, and it appears in no score, no table and no figure.
#
# The count compared against --min_cells is the one every downstream method
# actually sees, which is what is left *after* downsampling -- min(identity
# size, --n_cells) -- and not the identity's size in the input object. That is
# the same number --hdwgcna_min_cells is measured against, one step further
# down, so the two thresholds mean the same thing.
identity_sizes <- table(as.character(seuratObj@meta.data[[column]]))
retained <- pmin(as.integer(identity_sizes), n_cells)
names(retained) <- names(identity_sizes)

keep_identities <- names(retained)[retained >= min_cells]
dropped_identities <- setdiff(names(retained), keep_identities)

if (length(dropped_identities) > 0) {
  message("Below --min_cells (", min_cells, "), dropped: ",
          paste0(dropped_identities, " (", retained[dropped_identities],
                 " cell(s))", collapse = ", "), ".")
}

if (length(keep_identities) == 0) {
  stop("No identity in '", column, "' reaches --min_cells (", min_cells,
       "); the largest has ", max(retained), " cell(s) after downsampling",
       if (n_cells < min_cells) {
         paste0(", and --n_cells (", n_cells,
                ") is itself below --min_cells, so nothing could pass")
       } else {
         ""
       },
       ".")
}

message("Keeping ", length(keep_identities), " identity/identities: ",
        paste0(keep_identities, " (", retained[keep_identities],
               " cell(s))", collapse = ", "), ".")

# Downsample cells by celltype
set.seed(seed)
downsampled_cells <- seuratObj@meta.data %>% tibble::rowid_to_column("id_cell") %>%
  filter(!!sym(column) %in% keep_identities) %>%
  group_by(!!sym(column)) %>%
  slice_sample(n = n_cells) %>%
  pull(id_cell)

ncells <- length(downsampled_cells)
seurat_downsample <- seuratObj[, downsampled_cells]

# Cells per identity, for the report's cell count figure: how many the input
# object had, and how many of them every network was built from. An identity
# dropped for --min_cells is listed with 0 used, so the figure shows what was
# left out as well as what was kept. The gene columns are added, and the table
# written, once gene4use is known further down.
used_sizes <- table(as.character(seurat_downsample@meta.data[[column]]))
cell_counts <- data.frame(
  identity = names(identity_sizes),
  n_input  = as.integer(identity_sizes),
  n_used   = as.integer(ifelse(names(identity_sizes) %in% names(used_sizes),
                               used_sizes[names(identity_sizes)], 0L)),
  status   = ifelse(names(identity_sizes) %in% keep_identities,
                    "kept", "dropped (below --min_cells)"),
  stringsAsFactors = FALSE
)

# Target QC ----------------------------------------------------------------
# A gene the object does not carry, or one with no counts in any cell kept
# above, has no edges in any network and nothing for a knockout to move, and
# every method downstream fails on it in its own way (CreateScRank refuses it,
# rank_celltype finds it missing from the network, scTenifoldKnk divides by
# zero). So it is taken out here, once, for the whole run: the targets that
# pass are written to targets_qc.txt, which every later step reads instead of
# --target, and the ones that do not are listed in target_qc.tsv for the
# report. A combined target loses only its failing gene -- "A;B" with B absent
# is the same knockout as "A" -- and is dropped only when none of it is left.
gene_counts <- Matrix::rowSums(seurat_downsample[["RNA"]]@counts)

qc_reason <- function(gene) {
  if (!gene %in% names(gene_counts)) {
    "absent from the expression profile"
  } else if (gene_counts[[gene]] == 0) {
    "zero counts in every retained cell"
  } else {
    NA_character_
  }
}

qc_rows <- list()
passing <- character()

for (i in seq_along(target_lines)) {
  genes   <- target_genes[[i]]
  reasons <- vapply(genes, qc_reason, character(1))
  kept    <- genes[is.na(reasons)]
  failed  <- genes[!is.na(reasons)]

  if (length(kept) > 0) {
    passing <- c(passing, paste(kept, collapse = ";"))
  }

  if (length(failed) > 0) {
    qc_rows[[length(qc_rows) + 1]] <- data.frame(
      target = target_lines[i],
      gene   = failed,
      reason = unname(reasons[failed]),
      action = if (length(kept) > 0) {
        paste0("dropped from the target; analysed as ", paste(kept, collapse = ";"))
      } else {
        "target not analysed"
      },
      stringsAsFactors = FALSE
    )
  }
}

target_qc <- if (length(qc_rows) > 0) {
  do.call(rbind, qc_rows)
} else {
  data.frame(target = character(), gene = character(), reason = character(),
             action = character(), stringsAsFactors = FALSE)
}

if (nrow(target_qc) > 0) {
  message("Not analysed due to QC checking: ",
          paste(unique(paste0(target_qc$gene, " (", target_qc$reason, ")")), collapse = ", "), ".")
}

write.table(target_qc, "target_qc.tsv", quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)

passing <- unique(passing)

if (length(passing) == 0) {
  stop("No target passed QC checking (see target_qc.tsv); nothing to analyse.")
}

writeLines(passing, "targets_qc.txt")

# The first passing target is what CreateScRank is handed downstream; all of
# them go into gene4use.
target  <- strsplit(passing[1], split = ";")[[1]]
targets <- unique(unlist(strsplit(passing, split = ";")))

# Gene set carried downstream. This reproduces the feature selection that
# scRank::CreateScRank does internally (R/method.R), so no scRank object has to
# be built here: highly variable genes, every TF and every drug target known for
# the species, and the requested targets.
if (!species %in% c("human", "mouse")) {
  stop("species must be 'human' or 'mouse', got: ", species)
}

# --n_hvg variable features. Those the object already carries are reused when
# there are at least that many -- the first n_hvg of them -- and computed
# otherwise, since reusing a shorter list would quietly hand back fewer genes
# than were asked for. nfeatures cannot exceed the number of genes present.
stored_hvg <- intersect(VariableFeatures(seurat_downsample), rownames(seurat_downsample))
hvg <- if (length(stored_hvg) >= n_hvg) {
  stored_hvg
} else {
  VariableFeatures(FindVariableFeatures(seurat_downsample,
                                        selection.method = "vst",
                                        nfeatures = min(n_hvg, nrow(seurat_downsample)),
                                        verbose = FALSE))
}
hvg <- head(hvg, n_hvg)
message("Using ", length(hvg), " highly variable gene(s) (--n_hvg ", n_hvg, ").")

utile_database <- scRank::utile_database
tf_gene <- utile_database$Gene_TF[[species]]$Symbol
drug_gene <- if (species == "human") {
  utile_database$Drug_Target$human$Symbol
} else {
  utile_database$Drug_Target$mouse$mousegene
}

genes_4_use <- unique(c(target, hvg, tf_gene, drug_gene))

# Drop the clone-named lncRNAs (RP11-..., RP5-...), which scRank's own
# selection removes alongside the mitochondrial and ribosomal genes already
# gone from the object. Indexed with a logical, so nothing matching removes
# nothing.
genes_4_use <- genes_4_use[!grepl("^RP[0-9]+-", toupper(genes_4_use))]

# Targets are added back after that filter so a target is never dropped for
# its name, and anything missing from the object is then dropped.
genes_4_use <- unique(c(genes_4_use, targets))
genes_4_use <- genes_4_use[genes_4_use %in% rownames(seurat_downsample)]

# Genes per identity, next to its cells in the report. Every network is
# offered the same gene4use, but a gene with no counts in an identity's
# retained cells can have no edges in that identity's network: scTenifoldKnk
# drops it outright, hdWGCNA pads it with zeros, and GENIE3 and scRank give it
# nothing to regress on. So genes_expressed is what each network is really
# built on, and gene4use minus it is what that identity filters out. Dropped
# identities have no network and are left NA.
retained_ident <- as.character(seurat_downsample@meta.data[[column]])
g4u_counts <- seurat_downsample[["RNA"]]@counts[genes_4_use, , drop = FALSE]
genes_expressed <- vapply(cell_counts$identity, function(id) {
  cells <- retained_ident == id
  if (!any(cells)) return(NA_integer_)
  sum(Matrix::rowSums(g4u_counts[, cells, drop = FALSE]) > 0)
}, integer(1))

cell_counts$genes_total     <- nrow(seurat_downsample)
cell_counts$genes_gene4use  <- ifelse(is.na(genes_expressed), NA_integer_, length(genes_4_use))
cell_counts$genes_expressed <- unname(genes_expressed)
cell_counts$genes_filtered  <- cell_counts$genes_gene4use - cell_counts$genes_expressed

# Genes in each identity's scTenifoldKnk network, by the same rule
# sctenifoldknk_build.R applies: every gene of the object detected in more
# than --sctknk_min_pct of the identity's cells, plus the target genes with any
# count there, which are kept below that threshold. The mitochondrial and
# ribosomal protein genes are already gone from the object.
all_counts <- seurat_downsample[["RNA"]]@counts
cell_counts$genes_sctknk <- unname(vapply(cell_counts$identity, function(id) {
  cells <- retained_ident == id
  if (!any(cells)) return(NA_integer_)
  sub <- all_counts[, cells, drop = FALSE]
  detected <- Matrix::rowMeans(sub > 0) > sctknk_min_pct
  kept_targets <- intersect(targets, rownames(sub))
  kept_targets <- kept_targets[Matrix::rowSums(sub[kept_targets, , drop = FALSE]) > 0]
  length(union(rownames(sub)[detected], kept_targets))
}, integer(1)))

# How many mitochondrial/ribosomal protein genes were removed on load, so the
# report can say what genes_total is short of. The same for every identity.
cell_counts$genes_mt_rb <- ifelse(cell_counts$status == "kept", length(mt_rb_genes), NA_integer_)

write.table(cell_counts, "cell_counts.tsv", quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)

# Expression of every target gene in every kept identity, for the report's
# heatmap. One row per single gene: a ';'-joined target is split into its
# genes, each shown on its own. Only the retained cells count, since those are
# what the networks were built from. Log-normalised here from the counts
# (per 10,000, then log1p) rather than read from the data layer, which holds
# raw counts in an object that was never normalised; the library size is the
# cell's total over every gene, as NormalizeData() would take it -- taken
# before the mitochondrial and ribosomal genes were removed, so the values do
# not depend on that removal.
target_counts <- seurat_downsample[["RNA"]]@counts[targets, , drop = FALSE]
lib_size <- seurat_downsample$netperturb_lib_size
target_lognorm <- log1p(t(t(as.matrix(target_counts)) / pmax(lib_size, 1)) * 1e4)

target_expression <- do.call(rbind, lapply(keep_identities, function(id) {
  cells <- retained_ident == id
  data.frame(
    identity       = id,
    gene           = targets,
    avg_expression = unname(rowMeans(target_lognorm[, cells, drop = FALSE])),
    pct_expressing = unname(100 * rowMeans(as.matrix(target_counts[, cells, drop = FALSE]) > 0)),
    stringsAsFactors = FALSE
  )
}))

# Data quality ---------------------------------------------------------------
# Diagnostics for the report, not filters: each one is a reason an identity's
# network, and so every score from it, might be skewed. All are taken over the
# retained cells, since those are what the networks were built from.
#
# Per target and identity, how much the target varies beyond what its mean
# predicts: Seurat's vst standardised variance, over every gene of the
# identity. Around 1 is what sampling noise alone gives, so a target near it
# has no structure for a network to connect it to, however well it is detected.
# The percentile places it among the identity's genes.
# Seurat's vst standardised variance, computed here from the sparse counts
# rather than through FindVariableFeatures()/HVFInfo(), whose arguments differ
# between the Seurat v4 in the container and v5: a loess fit of log10 variance
# on log10 mean gives each gene's expected variance, values are standardised
# by it and clipped above at sqrt(cells), and their variance is returned. Zeros
# are handled in closed form, so the matrix is never made dense.
vst_variance <- function(counts) {
  n  <- ncol(counts)
  mu <- Matrix::rowMeans(counts)
  vr <- (Matrix::rowSums(counts^2) - n * mu^2) / (n - 1)
  ok <- vr > 0
  out <- setNames(rep(NA_real_, nrow(counts)), rownames(counts))
  if (sum(ok) < 10) return(out)
  fit <- loess(log10(vr[ok]) ~ log10(mu[ok]), span = 0.3)
  sd_exp <- sqrt(10^fit$fitted)
  clip <- sqrt(n)
  trip <- Matrix::summary(as(counts[ok, , drop = FALSE], "dgCMatrix"))
  z <- pmin((trip$x - mu[ok][trip$i]) / sd_exp[trip$i], clip)
  nz <- tabulate(trip$i, nbins = sum(ok))
  z0 <- pmin(-mu[ok] / sd_exp, clip)
  s1 <- (n - nz) * z0 + vapply(split(z, factor(trip$i, levels = seq_len(sum(ok)))), sum, numeric(1))
  s2 <- (n - nz) * z0^2 + vapply(split(z^2, factor(trip$i, levels = seq_len(sum(ok)))), sum, numeric(1))
  out[ok] <- (s2 - s1^2 / n) / (n - 1)
  out
}

target_var <- do.call(rbind, lapply(keep_identities, function(id) {
  v <- vst_variance(all_counts[, retained_ident == id, drop = FALSE])
  data.frame(identity = id, gene = targets,
             var_standardized = unname(v[targets]),
             var_percentile = unname(vapply(targets, function(g) {
               if (is.na(v[g])) NA_real_ else 100 * mean(v <= v[g], na.rm = TRUE)
             }, numeric(1))),
             stringsAsFactors = FALSE)
}))
# Joined by position rather than merge(), which would reorder the rows the
# report reads targets in the order they were given.
m <- match(paste(target_expression$identity, target_expression$gene),
           paste(target_var$identity, target_var$gene))
target_expression$var_standardized <- target_var$var_standardized[m]
target_expression$var_percentile   <- target_var$var_percentile[m]

write.table(target_expression, "target_expression.tsv", quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)

# Per identity:
#   median_counts / median_genes: sequencing depth and genes detected per cell.
#   median_mt_rb_frac: share of each cell's counts in the mitochondrial and
#     ribosomal genes removed on load; high values point at stressed or
#     damaged cells, and at how much normalisation leaned on those genes.
#   sparsity_gene4use: share of zeros in the matrix the rank-score networks are
#     fitted on; mostly-zero genes give unstable regressions and correlations.
#   pc1_var_frac / top5_var_frac: variance carried by the leading principal
#     components of that matrix. A large share in one component means one
#     programme, or a mixture of cell states, dominates the identity, and its
#     co-expression is then about that rather than regulation.
#   pc1_depth_rho, max_pc_depth_rho (and which PC): Spearman correlation of
#     component scores with log sequencing depth. When the main axis of
#     variation is depth, co-expression edges mostly are too.
#   with --batch: n_batches, largest_batch_frac, batch_entropy (0 = one batch,
#     1 = evenly split) and batch_r2, the share of the top 10 components'
#     variance explained by batch -- correlation that comes from donors
#     differing, not from regulation within them.
log_depth <- log10(pmax(seurat_downsample$netperturb_lib_size, 1))
lib_after <- Matrix::colSums(all_counts)
identity_qc <- do.call(rbind, lapply(keep_identities, function(id) {
  cells <- retained_ident == id
  cnt   <- all_counts[genes_4_use, cells, drop = FALSE]
  lib   <- seurat_downsample$netperturb_lib_size[cells]
  row <- data.frame(
    identity          = id,
    identity_id       = gsub("[^A-Za-z0-9_\\-]", "_", id),
    n_cells           = sum(cells),
    median_counts     = median(lib),
    median_genes      = median(Matrix::colSums(all_counts[, cells, drop = FALSE] > 0)),
    median_mt_rb_frac = median(1 - lib_after[cells] / pmax(lib, 1)),
    sparsity_gene4use = 1 - Matrix::nnzero(cnt) / length(cnt),
    pc1_var_frac = NA_real_, top5_var_frac = NA_real_,
    pc1_depth_rho = NA_real_, max_pc_depth_rho = NA_real_, max_pc_depth = NA_integer_,
    n_batches = NA_integer_, largest_batch_frac = NA_real_,
    batch_entropy = NA_real_, batch_r2 = NA_real_,
    stringsAsFactors = FALSE
  )
  # Log-normalised by the pre-removal library size, as target_expression is,
  # then scaled per gene; genes with no variance here carry no information.
  x <- log1p(t(t(as.matrix(cnt)) / pmax(lib, 1)) * 1e4)
  x <- x[apply(x, 1, var) > 0, , drop = FALSE]
  n_pc <- min(10, nrow(x) - 1, ncol(x) - 1)
  if (n_pc >= 2) {
    set.seed(seed)
    pca <- tryCatch(irlba::prcomp_irlba(t(x), n = n_pc, center = TRUE, scale. = TRUE),
                    error = function(e) NULL)
    if (!is.null(pca)) {
      total <- nrow(x)  # scaled genes each carry unit variance
      pc_var <- pca$sdev^2
      rho <- apply(pca$x, 2, function(s) suppressWarnings(cor(s, log_depth[cells], method = "spearman")))
      row$pc1_var_frac     <- pc_var[1] / total
      row$top5_var_frac    <- sum(head(pc_var, 5)) / total
      row$pc1_depth_rho    <- rho[1]
      row$max_pc_depth     <- which.max(abs(head(rho, 5)))
      row$max_pc_depth_rho <- rho[row$max_pc_depth]
      if (!is.na(batch)) {
        b <- factor(seurat_downsample@meta.data[[batch]][cells])
        p <- as.numeric(table(b)) / length(b)
        p <- p[p > 0]
        row$n_batches          <- length(p)
        row$largest_batch_frac <- max(p)
        row$batch_entropy      <- if (length(p) > 1) -sum(p * log(p)) / log(length(p)) else 0
        if (length(p) > 1) {
          r2 <- apply(pca$x, 2, function(s) summary(lm(s ~ b))$r.squared)
          row$batch_r2 <- sum(r2 * pc_var) / sum(pc_var)
        }
      }
    }
  }
  row
}))

write.table(identity_qc, "identity_qc.tsv", quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)

split_obj <- SplitObject(seurat_downsample, split.by = column)

# Create Seurat split objects
sc_obj <- lapply(split_obj, function(seuobj){
  obj <- seuobj
  obj@misc$gene4use <- genes_4_use
  # the QC-passing target genes, which SCTENIFOLDKNK_BUILD keeps in its
  # network even below --sctknk_min_pct
  obj@misc$targets <- targets
  return(obj)
})


clean_name <- function(name) {
  gsub("[^A-Za-z0-9_\\-]", "_", name)  # Replace any non-safe character with "_"
}

# UMAP of the cells that survive downsampling, coloured by the identity column
# the run scores on. An embedding the object already carries is reused, so the
# figure matches whatever has been published for this dataset; one is computed
# only when the object has none. This is a QC figure, so a failure to draw it
# must not sink a run that is otherwise fine: it is guarded, and a placeholder
# carrying the reason is written instead.
umap_file <- paste0("umap_", clean_name(column), ".png")

build_umap <- function(obj) {
  reductions <- Reductions(obj)
  embedding <- reductions[tolower(reductions) %in% c("umap", "tsne")]

  if (length(embedding) == 0) {
    message("No UMAP/t-SNE reduction found; computing a UMAP for the plot.")
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, verbose = FALSE)
    obj <- ScaleData(obj, verbose = FALSE)
    # npcs cannot exceed either dimension of the matrix being decomposed, and
    # the downsampled object can be small on both.
    npcs <- max(2, min(30, ncol(obj) - 1, nrow(obj) - 1))
    obj <- RunPCA(obj, npcs = npcs, seed.use = seed, verbose = FALSE)
    obj <- RunUMAP(obj, dims = seq_len(npcs), seed.use = seed, verbose = FALSE)
    embedding <- "umap"
  }

  DimPlot(obj,
          reduction = embedding[1],
          group.by  = column,
          label     = TRUE,
          repel     = TRUE) +
    labs(
      title    = sprintf("Cells retained after downsampling (n = %d)", ncol(obj)),
      subtitle = sprintf("coloured by '%s'", column)
    ) +
    theme(plot.title = element_text(face = "bold"))
}

umap_plot <- tryCatch(
  build_umap(seurat_downsample),
  error = function(e) {
    message("UMAP plot failed: ", conditionMessage(e))
    ggplot() +
      annotate("text", x = 0, y = 0, size = 5,
               label = paste0("UMAP unavailable\n", conditionMessage(e))) +
      theme_void()
  }
)

ggsave(umap_file, umap_plot, width = 8, height = 6, dpi = 150, bg = "white")

# Save each object with a cleaned file name
invisible(lapply(names(sc_obj), function(name) {
  file_name <- paste0(clean_name(name), ".RDS")
  saveRDS(sc_obj[[name]], file = file_name)
}))

