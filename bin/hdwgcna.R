#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(WGCNA)
})

args <- commandArgs(trailingOnly = TRUE)

seuratObj <- args[1]
column    <- args[2]
n_cores   <- as.integer(args[3])
cut_ratio <- as.numeric(args[4])
min_cells <- as.integer(args[5])

# Metacell aggregation settings
k          <- 25
max_shared <- 10

# One soft power for every cell identity, see the note at TestSoftPowers below.
# 6 is WGCNA's standing recommendation for an unsigned network on few samples.
soft_power <- 6

allowWGCNAThreads(nThreads = n_cores)

cell_type <- sub("\\.RDS$", "", basename(seuratObj))

sc_obj <- readRDS(seuratObj)
n_cells <- ncol(sc_obj)

# Cell identities that cannot yield a usable network are skipped instead of
# failing the run, and are simply absent from the final ranking.
drop_celltype <- function(reason) {
  message("Dropping ", cell_type, ": ", reason)
  quit(save = "no", status = 0)
}

if (n_cells < min_cells) {
  drop_celltype(paste0(n_cells, " cells available, min_cells is ", min_cells))
}

# gene4use is set by DOWNSAMPLE and already contains the requested targets. It
# is the gene universe every cell type network must share, so scRank can align
# them against each other.
genes_4_use <- sc_obj@misc$gene4use
genes_4_use <- intersect(genes_4_use, rownames(sc_obj))

if (length(genes_4_use) > 6000) {
  message("Constructing a ", length(genes_4_use), " x ", length(genes_4_use),
          " TOM in a single block, this needs a lot of memory")
}

# hdWGCNA fails on cell identities that are too small or too homogeneous to
# aggregate into metacells, which is expected rather than exceptional here.
built <- tryCatch({
  sc_obj <- SetupForWGCNA(
    sc_obj,
    gene_select = "custom",
    features    = genes_4_use,
    wgcna_name  = "network"
  )

  sc_obj <- MetacellsByGroups(
    seurat_obj  = sc_obj,
    group.by    = column,
    k           = k,
    min_cells   = min_cells,
    max_shared  = max_shared,
    ident.group = column
  )
  sc_obj <- NormalizeMetacells(sc_obj)

  sc_obj <- SetDatExpr(
    sc_obj,
    group_name = as.character(sc_obj@meta.data[[column]][1]),
    group.by   = column,
    assay      = "RNA",
    layer      = "data"
  )

  # An unsigned network is used on purpose. A signed adjacency is
  # ((1 + cor) / 2)^power, which drives anti-correlated pairs towards zero, so
  # repressive edges are cut before they ever reach scRank and the sign recovery
  # below has nothing to recover. Unsigned keeps |cor|^power, letting magnitude
  # carry the strength of the association and the recovered sign carry its
  # direction, which is the same semantics as scRank's signed regression betas.
  # The power is deliberately NOT taken from the estimate. powerEstimate is the
  # first power whose scale-free fit clears RsquaredCut, so it lands somewhere
  # different for every cell identity and is NA whenever the fit never clears it
  # at all. A per-cell-type exponent means a per-cell-type |cor|^power
  # transform, and the resulting edge weights are then not on a common scale.
  # scRank compares perturbation scores across cell identities, so the exponent
  # has to be the same everywhere. The estimate is still computed and logged so
  # the per-cell-type fit stays auditable.
  sc_obj <- TestSoftPowers(sc_obj, networkType = "unsigned")
  power_est <- GetActiveWGCNA(sc_obj)$sft$powerEstimate
  message("Soft power: using fixed ", soft_power, ", per-cell-type estimate was ",
          if (is.null(power_est) || is.na(power_est) || is.infinite(power_est)) "NA" else power_est)

  sc_obj <- ConstructNetwork(
    sc_obj,
    soft_power    = soft_power,
    setDatExpr    = FALSE,
    overwrite_tom = TRUE,
    tom_name      = "network",
    minModuleSize = 20,
    networkType   = "unsigned",
    TOMType       = "unsigned",
    maxBlockSize  = length(genes_4_use)
  )
  TRUE
}, error = function(e) {
  drop_celltype(paste0("hdWGCNA failed: ", conditionMessage(e)))
})

TOM <- as.matrix(GetTOM(sc_obj))
datExpr <- GetDatExpr(sc_obj)

common <- intersect(rownames(TOM), colnames(datExpr))
if (length(common) < 2) {
  drop_celltype("fewer than two genes shared between the TOM and the metacell matrix")
}

# --- 1. Sign recovery -------------------------------------------------------
# A TOM is always positive, so activation and repression collapse onto each
# other. scRank reads edge sign directly (agonist mode saturates positive
# out-edges, and the score weights edges by magnitude), so the sign is put back
# while keeping the TOM magnitude, which is what the shared-neighbour denoising
# actually buys us over a plain correlation.
cormat <- WGCNA::cor(datExpr[, common], use = "pairwise.complete.obs")
cormat[is.na(cormat)] <- 0

weight <- TOM[common, common] * sign(cormat)
weight[is.na(weight)] <- 0

# --- 2. Pad to the shared gene universe -------------------------------------
# hdWGCNA drops genes during its own QC, which would leave each cell type with a
# different gene set. scRank's .align_net refuses networks whose features are
# not identical, so missing genes come back as zero rows and columns.
full <- matrix(
  0,
  nrow = length(genes_4_use),
  ncol = length(genes_4_use),
  dimnames = list(genes_4_use, genes_4_use)
)
keep <- intersect(common, genes_4_use)
full[keep, keep] <- weight[keep, keep]

n_padded <- length(genes_4_use) - length(keep)
if (n_padded > 0) {
  message(n_padded, " of ", length(genes_4_use),
          " genes were not in the hdWGCNA network and were padded with zeros")
}

# --- 3. Sparsify ------------------------------------------------------------
# A TOM is fully dense. scRank derives degree and entropy from the count of
# non-zero edges and sums its marker distance over the target's neighbourhood,
# so a dense matrix makes every gene a neighbour of every other one.
#
# The quantile is taken over the genes hdWGCNA actually kept rather than over
# the padded matrix. How much padding step 2 added varies a lot between cell
# identities, so a quantile spanning those zeros makes cut_ratio buy a
# different amount of sparsification per cell type. Worse, once the zero
# fraction exceeds cut_ratio the threshold is 0, the comparison below is never
# true, and the network is not sparsified at all.
diag(full) <- 0
sub <- full[keep, keep, drop = FALSE]

nz <- sub[sub != 0]
if (length(nz) == 0) {
  drop_celltype("no non-zero edges in the hdWGCNA network")
}

threshold <- quantile(abs(nz), cut_ratio, na.rm = TRUE)
sub[abs(sub) < threshold] <- 0
full[keep, keep] <- sub

# --- 4. Rank-normalise the surviving edges ----------------------------------
# scRank's manifold alignment offsets the network by 1 before building the
# Laplacian and propagates the perturbation multiplicatively along paths, so
# the score is driven by the typical edge weight, not the largest one. Scaling
# by max |w| pins the maximum to 1 but leaves the rest of the distribution
# free, and that distribution tracks metacell count: a cell identity with few
# metacells has correlations inflated by small-sample noise, which survive
# |cor|^power far better than the true, more modest correlations of a
# well-sampled identity. In practice that left the median surviving weight two
# orders of magnitude apart across cell types and made the perturbation score
# an almost perfectly monotone function of how many metacells the identity
# could produce.
#
# Ranking gives every cell type the same edge-weight distribution while keeping
# the ordering and the sign within each network, which is all scRank reads off
# a network, so the cross-cell-type comparison is no longer confounded by
# sample size.
i <- which(full != 0)
if (length(i) == 0) {
  drop_celltype("no non-zero edges left after sparsification")
}
full[i] <- sign(full[i]) * rank(abs(full[i])) / length(i)

n_metacells <- ncol(GetMetacellObject(sc_obj))

# Save the object
saveRDS(full, file = paste0(cell_type, "_weight_hdWGCNA_", n_metacells, ".rds"))
