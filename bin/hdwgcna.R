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
  sc_obj <- TestSoftPowers(sc_obj, networkType = "unsigned")
  power_est <- GetActiveWGCNA(sc_obj)$sft$powerEstimate
  final_power <- if (is.null(power_est) || is.na(power_est) || is.infinite(power_est)) 6 else power_est

  sc_obj <- ConstructNetwork(
    sc_obj,
    soft_power    = final_power,
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
# so a dense matrix makes every gene a neighbour of every other one. The
# quantile spans the whole matrix, matching scRank's own Constr_net.
diag(full) <- 0
abs_full <- abs(full)
threshold <- quantile(abs_full, cut_ratio, na.rm = TRUE)
full[abs_full < threshold] <- 0

# --- 4. Scale to max |w| = 1 ------------------------------------------------
# scRank's manifold alignment offsets the network by 1 before building the
# Laplacian, so all the signal sits in the deviation from 1. Raw TOM values are
# orders of magnitude smaller and would collapse into numerical noise.
max_weight <- max(abs(full))
if (max_weight == 0) {
  drop_celltype("no non-zero edges left after sparsification")
}
full <- full / max_weight

n_metacells <- ncol(GetMetacellObject(sc_obj))

# Save the object
saveRDS(full, file = paste0(cell_type, "_weight_hdWGCNA_", n_metacells, ".rds"))
