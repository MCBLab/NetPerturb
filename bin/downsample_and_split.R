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

if (!column %in% colnames(seuratObj@meta.data)) {
  stop("--column '", column, "' is not a metadata column of this object. ",
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
stored_hvg <- VariableFeatures(seurat_downsample)
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

# Drop mitochondrial and ribosomal genes. The match is guarded because `-x` on
# an empty index vector would empty the whole set instead of removing nothing.
mt_rb <- grep("^RP[[:digit:]]+|^RPL|^RPS|^MT-", toupper(genes_4_use))
if (length(mt_rb) > 0) {
  genes_4_use <- genes_4_use[-mt_rb]
}

# Targets are added back after that filter so a target is never dropped for
# looking ribosomal, and anything missing from the object is then dropped.
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
# count there, which are kept below that threshold.
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

write.table(cell_counts, "cell_counts.tsv", quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)

# Expression of every target gene in every kept identity, for the report's
# heatmap. One row per single gene: a ';'-joined target is split into its
# genes, each shown on its own. Only the retained cells count, since those are
# what the networks were built from. Log-normalised here from the counts
# (per 10,000, then log1p) rather than read from the data layer, which holds
# raw counts in an object that was never normalised; the library size is the
# cell's total over every gene, as NormalizeData() would take it.
target_counts <- seurat_downsample[["RNA"]]@counts[targets, , drop = FALSE]
lib_size <- Matrix::colSums(seurat_downsample[["RNA"]]@counts)
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
write.table(target_expression, "target_expression.tsv", quote = FALSE, sep = "\t",
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

