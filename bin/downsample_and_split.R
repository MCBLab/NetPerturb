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

targets <- readLines(targets)
target <- strsplit(targets[1], split = ";")[[1]]

targets <- unlist(strsplit(targets, split = ";"))

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

if (seuratObj == 'AML_object.rda') {
    load(seuratObj)
    seuratObj <- seuratObj[c(VariableFeatures(seuratObj)[1:200], target),]
} else {
    seuratObj <- readRDS(seuratObj)
}

seuratObj <- as_v3_assay(seuratObj)

# Downsample cells by celltype
downsampled_cells <- seuratObj@meta.data %>% tibble::rowid_to_column("id_cell") %>%
  group_by(!!sym(column)) %>%
  slice_sample(n = n_cells) %>%
  pull(id_cell)

ncells <- length(downsampled_cells)
seurat_downsample <- seuratObj[, downsampled_cells]

# Gene set carried downstream. This reproduces the feature selection that
# scRank::CreateScRank does internally (R/method.R), so no scRank object has to
# be built here: highly variable genes, every TF and every drug target known for
# the species, and the requested targets.
if (!species %in% c("human", "mouse")) {
  stop("species must be 'human' or 'mouse', got: ", species)
}

# Variable features the object already carries are reused; they are computed
# only when it has none. nfeatures cannot exceed the number of genes present.
hvg <- if (length(VariableFeatures(seurat_downsample)) > 0) {
  VariableFeatures(seurat_downsample)
} else {
  VariableFeatures(FindVariableFeatures(seurat_downsample,
                                        selection.method = "vst",
                                        nfeatures = min(2000, nrow(seurat_downsample)),
                                        verbose = FALSE))
}
hvg <- head(hvg, 2000)

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

split_obj <- SplitObject(seurat_downsample, split.by = column)

# Create Seurat split objects
sc_obj <- lapply(split_obj, function(seuobj){
  obj <- seuobj
  obj@misc$gene4use <- genes_4_use
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
    obj <- RunPCA(obj, npcs = npcs, verbose = FALSE)
    obj <- RunUMAP(obj, dims = seq_len(npcs), verbose = FALSE)
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

