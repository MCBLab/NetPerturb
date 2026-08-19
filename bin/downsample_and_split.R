#!/usr/bin/env Rscript
library(Seurat)
library(scRank)
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

if (seuratObj == 'AML_object.rda') {
    load(seuratObj)
    seuratObj <- seuratObj[c(VariableFeatures(seuratObj)[1:200], target),]
} else {
    seuratObj <- readRDS(seuratObj)
}

# Downsample cells by celltype
downsampled_cells <- seuratObj@meta.data %>% tibble::rowid_to_column("id_cell") %>%
  group_by(!!sym(column)) %>%
  slice_sample(n = n_cells) %>%
  pull(id_cell)

ncells <- length(downsampled_cells)
seurat_downsample <- seuratObj[, downsampled_cells]

non_targets <- targets[!targets %in% rownames(seuratObj)]

obj <- CreateScRank(input = seurat_downsample,
                    species = species, 
                    cell_type = column,
                    target = target)

genes_4_use <- unique(c(obj@para$gene4use, targets))
genes_4_use <- setdiff(genes_4_use, non_targets)

split_obj <- SplitObject(seurat_downsample, split.by = column)

# Create scRank object
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

