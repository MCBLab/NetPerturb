#!/usr/bin/env Rscript
library(Seurat)
library(scRank)
library(dplyr)

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
set.seed(42)
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

# Save each object with a cleaned file name
invisible(lapply(names(sc_obj), function(name) {
  file_name <- paste0(clean_name(name), ".RDS")
  saveRDS(sc_obj[[name]], file = file_name)
}))

