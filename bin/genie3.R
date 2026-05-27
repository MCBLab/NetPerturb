#!/usr/bin/env Rscript
library(Seurat)
library(dplyr)
library(GENIE3)

args <- commandArgs(trailingOnly = TRUE)

seuratObj <- args[1]
n_cores <- args[2]
n_cores <- as.integer(n_cores)

cell_type <- sub(".RDS", "", seuratObj)

sc_obj <- readRDS(seuratObj)

weight <- GENIE3(as.matrix(sc_obj[sc_obj@misc$gene4use]@assays$RNA$data), nCores = n_cores)
weight <- weight[colnames(weight),]

n_cells <- dim(sc_obj)[2]

# Save the object
saveRDS(weight, file = paste0(cell_type, "_weight_GENIE3_", n_cells, ".rds"))
