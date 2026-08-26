#!/usr/bin/env Rscript
library(Seurat)
library(dplyr)
library(GENIE3)

args <- commandArgs(trailingOnly = TRUE)

seuratObj <- args[1]
n_cores <- args[2]
n_cores <- as.integer(n_cores)
n_trees <- if (length(args) >= 3) as.integer(args[3]) else 1000L

cell_type <- sub(".RDS", "", seuratObj)

sc_obj <- readRDS(seuratObj)

# GENIE3 não tem parâmetro de seed próprio; com nCores > 1 ela usa
# doRNG::`%dorng%`, que deriva os streams de RNG por worker do estado do
# RNG no momento da chamada — set.seed() aqui cobre tanto o caminho
# single-core quanto o multi-core.
set.seed(42)
weight <- GENIE3(as.matrix(sc_obj[sc_obj@misc$gene4use]@assays$RNA$data), nCores = n_cores, nTrees = n_trees)
weight <- weight[colnames(weight),]

n_cells <- dim(sc_obj)[2]

# Save the object
saveRDS(weight, file = paste0(cell_type, "_weight_GENIE3_", n_cells, ".rds"))
