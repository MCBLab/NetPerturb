#!/usr/bin/env Rscript
library(Seurat)
library(dplyr)
library(scTenifoldNet)

args <- commandArgs(trailingOnly = TRUE)

seuratObj <- args[1]
n_cores <- args[2]
n_cores <- as.integer(n_cores)

cell_type <- sub(".RDS", "", seuratObj)

sc_obj <- readRDS(seuratObj)

mat <- as.matrix(sc_obj[sc_obj@misc$gene4use]@assays$RNA$counts)
mat <- mat[,colSums(mat) > 30]

n_cells <- dim(sc_obj)[2]

X <- cpmNormalization(mat)

n_cells <- min(1000, dim(sc_obj)[2])

xList <- makeNetworks(X = X, nCells = n_cells, nNet = 30,
                        nComp = 5, scaleScores = TRUE,
                        symmetric = FALSE, q = 0.9,
                        nCores = n_cores)

tensorOut <- tensorDecomposition(xList = xList, K = 5,
                                 nDecimal = 2, maxIter = 2e3,
                                 maxError = 1e-6)

weight <- as.matrix(tensorOut$X)

# Save the object
saveRDS(weight, file = paste0(cell_type, "_weight_sctenifoldnet_", n_cells, ".rds"))
