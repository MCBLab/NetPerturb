#!/usr/bin/env Rscript
library(Seurat)
library(dplyr)
library(GENIE3)

args <- commandArgs(trailingOnly = TRUE)

seuratObj <- args[1]
n_cores <- args[2]
n_cores <- as.integer(n_cores)

# --seed for the random forests. GENIE3 runs its workers under doRNG, so this
# one set.seed() gives the same network whatever --n_cores is. Missing or
# unparseable falls back to 1, the pipeline default.
seed <- suppressWarnings(as.integer(if (length(args) >= 3) args[3] else NA))
if (is.na(seed)) seed <- 1L

cell_type <- sub(".RDS", "", seuratObj)

sc_obj <- readRDS(seuratObj)

set.seed(seed)
weight <- GENIE3(as.matrix(sc_obj[sc_obj@misc$gene4use]@assays$RNA$data), nCores = n_cores)
weight <- weight[colnames(weight),]

n_cells <- dim(sc_obj)[2]

# Save the object
saveRDS(weight, file = paste0(cell_type, "_weight_GENIE3_", n_cells, ".rds"))
