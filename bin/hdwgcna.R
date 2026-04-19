#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(WGCNA)
})

args <- commandArgs(trailingOnly = TRUE)
input_rds <- args[1]
column <- args[2]
n_cores <- as.numeric(args[3])
targets_file <- args[4]

allowWGCNAThreads(nThreads = n_cores)

cell_type <- sub("\\.RDS$|\\.rds$", "", basename(input_rds))
seuratObj <- readRDS(input_rds)
n_cells <- ncol(seuratObj)

if (n_cells < 150) {
  quit(save = "no", status = 0)
}

gene_4_use <- seuratObj@misc$gene4use

if (file.exists(targets_file)) {
  targets <- readLines(targets_file)
  targets <- targets[targets != ""]
  gene_4_use <- unique(c(gene_4_use, intersect(targets, rownames(seuratObj))))
}

nome_arquivo <- paste0(cell_type, "_weight_hdWGCNA_", n_cells, ".rds")

seuratObj <- SetupForWGCNA(
  seuratObj,
  gene_select = "custom",
  features = gene_4_use,
  wgcna_name = "network"
)

seuratObj <- MetacellsByGroups(
  seurat_obj = seuratObj,
  group.by = column,    
  k = 25,
  min_cells = 100, 
  max_shared = 10,
  ident.group = column  
)
seuratObj <- NormalizeMetacells(seuratObj)

seuratObj <- SetDatExpr(
  seuratObj,
  group_name = unique(seuratObj@meta.data[[column]])[1],
  group.by = column,
  assay = 'RNA',
  layer = 'data' 
)

seuratObj <- TestSoftPowers(seuratObj, networkType = 'signed')
power_est <- GetActiveWGCNA(seuratObj)$sft$powerEstimate
final_power <- if(is.null(power_est) || is.na(power_est) || is.infinite(power_est)) 9 else power_est

seuratObj <- ConstructNetwork(
  seuratObj, 
  soft_power = final_power,
  setDatExpr = FALSE, 
  overwrite_tom = TRUE, 
  tom_name = "network", 
  minModuleSize = 20
)
n_metacells <- ncol(GetMetacellObject(seuratObj))
nome_arquivo <- paste0(cell_type, "_weight_hdWGCNA_", n_metacells, ".rds")

weight <- as.matrix(GetTOM(seuratObj))
weight <- weight[colnames(weight), ]

saveRDS(weight, file = nome_arquivo)