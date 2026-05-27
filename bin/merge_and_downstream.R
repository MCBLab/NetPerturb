#!/usr/bin/env Rscript
library(Seurat)
library(scRank)
library(dplyr)
library(GENIE3)


args <- commandArgs(trailingOnly = TRUE)

seuratObj <- args[1]
target <- args[2]
species <- args[3]
column <- args[4]
rds_files <- args[5:length(args)]

cell_types <- sub("_weight.*", "", basename(rds_files))

target_rank <- strsplit(target[1], split = ";")[[1]]

sc_objs <- lapply(rds_files, readRDS)

if (seuratObj == 'AML_object.rda') {
  load(seuratObj)
  seuratObj <- seuratObj[c(VariableFeatures(seuratObj)[1:200], target),]
} else {
  seuratObj <- readRDS(seuratObj)
}

obj <- CreateScRank(input = seuratObj,
                    species = species, 
                    cell_type = column,
                    target = target_rank)

obj@net <- sc_objs
names(obj@net) <- cell_types

obj@para$ct.keep = names(obj@net)

saveRDS(obj, "merged_obj.RDS")

all_ranks <- data.frame()

for (target_sc in target) {
  message("Processing target: ", target_sc)
  
  # Set the target
  obj@para$target <- strsplit(target_sc, split = ";")[[1]]
  
  # Try running rank_celltype
  tryCatch({
    obj <- rank_celltype(obj, n.core = 4)
    
    # Extract data and convert to long format
    perb_scores <- obj@cell_type_rank$perb_score
    df_long <- data.frame(
      cell_type = cell_types,
      target = target_sc,
      perb_score = as.numeric(perb_scores)
    )
    
    # Append to the main data frame
    all_ranks <- rbind(all_ranks, df_long)
    
    message("Finished: ", target_sc)
    
  }, error = function(e) {
    message("Failed for target ", target_sc, ": ", e$message)
  })
}

# Save results (even if partial)
write.table(
  all_ranks, 
  paste0("perbscore_all_targets.", gsub(";", "_", target), ".txt"), 
  quote = FALSE, 
  row.names = FALSE, 
  col.names = TRUE, 
  sep = "\t"
)
