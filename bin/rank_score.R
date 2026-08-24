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
binding <- args[5]
top_n <- args[6]
rds_files <- args[7:length(args)]

cell_types <- sub("_weight.*", "", basename(rds_files))

target_rank <- strsplit(target[1], split = ";")[[1]]
target_id <- gsub("[^A-Za-z0-9_.-]+", "_", paste(target_rank, collapse = "_"))

#add a print statement to check the target variable before processing the targets
print(paste("Processing targets:", paste(target_rank, collapse = ", ")))

sc_objs <- lapply(rds_files, readRDS)

if (seuratObj == 'AML_object.rda') {
  load(seuratObj)
  seuratObj <- seuratObj[c(VariableFeatures(seuratObj)[1:200], target_rank[1]),]
} else {
  seuratObj <- readRDS(seuratObj)
}

obj <- CreateScRank(input = seuratObj,
                    species = species, 
                    cell_type = column,
                    target = target_rank[1])

obj@net <- sc_objs
names(obj@net) <- cell_types

obj@para$ct.keep = names(obj@net)

# saveRDS(obj, paste0("merged_obj.", target_id, ".RDS"))

# Top connections per cell type -------------------------------------------
# scRank's init_mod() reads a target's neighbours straight off its row in the
# cell type network -- colnames(net)[abs(net[target, ]) > 0] -- and then goes
# on to modularise a single cell type. Only that neighbour row is wanted here,
# for every cell type rather than one, so the row is read directly. This runs
# before rank_celltype so the connections survive a target whose ranking
# fails; rank_celltype does not touch obj@net anyway.
top_n <- suppressWarnings(as.integer(top_n))
if (is.na(top_n) || top_n < 1) {
  top_n <- 15
}

connections <- list()

for (ct in cell_types) {
  net <- obj@net[[ct]]

  if (is.null(net) || is.null(rownames(net)) || nrow(net) == 0) {
    message("No network for ", ct, "; no connections recorded.")
    next
  }

  for (gene in target_rank) {
    if (!(gene %in% rownames(net))) {
      message("Target ", gene, " is absent from the ", ct, " network.")
      next
    }

    weights <- as.numeric(net[gene, ])
    names(weights) <- colnames(net)

    # A gene is not its own connection, and a zero weight is the absence of an
    # edge rather than a weak one.
    weights <- weights[names(weights) != gene]
    weights <- weights[is.finite(weights) & weights != 0]

    if (length(weights) == 0) {
      message("Target ", gene, " has no non-zero edge in ", ct, ".")
      next
    }

    # Ranked on magnitude: a strong repressive edge matters as much as a strong
    # activating one, and the direction is kept in the weight itself.
    weights <- weights[order(abs(weights), decreasing = TRUE)]
    weights <- weights[seq_len(min(top_n, length(weights)))]

    connections[[length(connections) + 1]] <- data.frame(
      cell_type   = ct,
      target      = target[1],
      binding     = binding,
      target_gene = gene,
      partner     = names(weights),
      weight      = as.numeric(weights),
      rank        = seq_along(weights),
      stringsAsFactors = FALSE
    )
  }
}

top_connections <- if (length(connections) > 0) {
  do.call(rbind, connections)
} else {
  data.frame(
    cell_type = character(), target = character(), binding = character(),
    target_gene = character(), partner = character(), weight = numeric(),
    rank = integer(), stringsAsFactors = FALSE
  )
}

write.table(
  top_connections,
  paste0("top_connections.", target_id, ".txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t"
)

all_ranks <- data.frame()

for (target_sc in target) {
  message("Target ", target_sc, " found. Proceeding with rank_celltype.") 
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
      binding = binding,
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
  paste0("perbscore_all_targets.", target_id, ".txt"), 
  quote = FALSE, 
  row.names = FALSE, 
  col.names = TRUE, 
  sep = "\t"
)
