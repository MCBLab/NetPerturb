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
assay <- args[7]
rds_files <- args[8:length(args)]

cell_types <- sub("_weight.*", "", basename(rds_files))

target_rank <- strsplit(target[1], split = ";")[[1]]
target_id <- gsub("[^A-Za-z0-9_.-]+", "_", paste(target_rank, collapse = "_"))

# Both output tables, empty but named, for the paths that have nothing to write.
# Explicit columns for the reason the all_ranks initialiser below gives: a bare
# data.frame() writes no header, and MERGE builds its header from whichever
# table it reads first.
empty_ranks <- data.frame(
  cell_type = character(), target = character(), binding = character(),
  perb_score = numeric(), stringsAsFactors = FALSE
)

empty_connections <- data.frame(
  cell_type = character(), target = character(), binding = character(),
  target_gene = character(), partner = character(), weight = numeric(),
  rank = integer(), stringsAsFactors = FALSE
)

write_table <- function(x, prefix) {
  write.table(x, paste0(prefix, ".", target_id, ".txt"), quote = FALSE,
              row.names = FALSE, col.names = TRUE, sep = "\t")
}

# Leaves both tables headered and empty and exits clean. Unlike the
# rank_celltype failures at the bottom of this script -- a worker dying in
# scRank's own parallel backend, which is worth failing loudly so Nextflow
# retries rather than caches -- everything that calls this is deterministic:
# the gene is not in the object and will not be there on a retry either. So
# the target is skipped, MERGE concatenates two header-only tables into
# nothing, and the rest of the run carries on without it.
skip_target <- function(...) {
  message(...)
  write_table(empty_connections, "top_connections")
  write_table(empty_ranks, "perbscore_all_targets")
  quit(save = "no", status = 0)
}

#add a print statement to check the target variable before processing the targets
print(paste("Processing targets:", paste(target_rank, collapse = ", ")))

sc_objs <- lapply(rds_files, readRDS)

# The container pairs Seurat v4 with SeuratObject v5, so an object saved with a
# v5 assay is invisible to every Seurat v4 entry point (see the same fix in
# downsample_and_split.R). This script re-reads the original object directly
# (rather than the already-converted per-celltype splits DOWNSAMPLE wrote), so
# it needs its own conversion before CreateScRank touches it.
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
  Key(obj[[assay]]) <- key
  obj
}

# scRank, GENIE3, hdWGCNA and scTenifoldKnk all read the assay named "RNA"
# (scRank hardcodes GetAssayData(assay = "RNA")), so an object whose counts sit
# under another name -- e.g. "originalexp" from a SingleCellExperiment
# conversion -- fails in all of them. The --assay chosen is made the default,
# converted like any other, and renamed to "RNA", replacing whatever assay
# already went by that name, so everything downstream reads the right counts.
select_assay <- function(obj, assay) {
  if (!assay %in% Assays(obj)) {
    stop("--assay '", assay, "' is not an assay of this object. ",
         "Available: ", paste(Assays(obj), collapse = ", "))
  }
  DefaultAssay(obj) <- assay
  obj <- as_v3_assay(obj)
  if (assay != "RNA") {
    message("Using assay '", assay, "' as 'RNA'.")
    if ("RNA" %in% Assays(obj)) {
      obj[["RNA"]] <- NULL
    }
    obj <- do.call(RenameAssays, c(list(object = obj), setNames(list("RNA"), assay)))
  }
  # An object converted from AnnData often carries X alone, which lands in
  # `data` and leaves `counts` empty -- and scRank and scTenifoldKnk read
  # counts. Integer values there are raw counts filed in the wrong slot, so
  # they are moved over and normalised; anything else is already transformed
  # and no counts can be recovered from it.
  if (length(obj[["RNA"]]@counts) == 0) {
    x <- obj[["RNA"]]@data
    if (length(x) == 0 || !all(x@x == round(x@x))) {
      stop("Assay '", assay, "' has no counts, and its data layer is not raw ",
           "counts either; pass --assay pointing at an assay with raw counts.")
    }
    message("Assay '", assay, "' has no counts; using its integer data layer ",
            "as counts and log-normalising it.")
    obj[["RNA"]] <- CreateAssayObject(counts = x)
    obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
  }
  obj
}

if (seuratObj == 'AML_object.rda') {
  load(seuratObj)
  # Every target gene the object actually carries, not just the first: the
  # profile check below reads rownames(), so a gene dropped here would look
  # unexpressed rather than un-subsetted. intersect() also keeps a target that
  # is missing outright from erroring the subset itself.
  seuratObj <- seuratObj[unique(c(VariableFeatures(seuratObj)[1:200],
                                  intersect(target_rank, rownames(seuratObj)))), ]
} else {
  seuratObj <- readRDS(seuratObj)
}

seuratObj <- select_assay(seuratObj, assay)

# CreateScRank() validates the target it is handed against the expression
# profile and stops with "Please check if the target gene is in the gene
# expression profile." when it is not there, which takes the whole task down
# before either table is written. A gene the object does not carry cannot be
# in a network built from that same object either, so it has no edges to
# remove and contributes nothing to a perturbation score: dropping it and
# scoring the rest of the target gives the same answer as scoring all of it
# would have. That is already how the top-connections loop below treats a gene
# missing from a network, and it is why the target keeps its original name --
# "A;B" with B unexpressed is the same knockout as "A".
usable <- target_rank[target_rank %in% rownames(seuratObj)]

if (length(usable) < length(target_rank)) {
  message("Not in the expression profile, dropped from this target: ",
          paste(setdiff(target_rank, usable), collapse = ", "), ".")
}

if (length(usable) == 0) {
  skip_target("No gene of '", target[1], "' is in the expression profile; ",
              "skipping this target.")
}

# The gene handed to CreateScRank is only the one it validates: obj@para$target
# is reset below, before anything is ranked. It still has to be a gene the
# function accepts, though, and rownames() is not quite the same test -- a gene
# can sit in the matrix with no counts anywhere in it -- so the candidates are
# tried in turn and the ones it refuses are dropped alongside the ones that
# were never there.
#
# Only that one error is caught. Anything else out of CreateScRank -- a version
# skew between Seurat and SeuratObject in the image, a malformed object -- is
# re-raised to kill the task, because it is not a fact about this target, and
# swallowing it would hand Nextflow an empty table to cache as a success. That
# is the failure the exit status at the bottom of this script exists to avoid,
# and a blanket tryCatch here would reintroduce it through the back door.
is_absent_target <- function(e) {
  msg <- conditionMessage(e)
  grepl("target gene", msg, fixed = TRUE) &&
    grepl("expression profile", msg, fixed = TRUE)
}

obj <- NULL
refused <- character()

for (candidate in usable) {
  obj <- tryCatch(
    CreateScRank(input = seuratObj,
                 species = species,
                 cell_type = column,
                 target = candidate),
    error = function(e) {
      if (!is_absent_target(e)) {
        stop(e)
      }
      message("CreateScRank does not find ", candidate,
              " in the expression profile.")
      NULL
    }
  )

  if (!is.null(obj)) {
    break
  }

  refused <- c(refused, candidate)
}

usable <- setdiff(usable, refused)

if (is.null(obj)) {
  skip_target("CreateScRank accepted no gene of '", target[1], "'; ",
              "skipping this target.")
}

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

  for (gene in usable) {
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
  empty_connections
}

write_table(top_connections, "top_connections")

# rank_celltype zeroes the target's row in every cell type network
# (dpGRN[target, ] <- 0), so a gene missing from any one of them stops it with
# "Drug target gene is not in the network". Dropping such a gene, as this used
# to, loses the target from the score table whenever no gene is left -- even
# for the cell types whose networks do carry it. A gene missing from a network
# has no edges there, which is what knocking it out would leave anyway, so it
# is added to that network as an isolated gene (a zero row and column) and the
# target is ranked in every cell type. The same reasoning is why a gene the
# expression profile lacks is dropped above without changing the score.
present_in <- lapply(obj@net, function(net) intersect(usable, rownames(net)))

# A gene in no network at all adds nothing anywhere, so it leaves the target
# rather than being padded into every network.
in_any_net <- usable %in% unlist(present_in)
if (!all(in_any_net)) {
  message("Not in any cell type network, dropped from ranking: ",
          paste(usable[!in_any_net], collapse = ", "), ".")
}
usable <- usable[in_any_net]

# A cell type whose network has none of the target's genes cannot be perturbed
# by it: its score is 0. When that is every cell type, there is nothing for
# rank_celltype to do, and the target still gets its rows.
no_target_gene <- vapply(present_in, function(g) length(g) == 0, logical(1))

zero_ranks <- function(cts) {
  data.frame(cell_type = cts, target = target[1], binding = binding,
             perb_score = rep(0, length(cts)), stringsAsFactors = FALSE)
}

if (length(usable) == 0) {
  message("No gene of '", target[1], "' is in any cell type network; ",
          "writing a perturbation score of 0 for every cell type.")
  write_table(zero_ranks(cell_types), "perbscore_all_targets")
  quit(save = "no", status = 0)
}

pad_net <- function(net, genes) {
  add <- setdiff(genes, rownames(net))
  if (length(add) == 0) {
    return(net)
  }
  m <- as.matrix(net)
  m <- rbind(m, matrix(0, length(add), ncol(m), dimnames = list(add, colnames(m))))
  add_cols <- setdiff(genes, colnames(m))
  m <- cbind(m, matrix(0, nrow(m), length(add_cols), dimnames = list(rownames(m), add_cols)))
  if (inherits(net, "Matrix")) Matrix::Matrix(m, sparse = TRUE) else m
}

for (ct in cell_types) {
  missing <- setdiff(usable, rownames(obj@net[[ct]]))
  if (length(missing) > 0) {
    message("Added to the ", ct, " network with no edges: ",
            paste(missing, collapse = ", "), ".")
    obj@net[[ct]] <- pad_net(obj@net[[ct]], usable)
  }
}

# Empty but named for the same reason empty_ranks is: a target that fails for
# every cell type still has to leave a properly headered table behind.
all_ranks <- empty_ranks

n_failed <- 0

# rank_celltype forks one worker per cell type (foreach/doParallel, i.e.
# mclapply), and each worker builds dense 2n x 2n Laplacians for its networks.
# When the kernel kills workers for memory, mclapply only warns ("scheduled
# cores ... did not deliver results") and hands back NULL for their cell types,
# and rank_celltype then dies further on with the unhelpful
# "seq_len(n): argument must be coercible to non-negative integer". A dead
# worker says nothing about the target, so the ranking is rerun one cell type
# at a time, which keeps a single Laplacian in memory at once. Any other error
# is left alone: it would fail the same way sequentially.
rank_celltype_safe <- function(obj, n.core = 4) {
  workers_died <- FALSE
  res <- tryCatch(
    withCallingHandlers(
      rank_celltype(obj, n.core = n.core),
      warning = function(w) {
        if (grepl("did not deliver", conditionMessage(w), fixed = TRUE)) {
          workers_died <<- TRUE
        }
      }
    ),
    error = function(e) {
      if (!workers_died) {
        stop(e)
      }
      NULL
    }
  )
  if (!is.null(res)) {
    return(res)
  }
  message("rank_celltype workers were killed (likely out of memory); ",
          "retrying with n.core = 1.")
  doParallel::stopImplicitCluster()
  foreach::registerDoSEQ()
  rank_celltype(obj, n.core = 1)
}

for (target_sc in target) {
  # usable rather than the label split on ';': a gene the expression profile
  # does not carry was dropped above, and handing it to rank_celltype anyway
  # is the failure this guard exists to avoid. The label is unchanged, since
  # an unexpressed gene has no edges and so no effect on the score.
  message("Target ", target_sc, " found. Proceeding with rank_celltype on: ",
          paste(usable, collapse = ", "))
  # Set the target
  obj@para$target <- usable

  # Try running rank_celltype
  tryCatch({
    obj <- rank_celltype_safe(obj)

    # Extract data and convert to long format
    perb_scores <- obj@cell_type_rank$perb_score
    df_long <- data.frame(
      cell_type = cell_types,
      target = target_sc,
      binding = binding,
      perb_score = as.numeric(perb_scores)
    )
    # Set outright rather than read off rank_celltype, which has only an
    # isolated padded gene to perturb in these cell types.
    df_long$perb_score[no_target_gene[df_long$cell_type]] <- 0

    # Append to the main data frame
    all_ranks <- rbind(all_ranks, df_long)

    message("Finished: ", target_sc)

  }, error = function(e) {
    message("Failed for target ", target_sc, ": ", e$message)
    n_failed <<- n_failed + 1
  })
}

# Save results (even if partial)
write_table(all_ranks, "perbscore_all_targets")

# rank_celltype's failures are usually a worker in its own parallel backend
# (mclapply) getting killed rather than a real absence of signal for that
# target. Exiting 0 here would let Nextflow cache an empty result as a
# success, silently dropping the target from every downstream table on
# every future -resume. Every target in this invocation failing is worth
# retrying instead of caching, so fail loudly; a partial result (at least
# one target succeeded) is still saved and still exits clean.
if (n_failed == length(target) && n_failed > 0) {
  quit(save = "no", status = 1)
}
