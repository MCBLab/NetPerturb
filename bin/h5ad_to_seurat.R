#!/usr/bin/env Rscript
# Converts an AnnData .h5ad (scanpy, CELLxGENE, ...) into the Seurat object the
# rest of the pipeline reads, so --obj can point at either. The file is read
# directly with hdf5r -- no Python, anndata or reticulate -- and written with
# the same Seurat and SeuratObject that DOWNSAMPLE and RANK_SCORE load it with.
#
# What is carried over:
#   counts     -> the RNA assay's counts. Chosen by --h5ad_counts (see
#                 pick_counts() below); always raw, integer counts by default.
#   data       -> the RNA assay's data: the object's own log-normalised matrix
#                 (X, or raw.X) when there is one covering the same genes, so a
#                 normalisation done in scanpy is kept, as it is for a Seurat
#                 object; LogNormalize over the counts otherwise.
#   obs        -> meta.data, categoricals as factors, so --column and --batch
#                 name obs columns.
#   var index  -> feature names, or the symbol column --h5ad_gene_column names.
#   obsm       -> reductions: X_umap, X_tsne and X_pca as umap, tsne and pca
#                 (which DOWNSAMPLE's figure and hdWGCNA pick up), any other
#                 cells x dimensions matrix under its own name.
# Everything else (uns, obsp, varm, ...) is left behind; no step reads it.
#
# Both the current on-disk format (anndata >= 0.8, every element tagged with
# an encoding-type) and the older one (anndata 0.7: categoricals as codes plus
# obs/__categories, sparse matrices tagged h5sparse_format) are read.

suppressPackageStartupMessages({
  library(hdf5r)
  library(Matrix)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)

h5ad_file <- args[1]
out_file  <- args[2]

# Missing, empty or "null" (what Nextflow passes for an unset parameter) means
# the default.
arg_or <- function(i, default) {
  if (length(args) >= i && nzchar(args[i]) && args[i] != "null") args[i] else default
}
# --h5ad_counts: "auto", "X", "raw" (adata.raw.X) or the name of a layer.
counts_arg  <- arg_or(3, "auto")
# --h5ad_gene_column: var column holding the gene symbols; NA = the var index,
# or a recognised symbol column when the index is Ensembl IDs.
gene_column <- arg_or(4, NA_character_)

# Layers searched, in this order, by --h5ad_counts auto before raw.X and X.
# These are the names scanpy's tutorials, scvi-tools and most atlases use for a
# copy of the raw counts taken before normalisation.
count_layers <- c("counts", "raw_counts", "count", "umi_counts")

# Symbol columns looked for, in this order, when the var index is Ensembl IDs
# and --h5ad_gene_column is not given. feature_name is CELLxGENE's.
symbol_columns <- c("feature_name", "gene_symbols", "gene_symbol", "gene_name",
                    "gene_names", "symbol", "symbols", "gene")

# A log-normalised matrix is kept as the data layer only below this maximum:
# log1p of counts scaled to 1e4 (or even 1e6) per cell stays far under it,
# while counts normalised without the log reach it in nearly every cell.
max_log_value <- 50

h5 <- H5File$new(h5ad_file, mode = "r")
on.exit(h5$close_all(), add = TRUE)

# Generic helpers ----------------------------------------------------------

# Whether a path exists, checked one level at a time: H5Lexists fails, rather
# than answering FALSE, when an intermediate group is missing.
has <- function(node, path) {
  for (part in strsplit(path, "/", fixed = TRUE)[[1]]) {
    if (!inherits(node, c("H5File", "H5Group")) || !node$exists(part)) return(FALSE)
    node <- node[[part]]
  }
  TRUE
}

attr_or <- function(node, name, default = NULL) {
  if (!name %in% h5attr_names(node)) return(default)
  # An empty array attribute (an obs with no columns writes column-order as
  # one) makes H5Aread fail on its null buffer rather than return nothing.
  tryCatch(h5attr(node, name), error = function(e) default)
}

encoding <- function(node) {
  enc <- attr_or(node, "encoding-type")
  if (!is.null(enc)) return(as.character(enc))
  old <- attr_or(node, "h5sparse_format")
  if (!is.null(old)) return(paste0(old, "_matrix"))
  if (inherits(node, "H5D")) "array" else ""
}

# A 1-d element: a plain or string array, or one of the groups anndata writes
# for categoricals and nullable columns. Returns an R vector.
read_vector <- function(node) {
  enc <- encoding(node)
  if (inherits(node, "H5D")) {
    if (length(node$dims) != 1) stop("not one-dimensional")
    return(node$read())
  }
  if (enc == "categorical") {
    cats  <- as.character(read_vector(node[["categories"]]))
    codes <- as.integer(node[["codes"]]$read()) + 1L
    codes[codes < 1L] <- NA_integer_  # -1 is a missing value
    return(factor(cats[codes], levels = cats,
                  ordered = isTRUE(as.logical(attr_or(node, "ordered", FALSE)))))
  }
  if (enc %in% c("nullable-integer", "nullable-boolean", "nullable-string-array")) {
    values <- node[["values"]]$read()
    if (node$exists("mask")) values[as.logical(node[["mask"]]$read())] <- NA
    return(values)
  }
  stop("unsupported encoding '", enc, "'")
}

# An obs/var data frame: its index and every column that reads as a vector of
# the right length. A column that does not (a nested or 2-d element) is
# skipped with a message rather than failing the conversion over metadata
# nothing may ever read.
read_dataframe <- function(group) {
  index_name <- as.character(attr_or(group, "_index", "_index"))
  index <- as.character(read_vector(group[[index_name]]))
  columns <- as.character(attr_or(group, "column-order", character()))
  if (length(columns) == 0) {
    columns <- setdiff(names(group), c(index_name, "__categories"))
  }
  df <- data.frame(row.names = seq_along(index))
  for (col in columns) {
    value <- tryCatch({
      node <- group[[col]]
      # anndata 0.7: integer codes, with the categories stored apart
      if (inherits(node, "H5D") && "categories" %in% h5attr_names(node) &&
          has(group, paste0("__categories/", col))) {
        cats  <- as.character(group[["__categories"]][[col]]$read())
        codes <- as.integer(node$read()) + 1L
        codes[codes < 1L] <- NA_integer_
        factor(cats[codes], levels = cats)
      } else {
        read_vector(node)
      }
    }, error = function(e) {
      message("Skipping column '", col, "' of ", group$get_obj_name(), ": ",
              conditionMessage(e), ".")
      NULL
    })
    if (is.null(value)) next
    if (length(value) != length(index)) {
      message("Skipping column '", col, "' of ", group$get_obj_name(),
              ": ", length(value), " values for ", length(index), " rows.")
      next
    }
    df[[col]] <- value
  }
  list(index = index, df = df)
}

# A cells x genes matrix as the genes x cells dgCMatrix Seurat stores.
#   CSR over cells is already CSC over genes, so it is taken as it is.
#   CSC over cells is built cells x genes and transposed.
#   Dense: hdf5r reports dimensions in R's order, i.e. reversed from anndata's,
#   so a dense X already reads as genes x cells. It is read a block of cells at
#   a time and made sparse per block, so a large dense X is never held twice.
# sparseMatrix() rather than new("dgCMatrix"): it sorts row indices within each
# column and sums duplicates, neither of which scipy guarantees.
read_matrix <- function(node) {
  enc <- encoding(node)
  if (enc %in% c("csr_matrix", "csc_matrix")) {
    shape <- as.integer(attr_or(node, "shape", attr_or(node, "h5sparse_shape")))
    x <- as.numeric(node[["data"]]$read())
    i <- as.integer(node[["indices"]]$read())
    p <- as.numeric(node[["indptr"]]$read())
    if (enc == "csr_matrix") {
      return(sparseMatrix(i = i, p = p, x = x, dims = rev(shape), index1 = FALSE,
                          repr = "C"))
    }
    return(t(sparseMatrix(i = i, p = p, x = x, dims = shape, index1 = FALSE,
                          repr = "C")))
  }
  if (inherits(node, "H5D") && length(node$dims) == 2) {
    dims <- node$dims  # genes x cells
    step <- max(1L, floor(5e7 / max(dims[1], 1)))
    starts <- seq(1L, dims[2], by = step)
    blocks <- lapply(starts, function(s) {
      cells <- s:min(s + step - 1L, dims[2])
      as(Matrix(node[, cells, drop = FALSE], sparse = TRUE), "CsparseMatrix")
    })
    return(as(do.call(cbind, blocks), "dgCMatrix"))
  }
  stop("unsupported matrix encoding '", enc, "'")
}

is_counts <- function(m) {
  length(m@x) > 0 && all(m@x >= 0) && all(m@x == round(m@x))
}

# Matrices the object carries, by the name --h5ad_counts uses for them, with
# where each one's genes are described. Read lazily and at most once.
matrix_sources <- list()
if (has(h5, "X")) {
  matrix_sources[["X"]] <- list(path = "X", var = "var")
}
if (has(h5, "raw/X")) {
  matrix_sources[["raw"]] <- list(path = "raw/X", var = "raw/var")
}
if (has(h5, "layers")) {
  for (layer in names(h5[["layers"]])) {
    if (!is.null(matrix_sources[[layer]])) {
      message("Layer '", layer, "' shares its name with adata.", layer,
              "; --h5ad_counts ", layer, " means adata.", layer, ".")
      next
    }
    matrix_sources[[layer]] <- list(path = paste0("layers/", layer), var = "var")
  }
}
describe <- function(name) {
  if (name %in% c("X", "raw")) c(X = "X", raw = "raw.X")[[name]] else paste0("layers['", name, "']")
}

matrix_cache <- list()
get_matrix <- function(name) {
  if (is.null(matrix_cache[[name]])) {
    matrix_cache[[name]] <<- read_matrix(h5[[matrix_sources[[name]]$path]])
  }
  matrix_cache[[name]]
}

var_cache <- list()
get_var <- function(name) {
  path <- matrix_sources[[name]]$var
  if (is.null(var_cache[[path]])) var_cache[[path]] <<- read_dataframe(h5[[path]])
  var_cache[[path]]
}

# Cells --------------------------------------------------------------------

obs <- read_dataframe(h5[["obs"]])
n_obs <- length(obs$index)
cells <- obs$index
if (anyDuplicated(cells)) {
  message(sum(duplicated(cells)), " duplicated cell name(s) made unique.")
  cells <- make.unique(cells)
}
rownames(obs$df) <- cells

# Counts -------------------------------------------------------------------
# auto takes the first of the usual count layers, raw.X and X that holds
# non-negative whole numbers, since scRank, scTenifoldKnk and the report's
# QC all read raw counts. A source named explicitly is taken as it is, with a
# warning if it does not look like counts.
pick_counts <- function(choice) {
  if (choice == "auto") {
    candidates <- c(intersect(count_layers, names(matrix_sources)),
                    intersect(c("raw", "X"), names(matrix_sources)))
    for (name in candidates) {
      if (is_counts(get_matrix(name))) {
        message("Counts: ", describe(name), " (--h5ad_counts auto).")
        return(name)
      }
      message(describe(name), " is not raw counts (not non-negative whole numbers).")
    }
    stop("No raw counts found: checked ",
         if (length(candidates)) paste(vapply(candidates, describe, ""), collapse = ", ") else "nothing",
         ". Layers in this file: ",
         if (has(h5, "layers") && length(names(h5[["layers"]]))) paste(names(h5[["layers"]]), collapse = ", ") else "none",
         ". The networks need raw counts; pass --h5ad_counts naming the layer that holds them.")
  }
  if (is.null(matrix_sources[[choice]])) {
    stop("--h5ad_counts '", choice, "' is not in this file. Available: ",
         paste(c("auto", names(matrix_sources)), collapse = ", "),
         " ('raw' is adata.raw.X; anything else is a layer).")
  }
  if (!is_counts(get_matrix(choice))) {
    warning(describe(choice), " (--h5ad_counts ", choice, ") does not hold ",
            "non-negative whole numbers; it is used as counts regardless.",
            call. = FALSE, immediate. = TRUE)
  }
  message("Counts: ", describe(choice), " (--h5ad_counts ", choice, ").")
  choice
}

counts_name <- pick_counts(counts_arg)

# X and the layers share one var, but raw.X keeps every gene X had before it
# was subset (to highly variable genes, typically), so a layer of counts taken
# after that subset covers only those genes -- the classic scanpy tutorial
# does exactly this. A gene outside the counts reaches no network, and a
# target outside them is dropped by target QC, so auto moves to raw.X when it
# holds counts for more genes, and otherwise the gap is warned about.
if (counts_name != "raw" && !is.null(matrix_sources[["raw"]])) {
  n_counts <- length(get_var(counts_name)$index)
  n_raw <- length(get_var("raw")$index)
  if (n_raw > n_counts) {
    if (counts_arg == "auto" && is_counts(get_matrix("raw"))) {
      message("raw.X also holds counts, for ", n_raw, " genes against ", n_counts,
              " in ", describe(counts_name), "; using raw.X.")
      counts_name <- "raw"
    } else {
      warning(describe(counts_name), " covers ", n_counts, " genes, but raw.X ",
              "has ", n_raw, ": the counts were likely copied after X was ",
              "subset (e.g. to highly variable genes). Only these ", n_counts,
              " genes reach the networks, and a target outside them is dropped ",
              "by target QC. raw.X is not raw counts, so it cannot stand in; ",
              "save the counts before subsetting to keep every gene.",
              call. = FALSE, immediate. = TRUE)
    }
  }
}

counts <- get_matrix(counts_name)
counts_var <- get_var(counts_name)

if (ncol(counts) != n_obs || nrow(counts) != length(counts_var$index)) {
  stop(describe(counts_name), " is ", ncol(counts), " cells x ", nrow(counts),
       " genes, but obs has ", n_obs, " cells and its var ",
       length(counts_var$index), " genes.")
}

# Gene names ---------------------------------------------------------------
# The var index, unless --h5ad_gene_column names a symbol column, or the index
# is Ensembl IDs and a usual symbol column is there to use instead: targets,
# transcription factors and drug targets are all matched by symbol, so a run on
# Ensembl IDs would find none of them. The column is looked up in the counts'
# own var first and then, by index, in the main var, since raw.var often
# carries the index alone.
gene_lookup <- function(column) {
  if (column %in% names(counts_var$df)) {
    return(as.character(counts_var$df[[column]]))
  }
  if (counts_name == "raw") {
    main <- read_dataframe(h5[["var"]])
    if (column %in% names(main$df)) {
      return(as.character(main$df[[column]])[match(counts_var$index, main$index)])
    }
  }
  NULL
}

looks_ensembl <- function(x) mean(grepl("^ENS[A-Z]*G[0-9]{6,}", x)) > 0.5

genes <- counts_var$index
gene_source <- "var index"
if (!is.na(gene_column)) {
  symbols <- gene_lookup(gene_column)
  if (is.null(symbols)) {
    stop("--h5ad_gene_column '", gene_column, "' is not a var column. Available: ",
         paste(names(counts_var$df), collapse = ", "), ".")
  }
  gene_source <- paste0("var['", gene_column, "']")
} else if (looks_ensembl(genes)) {
  symbols <- NULL
  for (column in symbol_columns) {
    candidate <- gene_lookup(column)
    if (!is.null(candidate) && !looks_ensembl(candidate[!is.na(candidate)])) {
      symbols <- candidate
      gene_source <- paste0("var['", column, "']")
      message("The var index is Ensembl IDs; using var['", column,
              "'] as gene symbols (set --h5ad_gene_column to choose another).")
      break
    }
  }
  if (is.null(symbols)) {
    message("The var index looks like Ensembl IDs and no symbol column was found; ",
            "targets, TFs and drug targets are matched by symbol, so set ",
            "--h5ad_gene_column to the var column holding them.")
  }
} else {
  symbols <- NULL
}
if (!is.null(symbols)) {
  missing <- is.na(symbols) | !nzchar(trimws(symbols)) | symbols == "nan"
  if (any(missing)) {
    message(sum(missing), " gene(s) with no symbol keep their var index name.")
  }
  genes <- ifelse(missing, counts_var$index, symbols)
}

# Seurat replaces underscores with dashes on its own; doing it first means two
# names that would only collide afterwards are told apart by make.unique too.
if (any(grepl("_", genes, fixed = TRUE))) {
  message(sum(grepl("_", genes, fixed = TRUE)), " gene name(s) with '_' ",
          "written with '-', as Seurat requires.")
  genes <- gsub("_", "-", genes, fixed = TRUE)
}
if (anyDuplicated(genes)) {
  message(sum(duplicated(genes)), " duplicated gene name(s) made unique ",
          "(the second copy of X becomes X.1).")
  genes <- make.unique(genes)
}

dimnames(counts) <- list(genes, cells)

# Normalised data ----------------------------------------------------------
# The object's own log-normalised values are kept when there are any, as a
# Seurat object's data layer is, since GENIE3 and hdWGCNA read data. X is
# tried before raw.X; either has to be non-integer (else it is just counts),
# non-negative (else it has been scaled), on a log scale, and cover every gene
# of the counts.
data <- NULL
data_source <- "LogNormalize over the counts"
for (name in setdiff(intersect(c("X", "raw"), names(matrix_sources)), counts_name)) {
  m <- tryCatch(get_matrix(name), error = function(e) NULL)
  if (is.null(m) || ncol(m) != n_obs || length(m@x) == 0) next
  why <- if (is_counts(m)) {
    "holds whole numbers"
  } else if (min(m@x) < 0) {
    "has negative values (scaled)"
  } else if (max(m@x) > max_log_value) {
    "is not log-transformed"
  } else if (!all(counts_var$index %in% get_var(name)$index)) {
    "does not cover every gene of the counts"
  } else {
    NA_character_
  }
  if (!is.na(why)) {
    message(describe(name), " not used as normalised data: it ", why, ".")
    next
  }
  data <- m[match(counts_var$index, get_var(name)$index), , drop = FALSE]
  dimnames(data) <- dimnames(counts)
  data_source <- describe(name)
  break
}
message("Normalised data: ", data_source, ".")

# The object ---------------------------------------------------------------
# A v3 Assay, which is what the Seurat v4 in the pipeline's image works with
# (see as_v3_assay() in downsample_and_split.R).
assay <- CreateAssayObject(counts = counts)
if (!is.null(data)) {
  assay <- SetAssayData(assay, layer = "data", new.data = data)
}
obj <- CreateSeuratObject(assay, assay = "RNA", meta.data = obs$df,
                          project = tools::file_path_sans_ext(basename(h5ad_file)))
if (is.null(data)) {
  obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
}

# Reductions ---------------------------------------------------------------
reduction_names <- c(X_umap = "umap", X_tsne = "tsne", X_pca = "pca")
reduction_keys  <- c(umap = "UMAP_", tsne = "tSNE_", pca = "PC_")
if (has(h5, "obsm")) {
  used_keys <- character()
  for (entry in names(h5[["obsm"]])) {
    emb <- tryCatch({
      node <- h5[["obsm"]][[entry]]
      if (!inherits(node, "H5D") || length(node$dims) != 2) stop("not a dense matrix")
      m <- t(node$read())  # hdf5r reads it dimensions x cells
      if (nrow(m) != n_obs || !is.numeric(m)) stop("not cells x dimensions")
      m
    }, error = function(e) {
      message("Skipping obsm['", entry, "']: ", conditionMessage(e), ".")
      NULL
    })
    if (is.null(emb)) next
    name <- if (entry %in% names(reduction_names)) reduction_names[[entry]] else {
      gsub("[^a-z0-9]+", "_", tolower(sub("^X_", "", entry)))
    }
    if (name %in% Reductions(obj) || name %in% Assays(obj)) next
    key <- if (name %in% names(reduction_keys)) reduction_keys[[name]] else {
      paste0(gsub("[^A-Za-z0-9]", "", name), "_")
    }
    if (!grepl("^[A-Za-z]", key)) key <- paste0("DR", key)
    if (key %in% used_keys) key <- paste0(sub("_$", "", key), length(used_keys), "_")
    used_keys <- c(used_keys, key)
    dimnames(emb) <- list(cells, paste0(key, seq_len(ncol(emb))))
    obj[[name]] <- CreateDimReducObject(embeddings = emb, key = key, assay = "RNA")
  }
}

obj@misc$h5ad <- list(file = basename(h5ad_file), counts = describe(counts_name),
                      data = data_source, genes = gene_source)

message("Converted ", basename(h5ad_file), ": ", ncol(obj), " cells x ",
        nrow(obj), " genes; gene names from ", gene_source, ".")
message("Metadata columns (for --column / --batch): ",
        paste(colnames(obj@meta.data), collapse = ", "), ".")
message("Reductions: ",
        if (length(Reductions(obj))) paste(Reductions(obj), collapse = ", ") else "none", ".")

saveRDS(obj, out_file)
