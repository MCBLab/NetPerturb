#!/usr/bin/env Rscript
suppressMessages(library(Seurat))

B <- "/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/aml_bonemarrowmap"
OUT <- "/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder"

merge_one <- function(study, rds_name) {
  obj <- readRDS(file.path(B, rds_name))
  pred <- read.csv(file.path(OUT, paste0("annotate_", study, "_predictions.csv")),
                    row.names = "barcode", stringsAsFactors = FALSE)
  pred <- pred[colnames(obj), , drop = FALSE]
  stopifnot(all(rownames(pred) == colnames(obj)))

  obj$scMalignantFinder_prediction <- pred$scMalignantFinder_prediction
  obj$scMalignantFinder_probability <- pred$malignancy_probability

  cat("\n==", study, "==\n")
  cat("original NormalvsMalignant:\n"); print(table(obj$NormalvsMalignant, useNA = "always"))
  cat("scMalignantFinder_prediction:\n"); print(table(obj$scMalignantFinder_prediction, useNA = "always"))
  cat("cross-tab (where original label exists):\n")
  has_label <- !is.na(obj$NormalvsMalignant)
  print(table(obj$NormalvsMalignant[has_label], obj$scMalignantFinder_prediction[has_label]))

  saveRDS(obj, file.path(B, rds_name))
  cat("saved", rds_name, "with 2 new meta.data columns\n")
}

merge_one("Naldini_2023", "filtered_Naldini_2023.seurat.rds")
merge_one("BeneytoCalabuig_2023", "filtered_BeneytoCalabuig_2023.seurat.rds")
merge_one("Lasry_2022", "filtered_Lasry_2022.seurat.rds")
cat("\ndone\n")
