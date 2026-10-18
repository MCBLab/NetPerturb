#!/usr/bin/env python3
import scipy.io
import anndata as ad
import pandas as pd

base = "/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder"

parts = []
for study in ["Naldini_2023", "BeneytoCalabuig_2023", "Lasry_2022"]:
    counts = scipy.io.mmread(f"{base}/train_{study}_counts.mtx").tocsr()  # genes x cells
    genes = [l.strip() for l in open(f"{base}/train_{study}_genes.tsv")]
    barcodes = [l.strip() for l in open(f"{base}/train_{study}_barcodes.tsv")]
    labels = pd.read_csv(f"{base}/train_{study}_labels.csv", index_col="barcode")
    labels = labels.loc[barcodes]

    a = ad.AnnData(
        X=counts.T.tocsr(),
        obs=pd.DataFrame({"Raw_annotation": labels["Raw_annotation"].values, "study": study}, index=barcodes),
        var=pd.DataFrame(index=genes),
    )
    a.var_names_make_unique()
    print(study, a.shape, a.obs["Raw_annotation"].value_counts().to_dict(), flush=True)
    parts.append(a)

train_adata = ad.concat(parts, join="inner")  # intersect gene sets across the two cohorts
train_adata.obs_names_make_unique()
print("combined train_adata:", train_adata.shape, flush=True)
print(train_adata.obs["Raw_annotation"].value_counts(), flush=True)
print(train_adata.obs.groupby(["study", "Raw_annotation"]).size(), flush=True)

train_adata.write_h5ad(f"{base}/combine_training_aml.h5ad")
print("saved combine_training_aml.h5ad", flush=True)
