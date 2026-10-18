#!/usr/bin/env python3
import sys
import scipy.io
import scipy.sparse
import anndata as ad
import pandas as pd

from scMalignantFinder import classifier

base = "/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder"

print("loading matrix...", flush=True)
counts = scipy.io.mmread(f"{base}/vangalen_counts.mtx").tocsr()  # genes x cells
genes = [l.strip() for l in open(f"{base}/vangalen_genes.tsv")]
barcodes = [l.strip() for l in open(f"{base}/vangalen_barcodes.tsv")]
print("counts shape (genes x cells):", counts.shape, "genes:", len(genes), "barcodes:", len(barcodes), flush=True)

adata = ad.AnnData(
    X=counts.T.tocsr(),  # cells x genes, as scanpy/AnnData convention
    obs=pd.DataFrame(index=barcodes),
    var=pd.DataFrame(index=genes),
)
print("adata:", adata, flush=True)

model = classifier.scMalignantFinder(
    test_input=adata,
    pretrain_dir=f"{base}/pretrained_model",
    norm_type=True,
    n_thread=4,
)
model.load()
result_adata = model.predict()

out = result_adata.obs[["scMalignantFinder_prediction", "malignancy_probability"]].copy()
out.index.name = "barcode"
out.to_csv(f"{base}/vangalen_scmalignantfinder_predictions.csv")
print("saved predictions:", out.shape, flush=True)
print(out["scMalignantFinder_prediction"].value_counts(), flush=True)
