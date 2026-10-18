#!/usr/bin/env python3
import sys
import scipy.io
import anndata as ad
import pandas as pd

from scMalignantFinder import classifier

base = "/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder"
model_method = sys.argv[1] if len(sys.argv) > 1 else "LogisticRegression"
out_suffix = sys.argv[2] if len(sys.argv) > 2 else model_method.lower()

print("loading vanGalen test matrix...", flush=True)
counts = scipy.io.mmread(f"{base}/vangalen_counts.mtx").tocsr()
genes = [l.strip() for l in open(f"{base}/vangalen_genes.tsv")]
barcodes = [l.strip() for l in open(f"{base}/vangalen_barcodes.tsv")]
test_adata = ad.AnnData(
    X=counts.T.tocsr(),
    obs=pd.DataFrame(index=barcodes),
    var=pd.DataFrame(index=genes),
)
print("test_adata:", test_adata.shape, flush=True)

model = classifier.scMalignantFinder(
    test_input=test_adata,
    pretrain_dir=None,
    train_h5ad_path=f"{base}/combine_training_aml.h5ad",
    feature_path=f"{base}/combined_tumor_up_down_degs.txt",
    model_method=model_method,
    norm_type=True,
    n_thread=4,
)
model.load()
result_adata = model.predict()

out = result_adata.obs[["scMalignantFinder_prediction", "malignancy_probability"]].copy()
out.index.name = "barcode"
out.to_csv(f"{base}/vangalen_scmalignantfinder_predictions_retrained_{out_suffix}.csv")
print(f"saved predictions ({model_method}):", out.shape, flush=True)
print(out["scMalignantFinder_prediction"].value_counts(), flush=True)
