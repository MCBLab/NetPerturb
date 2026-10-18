#!/usr/bin/env python3
import sys
import scipy.io
import anndata as ad
import pandas as pd

from scMalignantFinder import classifier

base = "/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder"
study = sys.argv[1]

print(f"loading {study} test matrix...", flush=True)
counts = scipy.io.mmread(f"{base}/annotate_{study}_counts.mtx").tocsr()
genes = [l.strip() for l in open(f"{base}/annotate_{study}_genes.tsv")]
barcodes = [l.strip() for l in open(f"{base}/annotate_{study}_barcodes.tsv")]
test_adata = ad.AnnData(X=counts.T.tocsr(), obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=genes))
print("test_adata:", test_adata.shape, flush=True)

model = classifier.scMalignantFinder(
    test_input=test_adata,
    pretrain_dir=None,
    train_h5ad_path=f"{base}/combine_training_aml.h5ad",
    feature_path=f"{base}/combined_tumor_up_down_degs.txt",
    model_method="LogisticRegression",
    norm_type=True,
    n_thread=4,
)
model.load()
result_adata = model.predict()

out = result_adata.obs[["scMalignantFinder_prediction", "malignancy_probability"]].copy()
out.index.name = "barcode"
out.to_csv(f"{base}/annotate_{study}_predictions.csv")
print(f"saved predictions ({study}):", out.shape, flush=True)
print(out["scMalignantFinder_prediction"].value_counts(), flush=True)
