#!/bin/bash
#SBATCH --partition=intel-128
#SBATCH --account=dmcoelho
#SBATCH --qos=preempt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --job-name=scmalignantfinder_retrain_amldegs
#SBATCH --chdir=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder
#SBATCH --output=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder/retrain_amldegs_%j.out
#SBATCH --error=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder/retrain_amldegs_%j.err

set -euo pipefail

module load softwares/python/3.10.5-gnu8
source /home/lgdqamorim/scratch/scmalignantfinder/venv/bin/activate
unset PYTHONPATH

echo "=== computing AML-specific DEGs (Malignant vs Normal, per cohort, >=2/3 consistent) ==="
python compute_aml_degs.py

echo "=== training + predicting (LogisticRegression, AML-specific features) ==="
python - <<'PY'
import sys
sys.path.insert(0, ".")
import scipy.io, anndata as ad, pandas as pd
from scMalignantFinder import classifier

base = "."
counts = scipy.io.mmread(f"{base}/vangalen_counts.mtx").tocsr()
genes = [l.strip() for l in open(f"{base}/vangalen_genes.tsv")]
barcodes = [l.strip() for l in open(f"{base}/vangalen_barcodes.tsv")]
test_adata = ad.AnnData(X=counts.T.tocsr(), obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=genes))

model = classifier.scMalignantFinder(
    test_input=test_adata,
    pretrain_dir=None,
    train_h5ad_path=f"{base}/combine_training_aml.h5ad",
    feature_path=f"{base}/aml_specific_degs.txt",
    model_method="LogisticRegression",
    norm_type=True,
    n_thread=4,
)
model.load()
result_adata = model.predict()
out = result_adata.obs[["scMalignantFinder_prediction", "malignancy_probability"]].copy()
out.index.name = "barcode"
out.to_csv(f"{base}/vangalen_scmalignantfinder_predictions_amldegs.csv")
print("saved predictions:", out.shape)
print(out["scMalignantFinder_prediction"].value_counts())
PY

echo "=== comparing ==="
python compare_predictions_generic.py vangalen_scmalignantfinder_predictions_amldegs.csv

echo "=== done ==="
