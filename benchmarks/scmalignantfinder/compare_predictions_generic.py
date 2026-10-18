#!/usr/bin/env python3
import sys
import pandas as pd
from sklearn.metrics import roc_auc_score

base = "/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder"
pred_file = sys.argv[1]

pred = pd.read_csv(pred_file, index_col="barcode")
ref = pd.read_csv(f"{base}/vangalen_prediction_refined.csv", index_col="barcode")

df = ref.join(pred, how="inner")
df["ref_label"] = df["PredictionRefined"].str.lower()
df["pred_label"] = df["scMalignantFinder_prediction"].str.lower()

core = df[df["ref_label"].isin(["malignant", "normal"])]
agree = (core["ref_label"] == core["pred_label"]).mean()
auc = roc_auc_score((core["ref_label"] == "malignant").astype(int), core["malignancy_probability"])

print(f"file: {pred_file}")
print(f"cells with clear vanGalen label: {len(core)} of {len(df)}")
print(f"Agreement (malignant/normal only): {agree:.4%}")
print(f"AUC-ROC: {auc:.4f}")
print("\nCrosstab:")
print(pd.crosstab(core["ref_label"], core["pred_label"]))
