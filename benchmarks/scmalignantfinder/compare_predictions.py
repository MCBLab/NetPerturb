#!/usr/bin/env python3
import pandas as pd

base = "/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder"

pred = pd.read_csv(f"{base}/vangalen_scmalignantfinder_predictions.csv", index_col="barcode")
ref = pd.read_csv(f"{base}/vangalen_prediction_refined.csv", index_col="barcode")

df = ref.join(pred, how="inner")
print("merged cells:", len(df), "of ref", len(ref), "and pred", len(pred))

# Normalise labels to the same vocabulary for comparison.
df["ref_label"] = df["PredictionRefined"].str.lower()
df["pred_label"] = df["scMalignantFinder_prediction"].str.lower()

print("\nFull crosstab (including 'unclear'):")
print(pd.crosstab(df["ref_label"], df["pred_label"]))

core = df[df["ref_label"].isin(["malignant", "normal"])]
agree = (core["ref_label"] == core["pred_label"]).mean()
print(f"\nCells with a clear vanGalen label (malignant/normal): {len(core)} of {len(df)}")
print(f"Agreement (malignant/normal only): {agree:.4%}")

print("\nCrosstab, clear-label cells only:")
print(pd.crosstab(core["ref_label"], core["pred_label"]))

unclear = df[df["ref_label"] == "unclear"]
if len(unclear):
    print(f"\n'unclear' cells ({len(unclear)}) scMalignantFinder called:")
    print(unclear["pred_label"].value_counts())

core.to_csv(f"{base}/vangalen_comparison_core.csv")
