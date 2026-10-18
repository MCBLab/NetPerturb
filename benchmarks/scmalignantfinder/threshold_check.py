#!/usr/bin/env python3
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, roc_curve

base = "/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder"

pred = pd.read_csv(f"{base}/vangalen_scmalignantfinder_predictions.csv", index_col="barcode")
ref = pd.read_csv(f"{base}/vangalen_prediction_refined.csv", index_col="barcode")

df = ref.join(pred, how="inner")
df["ref_label"] = df["PredictionRefined"].str.lower()
core = df[df["ref_label"].isin(["malignant", "normal"])].copy()
y_true = (core["ref_label"] == "malignant").astype(int)
y_score = core["malignancy_probability"]

auc = roc_auc_score(y_true, y_score)
print(f"n cells: {len(core)}")
print(f"AUC-ROC (malignancy_probability vs vanGalen label): {auc:.4f}")
print(f"probability summary by true label:")
print(core.groupby("ref_label")["malignancy_probability"].describe())

fpr, tpr, thr = roc_curve(y_true, y_score)
j = tpr - fpr
best_i = np.argmax(j)
print(f"\nBest Youden's-J threshold: {thr[best_i]:.4f}  (sensitivity={tpr[best_i]:.4f}, specificity={1-fpr[best_i]:.4f})")

for t in [0.5, thr[best_i], 0.3, 0.2, 0.1]:
    pred_bin = (y_score >= t).astype(int)
    acc = (pred_bin == y_true).mean()
    sens = (pred_bin[y_true == 1] == 1).mean()
    spec = (pred_bin[y_true == 0] == 0).mean()
    print(f"threshold={t:.4f}  accuracy={acc:.4%}  sensitivity(malignant)={sens:.4%}  specificity(normal)={spec:.4%}")
