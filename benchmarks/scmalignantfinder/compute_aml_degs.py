#!/usr/bin/env python3
import scanpy as sc
import pandas as pd
import numpy as np

base = "/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder"
TOP_N = 1500  # per direction, per cohort

adata = sc.read_h5ad(f"{base}/combine_training_aml.h5ad")
print("loaded:", adata.shape, flush=True)

per_study_top = {}
for study in adata.obs["study"].unique():
    sub = adata[adata.obs["study"] == study].copy()
    sc.pp.normalize_total(sub, target_sum=1e4)  # same normalization the classifier itself uses, no log
    sc.tl.rank_genes_groups(sub, groupby="Raw_annotation", groups=["Malignant"], reference="Normal", method="wilcoxon")
    res = sc.get.rank_genes_groups_df(sub, group="Malignant")
    res = res.sort_values("scores", key=lambda s: s.abs(), ascending=False)
    top = res.head(TOP_N * 2)  # both directions mixed in by |score|, generous cut
    up = top[top["logfoldchanges"] > 0]["names"].head(TOP_N).tolist()
    down = top[top["logfoldchanges"] < 0]["names"].head(TOP_N).tolist()
    per_study_top[study] = {"up": set(up), "down": set(down)}
    print(f"{study}: {len(up)} up, {len(down)} down (of {sub.n_obs} cells)", flush=True)

studies = list(per_study_top.keys())
all_genes = set()
for s in studies:
    all_genes |= per_study_top[s]["up"] | per_study_top[s]["down"]

consistent = []
for g in all_genes:
    up_votes = sum(g in per_study_top[s]["up"] for s in studies)
    down_votes = sum(g in per_study_top[s]["down"] for s in studies)
    if up_votes >= 2 or down_votes >= 2:
        consistent.append(g)

print(f"\ngenes consistent (same direction) in >=2 of {len(studies)} cohorts: {len(consistent)}", flush=True)

with open(f"{base}/aml_specific_degs.txt", "w") as f:
    for g in sorted(consistent):
        f.write(g + "\n")
print("saved aml_specific_degs.txt", flush=True)
