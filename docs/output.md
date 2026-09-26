# Output

Everything below is published under `--outdir` (default `results/`). Which folders appear depends on the tracks you ran:

```text
results/
├── downsample/          always: QC figure and tables
├── genie3/ | scrank/ | hdwgcna/
│                        --network: one network per identity
├── rank_scores/         --network: perturbation scores and top connections
├── sctknk/              --sctknk: knockout table, GSEA, optional plots
│   └── plots/           --sctknk_plot only
├── report/              --network: the HTML report
└── pipeline_info/       always: Nextflow execution report, timeline, trace, DAG
```

Intermediate files, such as the per-identity objects, per-target tables and the wild-type knockout networks, stay in Nextflow's `work/` directory. They are reused by `-resume` and can be deleted once a run is finished.

## Report

`report/netperturb_report.html` is a single self-contained file (a few MB) that opens in any browser, with no server or internet connection needed. It is rendered whenever `--network` runs and has these sections:

| Section | What it shows |
|---|---|
| **Overview** | The `--network` method used, and **Target QC**: target genes removed before analysis, with the reason. |
| **Perturbation score table** | Every identity × target score, searchable and filterable. Never cut. |
| **Perturbation score distribution** | Pooled log10 scores, with the `--score_quantile` cut drawn over them so you can check it falls between the two modes. |
| **Perturbation score heatmap** | Identities × targets, above the cut. |
| **Knockout networks** | `--sctknk`: one interactive view per knocked-out target and identity, with the target in the centre and up to `--sctknk_top_genes` genes it moved around it (FDR < 0.05). |
| **Knockout differentially-regulated genes** | `--sctknk`: every gene clearing FDR < 0.05, queryable. |
| **Knockout gene set enrichment** | `--gsea_gmt`: a dot plot of the strongest sets per combination (`NES`, leading-edge size, significance) over a table of every set with FDR < 0.25. |
| **Cell identity overview** | The UMAP, **cells per identity** (input vs. kept, dropped identities flagged, genes each network was built on) and **target expression per identity** (mean log-normalised expression heatmap). |
| **Pipeline run** | Command line, versions, settings, elapsed time and **running time per step** from the execution trace. |

A section whose track was not run shows a short "not run" note instead.

!!! note
    With `--sctknk` alone (no `--network`) no report is rendered. The knockout and enrichment tables are still published in `sctknk/`.

## downsample/

| File | Contents |
|---|---|
| `umap_<column>.png` | UMAP of the retained cells coloured by `--column`. An embedding stored on the object is reused. |
| `targets_qc.txt` | The targets that passed QC, in `--target` format. Every later step reads this file. |
| `target_qc.tsv` | Genes removed by QC: one row per gene, with its target, the reason and the action taken. Header only when everything passed. |
| `cell_counts.tsv` | Cells and genes per identity: input size, cells kept, status (`kept`/dropped), and the gene counts behind each network. |
| `target_expression.tsv` | Mean log-normalised expression and percentage of expressing cells for each target gene in each kept identity. |

```text
# target_qc.tsv
target      gene   reason                               action
Hhip        Hhip   zero counts in every retained cell   target not analysed
Stfa1;Mpo   Mpo    absent from the expression profile   dropped from the target; analysed as Stfa1
```

In `cell_counts.tsv`, `genes_gene4use` is the gene set offered to the scoring-track networks and `genes_expressed` is how many of those actually have counts in that identity. `genes_filtered` is the difference. `genes_sctknk` is the size of that identity's knockout network (see `--sctknk_min_pct`). `genes_total` counts the object's genes after `DOWNSAMPLE` removed the mitochondrial and ribosomal protein genes on load, and `genes_mt_rb` is how many it removed.

```text
# cell_counts.tsv
identity   n_input  n_used  status  genes_total  genes_gene4use  genes_expressed  genes_filtered  genes_sctknk
sensitive  812      500     kept    18000        3120            2874             246             6412
resistant  2304     500     kept    18000        3120            2951             169             6980
```

## genie3/, scrank/, hdwgcna/

One `<identity>_weight_<METHOD>_<n>.rds` per identity: the scRank object carrying that identity's inferred regulatory network. Only the folder for the chosen `--network` is created. hdWGCNA skips identities too small to aggregate into metacells, and those files are simply absent.

## rank_scores/

### `perbscore_all_targets.txt`

The main result of the scoring track: one row per identity × target.

```text
cell_type  target     binding     perb_score
sensitive  Stfa1;Mpo  antagonist  1.38837233985176e-06
resistant  Stfa1;Mpo  antagonist  3.86382149842044e-06
sensitive  Brd4       antagonist  1.5395821350262e-06
resistant  Brd4       antagonist  1.61320913209275e-06
```

`perb_score` measures how much the identity's network depends on the target. Higher means a stronger dependency. Compare scores within a run, or between runs that used the same `--network` method. An identity whose network lacks every gene of the target scores `0`.

### `top_connections_all_targets.txt`

Each target's strongest edges in each identity's network, up to `--top_connections` per target gene, ranked by absolute weight. The sign is kept: positive edges are activating and negative edges are repressive.

```text
cell_type  target  binding     target_gene  partner  weight   rank
sensitive  Brd4    antagonist  Brd4         Myc      0.8123   1
sensitive  Brd4    antagonist  Brd4         Ccnd1    -0.6640  2
resistant  Brd4    antagonist  Brd4         Myc      0.4471   1
```

## sctknk/

### `sctenifoldknk_all_targets.txt`

One row per identity × knocked-out target × gene, sorted by p-value within each combination. A combined target such as `Stfa1;Mpo` is one joint knockout with its own rows.

```text
cell_type  target  gene      distance  Z     FC    p.value   p.adj
resistant  Brd4    G17_Brd4  14.43217  7.09  41.1  1.41e-10  9.83e-08
```

| Column | Meaning |
|---|---|
| `distance` | How far the gene moved between the wild-type and knocked-out networks after manifold alignment. |
| `Z` | Box-Cox transformed, standardised distance. |
| `FC` | Squared distance over the mean squared distance of the genes not knocked out. This is a ratio to the average gene, not an expression fold change, and it **has no direction**. |
| `p.value`, `p.adj` | Chi-square test on `FC` and its Benjamini-Hochberg FDR. |

The knocked-out gene or genes are not in their own rows. They always move furthest by construction, so they say nothing about the effect of the knockout.

### `gsea_all_targets.txt`

With `--gsea_gmt`. One row per identity × knocked-out target × gene set, with fgsea's columns:

```text
cell_type  target  pathway                    pval      padj      log2err  ES     NES    size  leadingEdge
resistant  Brd4    GOBP_CHROMATIN_REMODELING  1e-50     4e-50     NA       1      3.126  30    G17_Brd4;G32_Brd4
sensitive  Brd4    GOBP_CHROMATIN_REMODELING  1.67e-11  3.33e-11  0.863    0.773  3.043  30    G5_Brd4;G1_Brd4
```

Genes are ranked by `log2(FC)`, so `NES` is signed. **Positive** means the set sits among the genes the knockout moved furthest. **Negative** means it sits among genes that stayed still while the rest of the network moved. Neither means up- or down-regulation. `size` is counted after intersecting the set with the knockout network's genes, so it is smaller than the set's nominal size. A `pval` of `1e-50` with `log2err` `NA` is fgsea's floor: the true value is even smaller.

### `plots/`

With `--sctknk_plot`: `<identity>_sctknk_plot_<target>.pdf`, scTenifoldKnk's `plotKO()` drawing of the knocked-out gene(s), the genes it moved (`p.adj < 0.05`) and the wild-type edges among them. Red edges are positive and blue edges are negative.

## pipeline_info/

Nextflow's standard run records, time-stamped per run: `execution_report_*.html` (resources per task), `execution_timeline_*.html`, `execution_trace_*.txt` and `pipeline_dag_*.html`. The report's *Pipeline run* section summarises the trace.
