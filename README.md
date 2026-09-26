<img src="docs/img/logo.svg" alt="NetPerturb logo" align="right" width="150">

# NetPerturb 🧬

📖 **Documentation: [mcblab.github.io/NetPerturb](https://mcblab.github.io/NetPerturb/)** (usage, output and FAQ)

**MCBLab/NetPerturb** is a scalable Nextflow pipeline designed to infer Gene Regulatory Networks (GRNs) and calculate single-cell expression ranking perturbation scores using the scRank algorithm. 

Previous GRN tools are often difficult to scale for large Single-Cell RNA-seq (scRNA-seq) datasets. In this context, `NetPerturb` was built to enable high-throughput perturbation scoring in a user-friendly, parallelized, and computationally effective way. The pipeline uses Singularity containers, making installation trivial and results highly reproducible across high-performance computing (HPC) environments.

## Pipeline Overview

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/images/netperturb_metro_dark.png">
  <img alt="NetPerturb metro map" src="docs/images/netperturb_metro_light.png">
</picture>

The three warm lines are the `--network` methods. They are mutually exclusive, so a
run rides exactly one of them from `--obj`/`--target` through to the HTML report:
they share downsampling, scoring and reporting, and diverge only at network
inference. The blue `sctknk` line is separate — `--sctknk` switches it on, and it
runs in parallel with whichever `--network` line was chosen, or on its own.

<details>
<summary>Regenerating this diagram</summary>

The map is defined in [`docs/netperturb_metro.mmd`](docs/netperturb_metro.mmd) and
rendered with [nf-metro](https://github.com/seqeralabs/nf-metro) (`pip install nf-metro`):

```bash
# FS bumps every text size and the label metrics that drive spacing, so the
# layout re-flows rather than just overprinting bigger glyphs.
FS=1.25

# Theme-aware SVG and an interactive pan/zoom page
nf-metro render docs/netperturb_metro.mmd -o docs/images/netperturb_metro.svg --font-scale $FS --embed-font --responsive
nf-metro render docs/netperturb_metro.mmd -o docs/images/netperturb_metro.html --format html --font-scale $FS --animate

# The baked light/dark PNGs used above (needs `pip install cairosvg`;
# --no-chrome-css bakes the colours, since rasterisers cannot resolve var())
for m in light dark; do
  nf-metro render docs/netperturb_metro.mmd -o /tmp/nm_$m.svg --mode $m --font-scale $FS --no-chrome-css --embed-font
  cairosvg /tmp/nm_$m.svg -s 2 -o docs/images/netperturb_metro_$m.png
done

# Print/poster assets: vector SVG and PDF, plus a 4x raster fallback
cp /tmp/nm_light.svg docs/images/netperturb_metro_poster.svg
cairosvg /tmp/nm_light.svg -f pdf -o docs/images/netperturb_metro_poster.pdf
cairosvg /tmp/nm_light.svg -s 4  -o docs/images/netperturb_metro_poster@4x.png
```

For print, use `docs/images/netperturb_metro_poster.pdf` or `.svg` — both are true
vector, so they stay sharp at any poster size. Raise `FS` above if the text still
reads small at your final dimensions; past roughly `1.4` the station and output
captions begin to collide.

Every station carries a `%%metro process:` mapping, so the map can also track a live
run. Serve it and point Nextflow's weblog at it:

```bash
nf-metro serve docs/netperturb_metro.mmd --port 8080
nextflow run main.nf -profile test,singularity -with-weblog http://localhost:8080/events
```

After changing the pipeline's processes, re-check the mappings still line up:

```bash
nextflow run main.nf -profile test -preview -with-dag dag.mmd
nf-metro check-mapping docs/netperturb_metro.mmd --dag dag.mmd
```

</details>

## Pipeline Summary

The workflow executes the following core modules:

### 1. Object Parsing and Downsampling (`DOWNSAMPLE`)
This is the initial step of the process. It ingests a fully processed Seurat object (`.rds`) and identifies the user-defined metadata column containing the cell identities (e.g., cell types or clones). To ensure statistical robustness and equitable GRN inference, it randomly downsamples the cells from each identity to a specified maximum number (`--n_cells`), balancing the computational load. Identities left with fewer than `--min_cells` cells after that are dropped here, before the split, so no network is ever inferred from one and it reaches no score, table or figure. It also writes a UMAP of the retained cells coloured by `--column` to `downsample/umap_<column>.png`, so the identities entering the analysis, and their relative sizes after downsampling, can be checked at a glance. An embedding already present on the object is reused; one is computed only if the object carries none. Drawing this figure is guarded, so a plotting failure leaves a placeholder image rather than stopping the run.

It also checks every gene named in `--target` against the retained cells. A gene absent from the object, or with zero counts in every retained cell, can have no edges in any network, so it is removed from the rest of the analysis here rather than failing a later step. A combined target loses only its failing gene (`Stfa1;Mpo` with `Mpo` absent is analysed as `Stfa1`) and is dropped only when none of it is left; the run aborts if no target passes. The passing targets go to `downsample/targets_qc.txt`, which every later step reads in place of `--target`. The excluded genes, with the reason and what happened to their target, go to `downsample/target_qc.tsv` and are listed in the report's Overview as "These genes were not analyzed due to QC checking".

No network in the pipeline is built on mitochondrial or ribosomal protein genes. They are highly expressed and tightly co-expressed, so they pull edges and alignment distance towards themselves whatever the question: on Kang 2018 PBMCs, before the knockout networks excluded them, ribosomal genes were 3% of those networks but 37% of their differentially-regulated genes, the same ones whatever the target. `DOWNSAMPLE` removes them from the object as soon as it is loaded, before downsampling, variable-gene selection or any QC, so no step downstream ever sees them: not `gene4use`, which GENIE3, scRank and hdWGCNA are built on, and not the knockout networks, which in `scTenifoldKnk()` itself would keep them. Each cell's total count is taken before the removal and used for the report's target expression, and subsetting genes leaves the `data` layer's values as they were, so an object normalised before it reached the pipeline keeps that normalisation. Steps that normalise on their own later, such as hdWGCNA's metacells and scTenifoldKnk's CPM, now do so over the remaining genes. The log and the report say how many genes went. Mitochondrial genes are `MT-`, or `MT.`/`MT_` as some converters write it (a bare `MT` prefix would take nuclear genes such as `MT2A` and `MTCH2`). Ribosomal genes are the cytosolic ribosomal proteins: `RPS`/`RPL` followed by a number and an optional A/X/Y variant, plus `RPLP0`–`RPLP2` and `RPSA`, so kinases such as `RPS6KA1` and paralogues such as `RPL22L1` stay in. Symbols are matched case-insensitively, so mouse `mt-`, `Rps` and `Rpl` count too. `gene4use` also drops the clone-named lncRNAs (`RP11-…`, `RP5-…`), as scRank's own selection does. A requested target is never removed for its name: it stays in the object, and in every network it has counts in.

### 2. Network Inference (`GENIE3`, `SCRANK`, `HDWGCNA`) and Knockout (`SCTENIFOLDKNK_BUILD`, `SCTENIFOLDKNK_KO`)
This is the heavy-lifting computational core. For each downsampled cellular identity, the pipeline infers a gene regulatory network using the method selected with `--network`: `genie3` runs [GENIE3](https://bioconductor.org/packages/release/bioc/html/GENIE3.html), `scrank` uses the scRank network strategy, and `hdwgcna` runs [hdWGCNA](https://smorabit.github.io/hdWGCNA/) on metacells. Each of these three returns regulatory interaction weights between genes for each cell state, which `RANK_SCORE` (below) turns into a `perb_score` for the requested target(s).

[scTenifoldKnk](https://github.com/cailab-tamu/scTenifoldKnk) is not one of those methods and is not selected with `--network`; it is a parallel track switched on with `--sctknk`. Rather than scoring how much a network leans on a target, it performs an actual in-silico knockout — it builds the wild-type network, zeros the target gene's outgoing edges, and compares the two networks by manifold alignment to get a genome-wide table of differentially-regulated genes with FDR. It runs in two steps. The wild-type network does not depend on the target, so `SCTENIFOLDKNK_BUILD` builds it once per cell identity — the expensive part: bootstrap networks and tensor decomposition. The knockout on it is target-specific, so `SCTENIFOLDKNK_KO` then runs the knockout, manifold alignment and differential regulation once per cell identity x target pair, in parallel, each reusing its identity's network. Both steps use scTenifoldKnk's own functions, arguments and seeds. The network is built on scTenifoldKnk's own gene set rather than `gene4use`: every gene detected in more than `--sctknk_min_pct` (5%) of the identity's cells, which is the gene filter `scTenifoldKnk()` applies itself. Mitochondrial and ribosomal protein genes are not among them, since `DOWNSAMPLE` removed them on load; `scTenifoldKnk()` would keep them. Two things differ from a single `scTenifoldKnk()` call. Its cell-level quality control is skipped, because `DOWNSAMPLE` has already chosen the cells. And a target detected in fewer cells than the threshold is kept in the network, where `scTenifoldKnk()` would refuse to knock it out, so every target gets a knockout wherever it has any counts. The differential regulation is scored with scTenifoldKnk 1.0.3's statistic, computed in the pipeline: each gene's `FC` is its squared alignment distance over the mean squared distance of the genes that were not knocked out. The current CRAN release, 1.1, averages over every gene, the knocked-out one included. That gene's distance is enormous by construction, so it inflates the mean for everyone else and leaves almost every gene at `p.adj` 1; on a 5,600-gene network one knockout went from 1 gene below FDR 0.05 to 12. Its output does not feed `RANK_SCORE` — it goes straight to its own merge step (`sctknk/sctenifoldknk_all_targets.txt`), which with the status table `sctknk/sctenifoldknk_status.txt` and the per-knockout summary is what this track publishes (with `--sctknk_plot`, the `plotKO()` PDFs go beside them in `sctknk/plots/`): the per-pair tables behind it are intermediate and stay in the work directory.

Three things sit on top of scTenifoldKnk's own statistic, because on its own it does not say what a target did. Whatever gene is knocked out of a network, the genes that move most in the alignment are the network's hubs, by an amount set by how many edges the knocked-out gene had, and `FC` only ranks a knockout's genes against each other. On Kang 2018 PBMCs every one of 25 targets returned the same two or three genes per cell type, and a target with no outgoing edge — which removes nothing — still returned "significant" genes from floating-point noise. So:

- **A knockout that removes nothing is skipped.** A target with no outgoing edge in a cell type's network (a gene kept below `--sctknk_min_pct` that got no edge, or one whose edges were all cut) is a no-op there. `SCTENIFOLDKNK_KO` writes no rows for it and says why on the status line, and the report labels it *no outgoing edges*.
- **A null model of random knockouts** (`--sctknk_null`, default 50). `SCTENIFOLDKNK_NULL` runs the same knockout, alignment and distance for random genes of each cell type's network, spread over the range of outgoing edge strength. Each gene of a target's knockout is then tested against what the same gene does under those knockouts: its log10 distance, centred on the knockout's median so the overall scale drops out, against a line fitted per gene over the random knockouts' log out-strength, in MAD units of the residuals, then standardised once more by their own median and MAD over the knockout's genes, since at the ends of the strength range a weak target's scores came out spread twice as wide as they should be. The table gains `z_null`, `p_null` and `p_null_adj`, a gene counts as moved when `p_null_adj` < 0.05 with `z_null` > 0, and the report, the summary and the enrichment ranking use them instead of `p.adj` and `log2FC`. Each knockout's `effect_pct` — the share of random knockouts that moved the network less — and `out_strength_pct` go on the status line. `--sctknk_null 0` switches the null off.
- **The tensor decomposition's rank is a parameter** (`--sctknk_td_k`, default 3, scTenifoldKnk's own). The CP decomposition returns a network of rank at most K, so at 3 every gene's outgoing edges are a mix of the same three patterns and every knockout removes a mix of the same three patterns — the hub pattern above is built in. On a test network, K = 10 made the knockouts of different targets move different genes (median Jaccard of their top-20 genes 0.05 against 0.54 at K = 3) while still tracking the edges each removed; 0 averages the bootstrap networks and skips the decomposition. The default stays at 3 until this is confirmed on a real dataset, so a run with the published setting is still the published method.

The knockout track takes the same targets `RANK_SCORE` does, line for line. A `;`-joined target is one joint knockout: the outgoing edges of all its genes are zeroed together in a single knocked-out network. `scTenifoldKnk()` itself only accepts one gene, but splitting the network build from the knockout makes a joint one possible. A `--target` file of `Brd4`, `Cstdc5`, `Stfa1;Mpo` therefore gives three knockouts per cell identity, with `Stfa1;Mpo` its own row set in the merged table. To get the single-gene knockouts as well, list `Stfa1` and `Mpo` on their own lines. A gene of a combined target that is not expressed in a cell identity is left out of that identity's knockout, and the rest is knocked out without it.

The two tracks are independent, so any combination of them can be run:

```bash
# a rank-score method on its own
nextflow run main.nf --network genie3 ...

# a rank-score method with scTenifoldKnk alongside it, in parallel
nextflow run main.nf --network genie3 --sctknk ...
nextflow run main.nf --network scrank --sctknk ...

# scTenifoldKnk on its own: DR gene table only, no RANK_SCORE/REPORT
nextflow run main.nf --sctknk ...
```

When both run, they share `DOWNSAMPLE` and the scTenifoldKnk results become sections 5 to 7 of that run's `REPORT` — the knockout network, the knockout gene table, and, with `--gsea_gmt`, the gene set enrichment over that table. Passing neither is an error.

The hdWGCNA module aggregates cells into metacells, builds an unsigned co-expression network and then adapts its topological overlap matrix (TOM) to what scRank expects from a network. Four things happen to the raw TOM:

1. **Sign recovery.** A TOM is always positive, so activation and repression are indistinguishable in it. Each edge is multiplied by the sign of the correlation between the same two genes across metacells, keeping the TOM magnitude but restoring its direction of effect. The network is built as `unsigned` for this reason: a `signed` adjacency already pushes anti-correlated pairs towards zero, so repressive edges would be cut before the sign could be recovered.
2. **Padding to a shared gene universe.** Genes dropped by hdWGCNA quality control return as zero rows and columns, so every cell identity is described over the same gene set. scRank refuses to align networks whose features differ.
3. **Sparsification.** A TOM is fully dense, which would make every gene a neighbour of every other one and flatten the degree and entropy terms scRank scores on. Edges whose absolute weight falls below the `--cut_ratio` quantile are cut.
4. **Rescaling.** Weights are divided by the largest absolute weight so that they span `[-1, 1]`, the range scRank's manifold alignment and its agonist mode both assume.

Cell identities that are too small to aggregate into metacells, or for which hdWGCNA otherwise fails, are skipped with a message in the log rather than failing the run. They are absent from the final ranking.

### 3. Perturbation Scoring (`RANK_SCORE`)
Using the list of target genes (`--target`) provided by the user, this module extracts the specific regulatory weight of the targets from the GENIE3 output. It calculates the perturbation score, which reflects how much the network relies on the specific target gene within that specific cell state.

A target gene the expression profile does not carry is dropped rather than fatal. `scRank::CreateScRank()` stops with *"Please check if the target gene is in the gene expression profile"* when handed one, which would take the whole task down before either table was written, so the genes are checked first and any that the object does not have — absent from it outright, or sitting in the matrix with no counts anywhere — are left out. Such a gene has no edges in a network built from that same object, so it contributes nothing to the score and a `;`-joined target keeps its name: `Stfa1;Mpo` with `Mpo` unexpressed is the same knockout as `Stfa1`. If *no* gene of a target survives, the target is skipped: both of its tables are written with headers and no rows, `MERGE` concatenates them into nothing, and the run carries on without it. Every other `CreateScRank` failure still aborts the task, since it is not a fact about the target and a cached empty table would hide it.

It also records the target's strongest connections. For each cell identity it reads the target's row of that identity's network — the same neighbourhood `scRank::init_mod()` uses to build a subnetwork, `abs(net[target, ]) > 0` — and keeps the `--top_connections` edges with the largest absolute weight, sign included. The result is a second table, one row per cell identity x target gene x partner.

### 4. Consolidate Results (`MERGE`)
This step collects the perturbation scores from all parallel `RANK_SCORE` tasks and merges them into a single, clean text file, ready for downstream visualization. The per-target connection tables are concatenated the same way, into `top_connections_all_targets.txt`.

### 5. Gene Set Enrichment (`GSEA_SCTENIFOLDKNK`)
Optional, and only on the knockout track: pass `--gsea_gmt` with a GMT file and the differentially-regulated genes are read as processes rather than as a list. For each cell identity x knocked-out gene combination, every gene is ranked by `log2FC` — the log2 of its squared manifold-alignment distance over the mean squared distance of the genes not knocked out, so still an ordering on how far the knockout moved it, but signed around the average gene — and that ranked list goes to [fgsea](https://bioconductor.org/packages/release/bioc/html/fgsea.html) against the supplied gene sets. This is the analysis [scTenifoldKnk's own paper](https://doi.org/10.1016/j.patter.2022.100434) runs on its output, and the ordering is the one it describes: distance itself, uncut, rather than the genes left after an FDR threshold.

Three things follow from that ranking and are worth knowing before reading the result. `log2FC` straddles zero, so fgsea runs two-sided (`scoreType = "std"`) and `NES` is signed: positive means a set clusters among the genes the knockout moved furthest, negative means it clusters among the ones that sat still while the rest of the network moved. Neither is up- or down-regulation — scTenifoldKnk compares network structure, not expression, so no direction of change exists in its output to recover. The knocked-out gene is not in its own ranking, since `SCTENIFOLDKNK_KO` leaves it out of the table: it sits at the top by construction and would otherwise hand a guaranteed hit to every set that annotates it. And the ranking covers the genes the knockout network was built on — those detected in more than `--sctknk_min_pct` of the identity's cells, typically several thousand — not the whole transcriptome. Sets are intersected with that universe before testing, so enrichment is relative to the genes this run actually modelled, and rarely detected genes never enter it. That is also why `--gsea_min_size` defaults to `10` rather than fgsea's conventional 15.

Gene sets are read from the file given and never fetched, so this runs on a compute node with no internet. The symbols have to match the data's own: MSigDB publishes mouse collections in MGI symbols and human ones in HGNC, and a human GMT against mouse data matches nothing — the run log names the symbols on each side and counts the overlap when that happens. Without `--gsea_gmt` the step does not run and the report says so. Results go to `sctknk/gsea_all_targets.txt`.

### 6. Report (`REPORT`)
Renders `perbscore_all_targets.txt` into a self-contained Quarto HTML report (`report/netperturb_report.html`): a searchable, filterable table of every cell type x target score, a heatmap of scores across all cell types and targets that were run, and the DOWNSAMPLE UMAP as a closing cell-identity overview. The overview names the `--network` method the run inferred its networks with, since scores are only comparable within one inference method. Pooled scores are bimodal, so the figures keep only the high mode, cut at `--score_quantile`; the table is always the full, uncut set. Sections 5 and 6 then go to the knockout track: an interactive plotly network per knocked-out gene and cell identity, the gene at the centre and the genes its knockout moved around it, node size following `log2FC`, with a menu to switch combination and the exact values on hover — followed by a queryable table of every differentially-regulated gene. Both are cut at FDR < 0.05, and the network is further capped at the `--sctknk_top_genes` strongest per combination while the table is not. Section 7 is the enrichment: a dot plot of the `--gsea_top_terms` strongest gene sets for the combination on display — one dot per set, `NES` on the x-axis about a zero line, leading-edge size as the dot's size and `-log(p.adj + 1)` as its shade, blue at the weak end of the cut through to red at the significant end — over a queryable table of every set clearing FDR < 0.25, set names shown cut at 50 characters with the whole name on hover — GSEA's own significance convention, looser than the 0.05 the knockout sections are cut at. A run without `--sctknk` shows a short "not run" note in place of each, and one without `--gsea_gmt` does the same for the enrichment alone. Every figure is embedded in the HTML, so the report is a single portable file, which the embedded plotly library makes a few megabytes.

After the knockout sections, a **Genes affected per knockout** heatmap counts, per cell type and target, the genes each knockout affected at FDR < 0.05, and marks knockouts that produced no table. The same numbers are published as `sctknk/sctenifoldknk_summary.txt`, which `MERGE_SCTENIFOLDKNK` writes even when no report runs.

A **Data quality** section then lists checks that can explain a skewed result, each with its reason: cells that are few, shallow or dominated by mito/ribo counts; leading principal components that track sequencing depth or batch (`--batch`); a network much denser or stronger than the other identities', which can shift every score in that identity; hdWGCNA networks over few metacells or with a weak scale-free fit; targets that are barely detected, vary no more than sampling noise, or sit at the bottom of their network; and knockouts that affect the same genes whatever the target. The numbers come from `downsample/identity_qc.tsv`, `downsample/target_expression.tsv` and `NETWORK_QC`'s `qc/network_qc.tsv` and `qc/target_network_qc.tsv`. The thresholds are rules of thumb, and nothing here filters the analysis.

## Quick Start
1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html) (`>=22.10.1`).
2. Install [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (highly recommended for full pipeline reproducibility).
3. Start running your analysis!

```bash
# Quick example
nextflow run netperturb/main.nf \
  -profile test,singularity

# Example with all parameters
nextflow run netperturb/main.nf \
  --obj /path/to/your/seurat_object.rds \
  --column clone_annotation \
  --species human \
  --n_cells 3000 \
  --binding antagonist \
  --n_cores 32 \
  --target /path/to/targets.txt \
  --network genie3 \
  --sctknk \
  --sctknk_top_genes 25 \
  --gsea_gmt /path/to/m5.go.bp.v2026.1.Mm.symbols.gmt \
  --gsea_top_terms 20 \
  --top_connections 15 \
  --score_quantile 0.75 \
  --outdir results \
  -profile singularity

``` 

### Inputs and References
NetPerturb requires the following main parameters:

`--obj`: Path to a fully processed Seurat object (`.rds` or compatible serialized object) containing normalized RNA assays and metadata annotations.

`--column`: Metadata column in the Seurat object that defines the cellular identities to compare, such as cell type, cluster, treatment group, clone, or phenotype.

`--species`: Species used by scRank when building and scoring regulatory networks. Accepted values depend on the underlying scRank annotation support, commonly `human` or `mouse`.

`--target`: Path to a text file containing the target genes to score. Each line should contain one target entry. If multiple genes should be evaluated together as one perturbation set, separate them with semicolons, for example `Stfa1;Mpo`. Max number is `two` targets at same time.

```sh
# Example
Brd4
Cstdc5
Stfa1;Mpo
```

`--assay`: Seurat assay holding the raw counts. Defaults to `RNA`. scRank, GENIE3, hdWGCNA and scTenifoldKnk all read an assay named `RNA`, so any other assay given here is renamed to `RNA` when the object is loaded, replacing an existing `RNA` assay if there is one. Use it for objects whose counts sit under another name, e.g. `--assay originalexp` for an object converted from a `SingleCellExperiment`.

`--network`: Network inference method driving the perturbation-scoring track. Supported values are `genie3`, `scrank` and `hdwgcna`. Optional only if `--sctknk` is given; otherwise the run has nothing to do and aborts.

`--sctknk`: Switches on the scTenifoldKnk track, which runs in parallel with whichever method `--network` selects and adds its differentially-regulated gene table as a section in that run's report. Defaults to `false`. With `--sctknk` and no `--network` it runs on its own, producing the table only — no `RANK_SCORE` or `REPORT`.

`--sctknk_plot`: Also draws scTenifoldKnk's `plotKO()` network for each cell identity x target knockout, written to `sctknk/plots/<identity>_sctknk_plot_<target>.pdf`. Defaults to `false`. The plot shows the knocked-out gene(s) together with the genes the knockout moved (`p.adj < 0.05`) and the wild-type edges among them, red positive and blue negative. It is drawn with `annotate = FALSE`, because annotation queries the Enrichr web service against human-only libraries, which a compute node without internet cannot reach. A plot that fails to draw is logged and skipped without affecting the table. Requires `--sctknk`.

`--sctknk_min_pct`: Fraction of an identity's cells a gene has to be detected in to enter that identity's scTenifoldKnk network. Defaults to `0.05`, the `qc_minPCT` gene filter `scTenifoldKnk()` applies itself, and is applied over every gene of the object, not `gene4use`. Target genes are kept below it as long as they have any count in the identity. Lowering it adds sparsely detected genes; the build cost grows roughly with the square of the gene count, so going much below the default can make `SCTENIFOLDKNK_BUILD` far slower and larger. The report's cells-per-identity section shows how many genes each identity's knockout network got. Requires `--sctknk`.

`--sctknk_null`: Random-gene knockouts `SCTENIFOLDKNK_NULL` runs per cell identity as the knockout track's null model. Defaults to `50`; `0` switches the null off and leaves scTenifoldKnk's own statistic only. Each costs one alignment, seconds to a couple of minutes on a large network, in one task per identity. The genes are drawn among those with an outgoing edge, never among the run's targets, spread over the range of outgoing edge strength with its two ends always in, so every target's strength lies inside the range the per-gene lines are fitted on. Fewer than 8 completed null knockouts means no calibration. Requires `--sctknk`.

`--sctknk_ndim`: Dimensions of the aligned manifold the knockout distances are measured in. Defaults to `2`, scTenifoldKnk's own. More dimensions carry more of a knockout's target-specific residue (on a test network the residual's correlation with the edges removed rose from 0.52 at 2 to 0.68 at 10) at a higher alignment cost (about 2.5x from 2 to 10). Requires `--sctknk`.

`--sctknk_td_k`: Components of the CP tensor decomposition that denoises the bootstrap networks into each identity's wild-type network. Defaults to `3`, scTenifoldKnk's own `td_K`. The result has rank at most K; see the knockout description above for what that does to the knockouts and why a larger value, or `0` to average the bootstrap networks instead, is worth trying. Requires `--sctknk`.

`--cut_ratio`: Quantile of absolute edge weight below which edges are cut, used by `--network hdwgcna`. Defaults to `0.95`, the same threshold scRank applies to its own networks, which keeps the strongest 5% of edges. Lower it to retain a denser network. Because a TOM has a different weight distribution than the regression coefficients scRank normally works with, this value is worth tuning on your data.

`--hdwgcna_min_cells`: Minimum number of cells an identity must have for `--network hdwgcna` to attempt metacell aggregation. Defaults to `150`. Identities below it are skipped.

`--top_connections`: Number of strongest edges `RANK_SCORE` keeps per target gene per cell identity, ranked on absolute weight. Defaults to `15`. These are written to `rank_scores/top_connections_all_targets.txt` for downstream use; the report no longer renders them.

`--sctknk_top_genes`: Number of differentially-regulated genes drawn around each knocked-out gene in the report's knockout network, taken from those clearing FDR < 0.05 and ranked on `log2FC`. Defaults to `25`. This caps the figure only — the table below it lists every significant gene.

`--gsea_gmt`: Path to a GMT file of gene sets for the enrichment step on the knockout track. No default — without it the step does not run and the report says so. Read from disk and never fetched, so it works on a node with no internet. The symbols must match `--species`: MSigDB publishes mouse collections as `.Mm.symbols.gmt` (MGI symbols) and human ones as `.Hs.symbols.gmt` (HGNC), so `m5.go.bp.*` or `m5.mpt.*` (the mammalian-phenotype sets the scTenifoldKnk paper uses) for mouse, `c5.go.bp.*` or `c2.cp.reactome.*` for human. Requires `--sctknk`; passed without it, the run warns and carries on.

`--gsea_min_size` / `--gsea_max_size`: Gene set size bounds, counted after intersecting each set with that combination's ranked list. Default to `10` and `500`. The lower bound sits below fgsea's conventional 15 on purpose: the ranked list is the knockout network's genes (detected in more than `--sctknk_min_pct` of cells) rather than the transcriptome, so sets arrive smaller than their nominal size and 15 would drop many of them before anything was tested.

`--gsea_top_terms`: Number of gene sets drawn per combination in the report's enrichment dot plot, taken from those clearing FDR < 0.25 and ranked on absolute `NES`, so the figure keeps the strongest sets at both ends of the ranking rather than only the positive end. Defaults to `20`. Caps the figure only — the table below it lists every significant set.

`--gsea_container`: Image carrying `fgsea` for the enrichment step, built from `container/gsea/Dockerfile`. Point it at a local `.sif` if the image has not been pushed.

`--score_quantile`: Quantile of the pooled log10 perturbation scores below which scores are cut from the report figures. Defaults to `0.75`, so the cut falls at q3 and the figures show the top quarter of scores. Pooled scores are bimodal, a low mode of near-zero values sitting well below the mode that carries the signal, and the report is only useful once the low one is gone. The distribution figure draws the cut over the full set of scores, so the line can be checked against where the two modes actually separate and this value tuned to land in the valley between them. The score table is never cut.

`--n_cells`: Maximum number of cells to keep per cellular identity during downsampling. If an identity has fewer cells than this value, the pipeline uses all available cells for that identity.

`--min_cells`: Fewest cells a cellular identity must contribute for `DOWNSAMPLE` to keep it. Defaults to `150`; `0` keeps every identity. The count it is compared against is what survives downsampling — `min(identity size, --n_cells)` — rather than the identity's size in the input object, since that is what every network method downstream actually receives, and it is the same number `--hdwgcna_min_cells` is measured against one step further down. An identity below it is dropped before the object is split, with a line in the log naming it and its size. Setting `--n_cells` below `--min_cells` drops every identity and aborts the run, which the error says explicitly.

`--n_hvg`: Number of highly variable genes `DOWNSAMPLE` puts into `gene4use`, the gene set every network method is built on, together with the species' transcription factors and drug targets and the requested targets. Defaults to `2000`. Variable features already stored on the object are reused when there are at least this many; otherwise they are computed with Seurat's `vst` on the downsampled cells. More genes mean larger networks and longer runs, most of all for GENIE3. It does not affect the scTenifoldKnk track, which picks its own genes with `--sctknk_min_pct`.

`--batch`: Metadata column naming each cell's donor, sample or batch. No default. Used only by the report's Data quality section, which then shows each identity's batch composition and how much of the variation in its leading principal components batch explains: correlation between genes that comes from donors differing rather than from regulation within them. It does not change any network or score.

`--seed`: Seed for every step that draws random numbers. Defaults to `1`. It fixes which cells `DOWNSAMPLE` keeps (and the UMAP it draws when the object has none), GENIE3's random forests (on any `--n_cores`, since GENIE3 seeds its workers through doRNG), hdWGCNA's metacell sampling and WGCNA's own `randomSeed`, scTenifoldKnk's bootstrap networks, tensor decomposition and manifold alignment, and fgsea's adaptive sampler. The default is the seed scTenifoldKnk hardcodes, so the knockout track gives the same table as `scTenifoldKnk()` run by hand. scRank's `Constr_net` also seeds itself with `1` internally and does not see `--seed`, so `--network scrank` networks do not change with it. Two runs with the same seed and inputs give the same results. Changing it changes which cells are kept, so every step reruns, even under `-resume`.

`--binding`: Perturbation mode passed to the scoring step, for example `antagonist` or `agonist`.

`--n_cores`: Number of CPU cores requested for parallelizable network inference and scoring steps.

`--outdir`: Directory where the final results and pipeline reports will be written. Defaults to `results`.

### Outputs
If successfully run, the workflow will generate its primary output in the specified --outdir:

rank_scores/perbscore_all_targets.txt: A consolidated table containing the cell identity (cell_type), the evaluated gene (target), and its final regulatory importance (perb_score).

```sh
# Example
cell_type	target	binding	perb_score
sensitive	Stfa1;Mpo	antagonist	1.38837233985176e-06
resistant	Stfa1;Mpo	antagonist	3.86382149842044e-06
sensitive	Brd4	antagonist	1.5395821350262e-06
resistant	Brd4	antagonist	1.61320913209275e-06
sensitive	Cstdc5	antagonist	1.30868421341405e-06
resistant	Cstdc5	antagonist	2.91128461301128e-06
```

downsample/target_qc.tsv: Target genes removed by `DOWNSAMPLE`'s QC check, one row per gene, with the requested target it came from, the reason (`absent from the expression profile` or `zero counts in every retained cell`) and the action taken. Header only when every gene passed. `downsample/targets_qc.txt` holds the targets that were actually analysed.

```sh
# Example
target	gene	reason	action
Hhip	Hhip	zero counts in every retained cell	target not analysed
Stfa1;Mpo	Mpo	absent from the expression profile	dropped from the target; analysed as Stfa1
```

rank_scores/top_connections_all_targets.txt: The strongest edges each target holds in each cell identity's network, up to `--top_connections` per target gene per identity, ranked on absolute weight. Published for downstream use; the report does not render it.

```sh
# Example
cell_type	target	binding	target_gene	partner	weight	rank
sensitive	Brd4	antagonist	Brd4	Myc	0.8123	1
sensitive	Brd4	antagonist	Brd4	Ccnd1	-0.6640	2
resistant	Brd4	antagonist	Brd4	Myc	0.4471	1
```

sctknk/sctenifoldknk_all_targets.txt: Written when `--sctknk` is passed, and the only table the knockout track publishes. One row per cell identity x knocked-out gene x differentially-regulated gene, with the manifold-alignment `distance`, `Z`-score, `FC`, and `p.value`/`p.adj`. A combined `--target` line appears under its own name, e.g. `Stfa1;Mpo`, as one joint knockout. `FC` is the gene's squared alignment distance over the mean squared distance of that run's genes, the knocked-out ones excluded, so it is a ratio to the average rather than a differential-expression fold change, and it carries no direction. The knocked-out gene or genes are left out of their own rows: they always move furthest by construction, so they say nothing about the knockout's effect. With `--sctknk_null` above 0 the table also carries `z_null`, `p_null` and `p_null_adj`: the gene tested against what it does under random knockouts of the same network (see the knockout description above).

sctknk/sctenifoldknk_status.txt: One line per cell identity x target: `status` (`ok`, `no outgoing edges`, `not expressed`, `no wild-type network` or `alignment failed`), the genes knocked out, the target's out-degree and out-strength in that identity's network with the share of null knockouts weaker than it (`out_strength_pct`), the knockout's median distance with the share of null knockouts that moved the network less (`effect_pct`), the number of null knockouts, the genes tested, and its hits under scTenifoldKnk's `p.adj` (`hits_raw`) and under the null (`hits_null`). A pair missing from this table produced nothing at all: it timed out on every attempt and was ignored.

sctknk/sctenifoldknk_summary.txt: One line per cell identity x target with the genes tested, the genes affected (`p_null_adj` < 0.05 with `z_null` > 0 when the null ran, `p.adj` < 0.05 otherwise), their share and the gene moved furthest.

```sh
# Example
cell_type	target	gene	distance	Z	FC	p.value	p.adj
resistant	Brd4	Brd4	16.90000	8.42	56.4	5.84e-14	8.76e-11
resistant	Brd4	G17_Brd4	14.43217	7.09	41.1	1.41e-10	9.83e-08
```

sctknk/gsea_all_targets.txt: Written when `--sctknk` is passed together with `--gsea_gmt`. One row per cell identity x knocked-out gene x gene set, carrying fgsea's own columns: `ES` and its size-normalised form `NES`, `pval` and the Benjamini-Hochberg `padj`, `size` — the set's size *after* intersection with that combination's ranked list, so smaller than the set's nominal size — and `leadingEdge`, the genes carrying the enrichment, `;`-joined. `NES` is signed: positive for a set concentrated among the genes the knockout moved furthest, negative for one concentrated among the genes it left alone.

```sh
# Example
cell_type	target	pathway	pval	padj	log2err	ES	NES	size	leadingEdge
resistant	Brd4	GOBP_CHROMATIN_REMODELING	1e-50	4e-50	NA	1	3.126	30	G17_Brd4;G32_Brd4;G24_Brd4
sensitive	Brd4	GOBP_CHROMATIN_REMODELING	1.67e-11	3.33e-11	0.863	0.773	3.043	30	G5_Brd4;G1_Brd4
```

A `pval` of `1e-50` with an `NA` `log2err` is fgsea's floor, not a failure: the true value is smaller than it estimates.

report/netperturb_report.html: A self-contained Quarto report, with a queryable table of every perturbation score, a cell type x target heatmap of the scores above the `--score_quantile` cut, then the knockout network and the knockout differentially-regulated gene table, both cut at FDR < 0.05, and the `--network` inference method the run used.

Other intermediate files (such as split matrices and raw GENIE3 weights) are temporarily stored in the work directory and can be retained or discarded based on standard Nextflow cache management.

### Development
[`docs/IMPLEMENTATION.md`](docs/IMPLEMENTATION.md) records how the pipeline was built, grouped into waves of work rather than individual commits.

#### Tests
The pipeline is covered by [nf-test](https://www.nf-test.com/). The default suite runs every process through its `stub` block, so it pulls no containers and executes no R, and finishes in about a minute:

```bash
nf-test test
```

It checks the wiring rather than the science: that `--network` selects exactly one inference process, that `--sctknk` adds its track without disturbing that choice and runs alone when `--network` is omitted, that each line of the target file becomes its own `RANK_SCORE` task, that `MERGE` gathers them under a single header, that an unsupported `--network` aborts before any task is launched, that `--gsea_gmt` adds the enrichment step to the knockout track while its absence routes a sentinel table to the report instead, and that every network module names its output so `rank_score.R` can still recover the cell identity from the file name.

The end to end run is opt-in, since it downloads the test object and pulls containers. It is excluded from the default suite and has its own config:

```bash
nf-test test -c nf-test.integration.config --tag integration --profile test,singularity
```

Note for machines running the uutils reimplementation of coreutils, the default on recent Ubuntu: Nextflow's task wrapper times tasks with `date +%s%3N`, and uutils ignores the `%3N` width modifier. The wrapper then aborts every task with `Unexpected: unbound variable` before the script runs. This affects any Nextflow pipeline on such a machine, not just this one. Installing GNU coreutils resolves it.

### Credits
NetPerturb is developed and maintained by the Marques-Coelho Bioinformatics Lab(MCBLab).

### Citations
If you use this pipeline in your research, please cite:

...
