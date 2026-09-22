# NetPerturb 🧬

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
This is the initial step of the process. It ingests a fully processed Seurat object (`.rds`) and identifies the user-defined metadata column containing the cell identities (e.g., cell types or clones). To ensure statistical robustness and equitable GRN inference, it randomly downsamples the cells from each identity to a specified maximum number (`--n_cells`), balancing the computational load. It also writes a UMAP of the retained cells coloured by `--column` to `downsample/umap_<column>.png`, so the identities entering the analysis, and their relative sizes after downsampling, can be checked at a glance. An embedding already present on the object is reused; one is computed only if the object carries none. Drawing this figure is guarded, so a plotting failure leaves a placeholder image rather than stopping the run.

### 2. Network Inference (`GENIE3`, `SCRANK`, `HDWGCNA`) and Knockout (`SCTENIFOLDKNK`)
This is the heavy-lifting computational core. For each downsampled cellular identity, the pipeline infers a gene regulatory network using the method selected with `--network`: `genie3` runs [GENIE3](https://bioconductor.org/packages/release/bioc/html/GENIE3.html), `scrank` uses the scRank network strategy, and `hdwgcna` runs [hdWGCNA](https://smorabit.github.io/hdWGCNA/) on metacells. Each of these three returns regulatory interaction weights between genes for each cell state, which `RANK_SCORE` (below) turns into a `perb_score` for the requested target(s).

[scTenifoldKnk](https://github.com/cailab-tamu/scTenifoldKnk) is not one of those methods and is not selected with `--network`; it is a parallel track switched on with `--sctknk`. Rather than scoring how much a network leans on a target, it performs an actual in-silico knockout — it builds the wild-type network, zeros the target gene's outgoing edges, and compares the two networks by manifold alignment to get a genome-wide table of differentially-regulated genes with FDR. Because the knockout is target-specific, it runs once per cell identity x gene pair rather than once per cell identity, and its output does not feed `RANK_SCORE` — it goes straight to its own merge step (`sctknk/sctenifoldknk_all_targets.txt`), which is the only file this track publishes: the per-pair tables behind it are intermediate and stay in the work directory.

scTenifoldKnk knocks out exactly one gene per run, so it cannot take a `;`-joined combined target the way `RANK_SCORE` does. Rather than skip those lines, the target file is flattened to the individual genes it names and each is knocked out separately, deduplicated across the whole file. A `--target` file of `Brd4`, `Cstdc5`, `Stfa1;Mpo` therefore gives four knockouts per cell identity — `Brd4`, `Cstdc5`, `Stfa1` and `Mpo`, each its own row set in the merged table. The scoring track is unaffected: `RANK_SCORE` still receives `Stfa1;Mpo` intact and scores it as one joint perturbation.

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

It also records the target's strongest connections. For each cell identity it reads the target's row of that identity's network — the same neighbourhood `scRank::init_mod()` uses to build a subnetwork, `abs(net[target, ]) > 0` — and keeps the `--top_connections` edges with the largest absolute weight, sign included. The result is a second table, one row per cell identity x target gene x partner.

### 4. Consolidate Results (`MERGE`)
This step collects the perturbation scores from all parallel `RANK_SCORE` tasks and merges them into a single, clean text file, ready for downstream visualization. The per-target connection tables are concatenated the same way, into `top_connections_all_targets.txt`.

### 5. Gene Set Enrichment (`GSEA_SCTENIFOLDKNK`)
Optional, and only on the knockout track: pass `--gsea_gmt` with a GMT file and the differentially-regulated genes are read as processes rather than as a list. For each cell identity x knocked-out gene combination, every gene is ranked by `log2FC` — the log2 of its squared manifold-alignment distance over the run's mean squared distance, so still an ordering on how far the knockout moved it, but signed around the average gene — and that ranked list goes to [fgsea](https://bioconductor.org/packages/release/bioc/html/fgsea.html) against the supplied gene sets. This is the analysis [scTenifoldKnk's own paper](https://doi.org/10.1016/j.patter.2022.100434) runs on its output, and the ordering is the one it describes: distance itself, uncut, rather than the genes left after an FDR threshold.

Three things follow from that ranking and are worth knowing before reading the result. `log2FC` straddles zero, so fgsea runs two-sided (`scoreType = "std"`) and `NES` is signed: positive means a set clusters among the genes the knockout moved furthest, negative means it clusters among the ones that sat still while the rest of the network moved. Neither is up- or down-regulation — scTenifoldKnk compares network structure, not expression, so no direction of change exists in its output to recover. The knocked-out gene is dropped from its own ranking, since it sits at the top by construction and would otherwise hand a guaranteed hit to every set that annotates it. And the ranking covers `gene4use` — the highly-variable, transcription-factor and drug-target set `DOWNSAMPLE` builds — not the whole transcriptome, so sets are intersected with that universe before testing and enrichment is relative to the genes this run actually modelled. That is also why `--gsea_min_size` defaults to `10` rather than fgsea's conventional 15.

Gene sets are read from the file given and never fetched, so this runs on a compute node with no internet. The symbols have to match the data's own: MSigDB publishes mouse collections in MGI symbols and human ones in HGNC, and a human GMT against mouse data matches nothing — the run log names the symbols on each side and counts the overlap when that happens. Without `--gsea_gmt` the step does not run and the report says so. Results go to `sctknk/gsea_all_targets.txt`.

### 6. Report (`REPORT`)
Renders `perbscore_all_targets.txt` into a self-contained Quarto HTML report (`report/netperturb_report.html`): a searchable, filterable table of every cell type x target score, a heatmap of scores across all cell types and targets that were run, and the DOWNSAMPLE UMAP as a closing cell-identity overview. The overview names the `--network` method the run inferred its networks with, since scores are only comparable within one inference method. Pooled scores are bimodal, so the figures keep only the high mode, cut at `--score_quantile`; the table is always the full, uncut set. Sections 5 and 6 then go to the knockout track: an interactive plotly network per knocked-out gene and cell identity, the gene at the centre and the genes its knockout moved around it, node size following `log2FC`, with a menu to switch combination and the exact values on hover — followed by a queryable table of every differentially-regulated gene. Both are cut at FDR < 0.05, and the network is further capped at the `--sctknk_top_genes` strongest per combination while the table is not. Section 7 is the enrichment: a dot plot of the `--gsea_top_terms` strongest gene sets for the combination on display — one dot per set, `NES` on the x-axis about a zero line, leading-edge size as the dot's size and `-log(p.adj + 1)` as its shade, blue at the weak end of the cut through to red at the significant end — over a queryable table of every set clearing FDR < 0.25, set names shown cut at 50 characters with the whole name on hover — GSEA's own significance convention, looser than the 0.05 the knockout sections are cut at. A run without `--sctknk` shows a short "not run" note in place of each, and one without `--gsea_gmt` does the same for the enrichment alone. Every figure is embedded in the HTML, so the report is a single portable file, which the embedded plotly library makes a few megabytes.

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

`--network`: Network inference method driving the perturbation-scoring track. Supported values are `genie3`, `scrank` and `hdwgcna`. Optional only if `--sctknk` is given; otherwise the run has nothing to do and aborts.

`--sctknk`: Switches on the scTenifoldKnk track, which runs in parallel with whichever method `--network` selects and adds its differentially-regulated gene table as a section in that run's report. Defaults to `false`. With `--sctknk` and no `--network` it runs on its own, producing the table only — no `RANK_SCORE` or `REPORT`.

`--cut_ratio`: Quantile of absolute edge weight below which edges are cut, used by `--network hdwgcna`. Defaults to `0.95`, the same threshold scRank applies to its own networks, which keeps the strongest 5% of edges. Lower it to retain a denser network. Because a TOM has a different weight distribution than the regression coefficients scRank normally works with, this value is worth tuning on your data.

`--hdwgcna_min_cells`: Minimum number of cells an identity must have for `--network hdwgcna` to attempt metacell aggregation. Defaults to `150`. Identities below it are skipped.

`--top_connections`: Number of strongest edges `RANK_SCORE` keeps per target gene per cell identity, ranked on absolute weight. Defaults to `15`. These are written to `rank_scores/top_connections_all_targets.txt` for downstream use; the report no longer renders them.

`--sctknk_top_genes`: Number of differentially-regulated genes drawn around each knocked-out gene in the report's knockout network, taken from those clearing FDR < 0.05 and ranked on `log2FC`. Defaults to `25`. This caps the figure only — the table below it lists every significant gene.

`--gsea_gmt`: Path to a GMT file of gene sets for the enrichment step on the knockout track. No default — without it the step does not run and the report says so. Read from disk and never fetched, so it works on a node with no internet. The symbols must match `--species`: MSigDB publishes mouse collections as `.Mm.symbols.gmt` (MGI symbols) and human ones as `.Hs.symbols.gmt` (HGNC), so `m5.go.bp.*` or `m5.mpt.*` (the mammalian-phenotype sets the scTenifoldKnk paper uses) for mouse, `c5.go.bp.*` or `c2.cp.reactome.*` for human. Requires `--sctknk`; passed without it, the run warns and carries on.

`--gsea_min_size` / `--gsea_max_size`: Gene set size bounds, counted after intersecting each set with that combination's ranked list. Default to `10` and `500`. The lower bound sits below fgsea's conventional 15 on purpose: the ranked list is `gene4use` rather than the transcriptome, so sets arrive smaller than their nominal size and 15 would drop most of a collection before anything was tested.

`--gsea_top_terms`: Number of gene sets drawn per combination in the report's enrichment dot plot, taken from those clearing FDR < 0.25 and ranked on absolute `NES`, so the figure keeps the strongest sets at both ends of the ranking rather than only the positive end. Defaults to `20`. Caps the figure only — the table below it lists every significant set.

`--gsea_container`: Image carrying `fgsea` for the enrichment step, built from `container/gsea/Dockerfile`. Point it at a local `.sif` if the image has not been pushed.

`--score_quantile`: Quantile of the pooled log10 perturbation scores below which scores are cut from the report figures. Defaults to `0.75`, so the cut falls at q3 and the figures show the top quarter of scores. Pooled scores are bimodal, a low mode of near-zero values sitting well below the mode that carries the signal, and the report is only useful once the low one is gone. The distribution figure draws the cut over the full set of scores, so the line can be checked against where the two modes actually separate and this value tuned to land in the valley between them. The score table is never cut.

`--n_cells`: Maximum number of cells to keep per cellular identity during downsampling. If an identity has fewer cells than this value, the pipeline uses all available cells for that identity.

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

rank_scores/top_connections_all_targets.txt: The strongest edges each target holds in each cell identity's network, up to `--top_connections` per target gene per identity, ranked on absolute weight. Published for downstream use; the report does not render it.

```sh
# Example
cell_type	target	binding	target_gene	partner	weight	rank
sensitive	Brd4	antagonist	Brd4	Myc	0.8123	1
sensitive	Brd4	antagonist	Brd4	Ccnd1	-0.6640	2
resistant	Brd4	antagonist	Brd4	Myc	0.4471	1
```

sctknk/sctenifoldknk_all_targets.txt: Written when `--sctknk` is passed, and the only file the knockout track publishes. One row per cell identity x knocked-out gene x differentially-regulated gene, with the manifold-alignment `distance`, `Z`-score, `FC`, and `p.value`/`p.adj`. Combined `--target` lines appear here as their individual genes, one knockout each. `FC` is the gene's squared alignment distance over the mean squared distance of that run, so it is a ratio to the average rather than a differential-expression fold change, and it carries no direction. The knocked-out gene appears in its own rows, normally at the top.

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
