# Usage

## Requirements

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) `>=22.10.1`.
- A container engine: [Singularity/Apptainer](https://docs.sylabs.io/) (recommended on HPC) or [Docker](https://docs.docker.com/engine/installation/). Select it with `-profile singularity` or `-profile docker`.
- Enough memory for the network steps. See [Resources](#resources) below.

!!! warning "scTenifoldKnk image"
    The knockout track's image is not published to a registry yet. Build it from `container/sctenifoldknk/Dockerfile` and point the two knockout processes at it before using `--sctknk`. The [FAQ](faq.md#how-do-i-provide-the-sctenifoldknk-container) shows how.

## Running the pipeline

```bash
nextflow run MCBLab/NetPerturb \
  -profile singularity \
  --obj seurat_object.rds \
  --column cell_type \
  --species mouse \
  --target targets.txt \
  --network genie3 \
  --outdir results
```

Nextflow writes its working files to `work/` in the launch directory and the published results to `--outdir`. Add `-resume` to a relaunched command to reuse every task that already finished. Only steps whose inputs or parameters changed are run again.

### Choosing what runs

The pipeline has two independent tracks. `--network` selects the **scoring track** and `--sctknk` switches on the **knockout track**. At least one of them is required:

```bash
# scoring track only: networks, perturbation scores, report
nextflow run MCBLab/NetPerturb --network genie3 ...

# both tracks in parallel, sharing DOWNSAMPLE; knockout sections added to the report
nextflow run MCBLab/NetPerturb --network scrank --sctknk ...

# knockout track only: differentially-regulated gene table (and GSEA), no report
nextflow run MCBLab/NetPerturb --sctknk ...
```

| `--network` | Method | Notes |
|---|---|---|
| `genie3` | [GENIE3](https://bioconductor.org/packages/release/bioc/html/GENIE3.html) random-forest regression | Scales with `--n_hvg`; the slowest for large gene sets. |
| `scrank` | scRank's own `Constr_net` | Seeds itself internally; does not follow `--seed`. |
| `hdwgcna` | [hdWGCNA](https://smorabit.github.io/hdWGCNA/) on metacells | The TOM is signed from correlations, sparsified at `--cut_ratio` and rescaled to `[-1, 1]`. Identities below `--hdwgcna_min_cells` are skipped. |

Scores are only comparable between runs that used the same `--network` method.

## Inputs

### Seurat object (`--obj`)

A processed Seurat object saved with `saveRDS()`. The pipeline needs:

- **raw counts** in an assay. It reads `RNA` by default; set `--assay` if yours is named differently, e.g. `originalexp` for an object converted from a `SingleCellExperiment`.
- a **metadata column** with the identities to compare, passed as `--column`.
- optionally, a UMAP embedding and variable features. Both are reused if present and computed otherwise.

### Target file (`--target`)

A plain text file with one target per line. Genes joined by `;` form a single combined target: scored as one joint perturbation, and knocked out together as one joint knockout. The largest combination supported is two genes.

```text
Brd4
Cstdc5
Stfa1;Mpo
```

The file above produces three `RANK_SCORE` tasks and three knockouts per identity. To also get `Stfa1` and `Mpo` individually, list them on their own lines.

Every gene is checked against the retained cells before any network is built. A gene absent from the object, or with no counts in any retained cell, is removed: a combined target loses only that gene, and a target left with no genes is dropped. The run stops only if no target survives. What was removed, and why, is listed in `downsample/target_qc.tsv` and in the report.

Symbols must match the object (MGI for mouse, HGNC for human), and `--species` must match both.

### Adding targets to a finished run (`--extra_target`)

Changing `--target` reruns everything from `DOWNSAMPLE` on, because the targets are part of the gene set the networks are built on. To add targets to a run whose networks are already built, keep `--target` as it was and list the new ones in a second file, in the same format, given as `--extra_target`, with `-resume`:

```bash
nextflow run main.nf <same options as before> --extra_target more_targets.txt -resume
```

The extras skip `DOWNSAMPLE` and every network step. `EXTRA_TARGET_QC` checks them against the networks each track already has, and only the ones those carry are scored by `RANK_SCORE` and knocked out by `SCTENIFOLDKNK_KO`; every earlier task, including the `--target` scores and knockouts, comes from the cache. An extra no network carries — not among the highly variable genes, transcription factors, drug targets and `--target` genes the rank-score networks were built on, or below `--sctknk_min_pct` in every identity for the knockout track — is not run, and the report's Target QC section and `qc/extra_target_qc_<track>.tsv` say so. To analyse it, move it to `--target` and run from the start.

## Parameters

### Required

| Parameter | Description |
|---|---|
| `--obj` | Path or URL to the Seurat object (`.rds`). |
| `--column` | Metadata column defining the cell identities to compare. |
| `--species` | `human` or `mouse`. Selects the transcription factor and drug target lists and the scRank annotation. |
| `--target` | Target file, as described above. |
| `--network` and/or `--sctknk` | At least one track. |

### Downsampling and gene selection

| Parameter | Default | Description |
|---|---|---|
| `--assay` | `RNA` | Assay holding raw counts. Any other assay is renamed to `RNA` on load. |
| `--n_cells` | — | Maximum cells kept per identity. Smaller identities keep every cell. |
| `--min_cells` | `150` | Identities with fewer cells after downsampling are dropped before any network is built. `0` keeps all. |
| `--n_hvg` | `2000` | Highly variable genes added to `gene4use`, the gene set of the scoring-track networks. |
| `--batch` | none | Metadata column naming each cell's donor, sample or batch. Used only by the report's *Data quality* section: batch composition of each identity and the variance batch explains. |
| `--seed` | `1` | Seed for every random step: cell sampling, GENIE3, hdWGCNA/WGCNA, scTenifoldKnk and fgsea. |
| `--cell_subset` | — | Comma-separated `--column` values the run covers, e.g. `HSC,Prog,Mono`; every kept identity by default. Narrows both tracks' networks, scores and knockouts; `DOWNSAMPLE` and its tables still cover every identity. Names matching nothing are warned about; none matching stops the run. |
| `--extra_target` | — | Targets to add to a run whose networks are already built, in `--target`'s format; see [above](#adding-targets-to-a-finished-run-extra_target). Only those the networks carry are run. |

### Scoring track

| Parameter | Default | Description |
|---|---|---|
| `--network` | — | `genie3`, `scrank` or `hdwgcna`. |
| `--binding` | `antagonist` | Perturbation mode for scRank scoring: `antagonist` or `agonist`. |
| `--top_connections` | `15` | Strongest edges recorded per target gene per identity. |
| `--cut_ratio` | `0.95` | hdWGCNA only: quantile of absolute edge weight below which edges are cut. |
| `--hdwgcna_min_cells` | `150` | hdWGCNA only: identities below this are not aggregated into metacells. |
| `--score_quantile` | `0.75` | Report only: pooled log10 scores below this quantile are left out of the figures. Tables are never cut. |

### Knockout track

| Parameter | Default | Description |
|---|---|---|
| `--sctknk` | `false` | Switch on the scTenifoldKnk knockout track. |
| `--sctknk_min_pct` | `0.05` | A gene must be detected in more than this fraction of an identity's cells to enter its knockout network. Targets are always kept. Build cost grows with the square of the gene count. |
| `--sctknk_null` | `50` | Random-gene knockouts per identity used as the knockout track's null model; each gene of a target's knockout is tested against what it does under them (`z_null`, `p_null_adj`). `0` leaves scTenifoldKnk's own statistic only. |
| `--sctknk_ndim` | `2` | Dimensions of the aligned manifold (scTenifoldKnk's default). More carry more target-specific signal at a higher alignment cost. |
| `--sctknk_td_k` | `3` | Rank of the tensor decomposition that builds the wild-type network (scTenifoldKnk's `td_K`). At 3 every knockout moves the same hub genes; 10 separated targets on a test network; `0` averages the bootstrap networks instead. |
| `--sctknk_plot` | `false` | Also draw scTenifoldKnk's `plotKO()` network per identity × target, as PDFs under `sctknk/plots/`. |
| `--sctknk_top_genes` | `25` | Report only: genes drawn around each knocked-out target in the knockout network figure. |

### Gene set enrichment

| Parameter | Default | Description |
|---|---|---|
| `--gsea_gmt` | — | GMT file of gene sets. Without it no enrichment runs. Requires `--sctknk`. |
| `--gsea_min_size` | `10` | Smallest gene set tested, counted after intersecting with the ranked genes. |
| `--gsea_max_size` | `500` | Largest gene set tested. |
| `--gsea_top_terms` | `20` | Report only: gene sets drawn per combination in the enrichment dot plot. |
| `--gsea_container` | `diegomscoelho/netperturb-gsea:1.0` | Image carrying fgsea. May be a local `.sif`. |

Use gene sets in the same symbols as your data. MSigDB publishes mouse collections as `*.Mm.symbols.gmt` (e.g. `m5.go.bp`, or `m5.mpt`, the phenotype sets the scTenifoldKnk paper uses) and human ones as `*.Hs.symbols.gmt` (e.g. `c5.go.bp`, `c2.cp.reactome`). The file is read from disk, so compute nodes need no internet access.

### General

| Parameter | Default | Description |
|---|---|---|
| `--n_cores` | — | CPUs per network and scoring task. Network tasks request exactly this many. |
| `--outdir` | `results` | Where results are published. |
| `--tracedir` | `<outdir>/pipeline_info` | Nextflow's execution report, timeline, trace and DAG. |

## Profiles

| Profile | Purpose |
|---|---|
| `singularity`, `docker`, `podman`, `shifter`, `charliecloud` | Container engine. |
| `arm` | Docker on ARM machines (runs the amd64 images under emulation). |
| `test` | Small mouse AML object with `--network genie3`, 500 cells per identity, 2 cores. |
| `test_ocasio`, `test_vangalen` | Larger public datasets (human AML bone marrow for van Galen, scored with hdWGCNA). They need more memory. |

Combine a data profile with an engine profile: `-profile test,singularity`.

## Resources

Requests per task, from `conf/modules.config`:

| Process | Memory | CPUs |
|---|---|---|
| `DOWNSAMPLE` | 16 GB | 1 |
| `GENIE3` | 16 GB | `--n_cores` |
| `SCRANK` | 64 GB | `--n_cores` |
| `HDWGCNA` | 24 GB | `--n_cores` |
| `RANK_SCORE` | 64 GB | 4 |
| `SCTENIFOLDKNK_BUILD` | 24 GB | `--n_cores` |
| `SCTENIFOLDKNK_KO` | 16 GB | `--n_cores` |
| `GSEA_SCTENIFOLDKNK` | 8 GB | 2 |
| `REPORT` | 4 GB | 1 |

To change them, or to add an executor such as SLURM, pass your own config with `-c`:

```groovy
// custom.config
process {
    executor = 'slurm'
    queue    = 'long'
    withName: SCRANK              { memory = 128.GB }
    withName: SCTENIFOLDKNK_BUILD { memory = 48.GB; time = 24.h }
}
```

```bash
nextflow run MCBLab/NetPerturb -profile singularity -c custom.config ...
```

## Reproducibility

Two runs with the same inputs, parameters and `--seed` give the same results. The seed controls which cells are kept, so changing it re-runs every step, even under `-resume`. The one exception is `--network scrank`: scRank seeds its own network construction with `1` and ignores `--seed`.

The report records the exact command line, Nextflow and pipeline versions, container engine and every setting that shapes the result.

## Testing

The pipeline is covered by [nf-test](https://www.nf-test.com/). The default suite runs every process through its stub, so it pulls no containers and runs no R code. It finishes in a few minutes:

```bash
nf-test test
```

The end-to-end run downloads the test object and pulls every container, so it is opt-in:

```bash
nf-test test -c nf-test.integration.config --tag integration --profile test,singularity
```
