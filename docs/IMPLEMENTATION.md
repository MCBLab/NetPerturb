# Implementation History

How NetPerturb was built, grouped into waves of work rather than individual commits. Each wave is one coherent unit of change, usually a pull request or a short run of commits that only makes sense together. Dates are when the work was done; the PR column says how it reached `main`, and names the branch instead when it has not landed yet. Most recent first.

| Wave | Theme | Done | PR |
|---|---|---|---|
| 21 | Gene set enrichment on the knockout table | Sep 2026 | branch `sctknk` |
| 20 | The knockout track running in parallel | Sep 2026 | branch `sctknk` |
| 19 | scTenifoldKnk replaces scTenifoldNet | Sep 2026 | branch `sctknk` |
| 18 | Running it on real data | Sep 2026 | #17 |
| 17 | Bigger test data, and DOWNSAMPLE without scRank | Aug 2026 | #17 |
| 16 | The strongest edges a target holds | Aug 2026 | #17 |
| 15 | What the report actually shows | Aug 2026 | #17 |
| 14 | The metro map | Aug 2026 | #17 |
| 13 | REPORT task | Aug 2026 | #17 |
| 12 | nf-test suite | Aug 2026 | #16 |
| 11 | hdWGCNA, second attempt | Aug 2026 | #16 |
| 10 | Rename to NetPerturb | Jul 2026 | — |
| 9 | Targets scored in parallel | Jun 2026 | #12 |
| 8 | Multi-target support | May 2026 | — |
| 7 | scRank's own network builder | May 2026 | #9 |
| 6 | hdWGCNA, first attempt, reverted | Apr 2026 | #8, #10 |
| 5 | scTenifoldNet as a second network method | Apr 2026 | #7 |
| 4 | First README | Mar 2026 | #6 |
| 3 | Shared gene universe across cell types | May–Sep 2025 | #3, #4 |
| 2 | Runnable test profile | May 2025 | #2 |
| 1 | Prototype pipeline | Dec 2024 | — |

---

## Wave 21 — Gene set enrichment on the knockout table

**Sep 2026** · branch `sctknk`, not yet on `main`

The knockout track answered "which genes did this move" and stopped there. `GSEA_SCTENIFOLDKNK` takes the answer one step further: for every cell type and knocked-out gene, the genes are ranked by the manifold-alignment distance and that ranked list goes to `fgsea`. This is what scTenifoldKnk's own paper does with its output — *"genes were sorted according to the value of the distance to produce a ranked gene list, which was used as input of the gene set enrichment analysis"* — and the ranking matters: it is the distance itself, uncut, not the genes surviving an FDR threshold.

Three decisions in the ranking are worth recording, because each one is a departure from a default that would have been wrong.

**One-sided scoring.** A distance is never negative, so under fgsea's default `scoreType = "std"` the package warns that every value is positive and then scores a depleted tail that cannot exist. `"pos"` asks the question actually being asked. The consequence reaches the reader, so the report says it outright: `NES` is always positive here, and a set the knockout left alone fails to enrich rather than scoring negative.

**The knocked-out gene is dropped from its own ranking.** It is first by construction — zeroing its edges is what the distance measures — so leaving it in hands a guaranteed top-of-list hit to every gene set that annotates it, which is exactly the set a reader would most want to believe. This is the one place the ranking departs from the paper's literal "sort all genes", and it is the same call wave 20 already made when it dropped the gene from its own network ring and kept it in the table.

**The universe is `gene4use`.** Sets are intersected with the few thousand highly-variable, TF and drug-target genes `DOWNSAMPLE` selects, not the transcriptome, so they arrive smaller than their nominal size and enrichment is relative to the genes the run modelled. That is why `--gsea_min_size` defaults to 10 rather than fgsea's conventional 15, and why the report says so where someone reading a result will see it.

Gene sets come from a `--gsea_gmt` file and are never fetched, so the step runs on a node with no internet; without one it does not run and the report renders a note rather than an empty section. The commonest way it comes back empty is a species mismatch — MSigDB's mouse collections carry MGI symbols and its human ones HGNC — so the script prints the symbols on each side, counts the overlap, and says how many would match if case were folded. It does not fold: `Myc` to `MYC` is right and `Trp53` to `TP53` is not, and a half-working match is worse than none, because the half that lands looks like a real result.

`fgsea()` dispatches to an adaptive sampler, so the script seeds it and runs serial. Nextflow's cache hides that non-determinism on `-resume`, which is an argument for pinning it rather than against: without a seed nobody re-running by hand can reproduce the numbers in the published table. Serial also keeps fgsea's default `nproc = 0` from handing the work to every core on the node, the same oversubscription `conf/modules.config` already documents for `SCTENIFOLDKNK`.

The container is the first on this track to be referenced by a registry tag behind `params.gsea_container` rather than an absolute `.sif` path, which is what wave 19 recorded as the reason the knockout track runs on one machine only.

Still open, and unchanged by this wave: `--sctknk` on its own runs no `REPORT`, so in that mode the enrichment table is published with nothing rendering it — the same gap the DR table has. Closing it means `report.qmd` surviving without `perbscore_file`, which is the whole of its first four sections.

## Wave 20 — The knockout track running in parallel

**Sep 2026** · branch `sctknk`, not yet on `main`

Wave 19 hung scTenifoldKnk off `--network`, which forced a choice between a knockout and a perturbation score. This wave separates them. `--network` selects at most one rank-score method (`genie3`, `scrank`, `hdwgcna`) and runs it through `RANK_SCORE`/`MERGE`/`REPORT`; `--sctknk` is a boolean that switches the knockout track on beside it. Either runs alone, both run in parallel and share `DOWNSAMPLE`, and passing neither is an error rather than a run with nothing in it. `sctknk` is no longer a `--network` value, and a test asserts that.

**One knockout per gene.** scTenifoldKnk's `gKO` takes a single gene symbol, so wave 19 skipped `;`-joined targets outright. `main.nf` now flattens the target file to the individual genes it names, deduplicated across the whole file, and knocks each out on its own: `Stfa1;Mpo` yields a `Stfa1` result and an `Mpo` result. This is the same split `downsample_and_split.R` already does to build `gene4use`. The scoring track is untouched — `RANK_SCORE` still receives `Stfa1;Mpo` intact and scores it as one joint perturbation — so the two tracks deliberately read a combined line differently. The skip branch stays in `sctenifoldknk.R` for anyone running the script by hand. Tasks are now one per (cell type, gene) pair, and those per-pair tables stopped being published: there is one per pair and they are purely intermediate, so only the merged `sctknk/sctenifoldknk_all_targets.txt` is kept.

**Surviving real data.** Three guards went into `sctenifoldknk.R` after runs died inside the package. Cells with no counts across `gene4use` are dropped, because `cpmNormalization` is a bare `t(t(X)/colSums(X))` with no guard for a zero total: such a cell becomes a column of `NaN`, and `makeNetworks` then subsets with `NA` and dies before the first network is built — and `gene4use` is only a few hundred genes, so a transcriptome-wide healthy cell can easily carry zero counts across just those. Cell types left with fewer genes than `nc_nComp` or fewer than three cells are skipped, since `pcNet` requires `nc_nComp < nGenes`. A pair that still fails is skipped rather than aborted on, because one failure would otherwise take a whole sweep of pairs down with it. All three write the same header-only table, so a skip reaches the merge and the report looking exactly like a run that found nothing.

**The report changed shape.** The plotly star figure wave 16 built for target connections is now drawn for knockouts instead, and the connections section left the report entirely — that table is still published, just no longer rendered. Each view is one knocked-out gene in one cell type, with the genes its knockout moved on a ring around it: capped at `--sctknk_top_genes` (25) of those clearing FDR < 0.05, ranked on `log2FC`, a dropdown switching between combinations, and an uncapped queryable table of every significant gene below it. The knocked-out gene is dropped from the ring, where it would be a spoke back to itself, and kept in the table, where it reads as confirmation the knockout took.

What `log2FC` means here is worth recording, because it is not what the name suggests: `dRegulation()` reports `FC` as a gene's squared manifold-alignment distance over the mean squared distance of that run. It is a ratio to the run's average rather than a differential-expression fold change, and it carries no direction — a gene that moved a long way scores high whichever way it moved. That is why the figure sizes nodes by it but gives no edge a sign or a colour to read, and why the palette stays neutral.

The stub suite grew with the track: the workflow-level tests now also cover `--sctknk` running beside `--network`, running alone, `sctknk` no longer being a `--network` value, and the error when neither is passed. Twenty tests in the default suite, plus the opt-in end-to-end run.

The metro map caught up in the same wave, and one thing it had never shown came out in the process: the leg carrying the knockout table into the report was written as a Mermaid dashed `-.->`, and nf-metro's grammar takes only `-->`, `---` and `==>`. Every render since had dropped that leg with a warning, so the blue line stopped at its own table and never reached the report on the picture. It is an ordinary segment now — the arrow kind is discarded by the parser anyway, so there was no dashed style to keep — with the conditionality left to the comment beside it.

## Wave 19 — scTenifoldKnk replaces scTenifoldNet

**Sep 2026** · branch `sctknk`, not yet on `main`

scTenifoldNet exists to compare two conditions, a control matrix against a knockout one. This pipeline only ever fed it wild-type data and read a target's edge weight off the single network it built, the same way it treats GENIE3, scRank and hdWGCNA, so the comparison the package is built around never actually happened here.

scTenifoldKnk fits that single-matrix input natively: given a target gene it builds the wild-type network, zeros the gene's outgoing edges itself, and compares the two by manifold alignment, returning a genome-wide table of differentially-regulated genes with FDR. That is an actual in-silico knockout rather than a network built once and read from, and it needs no adaptation to work from one input matrix.

The cost is that it does not fit the `--network` branch. The other three methods build one target-agnostic network per cell type and let `RANK_SCORE` score every target against it afterwards; scTenifoldKnk's knockout is target-specific by construction, so it runs once per (cell type, target) pair and its table never enters `RANK_SCORE`. It got `MERGE_SCTENIFOLDKNK`, a concatenation step of its own, and a section in the report. `assets/NO_SCTKNK_TABLE`, a header-only sentinel, stands in when the track was not run, so `REPORT` always has a path to stage and `report.qmd` renders a "not run" note instead of a table — the same shape a real run that found nothing produces, so both take the same branch.

This wave wired it as `--network sctknk` to run it alone plus `--run_sctknk` to run it alongside another method: two switches for one track, replaced in wave 20 by the single `--sctknk` boolean.

There is no published image. The container is built locally from `container/sctenifoldknk/Dockerfile` (`satijalab/seurat:5.4.0`, scTenifoldNet from GitHub since scTenifoldKnk wraps its network construction, then scTenifoldKnk from CRAN) and the process points straight at a `.sif` cached under an absolute path on one machine, which is the main thing stopping this track running anywhere else as-is. `bin/sctenifoldnet.R`, its module, its test and its container recipe were deleted, and `test_ocasio` was repointed from `sctnet` to `genie3`.

## Wave 18 — Running it on real data

**Sep 2026** · PR #17

Six fixes from running the pipeline on real datasets rather than the test objects. They share a shape worth recording: each one failed silently, either producing an empty but plausible result or quietly not doing what the config said.

- **Seurat v5 assays.** The container pairs Seurat v4 with SeuratObject v5, and v4's entry points resolve assays with `FilterObjects(classes.keep = "Assay")`, which `Assay5` does not inherit from — so an object saved with a v5 assay is invisible to all of them. `as_v3_assay()` converts it once in `downsample_and_split.R`, and again in `rank_score.R`, which re-reads the original object rather than the converted splits. The conversion also has to go in beside the original under a temporary name and be renamed afterwards, since an `Assay5` cannot be replaced under its own name.
- **hdWGCNA without a PCA.** `MetacellsByGroups` builds its KNN graph from a `pca` reduction, which an object that never went through standard Seurat processing does not carry. One is now computed on a scratch copy, over the same gene universe the network is built on, so the object's own layers reach `NormalizeMetacells` untouched.
- **scRank's network lookup.** `obj@net` is keyed by the raw value of the identity column while `cell_type` comes from the already-sanitised file name, so any label with a space or a comma never matched and fell through to the zero-matrix fallback — for every affected cell type, without a word in the log. Each invocation holds exactly one cell type, so the network is now taken by position, with the length guarded because `Constr_net` can genuinely return an empty list for a population too small to build anything from.
- **`cpus = {$params.n_cores}`.** The stray `$` meant this was not a closure, so `GENIE3`, `SCTENIFOLDNET` and `SCRANK` silently ran with the default of one cpu. The local executor then had no reason to serialise them: every cell type's task launched at once, each spawning up to `n_cores` workers of its own, against 64 physical cores.
- **`RANK_SCORE` declared no cpus at all** while calling `rank_celltype(n.core = 4)`, so its forked workers stacked the same way. It now declares `cpus = 4` to match, with memory raised to 64 GB, which a target combining two genes needed once several targets' tasks overlapped.
- **Empty results cached as successes.** When `rank_celltype` threw for every target in an invocation — usually a worker killed under that contention, not an absence of signal — the script printed a message and exited 0 with a headerless empty table. Nextflow cached it as a success, and `MERGE`'s plain `head -n 1` could then pick that file as the header for the whole merged table. Columns are now written even when nothing succeeds, and the script exits non-zero when every target failed, so the task is retried rather than the gap cached.

The same PR caught the config up with newer Nextflow: `params.binding` moved out of the `params` block, where it collides with `Script.binding`; the trace timestamp became `params.trace_report_suffix` instead of a script-level `def`; and the removed `docker.userEmulation` became an explicit `runOptions`.

## Wave 17 — Bigger test data, and DOWNSAMPLE without scRank

**Aug 2026** · PR #17

`downsample_and_split.R` was building a whole scRank object with `CreateScRank` for one thing: the `gene4use` list hanging off it. It now reproduces that selection directly — highly variable genes capped at 2000, every TF and drug target scRank knows for the species, and the requested targets, minus mitochondrial and ribosomal genes, then intersected with what the object actually has. Targets are added back after the MT/RP filter so one can never be dropped for looking ribosomal, and the filter itself is guarded, since `-x` on an empty index vector empties the whole set instead of removing nothing. Variable features already on the object are reused rather than recomputed. scRank is still the scoring engine; it is just no longer a dependency of the first step.

`test_ocasio` moved to a Zenodo-hosted `.rds` and had its `--column` corrected to `annotation`. `test_vangalen` was added beside it: van Galen 2019 human AML bone marrow, 21 cell identities, scored with hdWGCNA. Its URL has to use the `ndownloader.figshare.com` host, because the `figshare.com/ndownloader/...` form answers 202 with an empty body and would stage a zero-byte object.

## Wave 16 — The strongest edges a target holds

**Aug 2026** · PR #17

`RANK_SCORE` gained a second output. For each cell type it reads the target's row of that cell type's network — `abs(net[target, ]) > 0`, the same neighbourhood `scRank::init_mod()` reads before modularising a single subnetwork — and keeps the `--top_connections` (15) edges with the largest absolute weight, sign included. Ranking on magnitude is deliberate: a strong repressive edge is as informative as a strong activating one, and the direction survives in the weight itself. A zero weight is treated as the absence of an edge rather than a weak one. This runs before `rank_celltype`, so the connections survive a target whose scoring fails.

`MERGE` now concatenates two one-header-per-target families instead of one, through a shared bash function, emitting `top_connections_all_targets.txt` beside the scores. The report drew these as a plotly star per target and cell type, which is what put `plotly` in the report image; `report.qmd` treats its absence as a missing figure rather than an error, so a container built before that section existed still renders the rest. Wave 20 kept the figure and pointed it at knockouts instead, leaving the table published but no longer rendered.

## Wave 15 — What the report actually shows

**Aug 2026** · PR #17

Three additions that turned the wave 13 draft into something readable.

**The bimodal cut.** Pooled perturbation scores come out in two modes: near-zero scores from genes the network barely leans on, and the mode that carries the signal. Imputing the low mode up to a bound was tried first, and only stacked it against that bound while leaving it in the picture, so it is cut instead — everything below the `--score_quantile` (0.75) quantile of the pooled log10 scores is dropped from the figures. The table is never cut. The distribution figure draws both panels, all scores and the retained mode, with the cut line over the full set so it can be checked against where the modes actually separate and the parameter moved into the valley between them. A second dashed line marks the strongest 2.5% of the run, fixed in the qmd rather than exposed as a parameter, and those are the scores the heatmap asterisks mark. Red and green being a poor pair for the commonest colour blindness, the two lines differ in dash as well as hue and the legend names both.

**Cell identities.** `DOWNSAMPLE` writes a UMAP of the cells that survive downsampling, coloured by `--column`, and it closes the report. An embedding already on the object is reused, so the figure matches whatever has been published for that dataset; one is computed only when the object carries none. It is a QC figure, so the whole thing is wrapped: a failure leaves a placeholder carrying the reason rather than sinking a run that is otherwise fine.

**Provenance.** The `--network` method is passed through to the report and named in the overview, since scores are only comparable within one inference method — the point wave 11 closed on. The report still renders when it is not supplied, for a qmd knitted by hand against a table someone already has.

## Wave 14 — The metro map

**Aug 2026** · PR #17

`docs/netperturb_metro.mmd` describes the pipeline as a transit map, one line per method, rendered with [nf-metro](https://github.com/seqeralabs/nf-metro) into the SVG, an interactive HTML page, light and dark PNGs, and a print PDF under `docs/images/`. The README picks the PNGs through a `<picture>` element so the map follows the reader's theme; they have to be baked per theme because rasterisers cannot resolve `var()`.

The map is not only a picture. Every station carries a `%%metro process:` mapping to the Nextflow process it stands for, so `nf-metro serve` plus `-with-weblog` overlays a live run onto it, and `nf-metro check-mapping` against a `-preview -with-dag` DAG catches the mappings drifting after the workflow changes. The regeneration commands live in a collapsed block in the README rather than here, since they are something to run rather than something to know.

## Wave 13 — REPORT task

**Aug 2026** · PR #17

Added `REPORT`, a fifth pipeline step that runs after `MERGE` and turns `perbscore_all_targets.txt` into a self-contained Quarto HTML report: a `DT` table filterable by cell type, target and binding, and a `ggplot2` heatmap of every target scored against every cell type. `bin/report.qmd` is a parameterized Quarto document (`perbscore_file` param) rather than an executable `bin/*.R` script, so it is passed into the process as an explicit `path` input and staged under its own name to avoid colliding with the file it's copied to before rendering.

The container started as `rocker/verse:4.4.1`, which already bundles Quarto and tidyverse, with `DT` installed at task runtime from RSPM's prebuilt binaries rather than baking a dedicated image. That did not survive the wave: `container/rquarto/Dockerfile` now builds `diegomscoelho/rquarto:1.5.54` from `rocker/r-ver:4.3.2` with a pinned Quarto CLI and every R package installed from Posit Package Manager binaries, failing the build loudly if any of them is missing. The process also points `HOME` and every XDG directory at the task work dir, because on HPC the container's own `$HOME` is usually read-only and neither Quarto nor the Deno runtime it embeds can create its cache there.

Scoped deliberately to only `perbscore_all_targets.txt` — no braak-stage or macro-cell-type grouping like the ad hoc `AD_scRank_2025` report this was modelled on, since NetPerturb's own output doesn't carry that structure yet.

## Wave 12 — nf-test suite

**Aug 2026** · `hdwgcna_network`

Every process gained a `stub` block, which lets the whole pipeline run without pulling a container or executing any R. On top of that sits a suite of fourteen tests: one per module asserting its output contract, and five at the workflow level asserting the shape of the run.

The contract worth naming is the file name. `rank_score.R` recovers a cell identity with `sub("_weight.*", "")` on the network file, so the `_weight_` separator is an undeclared coupling between four network scripts and the scoring script. Each network module now has a test that would fail if its naming drifted.

`MERGE` is tested for real rather than stubbed, since it is pure bash with no container.

The end to end run is kept out of the default suite. nf-test's `ignore` also blocks running an ignored file by path, so the opt-in run has its own `nf-test.integration.config`. A GitHub Actions workflow runs the stub suite and parses every R script on each push.

## Wave 11 — hdWGCNA, second attempt

**Aug 2026** · `hdwgcna_network`

Reintroduces `--network hdwgcna`, this time adapting the TOM to scRank's assumptions instead of handing it over raw. `bin/hdwgcna.R` applies four transformations after `ConstructNetwork`:

1. **Sign recovery** — each edge is multiplied by the sign of the correlation between the same two genes across metacells, keeping the TOM magnitude but restoring direction of effect. The network is deliberately built as `unsigned`, because a `signed` adjacency drives anti-correlated pairs towards zero and repressive edges would be cut before there is any sign left to recover.
2. **Padding** — genes dropped by hdWGCNA quality control return as zero rows and columns, so every identity is described over the same `gene4use` universe. scRank's `.align_net` refuses networks whose features differ.
3. **Sparsification** — edges below the `--cut_ratio` quantile of absolute weight are cut, since a dense TOM makes every gene a neighbour of every other one and flattens the degree and entropy terms the score is built on.
4. **Rescaling** — weights are divided by the largest absolute weight to span `[-1, 1]`.

Identities too small to aggregate into metacells are skipped with a message rather than failing the run, and are absent from the final ranking.

The normalisation is scoped to hdWGCNA only. `genie3.R` and `sctenifoldnet.R` are untouched, so `perb_score` values are not comparable across `--network` choices.

## Wave 10 — Rename to NetPerturb

**Jul 2026** · direct to `main`

`NF_scRank` became `NetPerturb` in the README. The `manifest` block in `nextflow.config` still carries the old name and homepage.

## Wave 9 — Targets scored in parallel

**Jun 2026** · PR #12 (`paralel-targets`)

Each target became its own task. Targets are read into a Nextflow channel in `main.nf` and fanned out, so the per-target loop inside the R script no longer carries the parallelism.

This split the old `merge_and_downstream` process in two: `RANK_SCORE`, which runs once per target and emits its own table, and `MERGE`, which concatenates them into `perbscore_all_targets.txt`. `--binding` was also added here, exposing scRank's agonist and antagonist perturbation modes, with `antagonist` as the default in `nextflow.config`.

## Wave 8 — Multi-target support

**May 2026** · direct to `main`

Scoring more than one target per run. The target file became a list, entries could combine genes with `;` to be perturbed together as one set, and `downsample_and_split.R` learned to fold all of them into `gene4use`. Targets were still scored in sequence inside a single process.

## Wave 7 — scRank's own network builder

**May 2026** · PR #9 (`scrank_net`)

Added `--network scrank`, which uses scRank's native `Constr_net` rather than an external inference tool, giving a baseline whose output is by definition in the format the scoring step expects. Includes a mock zero matrix for identities where `Constr_net` returns `NULL`, so one empty identity does not fail the run.

Merged after the wave 6 revert despite branching before it.

## Wave 6 — hdWGCNA, first attempt, reverted

**Apr 2026** · PR #8 (`test-hdwgcna`), reverted by PR #10 (`revert-8-test-hdwgcna`)

hdWGCNA was added as a fourth network module and reverted three days later. The module built a network and wrote `GetTOM()` straight to disk, in the same shape the other methods used.

The reason this did not hold up is worth recording, because it is the whole subject of wave 11: a topological overlap matrix does not satisfy the assumptions scRank makes about a network. It is dense where scRank expects roughly 5% of edges to survive, its values are orders of magnitude smaller than the `[-1, 1]` range the manifold alignment assumes, and it is unsigned, so activation and repression are indistinguishable. The revert also took `resources.config` and some `.gitignore` entries with it.

## Wave 5 — scTenifoldNet as a second network method

**Apr 2026** · PR #7 (`sctenifoldnet`)

The wave that turned a single-method pipeline into a multi-method one. `--network` was introduced, the module that had been called `SCRANK` was renamed `GENIE3` to reflect what it actually ran, and `SCTENIFOLDNET` was added beside it with its own container recipe under `container/sctenifoldnet/`.

Also in this wave: the `test_ocasio` profile for a larger and more realistic dataset than the AML test object, a fix for cells with zero counts, and the `bin/*.R` scripts made executable so they resolve on `PATH` inside the containers.

## Wave 4 — First README

**Mar 2026** · PR #6 (`docs-readme`)

Documentation structure, parameters and usage examples, modelled on the layout of `juliaapolonio/Causeway`.

## Wave 3 — Shared gene universe across cell types

**May–Sep 2025** · PR #3 (`genie-patch`), PR #4 (`gene4use`)

The correctness wave. Each cell identity was inferring its network over its own gene set, so the resulting matrices were not comparable and could not be aligned against each other. `gene4use` was computed once in `downsample_and_split.R` and reused by every identity, and the target genes were concatenated onto it correctly so they always survive into the network.

GENIE3 was also switched from `counts` to `data`, running inference on normalised expression instead of raw counts. `bin/heatmap_scrank.R` arrived in this wave as a scratch plotting script, including greying out outliers.

## Wave 2 — Runnable test profile

**May 2025** · PR #2 (`fix-test-profile`)

Made the pipeline start from a clean checkout without local data. Arguments became proper Nextflow paths so staging worked, `.rda` objects were accepted alongside `.rds` so the public AML test object could be used directly, `outdir` got a default, and a quoting bug in `nextflow.config` was fixed.

## Wave 1 — Prototype pipeline

**Dec 2024** · direct to `main`

The first working skeleton: ingest a Seurat object, split it by cell identity, infer a network per identity, score a target. Everything ran as one path with no choice of method, and the scoring step was a single `merge_and_downstream` process doing both the ranking and the consolidation.
