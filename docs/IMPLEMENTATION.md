# Implementation History

How NetPerturb was built, grouped into waves of work rather than individual commits. Each wave is one coherent unit of change, usually a pull request or a short run of commits that only makes sense together. Dates are the point the work landed on `main`. Most recent first.

| Wave | Theme | Landed | PR |
|---|---|---|---|
| 13 | REPORT task | Aug 2026 | — |
| 12 | nf-test suite | Aug 2026 | — |
| 11 | hdWGCNA, second attempt | Aug 2026 | — |
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

## Wave 13 — REPORT task

**Aug 2026**

Added `REPORT`, a fifth pipeline step that runs after `MERGE` and turns `perbscore_all_targets.txt` into a self-contained Quarto HTML report: a `DT` table filterable by cell type, target and binding, and a `ggplot2` heatmap of every target scored against every cell type. `bin/report.qmd` is a parameterized Quarto document (`perbscore_file` param) rather than an executable `bin/*.R` script, so it is passed into the process as an explicit `path` input and staged under its own name to avoid colliding with the file it's copied to before rendering.

The container is `rocker/verse:4.4.1`, which already bundles Quarto and tidyverse; `DT` is not part of that image and is installed at task runtime from RSPM's prebuilt binaries rather than baking a dedicated image, since the report is a single lightweight render step. Revisit this if render time or reproducibility becomes a concern.

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
