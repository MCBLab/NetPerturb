# Implementation History

How NetPerturb was built, grouped into waves of work rather than individual commits. Each wave is one coherent unit of change, usually a pull request or a short run of commits that only makes sense together. Dates are when the work was done; the PR column says how it reached `main`, and names the branch instead when it has not landed yet. Most recent first.

| Wave | Theme | Done | PR |
|---|---|---|---|
| 26 | A documentation website | Sep 2026 | branch `mkdocs` |
| 25 | The report describes the run | Sep 2026 | #21 |
| 24 | One seed for the whole run | Sep 2026 | #21 |
| 23 | The knockout split in two | Sep 2026 | #21 |
| 22 | Targets and identities checked up front | Sep 2026 | #21 |
| 21 | Gene set enrichment on the knockout table | Sep 2026 | #20 |
| 20 | The knockout track running in parallel | Sep 2026 | #20 |
| 19 | scTenifoldKnk replaces scTenifoldNet | Sep 2026 | #20 |
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

## Wave 26 — A documentation website

**Sep 2026** · branch `mkdocs`, not yet on `main`

The README had grown into the only manual: one long page mixing a pipeline description, every parameter, every output file and the development notes. The pipeline now also has a website at `https://mcblab.github.io/NetPerturb/`, laid out like the one [Causeway](https://juliaapolonio.github.io/Causeway/) uses: MkDocs with the readthedocs theme and four pages, **Home**, **Usage**, **Output** and **FAQ**, built from `docs/` by `mkdocs.yml` at the repository root.

**Content.** The site is a reorganisation of the README for someone running the pipeline, not a copy of it:

- *Home* has the overview, the metro map and a quick start.
- *Usage* turns the parameter prose into grouped tables (required, downsampling, scoring, knockout, enrichment, general) and adds profiles, per-process resources and a custom `-c` config example.
- *Output* walks the `--outdir` tree folder by folder, including the report's sections.
- *FAQ* is new. Its questions are the failure modes the pipeline already handles and explains in its code: missing targets and identities, an empty GSEA from a species mismatch, what a signed `NES` or a directionless `FC` means, the unpublished scTenifoldKnk image and how to point the processes at a local `.sif`, and the uutils `date` bug.

The reasoning behind design decisions stays here in `IMPLEMENTATION.md`, which is kept out of the site with `exclude_docs`, along with the metro source, the poster exports and the logo script.

**Theme.** The site follows the MCB Lab site (`mcblab.github.io`), which is Quarto with Bootswatch Sandstone. Sandstone's palette is mapped onto the readthedocs theme in `docs/css/docs.css`: a dark sidebar (`#3e3f3a`), a blue header (`#325d88`), a sand frame (`#dfd7ca`) around cream content (`#f8f5f0`), and the lab's font, Atkinson Hyperlegible, self-hosted in `docs/fonts/`. The home page borrows the lab site's card grid for the two tracks and the report. Readthedocs has no dark mode, so the metro map is shown in its light rendering only. Inline code does not wrap on desktop, so a flag like `--network` is never split at its hyphens; on phones it may wrap, so a long path cannot push the page sideways. Wide tab-separated examples scroll inside their own box instead of being clipped.

**Logo.** A hexagon sticker in the same palette: a small network whose red node reaches its neighbours through dashed edges (the perturbed gene and the edges a knockout removes), solid edges among the rest, and a cell beside it. `docs/img/make_logo.py` generates it. The title is drawn from Atkinson Hyperlegible's glyph outlines rather than as SVG `<text>`, so it renders the same in a browser, in cairosvg and on GitHub, where the font is not installed. The script also writes a text-free mark for the favicons, since the title cannot be read at 16–32 px. The README now shows the logo and links to the site.

**Deployment.** `.github/workflows/docs.yml` runs `mkdocs gh-deploy --strict` on pushes to `main` that touch `docs/` or `mkdocs.yml`, publishing to the `gh-pages` branch. Strict mode fails the build on a broken internal link or anchor, which matters because the FAQ and Usage link into each other's sections. GitHub Pages has to be set to serve that branch once, in the repository settings.

## Wave 25 — The report describes the run

**Sep 2026** · #21

Until this wave the report showed results and nothing about what produced them. A reader could not tell which targets had been quietly left out, how many cells and genes each network had actually been built on, whether a target was even expressed in a given identity, or how the run had been launched. Four sections close that gap, and each reads a small table the pipeline already had the facts for.

**Target QC.** Wave 22's check writes the genes it removed to `downsample/target_qc.tsv`. The Overview lists them as "These genes were not analyzed due to QC checking", with the reason and what happened to their target. A header-only file reads the same as a run where every gene passed.

**Cells per identity.** `DOWNSAMPLE` writes `cell_counts.tsv`: each identity's size in the input, what was kept, and the identities `--min_cells` dropped. The gene columns show what each network was really built on. Every method is offered the same `gene4use`, but a gene with no counts in an identity's retained cells has no edges there: scTenifoldKnk drops it, hdWGCNA pads it with zeros, and GENIE3 and scRank have nothing to regress on. So `genes_expressed` is the honest count and `genes_filtered` is what that identity lost. The figure puts cells and genes side by side as facets on one identity axis. Each panel draws what was available as a back bar and what the network got as a front bar over it. A dropped identity has no used cells, so its input bar carries the flag. The knockout networks choose their own genes (`genes_sctknk`, wave 23), so they get their own panel, shown only when the run built them. An older `cell_counts.tsv` without gene columns falls back to cells only.

**Target expression per identity.** This is a heatmap of each single target gene's mean expression in each kept identity. A `;`-joined target is split into its genes. It is log-normalised from the counts in `DOWNSAMPLE` (per 10,000, then `log1p`, with the library size taken over every gene the way `NormalizeData()` would), rather than read from the `data` layer, which holds raw counts in an object nobody normalised. Only the retained cells count, since those are what the networks saw. Genes keep the order they were given in, and identities are sorted as in the perturbation score heatmap so the two can be read against each other. A target that scores near zero in an identity where it is barely expressed is now visibly that, not a finding.

**Pipeline run.** `REPORT`'s container cannot see Nextflow's `workflow` object, so `main.nf` writes what the report needs as a key/value `run_info.tsv`: run name, session, start time, command line, profile, container engine, versions, revision, directories, and every setting that shapes the result. The start is also passed as an epoch, so elapsed time does not depend on the container's time zone, and the render time is printed in the run's own offset. **Running time per step** reads the execution trace. That file is appended to as each task finishes, so `main.nf` snapshots it only once everything `REPORT` depends on is done, which is every task but `REPORT` itself. It waits three seconds first, because rows are written asynchronously just after a task ends. The path it reads is the one Nextflow is actually writing (`-with-trace <file>` or a `trace.file` from `-c` is followed), with the `nextflow.config` default only as a fallback. A missing trace becomes a one-line note instead of an empty string: `collectFile` emits no file for empty content, and `REPORT` would have waited on it forever. Steps are process names without their tag, so `GENIE3 (3)` and `GENIE3 (7)` are one step. `realtime` is plotted rather than `duration`, which includes queueing. A task restored by `-resume` is `CACHED` and keeps the times of the run that made it. The axis unit is chosen from the longest task, so a run of seconds is not drawn in fractions of an hour.

Every new input is optional in `report.qmd`. A report rendered by hand, or against a run from before these tables existed, renders "not recorded" or omits the section instead of failing. `DOWNSAMPLE`'s `publishDir` pattern now also takes `target_qc.tsv`, `targets_qc.txt`, `cell_counts.tsv` and `target_expression.tsv`, and fixtures for each were added so the `REPORT` module tests render every section under stub.

## Wave 24 — One seed for the whole run

**Sep 2026** · #21

`--seed` (default `1`) now reaches every step that draws random numbers: which cells `DOWNSAMPLE` keeps and the UMAP it computes when the object has none, GENIE3's random forests (through doRNG, so any `--n_cores` gives the same result), hdWGCNA's metacell sampling and WGCNA's `randomSeed`, scTenifoldKnk's bootstrap networks, tensor decomposition and manifold alignment, and fgsea's sampler (seeded since wave 21, now from the same parameter). The default is `1` because scTenifoldKnk hardcodes that seed, so with defaults the knockout track matches `scTenifoldKnk()` run by hand. The exception is scRank's `Constr_net`, which seeds itself with `1` internally and ignores `--seed`, so `--network scrank` networks do not change with it. This is documented, not worked around. Changing the seed changes which cells are kept, so every step reruns, even under `-resume`.

## Wave 23 — The knockout split in two

**Sep 2026** · #21

**Build once, knock out many times.** `scTenifoldKnk()` builds the wild-type network and knocks one gene out of it in a single call. With one call per (cell type, gene), a run with *n* targets rebuilt the same network *n* times. The network does not depend on the target and is seeded, so every rebuild came out identical, and the build (bootstrap networks plus tensor decomposition) is the expensive half. `SCTENIFOLDKNK` became two processes. `SCTENIFOLDKNK_BUILD` runs once per cell type and saves the wild-type network as `<cell type>_sctknk_wt.rds`. `SCTENIFOLDKNK_KO` combines each network with each target and runs only the knockout, manifold alignment and differential regulation, in parallel. Both steps call scTenifoldKnk's own functions in its own order, with its defaults and seeding.

**Joint knockouts, reversing wave 20.** Once the knockout is a separate step, zeroing several genes' edges in one network is as simple as zeroing one. The `main.nf` flattening from wave 20 is gone. The knockout track now takes `DOWNSAMPLE`'s passing targets line for line, as `RANK_SCORE` does, and `Stfa1;Mpo` is one joint knockout under its own name in the merged table. To get the single-gene knockouts as well, list the genes on their own lines. A gene of a combined target that is not expressed in a cell type is left out of that cell type's knockout, and the rest is knocked out without it.

**Its own genes, not `gene4use`.** The network used to be built on `gene4use`, the HVG, TF and drug-target set `DOWNSAMPLE` picks for the rank-score methods. `scTenifoldKnk()` instead keeps every gene detected in more than `qc_minPCT` (5%) of cells. Most of `gene4use` falls below that in a given cell type, and most of what clears it is not in `gene4use`, so the knockout was being modelled on sparse genes while missing the well-measured ones. The knockout track never feeds `RANK_SCORE`, so nothing required a shared universe. `SCTENIFOLDKNK_BUILD` now applies scTenifoldKnk's own filter over every gene in the object, exposed as `--sctknk_min_pct`. It makes one deliberate departure: targets are kept below the threshold, where `scTenifoldKnk()` would refuse to knock them out, so every target gets a knockout wherever it has any counts. Cell-level QC is still skipped, since `DOWNSAMPLE` has already chosen the cells. `DOWNSAMPLE` counts the resulting genes per identity (`genes_sctknk`) by the same rule, for the report. This also moved the GSEA universe: enrichment is now relative to the knockout network's several thousand genes, which is why the `--gsea_min_size` rationale in the README was rewritten.

**Memory.** The tensor decomposition's memory grows with the square of the gene count. A 7,000-gene cell type peaked near 15 GB, so `SCTENIFOLDKNK_BUILD` asks for 24 GB. The alignment in `SCTENIFOLDKNK_KO` is an eigendecomposition of a 2n x 2n matrix and keeps 16 GB. The image caps R's heap at 8 GB (`R_MAX_VSIZE` in its Dockerfile), below what the build needs, so both modules export `R_MAX_VSIZE` from `task.memory` and R gets what the task was given.

**The differential-regulation statistic is scTenifoldKnk 1.0.3's.** `FC` is a gene's squared alignment distance over the mean squared distance, tested against a chi-square with one degree of freedom. In 1.0.3 that mean left out the knocked-out genes. 1.1, the current CRAN release, dropped the `gKO` argument and averages over every gene, and the knocked-out gene's distance is enormous by construction. On a 5,600-gene network one knocked-out gene held about 98% of the summed `FC`, which shrank every other gene's `FC` about fifty-fold and left one gene below FDR 0.05. `sctenifoldknk_ko.R` now computes the 1.0.3 statistic itself as `d_regulation()`, excluding the knocked-out genes by name and leaving them out of the FDR correction, since they are not among the genes tested. The same knockout then gave 12 genes below FDR 0.05. As a consequence, the knocked-out genes are now dropped in `SCTENIFOLDKNK_KO` rather than downstream. The equivalent filters in `report.qmd` and `gsea_sctenifoldknk.R` remain only for tables from older runs.

**`plotKO()` figures.** `--sctknk_plot` also draws scTenifoldKnk's own network plot per (cell type, target) into `sctknk/plots/`. It is drawn with `annotate = FALSE`, because annotation queries Enrichr against human-only libraries, which a node without internet cannot reach. A plot that fails to draw is logged and skipped without touching the table. It is off by default because it adds one file per pair. That `publishDir` is also switched off when the flag is off: Nextflow creates a publish directory even when its pattern matches nothing, and every `--sctknk` run was leaving an empty `sctknk/plots/` behind. Two workflow tests that count what `sctknk/` holds caught this.

## Wave 22 — Targets and identities checked up front

**Sep 2026** · #21

Runs on new datasets kept dying deep inside a network or scoring step, over a fact that was knowable at the start: a target the object did not have, an identity too small to model, or counts under an assay other than `RNA`. This wave moves those checks into `DOWNSAMPLE`, so later steps only receive inputs they can use, and makes `RANK_SCORE` tolerate the cases that can only be discovered once a network exists.

**Small identities are dropped before the split.** With `--min_cells` (default 150; `0` keeps everything), an identity with fewer cells is removed before it is written out, so no network is inferred from it and it appears in no table or figure. It is compared against `min(identity size, --n_cells)`, what survives downsampling, because that is what every network receives and what `--hdwgcna_min_cells` is measured against. Setting `--n_cells` below `--min_cells` drops every identity, and the error says so explicitly.

**Targets are checked against the retained cells.** A gene absent from the object, or with zero counts in every retained cell, has no edges in any network. It is now removed in `DOWNSAMPLE` instead of failing a later step. A combined target loses only its failing gene (`Stfa1;Mpo` with `Mpo` absent is analysed as `Stfa1`) and is dropped only when none of it is left. The run aborts if no target passes. The survivors go to `targets_qc.txt`, which `main.nf` reads in place of `--target` for both tracks and for `SCRANK`. The removed genes go to `target_qc.tsv` with a reason and an action, which the report shows (wave 25).

**`--assay` and `--n_hvg`.** Every network method reads an assay named `RNA`, so `--assay` names the one holding the counts and it is renamed to `RNA` on load, replacing any existing `RNA` assay. This is for objects converted from a `SingleCellExperiment`, where the counts sit under `originalexp`. `--n_hvg` (default 2000) sets how many highly variable genes go into `gene4use`. Variable features already stored on the object are reused when there are at least that many; otherwise they are computed with `vst` on the downsampled cells.

**`RANK_SCORE` stops failing on missing genes.** Three failures came from scRank rather than from anything wrong with the target:

- `CreateScRank()` stops when a target gene is not in the expression profile. The genes are now checked first and absent or all-zero ones are left out. Such a gene has no edges, so it changes nothing about the score, and a combined target keeps its name. If no gene of a target survives, both of its tables are written with headers and no rows, and the run carries on. Every other `CreateScRank` error still aborts, because a cached empty table would hide a real failure.
- `rank_celltype()` zeroes the target's row in every cell type's network, so a gene missing from even one network stopped it with "Drug target gene is not in the network". Dropping such genes, the first fix, lost the whole target from the score table whenever nothing was left, including for the cell types whose networks did carry it. A gene missing from a network has no edges there, which is exactly what knocking it out would leave, so it is now added to that network as an isolated gene (a zero row and column) and the target is ranked in every cell type. A cell type with none of the target's genes scores 0 outright. A gene in no network at all leaves the target, and a target with no gene in any network gets a score of 0 for every cell type instead of no rows.
- `rank_celltype()` forks one worker per cell type, each building dense 2n x 2n Laplacians. When the kernel kills workers for memory, `mclapply` only warns ("did not deliver results") and the ranking dies later with an unrelated-looking `seq_len` error. `rank_celltype_safe()` watches for that warning and reruns serially, holding one Laplacian at a time. Any other error is re-raised, since it would fail the same way serially.

Target lists for the datasets these runs were made on were added under `testdata/`: Kang, and expanded van Galen and Ocasio lists.

## Wave 21 — Gene set enrichment on the knockout table

**Sep 2026** · #20

The knockout track answered "which genes did this move" and stopped there. `GSEA_SCTENIFOLDKNK` takes the answer one step further: for every cell type and knocked-out gene, the genes are ranked by `log2FC` — the log2 of the squared manifold-alignment distance over the run's mean squared distance — and that ranked list goes to `fgsea`. This is what scTenifoldKnk's own paper does with its output — *"genes were sorted according to the value of the distance to produce a ranked gene list, which was used as input of the gene set enrichment analysis"* — and the ordering matters: it is distance itself, uncut, not the genes surviving an FDR threshold. `log2FC` is a monotone transform of that distance, so it produces the paper's order while placing the average gene at zero.

Three decisions in the ranking are worth recording, because each one is a departure from a default that would have been wrong.

**Two-sided scoring, and what the sign is not.** The ranking statistic is `log2FC` rather than the raw distance for one reason: a distance is never negative, so it gives fgsea a one-tailed list, and the whole enrichment collapses onto "is this set among the genes that moved". `log2FC` is the same ordering centred on the run's average gene, so it straddles zero, `scoreType` is fgsea's default `"std"`, and `NES` comes back signed — a set can now be reported as concentrated at the quiet end of the ranking as well as the moved end. The sign wants reading carefully and both the report and the module docstring say so: scTenifoldKnk compares network structure, not expression, so a negative `NES` means "this set sat still while the rest of the network moved", never "this set went down". (Before Sep 2026 this ran one-sided with `scoreType = "pos"` and `NES` was positive by construction.)

**The knocked-out gene is dropped from its own ranking.** It is first by construction — zeroing its edges is what the distance measures — so leaving it in hands a guaranteed top-of-list hit to every gene set that annotates it, which is exactly the set a reader would most want to believe. This is the one place the ranking departs from the paper's literal "sort all genes", and it is the same call wave 20 already made when it dropped the gene from its own network ring and kept it in the table.

**The universe is `gene4use`.** Sets are intersected with the few thousand highly-variable, TF and drug-target genes `DOWNSAMPLE` selects, not the transcriptome, so they arrive smaller than their nominal size and enrichment is relative to the genes the run modelled. That is why `--gsea_min_size` defaults to 10 rather than fgsea's conventional 15, and why the report says so where someone reading a result will see it. (Superseded in wave 23: the knockout network now picks its own genes, so the universe became every gene detected in more than `--sctknk_min_pct` of a cell type's cells.)

**Section 7 is a dot plot, not a bar.** A bar of `NES` per set was the first shape, and it stopped working the moment the score went signed: bars growing left and right of a zero line spend most of the figure on the axis and read as a diverging quantity, which `NES` is not — the sign is which end of the ranking a set sits at, not a direction of change. The dot plot puts one mark per set at its `NES`, sizes it by how many genes of it are in the leading edge and shades it by significance, so significance, breadth and position are three channels rather than one length plus a hover. Three consequences worth recording. The size, colour and x scales are computed once across every combination rather than per view, because the dropdown swaps traces in place and a scale that moved with it would make two identical-looking dots mean different things one click apart. The colour statistic is `-log(p.adj + 1)`, not `-log10(p.adj)`: the pseudocount bounds it in `[-log(2), 0)`, which is what makes `fgsea`'s 1e-50 p-value floor harmless — under `-log10` a single floored hit stretched the ramp to 50 and left every other dot the same pale shade, which an earlier 95th-percentile cap existed only to paper over. The trade is deliberate and worth knowing when reading a figure: `log` is close to linear near 1, so the statistic separates the marginal sets from the rest and treats everything below `p.adj` of about 0.01 as equally significant. Its range is therefore derived rather than measured — `[-log(1 + gsea_fdr), 0]`, the whole span the FDR cut can produce — so a shade means the same thing in a report with three significant sets as in one with three hundred. The ramp runs blue through a neutral grey to red, the convention an enrichment dot plot is read with, rather than the one-hue sequential ramp the rest of the report uses; the grey midpoint is taken a few steps down from the surface, since the near-white one the palette documents would put a dot where a reader sees nothing.

**The table below it lost three columns, and cuts a fourth.** `ES`, `pval` and `leadingEdge` are gone from the rendered table — `NES` is the size-normalised form of `ES` and the only one comparable across sets, `padj` is what the section's cut is made on and `pval` sorts identically within a combination, and the `;`-joined leading edge was a wide column of gene names that needed a JavaScript truncating renderer to fit. What the figure actually uses from it, `leading_edge_n`, stays. All three columns are still in `sctknk/gsea_all_targets.txt`, which is the file anyone doing something further with the enrichment reads anyway. `pathway` is rendered cut at 50 characters -- MSigDB names run past 120, which widens the column far enough to push `p.adj` off a laptop screen. The cut is in the DT renderer rather than in the data, so the whole name stays in the cell: the column filter still matches the part that is not shown, and the tooltip shows it in full. 50 is also what the figure's axis labels use, so a name is not cut at two different places in one section.

Gene sets come from a `--gsea_gmt` file and are never fetched, so the step runs on a node with no internet; without one it does not run and the report renders a note rather than an empty section. The commonest way it comes back empty is a species mismatch — MSigDB's mouse collections carry MGI symbols and its human ones HGNC — so the script prints the symbols on each side, counts the overlap, and says how many would match if case were folded. It does not fold: `Myc` to `MYC` is right and `Trp53` to `TP53` is not, and a half-working match is worse than none, because the half that lands looks like a real result.

`fgsea()` dispatches to an adaptive sampler, so the script seeds it and runs serial. Nextflow's cache hides that non-determinism on `-resume`, which is an argument for pinning it rather than against: without a seed nobody re-running by hand can reproduce the numbers in the published table. Serial also keeps fgsea's default `nproc = 0` from handing the work to every core on the node, the same oversubscription `conf/modules.config` already documents for `SCTENIFOLDKNK`.

The container is the first on this track to be referenced by a registry tag behind `params.gsea_container` rather than an absolute `.sif` path, which is what wave 19 recorded as the reason the knockout track runs on one machine only.

Still open, and unchanged by this wave: `--sctknk` on its own runs no `REPORT`, so in that mode the enrichment table is published with nothing rendering it — the same gap the DR table has. Closing it means `report.qmd` surviving without `perbscore_file`, which is the whole of its first four sections.

## Wave 20 — The knockout track running in parallel

**Sep 2026** · #20

Wave 19 hung scTenifoldKnk off `--network`, which forced a choice between a knockout and a perturbation score. This wave separates them. `--network` selects at most one rank-score method (`genie3`, `scrank`, `hdwgcna`) and runs it through `RANK_SCORE`/`MERGE`/`REPORT`; `--sctknk` is a boolean that switches the knockout track on beside it. Either runs alone, both run in parallel and share `DOWNSAMPLE`, and passing neither is an error rather than a run with nothing in it. `sctknk` is no longer a `--network` value, and a test asserts that.

**One knockout per gene.** scTenifoldKnk's `gKO` takes a single gene symbol, so wave 19 skipped `;`-joined targets outright. `main.nf` now flattens the target file to the individual genes it names, deduplicated across the whole file, and knocks each out on its own: `Stfa1;Mpo` yields a `Stfa1` result and an `Mpo` result. This is the same split `downsample_and_split.R` already does to build `gene4use`. The scoring track is untouched — `RANK_SCORE` still receives `Stfa1;Mpo` intact and scores it as one joint perturbation — so the two tracks deliberately read a combined line differently. The skip branch stays in `sctenifoldknk.R` for anyone running the script by hand. Tasks are now one per (cell type, gene) pair, and those per-pair tables stopped being published: there is one per pair and they are purely intermediate, so only the merged `sctknk/sctenifoldknk_all_targets.txt` is kept. (Reversed in wave 23: once the network build was split from the knockout, a `;`-joined line became one joint knockout again, and the two tracks read a combined line the same way.)

**Surviving real data.** Three guards went into `sctenifoldknk.R` after runs died inside the package. Cells with no counts across `gene4use` are dropped, because `cpmNormalization` is a bare `t(t(X)/colSums(X))` with no guard for a zero total: such a cell becomes a column of `NaN`, and `makeNetworks` then subsets with `NA` and dies before the first network is built — and `gene4use` is only a few hundred genes, so a transcriptome-wide healthy cell can easily carry zero counts across just those. Cell types left with fewer genes than `nc_nComp` or fewer than three cells are skipped, since `pcNet` requires `nc_nComp < nGenes`. A pair that still fails is skipped rather than aborted on, because one failure would otherwise take a whole sweep of pairs down with it. All three write the same header-only table, so a skip reaches the merge and the report looking exactly like a run that found nothing.

**The report changed shape.** The plotly star figure wave 16 built for target connections is now drawn for knockouts instead, and the connections section left the report entirely — that table is still published, just no longer rendered. Each view is one knocked-out gene in one cell type, with the genes its knockout moved on a ring around it: capped at `--sctknk_top_genes` (25) of those clearing FDR < 0.05, ranked on `log2FC`, a dropdown switching between combinations, and an uncapped queryable table of every significant gene below it. The knocked-out gene is dropped from the ring, where it would be a spoke back to itself, and kept in the table, where it reads as confirmation the knockout took.

What `log2FC` means here is worth recording, because it is not what the name suggests: `dRegulation()` reports `FC` as a gene's squared manifold-alignment distance over the mean squared distance of that run. It is a ratio to the run's average rather than a differential-expression fold change, and it carries no direction — a gene that moved a long way scores high whichever way it moved. That is why the figure sizes nodes by it but gives no edge a sign or a colour to read, and why the palette stays neutral.

The stub suite grew with the track: the workflow-level tests now also cover `--sctknk` running beside `--network`, running alone, `sctknk` no longer being a `--network` value, and the error when neither is passed. Twenty tests in the default suite, plus the opt-in end-to-end run.

The metro map caught up in the same wave, and one thing it had never shown came out in the process: the leg carrying the knockout table into the report was written as a Mermaid dashed `-.->`, and nf-metro's grammar takes only `-->`, `---` and `==>`. Every render since had dropped that leg with a warning, so the blue line stopped at its own table and never reached the report on the picture. It is an ordinary segment now — the arrow kind is discarded by the parser anyway, so there was no dashed style to keep — with the conditionality left to the comment beside it.

## Wave 19 — scTenifoldKnk replaces scTenifoldNet

**Sep 2026** · #20

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
