# FAQ

## Choosing an analysis

### What is the difference between the perturbation score and the knockout?

They answer different questions about the same target.

- The **perturbation score** (`--network`) is one number per identity × target: how much that identity's regulatory network depends on the target. Use it to rank identities, e.g. which clones are most sensitive to inhibiting a drug target.
- The **knockout** (`--sctknk`) removes the target's outgoing edges from the network and measures how far every other gene moves. You get a genome-wide table per identity × target. Use it to see *which genes and processes* depend on the target.

The two tracks can run in the same command, share downsampling, and end up in the same report.

### Which `--network` method should I use?

- `genie3` is the default choice. It is a well-established regression-based method and fully reproducible with `--seed`. It is the slowest when `--n_hvg` is large.
- `scrank` uses scRank's own network builder. It is the natural choice if you want to match scRank results, but it ignores `--seed`.
- `hdwgcna` builds co-expression networks on metacells, which is more robust for sparse data and large identities. Tune `--cut_ratio` for your data.

Scores are only comparable between runs that used the same method.

### Can I perturb several genes at once?

Yes. Join them with `;` on one line of the target file, e.g. `Stfa1;Mpo`. The scoring track scores them as one joint perturbation, and the knockout track zeroes all of them in a single knocked-out network. To get the single-gene results as well, list each gene on its own line too.

## Missing results

### Why is one of my targets missing from the results?

It did not pass the target QC in `DOWNSAMPLE`. A gene that is absent from the object, or has no counts in any retained cell, cannot have edges in any network, so it is removed before any network is built. Check `downsample/target_qc.tsv` or the *Target QC* part of the report's Overview, which list every removed gene and the reason.

Common causes are a symbol in the wrong case or species (`TP53` vs `Trp53`), an alias instead of the official symbol, or a gene filtered out of the object during preprocessing.

### Why is one of my cell identities missing?

It had fewer than `--min_cells` cells after downsampling (default 150). This is logged and flagged in the report's *Cells per identity* figure. Lower `--min_cells` to keep it, but networks inferred from very few cells are unreliable. With `--network hdwgcna`, identities below `--hdwgcna_min_cells` are also skipped, because they are too small to aggregate into metacells.

### Why did `--sctknk` alone not produce a report?

The report is built around the perturbation score table, so it is only rendered when `--network` runs. With `--sctknk` alone, the knockout table and the enrichment are still published in `sctknk/`. Add a `--network` method to get them in a report.

### My knockout found very few significant genes. Is something wrong?

Not necessarily. A target with few outgoing edges in an identity's network moves little when knocked out, so few genes clear FDR < 0.05. Check the target's expression in the report's *Target expression per identity* heatmap. A barely expressed target has little to knock out. Genes below `--sctknk_min_pct` are also not in the knockout network at all; lowering it adds sparser genes, at a steep memory cost.

### GSEA returned no gene sets.

Nearly always the gene set symbols do not match the data's. MSigDB's mouse collections (`*.Mm.symbols.gmt`) use MGI symbols and its human ones (`*.Hs.symbols.gmt`) use HGNC. A human GMT against mouse data matches almost nothing. The task log prints example symbols from both sides and the overlap between them. Also check that `--gsea_min_size` is not too high, because sets are counted after intersecting with the knockout network's genes.

## Interpreting results

### Does a high `FC` in the knockout table mean the gene went up?

No. scTenifoldKnk compares network *structure*, not expression. `FC` is a gene's squared alignment distance divided by the average gene's, so it only says how far the gene moved, not in which direction.

### What does a negative `NES` mean?

Genes are ranked by how far the knockout moved them, centred on the average gene. A positive `NES` means the set is concentrated among the genes that moved most. A negative `NES` means the set sat still while the rest of the network moved. Neither means up- or down-regulation.

### The report's heatmap shows only a few scores.

The figures keep only the high mode of the pooled scores, cut at `--score_quantile` (default `0.75`, the top quarter). The distribution figure draws the cut, so you can check it falls between the two modes and adjust it. The score table and `rank_scores/perbscore_all_targets.txt` are never cut.

## Running the pipeline

### How do I provide the scTenifoldKnk container?

The knockout image is not on a registry yet. Build it from the repository and point both knockout processes at it:

```bash
docker build -t netperturb-sctenifoldknk:1.0 container/sctenifoldknk
# for Singularity/Apptainer:
singularity build sctenifoldknk-v1.0.sif docker-daemon://netperturb-sctenifoldknk:1.0
```

```groovy
// sctknk.config
process {
    withName: 'SCTENIFOLDKNK_BUILD|SCTENIFOLDKNK_KO' {
        container = '/path/to/sctenifoldknk-v1.0.sif'
    }
}
```

```bash
nextflow run MCBLab/NetPerturb -c sctknk.config --sctknk ...
```

The enrichment image is set with `--gsea_container`, which also accepts a local `.sif`.

### My counts are not in the `RNA` assay.

Pass the assay name with `--assay`, e.g. `--assay originalexp` for objects converted from a `SingleCellExperiment`. It is renamed to `RNA` when loaded, because every network method reads that name.

### A task ran out of memory.

Raise that process's memory with a custom config passed with `-c` (see [Usage → Resources](usage.md#resources)). The largest consumers are:

- `SCTENIFOLDKNK_BUILD`: memory grows with the square of the knockout network's gene count. Raising `--sctknk_min_pct` shrinks it.
- `GENIE3`, `SCRANK` and `HDWGCNA`: these grow with `--n_hvg` and `--n_cells`.
- `RANK_SCORE`: if its parallel workers are killed for memory, it retries automatically one identity at a time.

### Can I run on compute nodes without internet?

Yes, once the inputs and images are local. Pass `--obj` and `--gsea_gmt` as local paths and pre-pull the containers, e.g. by setting `NXF_SINGULARITY_CACHEDIR` and running once on a login node. Nothing else is downloaded during the run. scTenifoldKnk's `plotKO()` is drawn without its online annotation for this reason.

### Why did changing `--seed` re-run everything under `-resume`?

The seed decides which cells `DOWNSAMPLE` keeps, and every later step depends on those cells. The same seed and inputs always give the same results, which is what makes `-resume` safe.

### Every task fails with `Unexpected: unbound variable`.

This happens on machines using the uutils reimplementation of coreutils, the default on recent Ubuntu. Nextflow's task wrapper times tasks with `date +%s%3N`, which uutils does not support. It affects every Nextflow pipeline, not only this one. Installing GNU coreutils fixes it.

### Where can I report a problem or ask a question?

Open an issue on [GitHub](https://github.com/MCBLab/NetPerturb/issues) with the command you ran, the `.nextflow.log` and the failing task's `.command.err` from its `work/` directory.
