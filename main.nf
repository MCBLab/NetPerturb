/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { GENIE3 } from "./modules/local/genie3/main.nf"
include { SCTENIFOLDKNK } from "./modules/local/sctenifoldknk/main.nf"
include { SCRANK } from "./modules/local/scrank/main.nf"
include { HDWGCNA } from "./modules/local/hdwgcna/main.nf"
include { DOWNSAMPLE } from "./modules/local/downsample_and_split/main.nf"
include { RANK_SCORE } from "./modules/local/rank_score/main.nf"
include { MERGE } from "./modules/local/merge/main.nf"
include { MERGE_SCTENIFOLDKNK } from "./modules/local/merge_sctenifoldknk/main.nf"
include { GSEA_SCTENIFOLDKNK } from "./modules/local/gsea_sctenifoldknk/main.nf"
include { REPORT } from "./modules/local/report/main.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    obj = file(params.obj)
    column = params.column
    species = params.species
    n_cells = params.n_cells
    n_cores = params.n_cores
    target = file(params.target)
    //create a list of targets from the input file, assuming one target per line
    target_list = target.readLines().collect { it.trim() }.findAll { it } // remove empty lines
    target_ch = Channel.fromList(target_list)

    // scTenifoldKnk knocks out exactly one gene per run, so it cannot take a
    // ';'-joined combined target the way RANK_SCORE does. Rather than skip
    // those lines, the file is flattened to the individual genes it names and
    // each is knocked out on its own, so "Stfa1;Mpo" yields a separate Stfa1
    // result and Mpo result. Deduplicated across the whole file, so a gene
    // that appears both alone and inside a combination is knocked out once.
    // This is the same split downsample_and_split.R does to build gene4use.
    sctknk_target_list = target_list
        .collectMany { line -> line.split(';').collect { gene -> gene.trim() } }
        .findAll { it }
        .unique()
    sctknk_target_ch = Channel.fromList(sctknk_target_list)
    network = params.network
    sctknk = params.sctknk

    // The two tracks are independent. --network picks at most one rank-score
    // method and runs it through RANK_SCORE/MERGE/REPORT; --sctknk adds the
    // scTenifoldKnk track beside it. Either runs on its own, both run in
    // parallel when both are asked for, but at least one has to be asked for.
    if( network != null && !(network in ['genie3', 'scrank', 'hdwgcna']) ) {
        error "Invalid --network '${network}'. Supported values: genie3, scrank or hdwgcna"
    }

    if( network == null && !sctknk ) {
        error "Nothing to run. Pass --network (genie3, scrank or hdwgcna), --sctknk, or both."
    }

    // Enrichment runs on the knockout table, so without that track there is
    // nothing for a GMT to enrich. Warned rather than fatal: every other part
    // of the run is still valid.
    if( params.gsea_gmt && !sctknk ) {
        log.warn "--gsea_gmt was given without --sctknk; there is no knockout table to enrich, so no GSEA will run."
    }

    DOWNSAMPLE( obj, target, column, species, n_cells, params.min_cells )

    DOWNSAMPLE.out.scrank_obj
    .flatten()
    .set { sc_obj }

    // scTenifoldKnk knocks out the target gene as part of building its
    // network, so unlike the rank-score methods it is target-specific by
    // construction and runs once per (cell type, gene) pair rather than
    // once per cell type. Its table never enters RANK_SCORE -- it goes to its
    // own merge, and from there into REPORT when --network is running too.
    if( sctknk ) {
        sc_obj
        .combine( sctknk_target_ch )
        .set { sctknk_input }

        SCTENIFOLDKNK( sctknk_input, n_cores )

        MERGE_SCTENIFOLDKNK( SCTENIFOLDKNK.out.dr_table.collect() )

        // Ranks each pair's genes by how far the knockout moved them and asks
        // which gene sets sit at the top of that ranking -- the analysis the
        // scTenifoldKnk paper runs on its own output. Gene sets come from a
        // file rather than a web service so the step runs on a node with no
        // internet, which means no --gsea_gmt is simply no enrichment.
        if( params.gsea_gmt ) {
            GSEA_SCTENIFOLDKNK(
                MERGE_SCTENIFOLDKNK.out.merged_dr_table,
                file(params.gsea_gmt),
                params.gsea_min_size,
                params.gsea_max_size
            )
        }
    }

    if( network ) {

        if( network == 'genie3' ) {
           GENIE3( sc_obj, n_cores )

            GENIE3.out.rank_obj
            .collect()
            .set { rank_cells  }
        }
        else if( network == 'scrank' ) {
            SCRANK( sc_obj, species, target, column, n_cores )

            SCRANK.out.rank_obj
            .collect()
            .set { rank_cells  }
        }
        else if( network == 'hdwgcna' ) {
            HDWGCNA( sc_obj, column, n_cores, params.cut_ratio, params.hdwgcna_min_cells )

            HDWGCNA.out.rank_obj
            .collect()
            .set { rank_cells  }
        }

        RANK_SCORE( obj, target_ch, species, column, params.binding, params.top_connections, rank_cells )

        MERGE( RANK_SCORE.out.rank_scores.collect(), RANK_SCORE.out.top_connections.collect() )

        // REPORT always takes a scTenifoldKnk table path; when it was not
        // requested this is a sentinel empty file report.qmd recognises and
        // renders as "not run for this session" rather than a real table.
        sctknk_table = sctknk
            ? MERGE_SCTENIFOLDKNK.out.merged_dr_table
            : file("${projectDir}/assets/NO_SCTKNK_TABLE")

        // Same sentinel arrangement for the enrichment table, which has two
        // ways of not existing: the knockout track was not run at all, or it
        // was run without a --gsea_gmt to enrich against.
        gsea_table = (sctknk && params.gsea_gmt)
            ? GSEA_SCTENIFOLDKNK.out.gsea_table
            : file("${projectDir}/assets/NO_GSEA_TABLE")

        REPORT(
            MERGE.out.merged_rank_scores,
            DOWNSAMPLE.out.umap,
            file("${projectDir}/bin/report.qmd"),
            network,
            params.score_quantile,
            sctknk_table,
            params.sctknk_top_genes,
            gsea_table,
            params.gsea_top_terms
        )
    }
}
