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
    network = params.network

    if( !(network in ['genie3', 'sctknk', 'scrank', 'hdwgcna']) ) {
        error "Invalid --network '${params.network}'. Supported values: genie3, sctknk, scrank or hdwgcna"
    }

    DOWNSAMPLE( obj, target, column, species, n_cells )

    DOWNSAMPLE.out.scrank_obj
    .flatten()
    .set { sc_obj }

    // scTenifoldKnk knocks out the target gene as part of building its
    // network, so unlike the other three methods it is target-specific by
    // construction and runs once per (cell type, target) pair rather than
    // once per cell type. --network sctknk runs it alone (table only, no
    // RANK_SCORE/REPORT); --run_sctknk runs it alongside whichever of the
    // other three methods --network selects, feeding its table into that
    // run's REPORT as well.
    want_sctknk = (network == 'sctknk') || params.run_sctknk

    if( want_sctknk ) {
        sc_obj
        .combine( target_ch )
        .set { sctknk_input }

        SCTENIFOLDKNK( sctknk_input, n_cores )

        MERGE_SCTENIFOLDKNK( SCTENIFOLDKNK.out.dr_table.collect() )
    }

    if( network != 'sctknk' ) {

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
        sctknk_table = want_sctknk
            ? MERGE_SCTENIFOLDKNK.out.merged_dr_table
            : file("${projectDir}/assets/NO_SCTKNK_TABLE")

        REPORT(
            MERGE.out.merged_rank_scores,
            MERGE.out.merged_top_connections,
            DOWNSAMPLE.out.umap,
            file("${projectDir}/bin/report.qmd"),
            network,
            params.score_quantile,
            sctknk_table
        )
    }
}
