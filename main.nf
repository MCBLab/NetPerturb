/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { GENIE3 } from "./modules/local/genie3/main.nf"
include { SCTENIFOLDNET } from "./modules/local/sctenifoldnet/main.nf"
include { SCRANK } from "./modules/local/scrank/main.nf"
include { DOWNSAMPLE } from "./modules/local/downsample_and_split/main.nf"
include { RANK_SCORE } from "./modules/local/rank_score/main.nf"
include { MERGE } from "./modules/local/merge/main.nf"
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
    n_trees = params.n_trees
    target = file(params.target)
    //create a list of targets from the input file, assuming one target per line
    target_list = target.readLines().collect { it.trim() }.findAll { it } // remove empty lines
    target_ch = Channel.fromList(target_list)
    network = params.network

    if( !(network in ['genie3', 'sctnet', 'scrank']) ) {
        error "Invalid --network '${params.network}'. Supported values: genie3 or sctnet"
    } 
	
    DOWNSAMPLE( obj, target, column, species, n_cells )

    DOWNSAMPLE.out.scrank_obj
    .flatten()
    .set { sc_obj }
    
    if( network == 'genie3' ) {
       GENIE3( sc_obj, n_cores, n_trees )

        GENIE3.out.rank_obj
        .collect()
        .set { rank_cells  }
    }
    else if( network == 'sctnet' ) {
        SCTENIFOLDNET( sc_obj, n_cores )

        SCTENIFOLDNET.out.rank_obj
        .collect()
        .set { rank_cells  }
    }
    else if( network == 'scrank' ) {
        SCRANK( sc_obj, species, target, column, n_cores )

        SCRANK.out.rank_obj
        .collect()
        .set { rank_cells  }
    }

    RANK_SCORE( obj, target_ch, species, column, params.binding, rank_cells )

    MERGE( RANK_SCORE.out.rank_scores.collect() )

    if( params.render_report ) {
        qmd_template = file(params.report_qmd)
        logo         = file(params.report_logo)

        // RANK_SCORE roda em paralelo (um por target) e cada instância reconstrói o
        // mesmo obj@net (target-independente, vem do mesmo rank_cells de entrada) —
        // só obj@para$target muda entre elas. Basta pegar UM merged_obj.RDS pro
        // REPORT, que já itera sobre todos os targets internamente.
        REPORT( RANK_SCORE.out.merged_obj.first(), qmd_template, logo, target_list.join(";"), params.report_n_cores )
    }
}
