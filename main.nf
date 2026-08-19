/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { GENIE3 } from "./modules/local/genie3/main.nf"
include { SCTENIFOLDNET } from "./modules/local/sctenifoldnet/main.nf"
include { SCRANK } from "./modules/local/scrank/main.nf"
include { HDWGCNA } from "./modules/local/hdwgcna/main.nf"
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
    target = file(params.target)
    //create a list of targets from the input file, assuming one target per line
    target_list = target.readLines().collect { it.trim() }.findAll { it } // remove empty lines
    target_ch = Channel.fromList(target_list)
    network = params.network

    if( !(network in ['genie3', 'sctnet', 'scrank', 'hdwgcna']) ) {
        error "Invalid --network '${params.network}'. Supported values: genie3, sctnet, scrank or hdwgcna"
    }
	
    DOWNSAMPLE( obj, target, column, species, n_cells )

    DOWNSAMPLE.out.scrank_obj
    .flatten()
    .set { sc_obj }
    
    if( network == 'genie3' ) {
       GENIE3( sc_obj, n_cores )

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
    else if( network == 'hdwgcna' ) {
        HDWGCNA( sc_obj, column, n_cores, params.cut_ratio, params.hdwgcna_min_cells )

        HDWGCNA.out.rank_obj
        .collect()
        .set { rank_cells  }
    }

    RANK_SCORE( obj, target_ch, species, column, params.binding, rank_cells ) 

    MERGE( RANK_SCORE.out.rank_scores.collect() )

    REPORT( MERGE.out.merged_rank_scores, DOWNSAMPLE.out.umap, file("${projectDir}/bin/report.qmd") )
}
