/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { GENIE3 } from "./modules/local/genie3/main.nf"
include { SCTENIFOLDNET } from "./modules/local/sctenifoldnet/main.nf"
include { SCRANK } from "./modules/local/scrank/main.nf"
include { DOWNSAMPLE } from "./modules/local/downsample_and_split/main.nf"
include { MERGE } from "./modules/local/merge_and_downstream/main.nf"

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
    network = params.network

    if( !(network in ['genie3', 'sctnet', 'scrank']) ) {
        error "Invalid --network '${params.network}'. Supported values: genie3 or sctnet"
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

    MERGE( obj, target_list.join(";"), species, column, params.binding, rank_cells ) 
}
