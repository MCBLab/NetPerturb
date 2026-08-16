process GENIE3 {
 """
  Generates WGN based on scRank obj
 """

  label "r_genie3"

  container "${ workflow.containerEngine == 'singularity' ? 'docker://juliaapolonio/scrank:latest':
            'docker.io/juliaapolonio/scrank:latest' }"

  input:
    path scobj
    val n_cores
    val n_trees

  output:
    path "*.rds", emit: rank_obj

  when:
  task.ext.when == null || task.ext.when

  script:
    """
   #!/bin/bash
   genie3.R ${scobj} ${n_cores} ${n_trees}
    
  stub:
    """
    touch ${scobj.baseName}_weight_GENIE3_100.rds
    """
}
