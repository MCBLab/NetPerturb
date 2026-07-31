process SCRANK {
  """
  Generates WGN based on scRank obj
  """

  label "r_scrank"

  container "${ workflow.containerEngine == 'singularity' ? 'docker://juliaapolonio/scrank:latest':
            'docker.io/juliaapolonio/scrank:latest' }"

  input:
    path scobj
    val species
    path target
    val column
    val n_cores

  output:
    path "*.rds", emit: rank_obj

  when:
  task.ext.when == null || task.ext.when  

  script:
    """
    #!/bin/bash
    scrank.R ${scobj} ${species} ${target} ${column} ${n_cores}
    """

  stub:
    """
    touch ${scobj.baseName}_weight_SCRANK_100.rds
    """
}
