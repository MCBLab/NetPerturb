process MERGE {
  """
  Merge the networks and rank for each target
  """

  label "r_scrank"

  container "${ workflow.containerEngine == 'singularity' ? 'docker://juliaapolonio/scrank:latest':
            'docker.io/juliaapolonio/scrank:latest' }"

  input:
    path obj
    path target
    val species
    val column
    path 'rank_obj/*'

  output:
    path "merged_obj.RDS", emit: merged_obj
    path "perbscore_all_targets*.txt", emit: rank_scores

  when:
  task.ext.when == null || task.ext.when  

  script:
    """
    Rscript ${baseDir}/bin/merge_and_downstream.R ${obj} ${target} ${species} ${column} rank_obj/*
    """
}

