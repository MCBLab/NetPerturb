process RANK_SCORE {
  """
  Merge the networks and rank for each target
  """

  label "r_scrank"
  tag "$target"

  container "${ workflow.containerEngine == 'singularity' ? 'docker://juliaapolonio/scrank:latest':
            'docker.io/juliaapolonio/scrank:latest' }"

  input:
    path obj
    val target
    val species
    val column
    val binding
    path(rank_obj)

  output:
    path "perbscore_all_targets*.txt", emit: rank_scores

  when:
  task.ext.when == null || task.ext.when  

  script:
    """
   #!/bin/bash
    rank_score.R ${obj} "${target}" ${species} ${column} ${binding} ${rank_obj}
    """

  stub:
    //mirrors the target_id sanitising rank_score.R applies to the output name
    def target_id = target.replaceAll(/[^A-Za-z0-9_.-]+/, "_")
    """
    printf 'cell_type\\ttarget\\tbinding\\tperb_score\\n' > perbscore_all_targets.${target_id}.txt
    printf 'sensitive\\t${target}\\t${binding}\\t1e-06\\n' >> perbscore_all_targets.${target_id}.txt
    printf 'resistant\\t${target}\\t${binding}\\t2e-06\\n' >> perbscore_all_targets.${target_id}.txt
    """
}
