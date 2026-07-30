process MERGE {
  """
  Merge per-target rank score tables into one file
  """

  label "r_merge"

  input:
    path rank_scores

  output:
    path "perbscore_all_targets.txt", emit: merged_rank_scores

  when:
    task.ext.when == null || task.ext.when

  script:
    """
    #!/bin/bash
    first=1

    for score_file in ${rank_scores}; do
        if [ "\$first" -eq 1 ]; then
            head -n 1 "\$score_file" > perbscore_all_targets.txt
            first=0
        fi

        tail -n +2 "\$score_file" >> perbscore_all_targets.txt
    done
    """
}
