process MERGE {
  """
  Merge the per-target rank score and top connection tables into one file each
  """

  label "r_merge"

  input:
    path rank_scores
    path top_connections

  output:
    path "perbscore_all_targets.txt", emit: merged_rank_scores
    path "top_connections_all_targets.txt", emit: merged_top_connections

  when:
    task.ext.when == null || task.ext.when

  script:
    """
    #!/bin/bash

    # Both tables are one-header-per-target files from the same RANK_SCORE
    # tasks, so they concatenate the same way.
    merge_tables() {
        local out=\$1
        shift
        local first=1

        for table in "\$@"; do
            if [ "\$first" -eq 1 ]; then
                head -n 1 "\$table" > "\$out"
                first=0
            fi

            tail -n +2 "\$table" >> "\$out"
        done
    }

    merge_tables perbscore_all_targets.txt ${rank_scores}
    merge_tables top_connections_all_targets.txt ${top_connections}
    """
}
