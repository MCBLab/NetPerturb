process MERGE_SCTENIFOLDKNK {
  """
  Concatenate the per-(cell type, target) differentially-regulated gene
  tables from SCTENIFOLDKNK into one file, the same way MERGE does for
  RANK_SCORE's tables -- kept as its own process since scTenifoldKnk has
  only the one table type, not RANK_SCORE's pair.
  """

  label "r_merge"

  input:
    path dr_tables

  output:
    path "sctenifoldknk_all_targets.txt", emit: merged_dr_table

  when:
    task.ext.when == null || task.ext.when

  script:
    """
    #!/bin/bash

    first=1
    for table in ${dr_tables}; do
      if [ "\$first" -eq 1 ]; then
        head -n 1 "\$table" > sctenifoldknk_all_targets.txt
        first=0
      fi
      tail -n +2 "\$table" >> sctenifoldknk_all_targets.txt
    done
    """

  stub:
    """
    touch sctenifoldknk_all_targets.txt
    """
}
