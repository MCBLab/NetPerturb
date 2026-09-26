process MERGE_SCTENIFOLDKNK {
  """
  Concatenate the per-(cell type, target) differentially-regulated gene
  tables from SCTENIFOLDKNK_KO into one file, the same way MERGE does for
  RANK_SCORE's tables -- kept as its own process since scTenifoldKnk has
  only the one table type, not RANK_SCORE's pair. Also writes a summary of
  it: how many genes each knockout affected, per cell type and target.
  """

  label "r_merge"

  input:
    path dr_tables

  output:
    path "sctenifoldknk_all_targets.txt", emit: merged_dr_table
    path "sctenifoldknk_summary.txt", emit: summary

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

    # One row per knockout (cell type x target): how many genes were tested,
    # how many it affected (differentially regulated, p.adj < 0.05), their
    # share, and the gene it moved furthest (largest FC). A knockout that
    # affected nothing still has its row, since all its tested genes are in
    # the table; one that never produced a table -- it failed, or timed out on
    # every attempt and was ignored -- has none. Columns are found by name, so
    # this does not depend on the table's column order.
    awk -F '\t' -v OFS='\t' '
      NR == 1 { for (i = 1; i <= NF; i++) col[\$i] = i; next }
      {
        key = \$col["cell_type"] OFS \$col["target"]
        if (!(key in tested)) order[++n] = key
        tested[key]++
        if (\$col["p.adj"] != "NA" && \$col["p.adj"] + 0 < 0.05) {
          affected[key]++
          if (!(key in best) || \$col["FC"] + 0 > best[key]) { best[key] = \$col["FC"] + 0; top[key] = \$col["gene"] }
        }
      }
      END {
        print "cell_type", "target", "genes_tested", "genes_affected", "share_affected", "top_gene"
        for (k = 1; k <= n; k++) {
          key = order[k]; a = affected[key] + 0
          print key, tested[key], a, sprintf("%.4f", a / tested[key]), (a > 0 ? top[key] : "NA")
        }
      }' sctenifoldknk_all_targets.txt > sctenifoldknk_summary.txt
    """

  stub:
    """
    touch sctenifoldknk_all_targets.txt
    printf 'cell_type\\ttarget\\tgenes_tested\\tgenes_affected\\tshare_affected\\ttop_gene\\n' > sctenifoldknk_summary.txt
    """
}
