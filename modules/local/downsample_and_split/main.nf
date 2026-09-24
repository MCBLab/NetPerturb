process DOWNSAMPLE {
  """
  Downsamples and splits the Seurat object into cell types
  """

  label 'process_medium'
  label "r_scrank"

  container "${ workflow.containerEngine == 'singularity' ? 'docker://juliaapolonio/scrank:latest':
            'docker.io/juliaapolonio/scrank:latest' }"

  input:
    path obj
    path target
    val column
    val species
    val n_cells
    val min_cells
    val assay
    val n_hvg
    val seed

  output:
    path "*.RDS", emit: scrank_obj
    path "*.png", emit: umap
    path "targets_qc.txt", emit: targets
    path "target_qc.tsv", emit: target_qc
    path "cell_counts.tsv", emit: cell_counts
    path "target_expression.tsv", emit: target_expression

  when:
  task.ext.when == null || task.ext.when  

  script:
    """
    #!/bin/bash

    downsample_and_split.R ${obj} ${target} ${column} ${species} ${n_cells} ${min_cells} ${assay} ${n_hvg} ${seed}
    """

  stub:
    """
    touch sensitive.RDS
    touch resistant.RDS
    touch umap.png
    grep -v '^[[:space:]]*\$' ${target} > targets_qc.txt
    printf 'target\\tgene\\treason\\taction\\n' > target_qc.tsv
    printf 'identity\\tn_input\\tn_used\\tstatus\\tgenes_total\\tgenes_gene4use\\tgenes_expressed\\tgenes_filtered\\n' > cell_counts.tsv
    printf 'sensitive\\t100\\t100\\tkept\\t500\\t200\\t180\\t20\\n' >> cell_counts.tsv
    printf 'resistant\\t100\\t100\\tkept\\t500\\t200\\t190\\t10\\n' >> cell_counts.tsv
    printf 'identity\\tgene\\tavg_expression\\tpct_expressing\\n' > target_expression.tsv
    printf 'sensitive\\tBrd4\\t1.2\\t60\\n' >> target_expression.tsv
    printf 'resistant\\tBrd4\\t0.8\\t45\\n' >> target_expression.tsv
    """
}

