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
    printf 'identity\\tn_input\\tn_used\\tstatus\\n' > cell_counts.tsv
    printf 'sensitive\\t100\\t100\\tkept\\n' >> cell_counts.tsv
    printf 'resistant\\t100\\t100\\tkept\\n' >> cell_counts.tsv
    """
}

