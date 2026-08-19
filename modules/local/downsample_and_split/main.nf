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

  output:
    path "*.RDS", emit: scrank_obj
    path "*.png", emit: umap

  when:
  task.ext.when == null || task.ext.when  

  script:
    """
    #!/bin/bash

    downsample_and_split.R ${obj} ${target} ${column} ${species} ${n_cells}
    """

  stub:
    """
    touch sensitive.RDS
    touch resistant.RDS
    touch umap.png
    """
}

