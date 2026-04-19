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
    path "seurat_downsample.rds", emit: seurat_downsample
    path "unique_cell_types.txt", emit: cell_types_list

  when:
  task.ext.when == null || task.ext.when  

  script:
    """
    #!/bin/bash
    downsample_and_split.R "${obj}" "${target}" "${column}" "${species}" "${n_cells}"
    """
}

