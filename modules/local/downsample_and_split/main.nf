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
    val target
    val column
    val species
    val n_cells

  output:
    path "*.RDS", emit: scrank_obj

  when:
  task.ext.when == null || task.ext.when  

  script:
    """
    Rscript ${baseDir}/bin/downsample_and_split.R ${obj} ${target} ${column} ${species} ${n_cells}
    """
}

