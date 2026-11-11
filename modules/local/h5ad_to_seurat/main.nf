process H5AD_TO_SEURAT {
  """
  Converts an AnnData .h5ad (scanpy, CELLxGENE, ...) given as --obj with
  --h5ad into a Seurat object, which then takes the .h5ad's place everywhere --obj is read:
  DOWNSAMPLE and RANK_SCORE. See bin/h5ad_to_seurat.R for what is carried over.
  """

  label 'process_medium'
  label "r_h5ad"
  tag "${h5ad.name}"

  // Built from container/h5ad/Dockerfile, on top of the scrank image so the
  // object is written by the Seurat that reads it.
  container "${ workflow.containerEngine == 'singularity' ? 'docker://diegomscoelho/netperturb-h5ad:1.0':
            'docker.io/diegomscoelho/netperturb-h5ad:1.0' }"

  input:
    path h5ad
    val counts
    val gene_column

  output:
    path "${h5ad.baseName}.rds", emit: obj

  when:
  task.ext.when == null || task.ext.when

  script:
    """
    #!/bin/bash

    h5ad_to_seurat.R ${h5ad} ${h5ad.baseName}.rds "${counts ?: 'auto'}" "${gene_column ?: 'null'}"
    """

  stub:
    """
    touch ${h5ad.baseName}.rds
    """
}
