process SCTENIFOLDNET {
  """
  Generates a network using sctenifoldnet from a Seurat object
  """

  label "r_sctenifoldnet"

  container "${ workflow.containerEngine == 'singularity' ? 'docker://diegomscoelho/sctenifoldnet:v1.2.1':
            'docker.io/diegomscoelho/sctenifoldnet:v1.2.1' }"

  input:
    path scobj
    val n_cores

  output:
    path "*.rds", emit: rank_obj

  when:
  task.ext.when == null || task.ext.when  

  script:
    """
    #!/bin/bash
    sctenifoldnet.R ${scobj} ${n_cores}
    """

  stub:
    """
    touch ${scobj.baseName}_weight_sctenifoldnet_100.rds
    """
}
