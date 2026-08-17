process HDWGCNA {
  """
  Generates a co-expression network using hdWGCNA from a Seurat object
  """

  label "r_hdwgcna"
  tag "${scobj.baseName}"

  container "${ workflow.containerEngine == 'singularity' ? 'docker://leoshow21/hdwgcna:v2':
            'docker.io/leoshow21/hdwgcna:v2' }"

  input:
    path scobj
    val column
    val n_cores
    val cut_ratio
    val min_cells

  output:
    path "*_weight_hdWGCNA_*.rds", emit: rank_obj, optional: true

  when:
  task.ext.when == null || task.ext.when

  script:
    """
    #!/bin/bash
    hdwgcna.R ${scobj} ${column} ${n_cores} ${cut_ratio} ${min_cells}
    """

  stub:
    """
    touch ${scobj.baseName}_weight_hdWGCNA_100.rds
    """
}
