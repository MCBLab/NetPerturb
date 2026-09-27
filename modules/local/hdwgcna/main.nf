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
    val seed

  output:
    path "*_weight_hdWGCNA_*.rds", emit: rank_obj, optional: true
    path "*_hdwgcna_qc.tsv", emit: qc, optional: true

  when:
  task.ext.when == null || task.ext.when

  script:
    """
    #!/bin/bash
    hdwgcna.R ${scobj} ${column} ${n_cores} ${cut_ratio} ${min_cells} ${seed}
    """

  stub:
    """
    touch ${scobj.baseName}_weight_hdWGCNA_100.rds
    printf 'identity\\tn_metacells\\tsoft_power\\tpower_estimate\\tsft_r2\\n${scobj.baseName}\\t100\\t6\\t6\\t0.85\\n' > ${scobj.baseName}_hdwgcna_qc.tsv
    """
}
