process SCTENIFOLDKNK_BUILD {
  """
  Builds scTenifoldKnk's wild-type gene regulatory network for one cell type,
  once, for every SCTENIFOLDKNK_KO knockout on that cell type to reuse. The
  network does not depend on the target, so building it per knockout only
  rebuilt the same network again.
  """

  tag "${scobj.baseName}"

  label "r_sctenifoldknk"

  // Same image as SCTENIFOLDKNK_KO -- see the note there.
  container "/home/lgdqamorim/scratch/singularity_images/sctenifoldknk-v1.0.sif"

  input:
    path scobj
    val n_cores
    val seed

  output:
    path "*_sctknk_wt.rds", emit: wt

  when:
  task.ext.when == null || task.ext.when

  script:
    """
    #!/bin/bash
    sctenifoldknk_build.R ${scobj} ${n_cores} ${seed}
    """

  stub:
    """
    touch ${scobj.baseName}_sctknk_wt.rds
    """
}
