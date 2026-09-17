process SCTENIFOLDKNK {
  """
  Runs a virtual knockout of the target gene(s) with scTenifoldKnk and
  reports the resulting genome-wide differentially-regulated (DR) gene
  table. Unlike GENIE3/SCRANK/HDWGCNA, the KO is target-specific by
  construction, so this runs once per (cell type, target) pair rather than
  once per cell type -- see the .combine() call in main.nf.
  """

  label "r_sctenifoldknk"

  // No published image yet -- built locally per container/sctenifoldknk/Dockerfile
  // (satijalab/seurat:5.4.0 + scTenifoldNet from GitHub + scTenifoldKnk from
  // CRAN) and cached at this path. Referenced directly rather than through
  // Nextflow's singularity pull cache, since there is no real docker:// tag
  // to key that cache off of.
  container "/home/lgdqamorim/scratch/singularity_images/sctenifoldknk-v1.0.sif"

  input:
    tuple path(scobj), val(target)
    val n_cores

  output:
    path "*.txt", emit: dr_table

  when:
  task.ext.when == null || task.ext.when

  script:
    """
    #!/bin/bash
    sctenifoldknk.R ${scobj} "${target}" ${n_cores}
    """

  stub:
    def target_id = target.replaceAll(/[^A-Za-z0-9_.-]+/, "_")
    """
    touch ${scobj.baseName}_sctenifoldknk_${target_id}.txt
    """
}
