process SCTENIFOLDKNK_NULL {
  """
  The knockout track's null model for one cell type: the knockout, alignment
  and distance SCTENIFOLDKNK_KO runs for a target, run for --sctknk_null
  random genes of the wild-type network instead. Every knockout of a network
  moves its hub genes most, by an amount set by the edges it removed, so
  SCTENIFOLDKNK_KO reads this to test each gene of a target's knockout
  against what that gene does under random knockouts of comparable strength.
  """

  tag "${wt.name.replace('_sctknk_wt.rds', '')}"

  label "r_sctenifoldknk"

  // Same image as SCTENIFOLDKNK_KO -- see the note there.
  container "/home/lgdqamorim/scratch/singularity_images/sctenifoldknk-v1.0.sif"

  input:
    path wt
    path targets
    val n_null
    val n_cores
    val seed
    val ndim

  output:
    path "*_sctknk_null.rds", emit: null_model

  when:
  task.ext.when == null || task.ext.when

  script:
    """
    #!/bin/bash
    export R_MAX_VSIZE=${ task.memory ? task.memory.toGiga() + 'Gb' : '100Gb' }
    sctenifoldknk_null.R ${wt} ${n_null} ${n_cores} ${seed} ${ndim} ${targets}
    """

  stub:
    """
    touch ${wt.name.replace('_sctknk_wt.rds', '')}_sctknk_null.rds
    """
}
