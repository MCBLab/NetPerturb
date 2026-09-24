process SCTENIFOLDKNK_KO {
  """
  Knocks one target out of the wild-type network SCTENIFOLDKNK_BUILD saved
  for a cell type and reports scTenifoldKnk's genome-wide differentially-
  regulated (DR) gene table. Unlike GENIE3/SCRANK/HDWGCNA, the KO is
  target-specific by construction, so this runs once per (cell type, target)
  pair -- see the .combine() call in main.nf. A ';'-joined target is one joint
  knockout of all its genes.
  """

  // One task per (cell type, target), so the pair goes in the task name the
  // way RANK_SCORE puts its target there -- otherwise a failure in the log is
  // just "SCTENIFOLDKNK_KO (5)" with no way to tell which pair it was.
  tag "${wt.name.replace('_sctknk_wt.rds', '')}:${target}"

  label "r_sctenifoldknk"

  // No published image yet -- built locally per container/sctenifoldknk/Dockerfile
  // (satijalab/seurat:5.4.0 + scTenifoldNet from GitHub + scTenifoldKnk from
  // CRAN) and cached at this path. Referenced directly rather than through
  // Nextflow's singularity pull cache, since there is no real docker:// tag
  // to key that cache off of.
  container "/home/lgdqamorim/scratch/singularity_images/sctenifoldknk-v1.0.sif"

  input:
    tuple path(wt), val(target)
    val n_cores
    val plot
    val seed

  output:
    path "*.txt", emit: dr_table
    // only written with --sctknk_plot
    path "*_sctknk_plot_*.pdf", optional: true, emit: plot

  when:
  task.ext.when == null || task.ext.when

  script:
    """
    #!/bin/bash
    sctenifoldknk_ko.R ${wt} "${target}" ${n_cores} ${plot} ${seed}
    """

  stub:
    def target_id = target.replaceAll(/[^A-Za-z0-9_.-]+/, "_")
    def cell_type = wt.name.replace("_sctknk_wt.rds", "")
    def plot_cmd = plot ? "touch ${cell_type}_sctknk_plot_${target_id}.pdf" : ""
    """
    touch ${cell_type}_sctenifoldknk_${target_id}.txt
    ${plot_cmd}
    """
}
