process GENIE3 {
 """
  Generates WGN based on scRank obj
 """
tag "Cell: ${scobj.baseName} | Target: ${target_gene}"

  label "r_genie3"

  container "${ workflow.containerEngine == 'singularity' ? 'docker://juliaapolonio/scrank:latest':
            'docker.io/juliaapolonio/scrank:latest' }"

  input:
    path scobj
    val n_cores
    val target_gene

  output:
    path "*.rds", emit: rank_obj

  when:
  task.ext.when == null || task.ext.when  

  script:
    """
  Rscript ${baseDir}/bin/genie3.R ${scobj} ${n_cores} ${target_gene}
    """
}
