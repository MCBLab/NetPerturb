process GSEA_SCTENIFOLDKNK {
  """
  Gene set enrichment over the knockout track's differentially-regulated gene
  table. Genes are ranked on log2FC -- the log2 of the squared manifold-
  alignment distance scTenifoldKnk scored them on over the run's mean squared
  distance, so furthest-moved first, signed around the average gene -- and that
  ranked list goes to fgsea, which is what the scTenifoldKnk paper does with
  it. The sign is what makes NES two-tailed: negative means a set sat still
  while the rest of the network moved, not that it went down. Runs once over
  the merged table rather than per pair, since fgsea is fast and the table is
  already collected by then.
  """

  label "r_gsea"

  // Built from container/gsea/Dockerfile. The tag is a param so a locally
  // built .sif can be substituted without editing this file -- SCTENIFOLDKNK_KO
  // hardcodes one and that is precisely why its track only runs on the machine
  // that path belongs to.
  // A registry tag gets the docker:// prefix under singularity; a path to a
  // locally built .sif is handed over as it is, which is what makes the param
  // a real escape hatch rather than just a renamed tag.
  container "${ workflow.containerEngine == 'singularity' && !params.gsea_container.endsWith('.sif')
            ? 'docker://' + params.gsea_container
            : params.gsea_container }"

  input:
    path dr_table
    path gmt
    val min_size
    val max_size
    val seed

  output:
    path "gsea_all_targets.txt", emit: gsea_table

  when:
    task.ext.when == null || task.ext.when

  script:
    """
    #!/bin/bash
    gsea_sctenifoldknk.R "${dr_table}" "${gmt}" ${min_size} ${max_size} ${seed}
    """

  stub:
    """
    printf 'cell_type\\ttarget\\tpathway\\tpval\\tpadj\\tlog2err\\tES\\tNES\\tsize\\tleadingEdge\\n' > gsea_all_targets.txt
    printf 'sensitive\\tBrd4\\tTEST_MOVED\\t1e-10\\t2e-09\\t0.86\\t0.77\\t3.04\\t30\\tGeneA;GeneB\\n' >> gsea_all_targets.txt
    """
}
