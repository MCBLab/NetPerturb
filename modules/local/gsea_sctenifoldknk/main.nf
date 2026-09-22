process GSEA_SCTENIFOLDKNK {
  """
  Gene set enrichment over the knockout track's differentially-regulated gene
  table. Genes are ranked by the manifold-alignment distance scTenifoldKnk
  scored them on -- furthest-moved first -- and that ranked list goes to fgsea,
  which is what the scTenifoldKnk paper does with it. Runs once over the merged
  table rather than per pair, since fgsea is fast and the table is already
  collected by then.
  """

  label "r_gsea"

  // Built from container/gsea/Dockerfile. The tag is a param so a locally
  // built .sif can be substituted without editing this file -- SCTENIFOLDKNK
  // hardcodes one and that is precisely why its track only runs on the machine
  // that path belongs to.
  container "${ workflow.containerEngine == 'singularity' ? "docker://${params.gsea_container}" :
            params.gsea_container }"

  input:
    path dr_table
    path gmt
    val min_size
    val max_size

  output:
    path "gsea_all_targets.txt", emit: gsea_table

  when:
    task.ext.when == null || task.ext.when

  script:
    """
    #!/bin/bash
    gsea_sctenifoldknk.R "${dr_table}" "${gmt}" ${min_size} ${max_size}
    """

  stub:
    """
    printf 'cell_type\\ttarget\\tpathway\\tpval\\tpadj\\tlog2err\\tES\\tNES\\tsize\\tleadingEdge\\n' > gsea_all_targets.txt
    printf 'sensitive\\tBrd4\\tTEST_MOVED\\t1e-10\\t2e-09\\t0.86\\t0.77\\t3.04\\t30\\tGeneA;GeneB\\n' >> gsea_all_targets.txt
    """
}
