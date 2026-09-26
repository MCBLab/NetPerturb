process NETWORK_QC {
  """
  Summarises the structure of every network the run built -- size, density,
  how concentrated its strength is -- and where each target gene sits in it,
  for the report's Data quality section.
  """

  label "r_scrank"

  // Needs only R and Matrix, which the DOWNSAMPLE image already carries.
  container "${ workflow.containerEngine == 'singularity' ? 'docker://juliaapolonio/scrank:latest':
            'docker.io/juliaapolonio/scrank:latest' }"

  input:
    path networks
    path targets
    val network

  output:
    path "network_qc.tsv", emit: network_qc
    path "target_network_qc.tsv", emit: target_network_qc

  when:
  task.ext.when == null || task.ext.when

  script:
    """
    #!/bin/bash
    network_qc.R ${network} ${targets} ${networks}
    """

  stub:
    """
    printf 'track\\tmethod\\tidentity\\tn_genes\\tn_edges\\tdensity\\tmean_abs_weight\\tisolated_frac\\thub_share_top1pct\\tstrength_gini\\tn_metacells\\tsoft_power\\tpower_estimate\\tsft_r2\\tmetacells_per_gene\\n' > network_qc.tsv
    printf 'track\\tmethod\\tidentity\\tgene\\tin_network\\tdegree\\tstrength\\tstrength_percentile\\n' > target_network_qc.tsv
    """
}
