process EXTRA_TARGET_QC {
  """
  Checks the --extra_target targets against the networks one track already
  built, and passes on the ones those networks carry. The extras never reach
  DOWNSAMPLE or a network step, so adding them to a finished run reuses every
  network it has; a gene no network was built on is listed for the report as
  not run instead. Included twice in main.nf, once per track.
  """

  label "r_scrank"
  tag "${track}"

  // Needs only R and Matrix, which the DOWNSAMPLE image already carries.
  container "${ workflow.containerEngine == 'singularity' ? 'docker://juliaapolonio/scrank:latest':
            'docker.io/juliaapolonio/scrank:latest' }"

  input:
    path networks
    path extra_target, stageAs: "extra_target.txt"
    path target, stageAs: "target.txt"
    val track
    val method

  output:
    path "extra_targets_${track}.txt", emit: targets
    path "extra_target_qc_${track}.tsv", emit: qc

  when:
  task.ext.when == null || task.ext.when

  script:
    """
    #!/bin/bash
    extra_target_qc.R ${track} ${method} extra_target.txt target.txt ${networks}
    """

  stub:
    """
    # a target naming "Absent" stands in for one no network carries
    awk 'NF && !/Absent/' extra_target.txt > extra_targets_${track}.txt
    printf 'track\\ttarget\\tgene\\tnetworks\\treason\\taction\\n' > extra_target_qc_${track}.tsv
    awk -v OFS='\\t' -v t=${track} 'NF && /Absent/ { print t, \$0, \$0, "0/2", "not in any network this run built", "not run" }' \\
      extra_target.txt >> extra_target_qc_${track}.tsv
    """
}
