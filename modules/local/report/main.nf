process REPORT {
  """
  Render the scRank Quarto dashboard for the current pipeline run
  """

  label "quarto_report"

  container params.report_container

  input:
    path merged_obj
    path perbscore
    path qmd_template, stageAs: 'template.qmd'
    path logo, stageAs: 'lab-logo.png'
    val target
    val n_cores

  output:
    path "nf_scrank_report.html", emit: html
    path "nf_scrank_report_files", optional: true, emit: files

  when:
  task.ext.when == null || task.ext.when

  script:
    """
    #!/bin/bash
    export XDG_CACHE_HOME="\$PWD/.cache"
    mkdir -p "\$XDG_CACHE_HOME"
    cp ${qmd_template} nf_scrank_report.qmd
    quarto render nf_scrank_report.qmd \\
      -P merged_obj:${merged_obj} \\
      -P perbscore:${perbscore} \\
      -P target:"${target}" \\
      -P n_cores:${n_cores}
    """
}
