process REPORT {
  """
  Render perbscore_all_targets.txt into a queryable Quarto HTML report
  """

  label "r_report"

  container "${ workflow.containerEngine == 'singularity' ? 'docker://rocker/verse:4.4.1':
            'docker.io/rocker/verse:4.4.1' }"

  input:
    path perbscore_table
    path report_template, stageAs: "report_template.qmd"

  output:
    path "netperturb_report.html", emit: report

  when:
    task.ext.when == null || task.ext.when

  script:
    """
    #!/bin/bash

    cp ${report_template} report.qmd
    quarto render report.qmd \
      -P perbscore_file:${perbscore_table} \
      --output netperturb_report.html
    """

  stub:
    """
    touch netperturb_report.html
    """
}
