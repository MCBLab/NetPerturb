process REPORT {
  """
  Render perbscore_all_targets.txt into a queryable Quarto HTML report
  """

  label "r_report"

  container "${ workflow.containerEngine == 'singularity' ? 'docker://diegomscoelho/rquarto:1.5.54':
            'docker.io/diegomscoelho/rquarto:1.5.54' }"

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

    # On HPC the container's \$HOME is often read-only, so quarto (and the Deno
    # runtime it embeds) cannot create \$HOME/.cache/quarto. Point every XDG dir
    # at the task work dir, which is always writable and cleaned up with it.
    export HOME="\$PWD"
    export XDG_CACHE_HOME="\$PWD/.cache"
    export XDG_DATA_HOME="\$PWD/.local/share"
    export XDG_CONFIG_HOME="\$PWD/.config"
    export XDG_RUNTIME_DIR="\$PWD/.run"
    mkdir -p "\$XDG_CACHE_HOME" "\$XDG_DATA_HOME" "\$XDG_CONFIG_HOME" "\$XDG_RUNTIME_DIR"

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
