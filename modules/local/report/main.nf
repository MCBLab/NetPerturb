process REPORT {
  """
  Render the scRank Quarto dashboard for the current pipeline run
  """

  label "quarto_report"

  container params.report_container

  containerOptions "--bind ${params.report_r_libs_user}:${params.report_r_libs_user}" +
    (params.report_extra_lib_path ? " --bind ${params.report_extra_lib_path}:${params.report_extra_lib_path}" : "")

  input:
    path merged_obj
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
    export R_LIBS_USER="${params.report_r_libs_user}"
    # Renviron.site do container só usa o default site-library/library pra
    # R_LIBS se essa variável ainda não estiver setada (\${R_LIBS-default});
    # setando aqui, nosso R_LIBS_USER passa a vir ANTES do site-library no
    # .libPaths(), evitando conflito de versão com pacotes já embutidos na
    # imagem (ex.: dplyr/rlang mais antigos que os do R_LIBS_USER).
    export R_LIBS="${params.report_r_libs_user}"
    # Pacotes compilados fora do container (ex.: numa env conda) podem linkar
    # contra libs de sistema que o container não tem na mesma versão (ex.:
    # libxml2.so.16, exigido pelo igraph). report_extra_lib_path permite
    # apontar pro diretório `lib/` da env de origem pra suprir isso.
    export LD_LIBRARY_PATH="${params.report_extra_lib_path}:\$LD_LIBRARY_PATH"
    export XDG_CACHE_HOME="\$PWD/.cache"
    mkdir -p "\$XDG_CACHE_HOME"
    cp ${qmd_template} nf_scrank_report.qmd
    quarto render nf_scrank_report.qmd \\
      -P merged_obj:${merged_obj} \\
      -P target:"${target}" \\
      -P n_cores:${n_cores}
    """
}
