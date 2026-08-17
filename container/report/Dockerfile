FROM rocker/geospatial:4.5.1

# Evita que o BLAS tente abrir threads demais durante a instalação de
# pacotes (estoura RLIMIT_NPROC em builds rodando em nós compartilhados).
ENV OPENBLAS_NUM_THREADS=1
ENV OMP_NUM_THREADS=1

# O site-library padrão da imagem base é dono de root e não é gravável em
# builds rootless (sem subuid/subgid configurado). Instala num diretório
# novo em vez disso, que é criado pelo próprio build e por isso gravável.
RUN mkdir -p /opt/R-extra-libs
ENV R_LIBS=/opt/R-extra-libs

# Bioconductor dependencies for the report (ComplexHeatmap for plot_net's
# heatmap mode, fgsea because scRank imports it)
RUN R -e "BiocManager::install(c( \
    'ComplexHeatmap', \
    'fgsea' \
), update=FALSE, ask=FALSE, lib='/opt/R-extra-libs')"

# CRAN dependencies for the report
RUN R -e "install.packages(c('DT', 'ggpubr'), repos='https://cloud.r-project.org', lib='/opt/R-extra-libs')"

# scRank from GitHub, pinned to the commit validated against this pipeline
RUN R -e "remotes::install_github('ZJUFanLab/scRank', ref='23ede783f6a4271453f7066e00cc362bd7250dd7', lib='/opt/R-extra-libs')"

WORKDIR /workspace

CMD ["R"]
