FROM satijalab/seurat:5.4.0

# Install Bioconductor dependencies
RUN R -e "BiocManager::install(c( \
    'BiocGenerics', \
    'SummarizedExperiment' \
), ask=FALSE)"

# scTenifoldKnk (Imports: scTenifoldNet >= 1.4) wraps scTenifoldNet's own
# network-construction functions (makeNetworks, tensorDecomposition,
# manifoldAlignment, cpmNormalization) internally, so scTenifoldNet is
# installed first, the same way the retired sctenifoldnet container did.
RUN R -e "devtools::install_github('cailab-tamu/scTenifoldNet')"

# Remaining CRAN dependencies (Matrix, MASS, cli, enrichR, igraph, reshape2
# come along automatically); locfdr is optional (dr_empiricalNull = TRUE
# only) and left out.
RUN R -e "install.packages('scTenifoldKnk', repos='https://cloud.r-project.org')"

# Optional: set parallelization defaults
ENV R_MAX_VSIZE=8Gb

WORKDIR /workspace

CMD ["R"]
