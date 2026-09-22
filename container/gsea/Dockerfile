# Same base as the report image: rocker/r-ver comes with Posit Package Manager
# configured, which serves precompiled Linux binaries for CRAN and Bioconductor
# alike, so nothing here compiles from source.
FROM rocker/r-ver:4.3.2

# fgsea needs a compiler toolchain only if PPM has no binary for this platform;
# these are the headers its dependencies link against either way.
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# fgsea brings BiocParallel and data.table along with it. The check mirrors the
# one in container/rquarto/Dockerfile: fail the build loudly rather than ship an
# image whose first task dies on a missing package.
RUN R -e "install.packages('BiocManager'); \
    BiocManager::install('fgsea', ask = FALSE, update = FALSE); \
    pkgs <- c('fgsea', 'BiocParallel'); \
    missing <- pkgs[!pkgs %in% rownames(installed.packages())]; \
    if (length(missing)) stop('failed to install: ', paste(missing, collapse=', '))"

WORKDIR /workspace

CMD ["R"]
