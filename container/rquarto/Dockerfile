# 1. Start with an official R base image
FROM rocker/r-ver:4.3.2

# 2. Install system dependencies needed for Quarto and R packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    gdebi-core \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libuv1-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libjpeg-dev \
    libtiff5-dev \
    && rm -rf /var/lib/apt/lists/*

# 3. Download and install Quarto CLI
ARG QUARTO_VERSION="1.5.54"
RUN curl -fsSL -o quarto.deb "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb" \
    && gdebi --n quarto.deb \
    && rm quarto.deb

# 4. Create an application directory
WORKDIR /app

# 5. Install required R packages (e.g., rmarkdown, knitr)
# Uses the repo already configured in rocker/r-ver (Posit Package Manager),
# which serves precompiled Linux binaries -- no source compilation needed.
RUN R -e "pkgs <- c('rmarkdown', 'knitr', 'DT', 'dplyr', 'tidyr', 'ggplot2'); \
    install.packages(pkgs, Ncpus = parallel::detectCores()); \
    missing <- pkgs[!pkgs %in% rownames(installed.packages())]; \
    if (length(missing)) stop('failed to install: ', paste(missing, collapse=', '))"

# 6. Render the Quarto document upon container execution
CMD ["quarto", "--help"]
