# 1. Start with an official R base image
FROM rocker/r-ver:4.3.2

# 2. Install system dependencies needed for Quarto and R packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    gdebi-core \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# 3. Download and install Quarto CLI
ARG QUARTO_VERSION="1.5.54"
RUN curl -o quarto.deb -L "https://github.com{QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb" \
    && gdebi --n quarto.deb \
    && rm quarto.deb

# 4. Create an application directory
WORKDIR /app

# 5. Install required R packages (e.g., rmarkdown, knitr)
RUN R -e "install.packages(c('rmarkdown', 'knitr', 'DT'), repos='https://r-project.org')"

# 7. Render the Quarto document upon container execution
CMD ["quarto", "render", "document.qmd", "--to", "html"]
