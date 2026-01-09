# syntax=docker/dockerfile:1
## Build with: docker build -f container/Dockerfile .

FROM rocker/r-ver:4.5.2

LABEL org.opencontainers.image.title="Multi-Gene Correlations to Signature" \
    org.opencontainers.image.description="Containerized R workflow for multi-gene to signature correlation analysis with publication-ready plots." \
    org.opencontainers.image.url="https://github.com/NIDAP-Community/multi-gene-correlations" \
    org.opencontainers.image.source="https://github.com/NIDAP-Community/multi-gene-correlations" \
    org.opencontainers.image.documentation="https://github.com/NIDAP-Community/multi-gene-correlations#readme"

# System dependencies for building common R packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install key R packages via apt to ensure availability at runtime
RUN apt-get update && apt-get install -y --no-install-recommends \
    r-cran-jsonlite \
    r-cran-ggplot2 \
    r-cran-dplyr \
    r-cran-tidyr \
    r-cran-stringr \
    r-cran-rcolorbrewer \
    r-cran-circlize \
    r-bioc-complexheatmap \
    && rm -rf /var/lib/apt/lists/*

# Install CRAN packages using littler's install2.r (more robust)
# Install CRAN packages from a vanilla mirror (complements apt-installed packages)
RUN R -q -e "options(repos=c(CRAN='https://cran.r-project.org')); install.packages(c('GetoptLong','GlobalOptions'))"

# Install Bioconductor ComplexHeatmap via BiocManager, with R-universe fallback (already apt-installed)
RUN R -q -e "options(repos=c(CRAN='https://cran.r-project.org')); install.packages('BiocManager'); tryCatch(BiocManager::install('ComplexHeatmap', update=FALSE, ask=FALSE), error=function(e) install.packages('ComplexHeatmap', repos=c('https://jokergoo.r-universe.dev','https://cran.r-project.org')))"

# Install ComplexHeatmap via R-universe (CRAN fallback provided)
# Optional fallback via R-universe if BiocManager fails
# RUN R -q -e "install.packages('ComplexHeatmap', repos=c('https://jokergoo.r-universe.dev','https://cloud.r-project.org'))"

# Install L2P and L2Psupp from provided tarballs
RUN R -q -e "install.packages('https://github.com/CCBR/l2p/raw/master/l2p_0.0-14.tar.gz', repos=NULL)" \
    && R -q -e "install.packages('https://github.com/CCBR/l2p/raw/master/l2psupp_0.0-14.tar.gz', repos=NULL)"

WORKDIR /app

# Copy source and CLI entrypoint
COPY Multi_Gene_Correlations_to_Signature.R /app/
COPY entrypoint.R /app/
COPY container/Dockerfile /app/Dockerfile
COPY README.md /app/

# Ensure R can see apt-installed libraries at runtime
ENV R_LIBS_SITE="/usr/lib/R/site-library:/usr/lib/R/library"
RUN echo "R_LIBS_SITE=/usr/lib/R/site-library:/usr/lib/R/library" >> /usr/local/lib/R/etc/Renviron

# Default entrypoint runs the R CLI. For debugging, override with --entrypoint bash
ENTRYPOINT ["Rscript", "/app/entrypoint.R"]

# Show CLI help when no arguments are provided
CMD ["--help"]
