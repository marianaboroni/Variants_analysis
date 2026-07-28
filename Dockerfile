FROM rocker/r-ver:4.3.3

ENV R_REPOS=https://cloud.r-project.org
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    pandoc \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libcairo2-dev \
    libxt-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages(c("data.table","yaml","jsonlite","digest","ggplot2","rmarkdown","knitr","remotes","testthat","httr2","BiocManager"), repos=Sys.getenv("R_REPOS"))' \
    && Rscript -e 'BiocManager::install(c("maftools","GenomicRanges","IRanges","Biostrings","rtracklayer"), ask=FALSE, update=FALSE)'

WORKDIR /opt/tumoronly
COPY . /opt/tumoronly
RUN R CMD INSTALL .

WORKDIR /work
CMD ["Rscript", "-e", "library(tumoronly); check_installation()"]
