FROM rocker/tidyverse:4.1.2

# Update apt-get and install other libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    libfftw3-dev \
    libbz2-dev \
    libgdal-dev \
    libgeos-dev \
    libglpk40 \
    liblzma-dev \
    libmagick++-dev \
    libproj-dev \
    libudunits2-dev \
    libxt-dev \
    python3-pip \
    python3-dev

# Install pyrefinebio v0.3.4
RUN pip3 install pyrefinebio==0.3.4

# Install some Bioconductor packages dependencies
RUN Rscript -e "options(warn = 2); BiocManager::install(c( \
    'DESeq2', \
    'EnsDb.Hsapiens.v86' , \
    'ensembldb', \
    'GenomicRanges', \
    'GenomeInfoDb', \
    'IRanges', \
    'limma', \
    'MAST', \
    'monocle', \
    'multtest', \
    'quantro', \
    'qvalue', \
    'rtracklayer', \
    'S4Vectors', \
    'SingleCellExperiment', \
    'SummarizedExperiment'), \
    update = FALSE, \
    version = 3.14)"

# R packages
RUN install2.r --error --deps TRUE \
    ape \
    binr \
    caret \
    cluster \
    corrplot \
    cowplot \
    data.table \
    devtools \
    doParallel \
    dplyr \
    e1071 \
    fastICA \
    flexclust \
    fpc \
    gdata \
    ggplot2 \
    ggupset \
    glmnet \
    gridExtra \
    here \
    Hmisc \
    huge \
    jsonlite \
    kernlab \
    locfit \
    optparse \
    plyr \
    ranger \
    RColorBrewer \
    reshape2 \
    scales \
    sdcMicro \
    Seurat \
    stringr \
    styler \
    viridis

# Threading issue with preprocessCore::normalize.quantiles
# https://support.bioconductor.org/p/122925/#124701
# https://github.com/bmbolstad/preprocessCore/issues/1#issuecomment-326756305
RUN Rscript -e "options(warn = 2); BiocManager::install( \
    'preprocessCore', \
    configure.args = '--disable-threading', \
    force = TRUE, \
    update = FALSE, \
    version = 3.14)"

# ref = 341eb77105e7efd2654b4f112578648584936e06 is latest greenelab/TDM commit (retrieved 2021-05-28)
RUN Rscript -e "options(warn = 2); remotes::install_github( \
    'greenelab/TDM', ref = 'b041807835d4076c5549356c86c44f087d713b1a')"

# ref = 08ed6b54e4efe5249107cb335cd8e169657cbc44 is wgmao/PLIER commit corresponding to v0.1.6 (retrieved 2021-11-09)
RUN Rscript -e "options(warn = 2); remotes::install_github( \
    'wgmao/PLIER', ref = '08ed6b54e4efe5249107cb335cd8e169657cbc44')"
