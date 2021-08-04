FROM rocker/tidyverse:3.6.3

# Update apt-get and install other libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
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
    optparse \
    plyr \
    ranger \
    RColorBrewer \
    reshape2 \
    scales \
    sdcMicro \
    stringr \
    styler \
    viridis

# R Bioconductor packages
RUN Rscript -e "options(warn = 2); BiocManager::install(c( \
    'EnsDb.Hsapiens.v86' , \
    'ensembldb', \
    'limma', \
    'quantro'), \
    update = FALSE)"

# Threading issue with preprocessCore::normalize.quantiles
# https://support.bioconductor.org/p/122925/#124701
# https://github.com/bmbolstad/preprocessCore/issues/1#issuecomment-326756305
RUN Rscript -e "options(warn = 2); BiocManager::install( \
    'preprocessCore', \
    configure.args = '--disable-threading', \
    update = FALSE)"

# ref = 341eb77105e7efd2654b4f112578648584936e06 is latest greenelab/TDM commit 2021-05-28 
RUN Rscript -e "options(warn = 2); remotes::install_github( \
    'greenelab/TDM', ref = '341eb77105e7efd2654b4f112578648584936e06')"
