FROM rocker/tidyverse:4.0.3

# Update apt-get and install other libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    apt-utils \
    awscli \
    bzip2 \
    curl \
    dialog \
    default-jdk \
    libbz2-dev \
    libglpk40 \
    liblzma-dev \
    libopenmpi-dev \
    libreadline-dev \
    libxt-dev \
    openmpi-bin \
    python3-pip \
    python3-dev \
    zlib1g

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
    kernlab \
    optparse \
    parallel \
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

RUN installGithub.r greenelab/TDM
