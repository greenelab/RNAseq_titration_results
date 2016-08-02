# Summary

The objective of this project is to identify the best practices for 
normalizing of microarray and RNA-seq data for machine learning applications 
when the training data consists of a mix of the two platforms.

It builds on [Thompson et al.](https://peerj.com/preprints/1460/).

# Scripts

Running classifier_repeat_wrapper.R will regenerate the breast cancer (BRCA)
subtype classifier results. To run one iteration of the BRCA subtype pipeline,
use run_experiments.R.

# Data

The Cancer Genome Atlas BRCA data used for these analyses
is [available at zenodo](https://zenodo.org/record/58862).
```
# To download data, run in top directory:
sh brca_data_download.sh
```

# Requirements

This requires the following R & Bioconductor packages be installed (see 
check_installs.R for confirmation of installation):

* doParallel
* parallel
* ggplot2
* reshape2
* Hmisc
* data.table
* scales
* sdcMicro
* flexclust
* fpc
* corrplot
* ape
* cluster
* plyr
* dplyr
* devtools
* quantro
* preprocessCore
* gridExtra
* huge
* caret
* limma
* glmnet
* e1071
* stringr
* gdata
* binr
* cowplot
* kernlab
* ranger

One github package (TDM) is required. To install, run:

    library(devtools)
    devtools::install_github("greenelab/TDM")

We have created an R script which will call 'require' on each of these packages
to make sure they are installed: check_installs.R
