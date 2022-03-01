#!/bin/sh
set -euo pipefail

# Usage: bash run_differential_expression_experiments.sh CANCER_TYPE SUBTYPE_VS_OTHERS SUBTYPE_VS_SUBTYPE SUBTYPE_VS_SUBTYPE_SMALL
# where CANCER_TYPE is one of BRCA or GBM
# SUBTYPE_VS_OTHER is the subtype you want to compare to all others (e.g. Basal)
# SUBTYPE_VS_SUBTYPE is the two subtypes you want to compare head to head (e.g. LumA,Her2) (comma-separated)
# SUBTYPE_VS_SUBTYPE_SMALL is the two subtypes you want to compare head to head when limiting the sample size (e.g. LumA,Her2) (comma-separated)

cancer_type=$1
subtype_vs_others=$2
subtype_vs_subtype=$3
subtype_vs_subtype_small=$4
ncores=$5

if [ $cancer_type != "BRCA" ] && [ $cancer_type != "GBM" ]; then
  echo Cancer type must be BRCA or GBM in run_differential_expression_experiments.sh [cancer_type]
  exit
fi

# Run differential expression scripts
Rscript 1A-detect_differentially_expressed_genes.R --cancer_type $cancer_type --subtype_vs_others $subtype_vs_others --subtype_vs_subtype $subtype_vs_subtype --ncores $ncores
Rscript 3A-small_n_differential_expression.R --cancer_type $cancer_type --subtype_vs_subtype $subtype_vs_subtype_small --ncores $ncores
