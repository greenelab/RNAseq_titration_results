#!/bin/sh
set -euo pipefail

# cancer type (either BRCA or GBM)
cancer_type=$1

if [ $cancer_type = "BRCA" ] || [ $cancer_type = "GBM" ]; then
  continue
else
  echo Cancer type must be BRCA or GBM in run_machine_learning_experiments.sh [cancer_type]
  exit
fi

# Run differential expression scripts
Rscript 1A-detect_differentially_expressed_genes.R --cancer_type $cancer_type --subtype_vs_others $2 --subtype_vs_subtype $3
Rscript 2A-plot_DE_results.R --cancer_type $cancer_type --subtype_vs_others $2 --subtype_vs_subtype $3
Rscript 3A-small_n_differential_expression.R --cancer_type $cancer_type --subtype_vs_subtype $3 --max_n $4
