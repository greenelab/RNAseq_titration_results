#!/bin/sh
set -euo pipefail

# cancer type (either BRCA or GBM)
cancer_type=$1

if [ $cancer_type != "BRCA" ] && [ $cancer_type != "GBM" ]; then
  echo Cancer type must be BRCA or GBM in run_machine_learning_experiments.sh [cancer_type]
  exit
fi

# Run ten repeats of the supervised analysis
Rscript classifier_repeat_wrapper.R --cancer_type $cancer_type --n_repeats 10

# Run the unsupervised analyses
Rscript 4-ica_pca_feature_reconstruction.R --cancer_type $cancer_type --n_components 50
Rscript 5-predict_subtype_reconstructed_data.R --cancer_type $cancer_type
Rscript 6-plot_recon_error_kappa.R --cancer_type $cancer_type
