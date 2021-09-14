#!/bin/sh
set -euo pipefail

# cancer type (either BRCA or GBM)
cancer_type=$1
predictor=$2

if [ $cancer_type != "BRCA" ] && [ $cancer_type != "GBM" ]; then
  echo Cancer type must be BRCA or GBM in run_machine_learning_experiments.sh [cancer_type] [predictor]
  exit
fi

if [ $predictor != "subtype" ] && [ $predictor != "TP53" ] && [ $predictor != "PIK3CA" ]; then
  echo Predictor must be subtype, TP53, or PIK3CA in run_machine_learning_experiments.sh [cancer_type] [predictor]
  exit
fi

# Run ten repeats of the supervised analysis
# if the predictor is a gene, also generate null models
if [ $predictor == "TP53" ] || [ $predictor == "PIK3CA" ]; then
  Rscript classifier_repeat_wrapper.R --cancer_type $cancer_type --predictor $predictor --n_repeats 10
  Rscript classifier_repeat_wrapper.R --cancer_type $cancer_type --predictor $predictor --n_repeats 10 --null_model
else
  Rscript classifier_repeat_wrapper.R --cancer_type $cancer_type --predictor $predictor --n_repeats 10
fi

# Run the unsupervised analyses
Rscript 4-ica_pca_feature_reconstruction.R --cancer_type $cancer_type --predictor $predictor --n_components 50
Rscript 5-predict_category_reconstructed_data.R --cancer_type $cancer_type --predictor $predictor
Rscript 6-plot_recon_error_kappa.R --cancer_type $cancer_type --predictor $predictor
