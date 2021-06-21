#!/bin/sh

# Run ten repeats of the supervised analysis
Rscript classifier_repeat_wrapper.R 1

# Run the unsupervised analyses 
Rscript 4-ica_pca_feature_reconstruction.R 50
Rscript 5-predict_subtype_reconstructed_data.R
Rscript 6-plot_recon_error_kappa.R
