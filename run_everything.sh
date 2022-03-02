#!/bin/sh

# GBM subtype ------------------------------------------------------------------

bash run_machine_learning_experiments.sh GBM subtype 7
bash run_differential_expression_experiments.sh GBM Proneural Classical,Mesenchymal Classical,Mesenchymal 7

Rscript visualize_expression.R --cancer_type GBM --predictor subtype

Rscript plots/scripts/0-plot_predictor_category_distributions.R \
  --cancer_type GBM \
  --predictor subtype \
  --output_directory plots/supplementary

Rscript plots/scripts/0-plot_predictor_category_distributions.R \
  --cancer_type GBM \
  --predictor subtype \
  --plot_all_seeds \
  --output_directory plots/supplementary

Rscript plots/scripts/3-plot_category_kappa.R \
  --cancer_type GBM \
  --predictor subtype \
  --output_directory plots/supplementary

Rscript plots/scripts/6-plot_recon_kappa.R \
  --cancer_type GBM \
  --predictor subtype \
  --output_directory plots/supplementary

Rscript plots/scripts/6-plot_recon_error.R \
  --cancer_type GBM \
  --predictor subtype \
  --output_directory plots/supplementary

Rscript plots/scripts/7-plot_plier_pathways.R \
  --cancer_type GBM \
  --predictor subtype \
  --output_directory plots/main

Rscript plots/scripts/2A-plot_DEGs.R \
  --cancer_type GBM \
  --subtype_vs_others Proneural \
  --subtype_vs_subtype Classical,Mesenchymal \
  --proportion_output_directory plots/supplementary \
  --overlap_output_directory plots/supplementary \
  --overlap_measure Jaccard

Rscript plots/scripts/3A-plot_small_n_differential_expression.R \
  --cancer_type GBM \
  --subtype_vs_subtype Classical,Mesenchymal \
  --output_directory plots/supplementary \
  --overlap_measure Jaccard,Spearman

# ------------------------------------------------------------------------------