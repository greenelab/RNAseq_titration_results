#!/bin/sh

# This script runs all analysis code and plotting scripts for the publication,
# including BRCA (subtype, TP53, PIK3CA) and GBM (subtype, TP53, PIK3CA).
# The script calls:
#  1. run_machine_learning_experiments.sh
#  2. run_differential_expression_experiments.sh (subtype only)
#  3. visualize_expression.R (subtype only)
#  4. plotting scripts, as appropriate

################################################################################
# BRCA
################################################################################

# BRCA subtype -----------------------------------------------------------------

bash run_machine_learning_experiments.sh BRCA subtype 7
bash run_differential_expression_experiments.sh BRCA Basal Her2,LumA Her2,LumA 7

Rscript visualize_expression.R --cancer_type BRCA --predictor subtype

Rscript plots/scripts/0-plot_predictor_category_distributions.R \
  --cancer_type BRCA \
  --predictor subtype \
  --output_directory plots/supplementary

Rscript plots/scripts/0-plot_predictor_category_distributions.R \
  --cancer_type BRCA \
  --predictor subtype \
  --plot_all_seeds \
  --output_directory plots/supplementary

Rscript plots/scripts/3-plot_category_kappa.R \
  --cancer_type BRCA \
  --predictor subtype \
  --output_directory plots/main

Rscript plots/scripts/6-plot_recon_kappa.R \
  --cancer_type BRCA \
  --predictor subtype \
  --output_directory plots/supplementary

Rscript plots/scripts/6-plot_recon_error.R \
  --cancer_type BRCA \
  --predictor subtype \
  --output_directory plots/supplementary

Rscript plots/scripts/7-plot_plier_pathways.R \
  --cancer_type BRCA \
  --predictor subtype \
  --output_directory plots/main

Rscript plots/scripts/2A-plot_DEGs.R \
  --cancer_type BRCA \
  --subtype_vs_others Basal \
  --proportion_output_directory plots/supplementary \
  --overlap_output_directory plots/supplementary \
  --overlap_measure Jaccard
  
Rscript plots/scripts/2A-plot_DEGs.R \
  --cancer_type BRCA \
  --subtype_vs_subtype Her2,LumA \
  --proportion_output_directory plots/supplementary \
  --overlap_output_directory plots/main \
  --overlap_measure Jaccard

Rscript plots/scripts/3A-plot_small_n_differential_expression.R \
  --cancer_type BRCA \
  --subtype_vs_subtype Her2,LumA \
  --output_directory plots/main \
  --overlap_measure Jaccard,Spearman

# ------------------------------------------------------------------------------

# BRCA TP53 --------------------------------------------------------------------

bash run_machine_learning_experiments.sh BRCA TP53 7

Rscript plots/scripts/0-plot_predictor_category_distributions.R \
  --cancer_type BRCA \
  --predictor TP53 \
  --output_directory plots/supplementary

Rscript plots/scripts/0-plot_predictor_category_distributions.R \
  --cancer_type BRCA \
  --predictor TP53 \
  --plot_all_seeds \
  --output_directory plots/supplementary

Rscript plots/scripts/3-plot_category_kappa.R \
  --cancer_type BRCA \
  --predictor TP53 \
  --output_directory plots/supplementary

# ------------------------------------------------------------------------------

# BRCA PIK3CA ------------------------------------------------------------------

bash run_machine_learning_experiments.sh BRCA PIK3CA 7

Rscript plots/scripts/0-plot_predictor_category_distributions.R \
  --cancer_type BRCA \
  --predictor PIK3CA \
  --output_directory plots/supplementary

Rscript plots/scripts/0-plot_predictor_category_distributions.R \
  --cancer_type BRCA \
  --predictor PIK3CA \
  --plot_all_seeds \
  --output_directory plots/supplementary

Rscript plots/scripts/3-plot_category_kappa.R \
  --cancer_type BRCA \
  --predictor PIK3CA \
  --output_directory plots/supplementary

# ------------------------------------------------------------------------------

################################################################################
# GBM
################################################################################

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

# GBM TP53 ---------------------------------------------------------------------

bash run_machine_learning_experiments.sh GBM TP53 7

Rscript plots/scripts/0-plot_predictor_category_distributions.R \
  --cancer_type GBM \
  --predictor TP53 \
  --output_directory plots/supplementary

Rscript plots/scripts/0-plot_predictor_category_distributions.R \
  --cancer_type GBM \
  --predictor TP53 \
  --plot_all_seeds \
  --output_directory plots/supplementary

Rscript plots/scripts/3-plot_category_kappa.R \
  --cancer_type GBM \
  --predictor TP53 \
  --output_directory plots/main

# ------------------------------------------------------------------------------

# GBM PIK3CA -------------------------------------------------------------------

bash run_machine_learning_experiments.sh GBM PIK3CA 7

Rscript plots/scripts/0-plot_predictor_category_distributions.R \
  --cancer_type GBM \
  --predictor PIK3CA \
  --output_directory plots/supplementary

Rscript plots/scripts/0-plot_predictor_category_distributions.R \
  --cancer_type GBM \
  --predictor PIK3CA \
  --plot_all_seeds \
  --output_directory plots/supplementary

Rscript plots/scripts/3-plot_category_kappa.R \
  --cancer_type GBM \
  --predictor PIK3CA \
  --output_directory plots/supplementary

# ------------------------------------------------------------------------------
