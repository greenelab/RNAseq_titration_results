#!/bin/sh

# This script runs all analysis code and plotting scripts for the publication,
# including BRCA only (subtype, TP53, PIK3CA).
# The script calls:
#  1. run_machine_learning_experiments.sh
#  2. run_differential_expression_experiments.sh (subtype only)
#  3. plots/scripts/visualize_expression.R (subtype only)
#  4. plotting scripts, as appropriate

################################################################################
# BRCA
################################################################################

# BRCA subtype -----------------------------------------------------------------

# run machine learning and DEG analysis scripts
bash run_machine_learning_experiments.sh BRCA subtype 7
bash run_differential_expression_experiments.sh BRCA Basal Her2,LumA Her2,LumA 7

# plot array vs. RNA-seq expression levels after normalization
Rscript plots/scripts/visualize_expression.R --cancer_type BRCA --predictor subtype

# plot difference in subtype prediction kappa between non-reconstructed and reconstructed data
Rscript plots/scripts/recon_kappa_difference.R --cancer_type BRCA --output_directory plots/supplementary

# stacked bar plot showing distribution of subtypes in train/test sets (one representative example)
Rscript plots/scripts/0-plot_predictor_category_distributions.R \
  --cancer_type BRCA \
  --predictor subtype \
  --output_directory plots/supplementary

# stacked bar plots showing distribution of subtypes in train/test sets (all seeds)
Rscript plots/scripts/0-plot_predictor_category_distributions.R \
  --cancer_type BRCA \
  --predictor subtype \
  --plot_all_seeds \
  --output_directory plots/supplementary

# violin + line plots showing kappa values from predictions on test data
Rscript plots/scripts/3-plot_category_kappa.R \
  --cancer_type BRCA \
  --predictor subtype \
  --output_directory plots/main

# violin + line plots showing kappa values from predictions on reconstructed test data
Rscript plots/scripts/6-plot_recon_kappa.R \
  --cancer_type BRCA \
  --predictor subtype \
  --output_directory plots/supplementary

# violin + line plots showing gene-level MASE values from reconstructed test data
Rscript plots/scripts/6-plot_recon_error.R \
  --cancer_type BRCA \
  --predictor subtype \
  --output_directory plots/supplementary

# violin plots showing proportion of pathways significant in PLIER analyses
Rscript plots/scripts/7-plot_plier_pathways.R \
  --cancer_type BRCA \
  --predictor subtype \
  --output_directory plots/main

# bar plot showing proportion of genes differentially expressed (Basal vs. Others)
# line plot showing overlap with silver standard DEGs (Basal vs. Others)
Rscript plots/scripts/1A-plot_DEGs.R \
  --cancer_type BRCA \
  --subtype_vs_others Basal \
  --proportion_output_directory plots/supplementary \
  --overlap_output_directory plots/supplementary \
  --overlap_measure Jaccard,Spearman
  
# bar plot showing proportion of genes differentially expressed (Her2 vs. LumA)
# line plot showing overlap with silver standard DEGs (Her2 vs. LumA)
Rscript plots/scripts/1A-plot_DEGs.R \
  --cancer_type BRCA \
  --subtype_vs_subtype Her2,LumA \
  --proportion_output_directory plots/supplementary \
  --overlap_output_directory plots/main \
  --overlap_measure Jaccard,Spearman

# line plot showing overlap with silver standard DEGs (Her2 vs. LumA) across small n values
Rscript plots/scripts/2A-plot_small_n_differential_expression.R \
  --cancer_type BRCA \
  --subtype_vs_subtype Her2,LumA \
  --output_directory plots/main \
  --overlap_measure Jaccard,Spearman

# ------------------------------------------------------------------------------

# BRCA TP53 --------------------------------------------------------------------

# run machine learning analysis scripts
bash run_machine_learning_experiments.sh BRCA TP53 7

# stacked bar plot showing distribution of subtypes in train/test sets (one representative example)
Rscript plots/scripts/0-plot_predictor_category_distributions.R \
  --cancer_type BRCA \
  --predictor TP53 \
  --output_directory plots/supplementary

# stacked bar plots showing distribution of subtypes in train/test sets (all seeds)
Rscript plots/scripts/0-plot_predictor_category_distributions.R \
  --cancer_type BRCA \
  --predictor TP53 \
  --plot_all_seeds \
  --output_directory plots/supplementary

# violin + line plots showing kappa values from predictions on test data
Rscript plots/scripts/3-plot_category_kappa.R \
  --cancer_type BRCA \
  --predictor TP53 \
  --null_model \
  --output_directory plots/supplementary

# ------------------------------------------------------------------------------

# BRCA PIK3CA ------------------------------------------------------------------

# run machine learning analysis scripts
bash run_machine_learning_experiments.sh BRCA PIK3CA 7

# stacked bar plot showing distribution of subtypes in train/test sets (one representative example)
Rscript plots/scripts/0-plot_predictor_category_distributions.R \
  --cancer_type BRCA \
  --predictor PIK3CA \
  --output_directory plots/supplementary

# stacked bar plots showing distribution of subtypes in train/test sets (all seeds)
Rscript plots/scripts/0-plot_predictor_category_distributions.R \
  --cancer_type BRCA \
  --predictor PIK3CA \
  --plot_all_seeds \
  --output_directory plots/supplementary

# violin + line plots showing kappa values from predictions on test data
Rscript plots/scripts/3-plot_category_kappa.R \
  --cancer_type BRCA \
  --predictor PIK3CA \
  --null_model \
  --output_directory plots/supplementary

# ------------------------------------------------------------------------------
