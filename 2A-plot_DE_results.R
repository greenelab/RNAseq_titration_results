# J. Taroni Feb 2016
# This script compares differential expression "silver standards" which are
# differential expression analysis results from standard pipelines (i.e.,
# log transformed 100% array data and RSEM counts 100% RNA-seq [processed with
# limma::voom]) with differential expression results from RNA-seq titrated data
# (0-100%) normalized various ways.
#
# USAGE: Rscript 2A-plot_DE_results.R

library(ggplot2)
source(file.path("util", "differential_expression_functions.R"))

#### plot Basal v. Other results -----------------------------------------------

# basal v. other fit objects
basal.fit.list <-
  readRDS(file.path(
    "results", "differential_expression",
    "BRCA_titration_differential_exp_eBayes_fits_BasalvOther.RDS"))

# plot proportion of genes that are diff expressed and get topTable(s)
basal.results <- PlotProportionDE(basal.fit.list)
ggsave(file.path("plots",
                 "differential_expr_proportion_ltFDR5perc_BasalvOther.pdf"),
       plot = basal.results$plot, width = 11, height = 4.25)

# plot jaccard similarity with silver standards
basal.jacc.plot <- PlotSilverStandardJaccard(basal.results$top.table.list,
                                             title = "Basal v. Other FDR < 5%")
ggsave(file.path("plots", "basal_v_other_jaccard_lineplot.pdf"),
       plot = basal.jacc.plot, width = 8.5, height = 4)


#### plot Her2 v. LumA results -------------------------------------------------

her2.fit.list <-
  readRDS(file.path(
    "results", "differential_expression",
    "BRCA_titration_differential_exp_eBayes_fits_Her2vLumA.RDS"))

# plot proportion of genes that are diff expressed and get topTable(s)
her2.results <- PlotProportionDE(her2.fit.list)
ggsave(file.path("plots",
                 "differential_expr_proportion_ltFDR5perc_Her2vLumA.pdf"),
       plot = her2.results$plot, width = 11, height = 4.25)

# plot jaccard similarity with silver standards
her2.jacc.plot <- PlotSilverStandardJaccard(her2.results$top.table.list,
                                            title = "Her2 v. LumA FDR < 5%")
ggsave(file.path("plots", "her2_v_lumA_jaccard_lineplot.pdf"),
       plot = her2.jacc.plot, width = 8.5, height = 4)

