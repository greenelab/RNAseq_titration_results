# J. Taroni Feb 2016
# This script compares differential expression "silver standards" which are
# differential expression analysis results from standard pipelines (i.e.,
# log transformed 100% array data and RSEM counts 100% RNA-seq [processed with
# limma::voom]) with differential expression results from RNA-seq titrated data
# (0-100%) normalized various ways.
#
# USAGE: Rscript 2A-plot_DE_results.R --cancer_type --subtype_vs_others --subtype_vs_subtype

option_list <- list(
  optparse::make_option("--cancer_type",
                        default = NULL,
                        help = "Cancer type"),
  optparse::make_option("--subtype_vs_others",
                        default = NULL,
                        help = "Subtype used for comparison against all others"),
  optparse::make_option("--subtype_vs_subtype",
                        default = NULL,
                        help = "Subtypes used in head-to-head comparison (comma-separated without space e.g. Type1,Type2)")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
source("util/option_functions.R")
check_options(opt)

# load libraries
suppressMessages(source("ggplot2"))
source(file.path("util", "differential_expression_functions.R"))

# set options
cancer_type <- opt$cancer_type
subtype_vs_others <- opt$subtype_vs_others
subtype_vs_subtype <- opt$subtype_vs_subtype
two_subtypes <- as.vector(stringr::str_split(subtype_vs_subtype, pattern = ",", simplify = TRUE))

#### plot Subtype v. Other results ---------------------------------------------

# subtype v. other fit objects
subtype_vs_others.fit.list <- readRDS(file.path("results", "differential_expression",
  paste0(cancer_type, "_titration_differential_exp_eBayes_fits_",
         subtype_vs_others, "vOther.RDS")))

# plot proportion of genes that are diff expressed and get topTable(s)
subtype_vs_others.results <- PlotProportionDE(subtype_vs_others.fit.list)
ggsave(file.path("plots",
                 paste0(cancer_type, "_differential_expr_proportion_ltFDR5perc_",
                        subtype_vs_others, "vOther.pdf")),
       plot = subtype_vs_others.results$plot, width = 11, height = 4.25)

# plot jaccard similarity with silver standards
subtype_vs_others.jacc.plot <- PlotSilverStandardJaccard(subtype_vs_others.results$top.table.list,
                                                         title = paste(subtype_vs_others, "v. Other FDR < 5%"))
ggsave(file.path("plots", paste0(subtype_vs_others, "_v_other_jaccard_lineplot.pdf")),
       plot = subtype_vs_others.jacc.plot, width = 8.5, height = 4)

#### plot Subtype v. Subtype results --------------------------------------------
subtypes_combination <- stringr::str_c(two_subtypes, collapse = "v")
subtypes_combination_nice <- stringr::str_c(two_subtypes, collapse = " vs. ")

last_subtype.fit.list <- readRDS(file.path(
  "results", "differential_expression",
  paste0(cancer_type, "_titration_differential_exp_eBayes_fits_",
         subtypes_combination, ".RDS")))

# plot proportion of genes that are diff expressed and get topTable(s)
last_subtype.results <- PlotProportionDE(last_subtype.fit.list)
ggsave(file.path("plots",
                 paste0(cancer_type, "_differential_expr_proportion_ltFDR5perc_",
                        subtypes_combination, ".pdf")),
       plot = last_subtype.results$plot, width = 11, height = 4.25)

# plot jaccard similarity with silver standards
last_subtype.jacc.plot <- PlotSilverStandardJaccard(last_subtype.results$top.table.list,
                                                    title = paste(subtypes_combination_nice, "FDR < 5%"))
ggsave(file.path("plots", paste0(subtypes_combination, "_jaccard_lineplot.pdf")),
                 plot = last_subtype.jacc.plot, width = 8.5, height = 4)

