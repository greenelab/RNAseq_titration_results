# J. Taroni Feb 2016, S. Foltz Feb 2022
# This script compares differential expression "silver standards" which are
# differential expression analysis results from standard pipelines (i.e.,
# log transformed 100% array data and RSEM counts 100% RNA-seq [processed with
# limma::voom]) with differential expression results from RNA-seq titrated data
# (0-100%) normalized various ways.
#
# Plot similarity of DEGs to silver standards (subtype vs. subtype only)
# 
# USAGE: Rscript 2A-plot_DEG_silver.subtype_vs_subtype.R --cancer_type --subtype_vs_subtype --supplementary

option_list <- list(
  optparse::make_option("--cancer_type",
                        default = NA_character_,
                        help = "Cancer type"),
  optparse::make_option("--subtype_vs_subtype",
                        default = NA_character_,
                        help = "Subtypes used in head-to-head comparison (comma-separated without space e.g. Type1,Type2)"),
  optparse::make_option("--supplementary",
                        action = "store_true",
                        default = FALSE,
                        help = "Save plot in plots/supplementary folder instead of plots/main")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
source(here::here("util/option_functions.R"))
check_options(opt)

# load libraries
suppressMessages(library(tidyverse))
source(here::here("util/color_blind_friendly_palette.R"))
source(here::here("util", "differential_expression_functions.R"))

# set options
cancer_type <- opt$cancer_type
subtype_vs_subtype <- opt$subtype_vs_subtype
supplementary <- opt$supplementary
two_subtypes <- as.vector(stringr::str_split(subtype_vs_subtype, pattern = ",", simplify = TRUE))
file_identifier <- str_c(cancer_type, "subtype", sep = "_") # we are only working with subtype models here

# define directories
plot.dir <- here::here("plots")
plot.main.dir <- file.path(plot.dir, "main")
plot.supp.dir <- file.path(plot.dir, "supplementary")
plot.data.dir <- file.path(plot.dir, "data")

output_directory <- file.path(ifelse(supplementary,
                                     plot.supp.dir,
                                     plot.main.dir))

#### plot Subtype v. Subtype results -------------------------------------------
subtypes_combination <- stringr::str_c(two_subtypes, collapse = "v")
subtypes_combination_nice <- stringr::str_c(two_subtypes, collapse = " vs. ")

last_subtype.silver_df <- read_tsv(
  file.path(plot.data.dir,
            paste0(file_identifier, "_titration_differential_exp_eBayes_fits_",
                   subtypes_combination, ".silver.tsv")),
  col_types = "dddcdc") %>%
  #col_types = "cdccd") %>%
  mutate(perc.seq = factor(perc.seq,
                           levels = seq(0, 100, 10))) %>%
  rename("Jaccard" = "jaccard",
         "Rand" = "rand",
         "Spearman" = "spearman",
         "Normalization" = "normalization",
         "Perc.Seq" = "perc.seq",
         "Platform" = "platform") %>%
  gather("Jaccard", "Rand", "Spearman", key = "measure", value = "value")

plot_obj_all3 <- PlotSilverStandardStats(
  last_subtype.silver_df,
  title = paste(cancer_type, subtypes_combination_nice, "FDR < 5%"))

plot_obj_just_jaccard <- PlotSilverStandardStats(
  last_subtype.silver_df %>%
    filter(measure == "Jaccard"),
  title = paste(cancer_type, subtypes_combination_nice, "FDR < 5%"),
  single_measure = TRUE)

ggsave(
  file.path(output_directory,
            paste0(file_identifier, "_silver_standard_similarity_lt5_",
                   subtypes_combination, ".pdf")),
  plot = plot_obj_all3,
  width = 7.25,
  height = 4)

ggsave(
  file.path(output_directory,
            paste0(file_identifier, "_silver_standard_similarity_lt5_",
                   subtypes_combination, ".jaccard_only.pdf")),
  plot = plot_obj_just_jaccard,
  width = 7.25,
  height = 3)