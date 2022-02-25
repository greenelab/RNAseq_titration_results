# J. Taroni Feb 2016, S. Foltz Feb 2022
# This script compares differential expression "silver standards" which are
# differential expression analysis results from standard pipelines (i.e.,
# log transformed 100% array data and RSEM counts 100% RNA-seq [processed with
# limma::voom]) with differential expression results from RNA-seq titrated data
# (0-100%) normalized various ways.
#
# Plot the proportion of genes that are differentially expressed between conditions
#
# USAGE: Rscript 2A-plot_DEG_proportions.R --cancer_type --subtype_vs_others --subtype_vs_subtype

option_list <- list(
  optparse::make_option("--cancer_type",
                        default = NA_character_,
                        help = "Cancer type"),
  optparse::make_option("--subtype_vs_others",
                        #default = NA_character_, not required
                        help = "Subtype used for comparison against all others."),
  optparse::make_option("--subtype_vs_subtype",
                        #default = NA_character_, not required
                        help = "Subtypes used in head-to-head comparison (comma-separated without space e.g. Type1,Type2)"),
  optparse::make_option("--supplementary",
                        action = "store_true",
                        default = FALSE,
                        help = "Save plot in plots/supplementary folder instead of plots/main")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
source(here::here("util/option_functions.R"))
check_options(opt)

# at least one of --subtype_vs_others or --subtype_vs_subtype should be given
if (any(c("subtype_vs_others", "subtype_vs_subtype") %in% names(opt))) {
  
  subtype_vs_others <- NA # first assume option is not provided
  subtype_vs_subtype <- NA # then update as available below
  
  if ("subtype_vs_others" %in% names(opt)) {
    subtype_vs_others <- opt$subtype_vs_others
  }
  
  if ("subtype_vs_subtype" %in% names(opt)) {
    subtype_vs_subtype <- opt$subtype_vs_subtype
  }
  
} else {
  message("  Errors: must include --subtype_vs_others and/or --subtype_vs_subtype in plots/scripts/2A-plot_DEG_proportions.R.\n")
  stop()
}

# load libraries
suppressMessages(library(tidyverse))
source(here::here("util/color_blind_friendly_palette.R"))
source(here::here("util", "differential_expression_functions.R"))

# set options
cancer_type <- opt$cancer_type
supplementary <- opt$supplementary
file_identifier <- str_c(cancer_type, "subtype", sep = "_") # we are only working with subtype models here

# define directories
plot.dir <- here::here("plots")
plot.main.dir <- file.path(plot.dir, "main")
plot.supp.dir <- file.path(plot.dir, "supplementary")
plot.data.dir <- file.path(plot.dir, "data")

output_directory <- file.path(ifelse(supplementary,
                                     plot.supp.dir,
                                     plot.main.dir))

#### functions -----------------------------------------------------------------

plot_DEG_and_save <- function(subtypes # c(subtype, Other) or c(subtype1, subtype2)
                              #plot_data_dir = plot_data_dir,
                              #file_identifier = file_identifier,
                              #output_directory = output_directory,
                              #cancer_type = cancer_type
  ){
  
  subtypes_path <- str_c(subtypes, collapse = "v")
  subtypes_nice <- str_c(subtypes, collapse = " vs. ")
  
  input_filename <- file.path(
    plot.data.dir,
    paste0(file_identifier, "_titration_differential_exp_eBayes_fits_",
           subtypes_path, ".propDE.tsv"))
  
  output_filename <- file.path(
    output_directory,
    paste0(file_identifier, "_differential_expr_proportion_ltFDR5perc_",
           subtypes_path, ".pdf"))
  
  propDEG_df <- read_tsv(input_filename,
                         col_types = "dcd") %>%
    mutate(perc.seq = factor(perc.seq,
                             levels = seq(0, 100, 10)))
  
  # plot proportion of genes that are diff expressed
  plot_obj <- PlotProportionDE(propDEG_df,
                               subtypes = subtypes_nice,
                               cancer_type = cancer_type)
  
  #return(plot_obj)
  ggsave(output_filename,
         plot = plot_obj,
         width = 7.25,
         height = 4)
  
}

#### plot Subtype v. Other results ---------------------------------------------

if (!is.na(subtype_vs_others)) {
  
  subtypes <- c(subtype_vs_others, "Other")
  plot_DEG_and_save(subtypes)

}

#### plot Subtype v. Subtype results -------------------------------------------

if (!is.na(subtype_vs_subtype)) {
  
  subtypes <- as.vector(
    stringr::str_split(subtype_vs_subtype, pattern = ",", simplify = TRUE))
  plot_DEG_and_save(subtypes)

}