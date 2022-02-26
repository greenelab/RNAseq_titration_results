# J. Taroni Feb 2016, S. Foltz Feb 2022
# This script compares differential expression "silver standards" which are
# differential expression analysis results from standard pipelines (i.e.,
# log transformed 100% array data and RSEM counts 100% RNA-seq [processed with
# limma::voom]) with differential expression results from RNA-seq titrated data
# (0-100%) normalized various ways.
#
# Plot the proportion of genes that are differentially expressed between conditions
# Plot similarity of DEGs to silver standards (subtype vs. others only)
#
# USAGE: Rscript 2A-plot_DEGs.R --cancer_type --subtype_vs_others --subtype_vs_subtype

option_list <- list(
  optparse::make_option("--cancer_type",
                        default = NA_character_,
                        help = "Cancer type"),
  optparse::make_option("--subtype_vs_others",
                        help = "Subtype used for comparison against all others."),
  optparse::make_option("--subtype_vs_subtype",
                        help = "Subtypes used in head-to-head comparison (comma-separated without space e.g. Type1,Type2)"),
  optparse::make_option("--proportion_output_directory",
                        help = "Output directory of DEG proportion plot. Include this option to plot DEG proportion plot."),
  optparse::make_option("--overlap_output_directory",
                        help = "Output directory of DEG overlap plot. Include this option to plot silver standard overlap plot."),
  optparse::make_option("--overlap_measure",
                        help = "Which overlap measures to include in silver standard overlap plot (comma-separated without space e.g. Jaccard,Rand,Spearman; must be one or more of Jaccard, Rand, Spearman)")
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

# at least one of --proportion_output_directory or --overlap_output_directory should be given
if (any(c("proportion_output_directory", "overlap_output_directory") %in% names(opt))) {
  
  proportion_output_directory <- NA # first assume option is not provided
  overlap_output_directory <- NA # then update as available below
  
  if ("proportion_output_directory" %in% names(opt)) {
    proportion_output_directory <- opt$proportion_output_directory
    plot_proportion <- TRUE

  }
  
  if ("overlap_output_directory" %in% names(opt)) {
    overlap_output_directory <- opt$overlap_output_directory
    plot_overlap <- TRUE
    
    # check that overlap measures requested are the ones present in data
    if ("overlap_measure" %in% names(opt)) {
      overlap_measures <- sort(stringr::str_split(opt$overlap_measure,
                                            pattern = ",", simplify = TRUE))
      
      if (!all(overlap_measures %in% c("Jaccard", "Rand", "Spearman"))) {
        message("  Errors: --overlap_measure must be one or more of Jaccard, Rand, Spearman in plots/scripts/2A-plot_DEG_proportions.R.\n")
        stop()  
      }
      
    } else {
      message("  Errors: must include --overlap_measure with --overlap_output_directory in plots/scripts/2A-plot_DEG_proportions.R.\n")
      stop()
    }
  }
  
} else {
  message("  Errors: must include --proportion_output_directory and/or --overlap_output_directory in plots/scripts/2A-plot_DEG_proportions.R.\n")
  stop()
}

# load libraries
suppressMessages(library(tidyverse))
source(here::here("util/color_blind_friendly_palette.R"))
source(here::here("util", "differential_expression_functions.R"))

# set options
cancer_type <- opt$cancer_type
file_identifier <- str_c(cancer_type, "subtype", sep = "_") # we are only working with subtype models here

# define directories
plot.dir <- here::here("plots")
plot.data.dir <- file.path(plot.dir, "data")

#### functions -----------------------------------------------------------------

plot_DEG_proportions <- function(subtypes){
  
  subtypes_path <- str_c(subtypes, collapse = "v")
  subtypes_nice <- str_c(subtypes, collapse = " vs. ")
  
  input_filename <- file.path(
    plot.data.dir,
    paste0(file_identifier, "_titration_differential_exp_eBayes_fits_",
           subtypes_path, ".propDE.tsv"))
  
  output_filename <- file.path(
    proportion_output_directory,
    paste0(file_identifier, "_differential_expr_proportion_lt5_",
           subtypes_path, ".pdf"))
  
  propDEG_df <- read_tsv(input_filename,
                         col_types = "dcd") %>%
    mutate(perc.seq = factor(perc.seq,
                             levels = seq(0, 100, 10)))
  
  # plot proportion of genes that are diff expressed
  plot_obj <- PlotProportionDE(propDEG_df,
                               subtypes = subtypes_nice,
                               cancer_type = cancer_type)
  
  ggsave(output_filename,
         plot = plot_obj,
         width = 7.25,
         height = 4)
  
}

plot_silver_overlap <- function(subtypes){
  
  subtypes_path <- str_c(subtypes, collapse = "v")
  subtypes_nice <- str_c(subtypes, collapse = " vs. ")
  measures_path <- str_c(overlap_measures, collapse = "_")
  
  input_filename <- file.path(
    plot.data.dir,
    paste0(file_identifier, "_titration_differential_exp_eBayes_fits_",
           subtypes_path, ".silver.tsv"))
  
  output_filename <- file.path(
    overlap_output_directory,
    paste0(file_identifier, "_silver_standard_similarity_lt5_",
           measures_path, "_", subtypes_path, ".pdf"))
  
  silver_df <- read_tsv(input_filename,
                        col_types = "cdccd") %>%
    mutate(Perc.Seq = factor(Perc.Seq,
                             levels = seq(0, 100, 10))) %>%
    filter(measure %in% overlap_measures)
  
  using_single_measure <- length(overlap_measures) == 1
  
  plot_obj <- PlotSilverStandardStats(
    silver_df,
    title = paste(cancer_type, subtypes_nice, " FDR < 5%"),
    single_measure = using_single_measure)

  ggsave(
    output_filename,
    plot = plot_obj,
    width = 7.25,
    height = c(3,4,5)[length(overlap_measures)]
  )
}

#### plot Subtype v. Other results ---------------------------------------------

if (!is.na(subtype_vs_others)) {
  
  subtypes <- c(subtype_vs_others, "Other")

  if (plot_proportion) {
    plot_DEG_proportions(subtypes)  
  }
  
  if (plot_overlap) {
    plot_silver_overlap(subtypes)  
  }
}

#### plot Subtype v. Subtype results -------------------------------------------

if (!is.na(subtype_vs_subtype)) {
  
  subtypes <- as.vector(
    stringr::str_split(subtype_vs_subtype, pattern = ",", simplify = TRUE))

  if (plot_proportion) {
    plot_DEG_proportions(subtypes)  
  }
  
  if (plot_overlap) {
    plot_silver_overlap(subtypes)  
  }

}
