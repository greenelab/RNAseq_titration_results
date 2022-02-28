# J. Taroni Feb 2016, S. Foltz Feb 2022
# With small n data, plot comparison of SOMETHING
#
# USAGE: Rscript 3A-plot_small_n_differential_expression.R --cancer_type --subtype_vs_others --subtype_vs_subtype --output_directory --overlap_measure

option_list <- list(
  optparse::make_option("--cancer_type",
                        default = NA_character_,
                        help = "Cancer type"),
  optparse::make_option("--subtype_vs_others",
                        help = "Subtype used for comparison against all others."),
  optparse::make_option("--subtype_vs_subtype",
                        help = "Subtypes used in head-to-head comparison (comma-separated without space e.g. Type1,Type2)"),
  optparse::make_option("--output_directory",
                        default = NA_character_,
                        help = "Output directory of DEG overlap plot (absolute or relative path)."),
  optparse::make_option("--overlap_measure",
                        default = NA_character_,
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
  message("  Errors: must include --subtype_vs_others and/or --subtype_vs_subtype in plots/scripts/3A-plot_small_n_differential_expression.R.\n")
  stop()
}

# check that overlap measures requested are the ones present in data
if ("overlap_measure" %in% names(opt)) {
  overlap_measures <- sort(stringr::str_split(opt$overlap_measure,
                                              pattern = ",", simplify = TRUE))
  
  if (!all(overlap_measures %in% c("Jaccard", "Rand", "Spearman"))) {
    message("  Errors: --overlap_measure must be one or more of Jaccard, Rand, Spearman in plots/scripts/3A-plot_small_n_differential_expression.R.\n")
    stop()  
  }
}

# load libraries
suppressMessages(library(tidyverse))
source(here::here("util/color_blind_friendly_palette.R"))

# set options
cancer_type <- opt$cancer_type
file_identifier <- str_c(cancer_type, "subtype", sep = "_") # we are only working with subtype models here

# define directories
plot.dir <- here::here("plots")
plot.data.dir <- file.path(plot.dir, "data")
output_directory <- opt$output_directory

#### functions -----------------------------------------------------------------

DataSummary <- function(x) {
  # This function is supplied to ggplot2::stat_summary in order to plot the
  # median value of a vector as a point and the "confidence interval on the
  # median" used in notched boxplots as a vertical line. See boxplot.stats for
  # more information.
  m <- median(x)
  conf <- boxplot.stats(x)$conf
  ymin <- min(conf)
  ymax <- max(conf)
  return(c(y = m, ymin = ymin, ymax = ymax))
}

plot_small_n <- function(subtypes){
  # This function creates a panel of line plots faceted by 
  # % RNA-seq, the measure (Jaccard, Rand, Spearman), and normalization method
  
  subtypes_path <- str_c(subtypes, collapse = "v")
  subtypes_nice <- str_c(subtypes, collapse = " vs. ")
  measures_path <- str_c(overlap_measures, collapse = "_")
  
  input_filename <- file.path(plot.data.dir,
                              paste0(file_identifier,
                                     "_small_n_",
                                     subtypes_path,
                                     "_results.tsv"))
  
  output_filename <- file.path(
    output_directory,
    paste0(file_identifier, "_small_n_",
           measures_path, "_", subtypes_path, ".pdf"))
  
  using_single_measure <- length(overlap_measures) == 1
  
  #stats_df <- read_tsv(input_filename,
  #                     col_types = "ccdddddd") %>%
  #  filter(seq_prop %in% str_c(c(30, 50, 70), "% RNA-seq"),
  #         !is.na(value),
  #         measure %in% overlap_measures) %>%
  #  mutate(no.samples = factor(no.samples))

  stats_df <- read_tsv(input_filename,
                       col_types = "ccdddddd") %>%
    rename("Jaccard" = "jaccard", "Rand" = "rand", "Spearman" = "spearman") %>%
    gather(key = "measure", value = "value", Jaccard, Rand, Spearman) %>%
    mutate(seq_prop = factor(str_c(seq_prop, "% RNA-seq"),
                             levels = str_c(seq(0, 100, 10), "% RNA-seq"))) %>%
    filter(seq_prop %in% str_c(c(30, 50, 70), "% RNA-seq"),
           !is.na(value),
           measure %in% overlap_measures) %>%
    mutate(no.samples = factor(no.samples))
  
  plot_obj <- ggplot(stats_df, aes(x = no.samples,
                                   y = value,
                                   color = platform)) +
    stat_summary(fun = median,
                 geom = "line",
                 aes(group = platform),
                 position = position_dodge(0.7),
                 show.legend = FALSE) +
    stat_summary(fun = median, # this makes the point size consistent with other plots
                 geom = "point",
                 aes(group = platform),
                 position = position_dodge(0.7),
                 show.legend = FALSE) +
    stat_summary(fun.data = DataSummary, # this adds the error bars without median points
                 geom = "linerange",
                 aes(group = platform),
                 position = position_dodge(0.7),
                 show.legend = FALSE) +
    expand_limits(y = c(0,1)) +
    scale_y_continuous(breaks = seq(0, 1, 0.25)) +
    theme_bw() +
    labs(x = "Number of Samples from Each Subtype",
         y = ifelse(using_single_measure,
                    unique(stats_df$measure),
                    "Measure of Similarity"),
         title = paste(cancer_type, subtypes_nice, "FDR < 10%")) +
    scale_colour_manual(values = cbPalette[c(2, 3)])
  
  if (using_single_measure) {
    plot_obj <- plot_obj +
      facet_grid(normalization ~ seq_prop)
  } else {
    plot_obj <- plot_obj +
      facet_grid(measure + normalization ~ seq_prop)
  }
  
  ggsave(filename = output_filename,
         plot = plot_obj,
         width = 7.25,
         height = c(3,4,5)[length(overlap_measures)])
}

#### plot Subtype v. Other results ---------------------------------------------

if (!is.na(subtype_vs_others)) {
  
  subtypes <- c(subtype_vs_others, "Other")
  plot_small_n(subtypes)
}

#### plot Subtype v. Subtype results -------------------------------------------

if (!is.na(subtype_vs_subtype)) {
  
  subtypes <- as.vector(
    stringr::str_split(subtype_vs_subtype, pattern = ",", simplify = TRUE))
  plot_small_n(subtypes)
}
