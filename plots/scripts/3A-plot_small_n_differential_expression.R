# J. Taroni Feb 2016, S. Foltz Feb 2022
# With small n data, plot comparison of SOMETHING
#
# USAGE: Rscript 3A-plot_small_n_differential_expression.R --cancer_type --subtype_vs_subtype --supplementary

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

plot_small_n <- function(df, single_measure = FALSE){
  # This function creates a panel of line plots faceted by 
  # % RNA-seq, the measure (Jaccard, Rand, Spearman), and normalization method
  # If single measure == TRUE, then only one measure should be present in df
  
  plot_obj <- ggplot(df, aes(x = no.samples,
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
    labs(x = "Number of samples",
         y = ifelse(single_measure,
                    unique(df$measure),
                    "Measure of similarity"),
         title = paste(cancer_type, subtypes_combination_nice, "FDR < 10%")) +
    scale_colour_manual(values = cbPalette[c(2, 3)])
  
  if (single_measure) {
    plot_obj <- plot_obj +
      facet_grid(normalization ~ seq_prop)
  } else {
    plot_obj <- plot_obj +
      facet_grid(measure + normalization ~ seq_prop)
  }
  
  return(plot_obj)
}

#### plot Subtype v. Subtype results -------------------------------------------
subtypes_combination <- stringr::str_c(two_subtypes, collapse = "v")
subtypes_combination_nice <- stringr::str_c(two_subtypes, collapse = " vs. ")

stats_df <- read_tsv(file.path(plot.data.dir,
                               paste0(file_identifier,
                                      "_small_n_",
                                      subtypes_combination,
                                      "_results.tsv")),
                     col_types = "ccdddddd") %>%
  filter(seq_prop %in% str_c(c(30, 50, 70), "% RNA-seq"),
         !is.na(value)) %>%
  mutate(no.samples = factor(no.samples))

# plot with plot_small_n() and save files

ggsave(filename = here::here(output_directory,
                             paste0(file_identifier, "_small_n_",
                                    subtypes_combination, ".pdf")),
       plot = plot_small_n(stats_df),
       width = 7.25,
       height = 6)

ggsave(filename = here::here(output_directory,
                             paste0(file_identifier, "_small_n_",
                                    subtypes_combination, ".jaccard_only.pdf")),
       plot = plot_small_n(stats_df %>%
                             filter(measure == "Jaccard"),
                           single_measure = TRUE),
       width = 7.25,
       height = 3)
