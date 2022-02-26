# S. Foltz Feb 2022
# This plots the rate of return for significant PLIER pathways
# for data coming from different normalization methods and titration levels.

option_list <- list(
  optparse::make_option("--cancer_type",
                        default = NA_character_,
                        help = "Cancer type"),
  optparse::make_option("--predictor",
                        default = NA_character_,
                        help = "Predictor used"),
  optparse::make_option("--output_directory",
                        default = NA_character_,
                        help = "Save plot to this directory")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
source(here::here("util/option_functions.R"))
check_options(opt)

# load libraries
suppressMessages(library(tidyverse))
source(here::here("util/color_blind_friendly_palette.R"))

# set options
cancer_type <- opt$cancer_type
predictor <- opt$predictor
output_directory <- opt$output_directory
file_identifier <- str_c(cancer_type, predictor, sep = "_")

# define directories
plot.dir <- here::here("plots")
plot.data.dir <- here::here("plots/data")

# define input file

plot_data_filename <- file.path(
  plot.data.dir,
  str_c(file_identifier, "_PLIER_jaccard.tsv")
)

# define output files

output_filename <- file.path(output_directory,
                             str_c(file_identifier,
                                   "_PLIER_jaccard.pdf"))

# Read in data
jaccard_df <- read_tsv(plot_data_filename,
                       col_types = "dddddcdcdd") %>%
  filter(pseq == 50) %>%
  mutate(sample_size = case_when(nmeth == "array_only" ~ "Single Platform",
                                 nmeth == "seq_only" ~ "Single Platform",
                                 TRUE ~ "Combine Array and RNA-seq (double sample size)"),
         nmeth = case_when(nmeth == "array_only" ~ "LOG\nArray",
                           nmeth == "seq_only" ~ "LOG\nRNA-seq",
                           TRUE ~ str_to_upper(nmeth)))
# Plot results
set.seed(1) # using jitter

plot_obj <- jaccard_df %>%
  ggplot(aes(x = nmeth,
             y = jaccard)) +
  geom_violin(draw_quantiles = .5,
              scale = "width") +
  geom_jitter(shape = 16,
              alpha = 0.5,
              height = 0,
              width = 0.1) +
  expand_limits(y = 0) +
  facet_grid(. ~ fct_rev(sample_size),
             scales = "free_x",
             space='free') +
  ggtitle(cancer_type) +
  xlab("Normalization Method") +
  ylab("Proportion of Pathways Significant") +
  theme_bw()

ggsave(output_filename,
       plot = plot_obj,
       height = 4, width = 7.25)
