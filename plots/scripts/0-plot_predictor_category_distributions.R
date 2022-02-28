# S. Foltz Feb 2022
# This plots the predictor category distribution for each seed

option_list <- list(
  optparse::make_option("--cancer_type",
                        default = NA_character_,
                        help = "Cancer type"),
  optparse::make_option("--predictor",
                        default = NA_character_,
                        help = "Predictor used"),
  optparse::make_option("--null_model",
                        action = "store_true",
                        default = FALSE,
                        help = "Use null model input data"),
  optparse::make_option("--plot_all_seeds",
                        action = "store_true",
                        default = FALSE,
                        help = "Plot all seeds instead of representative seed"),
  optparse::make_option("--output_directory",
                        default = NA_character_,
                        help = "Output directory for plot (absolute or relative path)")
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
null_model <- opt$null_model
plot_all_seeds <- opt$plot_all_seeds
file_identifier <- ifelse(null_model,
                          str_c(cancer_type, predictor, "null", sep = "_"),
                          str_c(cancer_type, predictor, sep = "_"))

# define directories
plot.dir <- here::here("plots")
plot.data.dir <- here::here("plots/data")
output_directory <- opt$output_directory

# list potential input files
input_files <- list.files(path = plot.data.dir,
                          pattern = paste0(file_identifier,
                                           ".dist_split_stacked_bar."),
                          full.names = TRUE)

# define input data
if (plot_all_seeds) { # read in all seed data to one data frame
  plot_df <- input_files %>%
    map(read_tsv,
        col_types = "cccd") %>%
    reduce(rbind)
} else { # default
  plot_df <- read_tsv(input_files[1],
                      col_types = "cccd")
  initial_seed <- plot_df %>%
    pull(initial_seed) %>%
    unique()
}

# define output file
if (plot_all_seeds) {
  category.distribution.plot <- file.path(output_directory,
                                          paste0(file_identifier,
                                                 ".dist_split_stacked_bar.",
                                                 "all_seeds",
                                                 ".pdf"))
} else { # default
  category.distribution.plot <- file.path(output_directory,
                                          paste0(file_identifier,
                                                 ".dist_split_stacked_bar.",
                                                 initial_seed,
                                                 ".pdf"))
}

# plot

plot_obj <- plot_df %>%
  ggplot(aes(x = fct_rev(split),
             fill = category)) +
  geom_bar() +
  scale_fill_manual(values = cbPalette) +
  labs(x = "Split",
       y = "Count",
       title = str_replace(file_identifier,
                           pattern = "_",
                           replacement = " "),
       fill = "Predictor") +
  theme_bw(base_size = 10)

if (plot_all_seeds) {
  plot_obj <- plot_obj +
    coord_flip() +
    facet_wrap(~ initial_seed,
               ncol = 5)
  ggsave(filename = category.distribution.plot,
         plot = plot_obj,
         height = 3,
         width = 7.25)
} else { # default
  ggsave(filename = category.distribution.plot,
         plot = plot_obj,
         height = 3,
         width = 3.5)
}
