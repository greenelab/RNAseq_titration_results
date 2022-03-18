# S. Foltz Feb 2022
# This plots reconstruction kappa values from PCA reconstruction

option_list <- list(
  optparse::make_option("--cancer_type",
                        default = NA_character_,
                        help = "Cancer type"),
  optparse::make_option("--predictor",
                        default = NA_character_,
                        help = "Predictor used"),
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
file_identifier <- str_c(cancer_type, predictor, sep = "_")

# define directories
plot.dir <- here::here("plots")
plot.data.dir <- here::here("plots/data")
output_directory <- opt$output_directory

# define input file
input_filename <- file.path(plot.data.dir,
                            paste0(file_identifier,
                                   "_kappa_reconstructed_data.tsv"))

# define output files
output_filename <- file.path(output_directory,
                             paste0(file_identifier,
                                    "_kappa_reconstructed.pdf"))

# read in data

plot_df <- readr::read_tsv(input_filename,
                           col_types = "dcccdcc") %>%
  mutate(Perc.seq = factor(Perc.seq,
                           levels = seq(0, 100, 10)))

# for each normalization method, plot kappa stats
plot_obj <- ggplot(plot_df,
                   aes(x = Perc.seq,
                       y = Kappa,
                       color = Platform,
                       fill = Platform)) +
  facet_grid(rows = vars(Classifier),
             cols = vars(Normalization)) +
  geom_violin(position = position_dodge(0.7),
              alpha = 0.25,
              show.legend = FALSE) +
  stat_summary(fun = median,
               geom = "line",
               aes(group = Platform),
               position = position_dodge(0.7)) +
  stat_summary(fun = median,
               geom = "point",
               aes(group = Platform),
               position = position_dodge(0.7),
               size = 1,
               shape = 16) +
  expand_limits(y = 1) +
  scale_x_discrete(labels = c("0", "", "", "", "",
                              "50", "", "", "", "",
                              "100")) + 
  labs(x = "% RNA-seq Samples in Training Data",
       color = "Test Data Platform",
       fill = "Test Data Platform",
       y = "Kappa",
       title = str_c("PCA reconstruction of",
                     cancer_type, predictor, sep = " ")) +
  theme_bw() +
  scale_colour_manual(values = cbPalette[2:3]) +
  theme(legend.position = "bottom")

ggsave(output_filename,
       plot = plot_obj,
       height = 5,
       width = 7.5)
