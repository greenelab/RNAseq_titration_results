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

median_df <- readr::read_tsv(input_filename,
                             col_types = "dcccdcc") %>%
  mutate(Perc.seq = factor(Perc.seq,
                           levels = seq(0, 100, 10))) %>%
  group_by(Perc.seq, Platform, Classifier, Normalization) %>%
  summarize(n_obs = n(),
            med = median(Kappa),
            IQR = quantile(Kappa, 0.75) - quantile(Kappa, 0.25),
            median_ci_upper = med + 1.58*IQR/sqrt(n_obs),
            median_ci_lower = med - 1.58*IQR/sqrt(n_obs),
            .groups = "drop")

kappa_df <- read_tsv(input_filename,
                     col_types = "dcccdcc") %>%
  mutate(Perc.Seq = factor(Perc.seq,
                           levels = seq(0, 100, 10)))

# for each normalization method, plot kappa stats
plot_obj <- ggplot(plot_df,
                   aes(x = Perc.seq,
                       y = med, # median
                       color = Platform,
                       fill = Platform)) +
  facet_grid(rows = vars(Classifier),
             cols = vars(Normalization)) +
  geom_errorbar(aes(x = Perc.seq,
                    ymin = median_ci_lower,
                    ymax = median_ci_upper),
                size = 0.25,
                width = 0.5,
                position = position_dodge(0.7)) +
  geom_line(aes(group = Platform),
            size = 0.5,
            position = position_dodge(0.7)) + 
  geom_point(shape = 16,
             size = 0.5,
             show.legend = FALSE,
             position = position_dodge(0.7)) +
  geom_point(data = kappa_df,
             aes(x = Perc.seq,
                 y = Kappa,
                 color = Platform,
                 fill = Platform),
             alpha = 0.5,
             size = 0.25,
             shape = 16,
             position = position_dodge(0.7),
             show.legend = FALSE) +
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
  theme(legend.position = "bottom",
        panel.grid.major = element_line(size = 0.25),
        panel.grid.minor = element_line(size = 0.25))

ggsave(output_filename,
       plot = plot_obj,
       height = 4,
       width = 7.25)
