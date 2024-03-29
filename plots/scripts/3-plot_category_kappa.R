# S. Foltz Feb 2022
# This plots kappa values from category prediction

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
                        help = "Use delta kappa input data"),
  optparse::make_option("--output_directory",
                        default = NA_character_,
                        help = "Output directory for plot (absolute or relative path)"),
  optparse::make_option("--include_seurat",
                        action = "store_true",
                        default = FALSE,
                        help = "Include Seurat results in plot (default: FALSE)")
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
file_identifier <- str_c(cancer_type, predictor, sep = "_")
include_seurat <- opt$include_seurat

# define directories
plot.dir <- here::here("plots")
plot.data.dir <- here::here("plots/data")
output_directory <- opt$output_directory

# define input file
input_filename <- ifelse(null_model,
                         file.path(plot.data.dir,
                                   paste0(file_identifier,
                                          "_train_3_models_delta_kappa.tsv")),
                         file.path(plot.data.dir,
                                   paste0(file_identifier,
                                          "_train_3_models_kappa.tsv")))

# define output files
output_filename <- file.path(output_directory,
                             ifelse(null_model,
                                    paste0(file_identifier,
                                           "_train_3_models_delta_kappa.pdf"),
                                    paste0(file_identifier,
                                           "_train_3_models_kappa.pdf")))

# read in data
median_df <- read_tsv(input_filename,
                      col_types = "dddddccc") %>%
  mutate(Perc.Seq = factor(Perc.Seq,
                           levels = seq(0, 100, 10))) %>%
  group_by(Perc.Seq, Platform, Classifier, Normalization) %>%
  summarize(n_obs = n(),
            med = median(Kappa),
            IQR = quantile(Kappa, 0.75) - quantile(Kappa, 0.25),
            median_ci_upper = med + 1.58*IQR/sqrt(n_obs),
            median_ci_lower = med - 1.58*IQR/sqrt(n_obs),
            .groups = "drop")

kappa_df <- read_tsv(input_filename,
                     col_types = "dddddccc") %>%
  mutate(Perc.Seq = factor(Perc.Seq,
                           levels = seq(0, 100, 10)))

# default behavior: exclude (!include) seurat results
if (!include_seurat) {
  median_df <- median_df %>%
    filter(Normalization != "SEURAT") %>%
    droplevels()
  kappa_df <- kappa_df %>%
    filter(Normalization != "SEURAT") %>%
    droplevels()
}

# plot

plot_obj <- ggplot(median_df,
                   aes(x = Perc.Seq,
                       y = med, # median
                       color = Platform,
                       fill = Platform)) +
  facet_grid(rows = vars(Classifier),
             cols = vars(Normalization)) +
  geom_errorbar(aes(x = Perc.Seq,
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
             aes(x = Perc.Seq,
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
       y = ifelse(null_model,
                  "Delta Kappa",
                  "Kappa"),
       title = str_c(cancer_type, predictor, sep = " ")) +
  theme_bw() +
  scale_fill_manual(values = cbPalette[2:3]) +
  scale_colour_manual(values = cbPalette[2:3]) +
  theme(legend.position = "bottom",
        panel.grid.major = element_line(size = 0.25),
        panel.grid.minor = element_line(size = 0.25),
        strip.text.y = element_text(size = 7))

ggsave(output_filename,
       plot = plot_obj,
       height = 4,
       width = 7.25)
