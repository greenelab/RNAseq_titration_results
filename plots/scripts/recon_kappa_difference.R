# S. Foltz Mar 2022
# This compares kappa values from category prediction with/out reconstruction

option_list <- list(
  optparse::make_option("--cancer_type",
                        default = NA_character_,
                        help = "Cancer type"),
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
predictor <- "subtype"
file_identifier <- str_c(cancer_type, predictor, sep = "_")

# define directories
plot.dir <- here::here("plots")
plot.data.dir <- here::here("plots/data")
output_directory <- opt$output_directory

# define input file
without_recon_input_filename <- file.path(plot.data.dir,
                                          paste0(file_identifier,
                                                 "_train_3_models_kappa.tsv"))
with_recon_input_filename <- file.path(plot.data.dir,
                                       paste0(file_identifier,
                                              "_kappa_reconstructed_data.tsv"))

# define output files
output_filename <- file.path(output_directory,
                             paste0(file_identifier,
                                           "_kappa_reconstruction_difference.pdf"))

# read in data
without_df <- read_tsv(without_recon_input_filename,
                       col_types = "ddccc") %>%
  mutate(Perc.Seq = factor(Perc.Seq,
                           levels = seq(0, 100, 10)))

with_df <- read_tsv(with_recon_input_filename,
                    col_types = "dcccdcc") %>%
  mutate(Perc.Seq = factor(Perc.seq,
                           levels = seq(0, 100, 10)))

# get data summary (median kappa at each setting)

without_summary_df <- without_df %>%
  group_by(Perc.Seq, Classifier, Normalization, Platform) %>%
  summarize(median_without = median(Kappa)) %>%
  ungroup()

with_summary_df <- with_df %>%
  filter(Reconstruction == "PCA",
         Measure == "kappa") %>%
  group_by(Perc.Seq, Classifier, Normalization, Platform) %>%
  summarize(median_with = median(Kappa)) %>%
  ungroup()

# combined data frames and calculate difference in median kappas

joint_df <- without_summary_df %>%
  left_join(with_summary_df,
            by = c("Perc.Seq", "Classifier", "Normalization", "Platform")) %>%
  mutate(kappa_difference = median_without - median_with) %>%
  filter(!is.na(kappa_difference))

# plot

plot_obj <- ggplot(joint_df,
                   aes(x = Perc.Seq,
                       y = kappa_difference,
                       color = Platform,
                       fill = Platform)) +
  facet_grid(rows = vars(Classifier),
             cols = vars(Normalization)) +
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
  scale_x_discrete(labels = c("0", "", "", "", "",
                              "50", "", "", "", "",
                              "100")) + 
  labs(x = "% RNA-seq Samples in Training Data",
       color = "Test Data Platform",
       fill = "Test Data Platform",
       y = "Difference in Kappa\n(No Reconstruction - Reconstruction)",
       title = str_c(cancer_type, predictor, "(reconstruction difference)", sep = " ")) +
  theme_bw() +
  scale_colour_manual(values = cbPalette[2:3]) +
  theme(legend.position = "bottom")

ggsave(output_filename,
       plot = plot_obj,
       height = 4,
       width = 7.25)
