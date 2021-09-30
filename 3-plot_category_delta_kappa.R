# J. Taroni Jul 2016 -- modified by S. Foltz Sep 2021
# Based on 3-plot_category_kappa.R, but updated to calculate delta kappa.
# The purpose of this script is to plot delta kappa statistics from category
# predictions on hold-out data. It assumes there is a regular and null model
# and calculates delta kappa == regular kappa - null kappa. It should be run
# from the command line through the classifier_repeat_wrapper.R script or
# USAGE: Rscript 3-plot_category_delta_kappa.R

option_list <- list(
  optparse::make_option("--cancer_type",
                        default = NA_character_,
                        help = "Cancer type"),
  optparse::make_option("--predictor",
                        default = NA_character_,
                        help = "Predictor used")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
source(here::here("util/option_functions.R"))
check_options(opt)

# load libraries
suppressMessages(source(here::here("load_packages.R")))
source(here::here("util", "color_blind_friendly_palette.R"))

# set options
cancer_type <- opt$cancer_type
predictor <- opt$predictor
regular_file_identifier <- str_c(cancer_type, predictor, sep = "_")
null_file_identifier <- str_c(cancer_type, predictor, "null", sep = "_")

# define directories
plot.dir <- here::here("plots")
res.dir <- here::here("results")

# list array and seq files from results directory
lf <- list.files(res.dir, full.names = TRUE)
regular_array.files <- lf[grepl(paste0(regular_file_identifier,
                                       "_train_3_models_array_kappa_"), lf)]
null_array.files <- lf[grepl(paste0(null_file_identifier,
                                    "_train_3_models_array_kappa_"), lf)]
regular_seq.files <- lf[grepl(paste0(regular_file_identifier,
                                     "_train_3_models_seq_kappa_"), lf)]
null_seq.files <- lf[grepl(paste0(null_file_identifier,
                                  "_train_3_models_seq_kappa_"), lf)]

# define output files
plot.file.lead <- paste0(regular_file_identifier, "_train_3_models_delta_kappa_")
summary.df.filename <- file.path(res.dir,
                                 paste0(regular_file_identifier,
                                        "_train_3_models_delta_kappa_summary_table.tsv"))

# check that we have ordered pairs of regular and null files for array and seq
regular_array_seeds <- stringr::str_sub(regular_array.files, -8, -5)
null_array_seeds <- stringr::str_sub(null_array.files, -8, -5))
regular_seq_seeds <- stringr::str_sub(regular_seq.files, -8, -5)
null_seq_seeds <- stringr::str_sub(null_seq.files, -8, -5))
if (!(all(regular_array_seeds == null_array_seeds) &
      all(regular_seq_seeds == null_seq_seeds) &
      all(regular_array_seeds == regular_seq_seeds))) {
  stop("Seeds do not match in delta kappa plotting script.")
}

#### read in data --------------------------------------------------------------

# read in the tables that contain the kappa statistics for predictions on test
# data
regular_array.list <- list()  # initialize list that will hold all array tables
null_array.list <- list()  # initialize list that will hold all array tables
regular_seq.list <- list()  # initialize list that will hold all the RNA-seq tables
null_seq.list <- list()  # initialize list that will hold all the RNA-seq tables
for (i in 1:length(regular_array.files)) {
  regular_array.list[[i]] <- fread(regular_array.files[i], data.table = F)
  null_array.list[[i]] <- fread(null_array.files[i], data.table = F)
  regular_seq.list[[i]] <- fread(regular_seq.files[i], data.table = F)
  null_array.list[[i]] <- fread(null_array.files[i], data.table = F)
}

# calculate delta kappa values
delta_kappa_array.list <- list()
delta_kappa_seq.list <- list()
for (i in 1:length(regular_array.files)) {
  delta_kappa_array.list[[i]] <- regular_array.list[[i]] %>%
    left_join(null_array.list[[i]],
              by = c("perc.seq", "classifier", "norm.method")) %>%
    mutate(delta_kappa = kappa.x - kappa.y) %>% # regular kappa - null kappa
    select(delta_kappa, perc.seq, classifier, norm.method)
  delta_kappa_seq.list[[i]] <- regular_seq.list[[i]] %>%
    left_join(null_seq.list[[i]],
              by = c("perc.seq", "classifier", "norm.method")) %>%
    mutate(delta_kappa = kappa.x - kappa.y) %>% # regular kappa - null kappa
    select(delta_kappa, perc.seq, classifier, norm.method)
}

# combine all tables from each platform into a data.frame
array.df <- data.table::rbindlist(delta_kappa_array.list)
seq.df <- data.table::rbindlist(delta_kapp_seq.list)

#### plot test set results -----------------------------------------------------
# bind all kappa stats together
test.df <- cbind(rbind(array.df, seq.df),
                 c(rep("Microarray", nrow(array.df)),
                   rep("RNA-seq", nrow(seq.df))))
colnames(test.df) <- c("Delta.Kappa", "Perc.Seq", "Classifier",
                       "Normalization", "Platform")

# order %seq to display 0-100
test.df$Perc.Seq <- factor(test.df$Perc.Seq, levels = seq(0, 100, 10))

# recode model types
cls.recode.str <-
  "'glmnet' = 'LASSO'; 'rf' = 'Random Forest'; 'svm' = 'Linear SVM'"
test.df$Classifier <- car::recode(test.df$Classifier,
                                  recodes = cls.recode.str)

# capitalize norm methods
test.df$Normalization <- as.factor(toupper(test.df$Normalization))

# plot performance of a classifier/model type on all normalization method in a
# single plot
test.df$Classifier <- as.factor(test.df$Classifier)
cls.methods <- unique(test.df$Classifier)
for (cls in cls.methods) {
  plot.nm <- file.path(plot.dir,
                       paste0(plot.file.lead, # includes file_identifier
                              stringr::str_replace_all(cls,
                                                       pattern = " ",
                                                       replacement = "_"),
                              "_VIOLIN_test.pdf"))
  ggplot(test.df[which(test.df$Classifier == cls), ],
         aes(x = Perc.Seq, y = Delta.Kappa, color = Platform, fill = Platform)) +
    facet_wrap(~ Normalization, ncol = 5) +
    geom_violin(colour = "black", position = position_dodge(0.8),
                alpha = 0.2) +
    stat_summary(fun = median, geom = "line", aes(group = Platform),
                 position = position_dodge(0.6)) +
    stat_summary(fun = median, geom = "point", aes(group = Platform),
                 position = position_dodge(0.7), size = 1) +
    ggtitle(paste0(cancer_type, ": ", cls)) +
    xlab("% RNA-seq samples") +
    ylab("Delta Kappa") +
    theme_bw() +
    scale_colour_manual(values = cbPalette[2:3]) +
    theme(text = element_text(size = 18)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
  ggsave(plot.nm, plot = last_plot(), height = 3.5, width = 15)
}

# get summary data.frame + write to file
summary.df <- test.df %>%
  dplyr::group_by(Classifier, Normalization, Platform, Perc.Seq) %>%
  dplyr::summarise(Median = median(Delta.Kappa),
                   Mean = mean(Delta.Kappa),
                   SD = sd(Delta.Kappa)) %>%
  dplyr::ungroup()
readr::write_tsv(summary.df,
                 summary.df.filename)
