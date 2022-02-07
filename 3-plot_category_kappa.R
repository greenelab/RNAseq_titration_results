# J. Taroni Jul 2016
# The purpose of this script is to plot Kappa statistics from category
# predictions on hold-out data. It should be run from the command line
# through the classifier_repeat_wrapper.R script or alternatively
# USAGE: Rscript 3-plot_category_kappa.R

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
                        help = "Use null model as baseline for plotting delta kappa")
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
null_model <- opt$null_model
file_identifier <- str_c(cancer_type, predictor, sep = "_")

# define directories
plot.dir <- here::here("plots")
plot.data.dir <- file.path(plot.dir, "data")
res.dir <- here::here("results")

# list array and seq files from results directory
lf <- list.files(res.dir, full.names = TRUE)
array.files <- lf[grepl(paste0(file_identifier,
                               "_train_3_models_array_kappa_"), lf)]
seq.files <- lf[grepl(paste0(file_identifier,
                             "_train_3_models_seq_kappa_"), lf)]
if (null_model) {
  null_array.files <- lf[grepl(paste0(file_identifier,
                                      "_null_train_3_models_array_kappa_"), lf)]
  null_seq.files <- lf[grepl(paste0(file_identifier,
                                    "_null_train_3_models_seq_kappa_"), lf)]
  
  # check that we have ordered pairs of regular and null files for array and seq
  array_seeds <- stringr::str_sub(array.files, -8, -5)
  null_array_seeds <- stringr::str_sub(null_array.files, -8, -5)
  seq_seeds <- stringr::str_sub(seq.files, -8, -5)
  null_seq_seeds <- stringr::str_sub(null_seq.files, -8, -5)
  if (!(all(array_seeds == null_array_seeds) &
        all(seq_seeds == null_seq_seeds))) {
    stop("Array or seq seeds do not match in delta kappa plotting script.")
  }
  
}

# define output files
plot.file.lead <- ifelse(null_model,
                         paste0(file_identifier, "_train_3_models_delta_kappa_"),
                         paste0(file_identifier, "_train_3_models_kappa_"))

test.df.filename <- ifelse(null_model,
                              file.path(plot.data.dir,
                                        paste0(file_identifier,
                                               "_train_3_models_delta_kappa.tsv")),
                              file.path(plot.data.dir,
                                        paste0(file_identifier,
                                               "_train_3_models_kappa.tsv")))

summary.df.filename <- ifelse(null_model,
                              file.path(plot.data.dir,
                                        paste0(file_identifier,
                                               "_train_3_models_delta_kappa_summary_table.tsv")),
                              file.path(plot.data.dir,
                                 paste0(file_identifier,
                                        "_train_3_models_kappa_summary_table.tsv")))

#### read in data --------------------------------------------------------------

# read in the tables that contain the kappa statistics for predictions on test
# data
array.list <- list()  # initialize list that will hold all array tables
seq.list <- list()  # initialize list that will hold all the RNA-seq tables
for (file_index in 1:length(array.files)) {
  array.list[[file_index]] <- fread(array.files[file_index], data.table = F)
  seq.list[[file_index]] <- fread(seq.files[file_index], data.table = F)
}

if (null_model) {
  null_array.list <- list() # initialize list that will hold null array tables
  null_seq.list <- list() # initialize list that will hold null RNA-seq tables
  for (null_file_index in 1:length(null_array.files)) {
    null_array.list[[null_file_index]] <- fread(null_array.files[null_file_index], data.table = F)
    null_seq.list[[null_file_index]] <- fread(null_seq.files[null_file_index], data.table = F)
  }
  
  # calculate delta kappa values
  delta_kappa_array.list <- list() # list for delta kappa array values
  delta_kappa_seq.list <- list() # list for delta kappa seq values
  for (pair_index in 1:length(array.files)) {
    delta_kappa_array.list[[pair_index]] <- array.list[[pair_index]] %>%
      left_join(null_array.list[[pair_index]],
                by = c("perc.seq", "classifier", "norm.method"),
                suffix = c(".true", ".null")) %>%
      mutate(delta_kappa = kappa.true - kappa.null) %>% # regular kappa - null kappa
      select(delta_kappa, perc.seq, classifier, norm.method)
    delta_kappa_seq.list[[pair_index]] <- seq.list[[pair_index]] %>%
      left_join(null_seq.list[[pair_index]],
                by = c("perc.seq", "classifier", "norm.method"),
                suffix = c(".true", ".null")) %>%
      mutate(delta_kappa = kappa.true - kappa.null) %>% # regular kappa - null kappa
      select(delta_kappa, perc.seq, classifier, norm.method)
  }
  
}

# combine all tables from each platform into a data.frame
# cannot use ifelse() because return value must be same dim as conditional test
if (null_model) {
  array.df <- data.table::rbindlist(delta_kappa_array.list)
  seq.df <- data.table::rbindlist(delta_kappa_seq.list)
} else {
  array.df <- data.table::rbindlist(array.list)
  seq.df <- data.table::rbindlist(seq.list)
}

#### plot test set results -----------------------------------------------------
# bind all kappa stats together
test.df <- cbind(rbind(array.df, seq.df),
                 c(rep("Microarray", nrow(array.df)),
                   rep("RNA-seq", nrow(seq.df))))

colnames(test.df) <- c("Kappa", "Perc.Seq", "Classifier",
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
test.df$Classifier <- as.factor(test.df$Classifier)

readr::write_tsv(test.df,
                 test.df.filename) # delta or not delta in file name

# get summary data.frame + write to file
summary.df <- test.df %>%
  dplyr::group_by(Classifier, Normalization, Platform, Perc.Seq) %>%
  dplyr::summarise(Median = median(Kappa),
                   Mean = mean(Kappa),
                   SD = sd(Kappa)) %>%
  dplyr::ungroup()

readr::write_tsv(summary.df,
                 summary.df.filename) # delta or not delta in file name

# plot performance of a classifier/model type on all normalization method in a
# single plot

cls.methods <- unique(test.df$Classifier)
for (cls in cls.methods) {
  plot.nm <- file.path(plot.dir,
                       paste0(plot.file.lead, # includes file_identifier and delta/not delta
                              stringr::str_replace_all(cls,
                                                       pattern = " ",
                                                       replacement = "_"),
                              "_VIOLIN_test.pdf"))
  
  plot.data.nm <- file.path(plot.data.dir,
                       paste0(plot.file.lead, # includes file_identifier and delta/not delta
                              stringr::str_replace_all(cls,
                                                       pattern = " ",
                                                       replacement = "_"),
                              "_VIOLIN_test.tsv"))
  
  plot_df <- test.df[which(test.df$Classifier == cls), ]
  
  ggplot(plot_df,
         aes(x = Perc.Seq, y = Kappa, color = Platform, fill = Platform)) +
    facet_wrap(~ Normalization, ncol = 5) +
    geom_violin(colour = "black", position = position_dodge(0.8),
                alpha = 0.2) +
    stat_summary(fun = median, geom = "line", aes(group = Platform),
                 position = position_dodge(0.6)) +
    stat_summary(fun = median, geom = "point", aes(group = Platform),
                 position = position_dodge(0.7), size = 1) +
    ggtitle(paste0(cancer_type, ": ", cls)) +
    xlab("% RNA-seq samples") +
    ylab(ifelse(null_model,
                "Delta Kappa",
                "Kappa")) +
    theme_bw() +
    scale_colour_manual(values = cbPalette[2:3]) +
    theme(text = element_text(size = 18)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
  ggsave(plot.nm, plot = last_plot(), height = 3.5, width = 15)
}
