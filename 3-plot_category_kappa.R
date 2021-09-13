# J. Taroni Jul 2016
# The purpose of this script is to plot Kappa statistics from category
# predictions on hold-out data. It should be run from the command line
# through the classifier_repeat_wrapper.R script or alternatively
# USAGE: Rscript 3-plot_category_kappa.R

option_list <- list(
  optparse::make_option("--cancer_type",
                        default = NULL,
                        help = "Cancer type"),
  optparse::make_option("--predictor",
                        default = NULL,
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
file_identifier <- str_c(cancer_type, predictor, sep = "_")

# define directories
plot.dir <- here::here("plots")
res.dir <- here::here("results")

# list array and seq files from results directory
lf <- list.files(res.dir, full.names = TRUE)
array.files <- lf[grepl(paste0(file_identifier,
                               "_train_3_models_array_kappa_"), lf)]
seq.files <- lf[grepl(paste0(file_identifier,
                             "_train_3_models_seq_kappa_"), lf)]

# define output files
plot.file.lead <- paste0(file_identifier, "_train_3_models_kappa_")
summary.df.filename <- file.path(res.dir,
                                 paste0(file_identifier,
                                        "_train_3_models_summary_table.tsv"))

#### read in data --------------------------------------------------------------

# read in the tables that contain the kappa statistics for predictions on test
# data
array.list <- list()  # initialize list that will hold all array tables
seq.list <- list()  # initialize list that will hold all the RNA-seq tables
for (i in 1:length(array.files)) {
  array.list[[i]] <- fread(array.files[i], data.table = F)
  seq.list[[i]] <- fread(seq.files[i], data.table = F)
}

# combine all tables from each platform into a data.frame
array.df <- data.table::rbindlist(array.list)
seq.df <- data.table::rbindlist(seq.list)
rm(list = c("array.list", "seq.list"))

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
    theme_bw() +
    scale_colour_manual(values = cbPalette[2:3]) +
    theme(text = element_text(size = 18)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
  ggsave(plot.nm, plot = last_plot(), height = 3.5, width = 15)
}

# get summary data.frame + write to file
summary.df <- test.df %>%
  dplyr::group_by(Classifier, Normalization, Platform, Perc.Seq) %>%
  dplyr::summarise(Median = median(Kappa),
                   Mean = mean(Kappa),
                   SD = sd(Kappa),
                   .groups = "drop")
readr::write_tsv(summary.df,
                 summary.df.filename)
