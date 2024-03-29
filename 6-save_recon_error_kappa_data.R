# J. Taroni Oct 2016
# This script plots reconstruction errors (MASE and RMSE) from
# 4-ica_pca_feature_reconstruction.R and the Kappa statistics associated with
# predictions on reconstructed data from 5-predict_category_reconstructed_data.R
# as violin plots, respectively.
# USAGE: Rscript 6-plot_recon_error_kappa.R --cancer_type --predictor --null_model

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
                        help = "Refer to models with permuted dependent variable (within subtype if predictor is a gene)")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
source(here::here("util/option_functions.R"))
check_options(opt)

# load libraries
suppressMessages(library(tidyverse))
source(here::here("util", "color_blind_friendly_palette.R"))

# set options
cancer_type <- opt$cancer_type
predictor <- opt$predictor
null_model <- opt$null_model
file_identifier <- ifelse(null_model,
                          str_c(cancer_type, predictor, "null", sep = "_"),
                          str_c(cancer_type, predictor, sep = "_"))

# define directories
plot.dir <- here::here("plots")
plot.data.dir <- file.path(plot.dir, "data")
rcn.res.dir <- here::here("results", "reconstructed_data")

# define input files
# pattern = "kappa" captures a downstream output file if this script is rerun
# pattern = "kappa_[0-9]+.tsv" captures the intended filenames including seeds between 1:10000
kappa.df.files <- list.files(rcn.res.dir,
                             pattern = paste0(file_identifier,
                                              "_prediction_reconstructed_data_kappa_[0-9]+.tsv"),
                                              full.names = TRUE)
error.files <- list.files(rcn.res.dir,
                          pattern = paste0(file_identifier,
                                           "_reconstruction_error"),
                          full.names = TRUE)

# define output files
kap.plot.file.lead <- file.path(plot.dir, paste0(file_identifier, "_kappa_reconstructed_data_"))
err.plot.file.lead <- file.path(plot.dir, paste0(file_identifier, "_reconstruction_error_"))
kap.plot.data.file <- file.path(plot.data.dir, paste0(file_identifier, "_kappa_reconstructed_data.tsv"))
err.plot.data.file <- file.path(plot.data.dir, paste0(file_identifier, "_reconstruction_error.tsv"))

#### plot kappa stats ----------------------------------------------------------

# read in kappa data.frames from each replicate and bind -- line plot with
# boxplot "confidence intervals"
kappa.df.list <- list()
fl.iter <- 1
for (fl in kappa.df.files) {
  kappa.df.list[[fl.iter]] <- data.table::fread(fl, data.table = FALSE)
  fl.iter <- fl.iter + 1
}
kappa.master.df <- as.data.frame(data.table::rbindlist(kappa.df.list))
rm(kappa.df.list)

# order Perc.seq so line plot displays 0-100
kappa.master.df$Perc.seq <- factor(kappa.master.df$Perc.seq,
                                   levels = seq(0, 100, 10))

# rename classifiers

cls.recode.str <-
  "'glmnet' = 'LASSO'; 'rf' = 'Random Forest'; 'svm' = 'Linear SVM'"
kappa.master.df$Classifier <- car::recode(kappa.master.df$Classifier,
                                          recodes = cls.recode.str)
kappa.master.df$Classifier <- as.factor(kappa.master.df$Classifier)

# get norm and reconstruction methods as factors
kappa.master.df$Normalization <- stringr::str_to_upper(kappa.master.df$Normalization)
kappa.master.df$Normalization <- as.factor(kappa.master.df$Normalization)
kappa.master.df$Reconstruction <- as.factor(kappa.master.df$Reconstruction)

# rename platforms
plt.recode.str <-
  "'array' = 'Microarray'; 'seq' = 'RNA-seq'"
kappa.master.df$Platform <- car::recode(kappa.master.df$Platform,
                                        recodes = plt.recode.str)
kappa.master.df$Platform <- as.factor(kappa.master.df$Platform)


write.table(kappa.master.df,
            file = kap.plot.data.file,
            quote = FALSE, sep = "\t", row.names = FALSE)

# get summary data.frame + write to file
kappa.summary.df <-
  kappa.master.df %>%
  dplyr::group_by(Classifier, Normalization, Platform, Perc.seq) %>%
  dplyr::summarise(Median = median(Kappa),
                   Mean = mean(Kappa),
                   SD = sd(Kappa),
                   .groups = "drop")
readr::write_tsv(kappa.summary.df,
                 file.path(rcn.res.dir,
                           paste0(file_identifier,
                                  "_kappa_reconstructed_data_summary_table.tsv")))

rm(kappa.master.df)

#### plot error measures -------------------------------------------------------

# read in error measure data.frames from each replicate and bind -- violin plot
error.df.list <- list()
for(fl.iter in seq_along(error.files)){
  error.df.list[[fl.iter]] <- data.table::fread(error.files[fl.iter],
                                                data.table = FALSE)
}
error.master.df <- as.data.frame(data.table::rbindlist(error.df.list))
rm(error.df.list)

# order perc.seq so plot displays 0-100
error.master.df$perc.seq <- factor(error.master.df$perc.seq,
                                   levels = seq(0, 100, by = 10))

# get norm and reconstruction methods as factors
error.master.df$norm.method <- stringr::str_to_upper(error.master.df$norm.method)
error.master.df$norm.method <- as.factor(error.master.df$norm.method)

# rename platforms -- same as above for kappa data.frame
error.master.df$platform <- car::recode(error.master.df$platform,
                                        recodes = plt.recode.str)
error.master.df$platform <- as.factor(error.master.df$platform)

# reconstruction method as factor
error.master.df$comp.method <- as.factor(error.master.df$comp.method)

# take the average of each genes error across replicates
error.mean.df <- error.master.df %>%
  dplyr::group_by(gene, perc.seq, norm.method, comp.method, platform) %>%
  dplyr::summarise(mean_mase = mean(MASE),
                   .groups = "drop")
rm(error.master.df)
colnames(error.mean.df) <- c("Gene", "Perc.seq", "Normalization",
                             "Method", "Platform", "Mean_Value")

write.table(error.mean.df,
            file = err.plot.data.file,
            quote = FALSE, sep = "\t", row.names = FALSE)
