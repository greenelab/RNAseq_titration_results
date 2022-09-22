# J. Taroni Oct 2016
# The purpose of this script is to perform category prediction
# (from 2-train_test_brca_category.R) on test/holdout data that has been
# reconstructed using the components from PCA on training data (the
# output of 4-ica_pca_feature_reconstruction.R). It outputs a list of
# confusionMatrix objects and a data.frame of Kappa statistics from these
# predictions.
# It should be run from the command line.
# USAGE: Rscript 5-predict_category_reconstructed_data.R --cancer_type --predictor --null_model

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
suppressMessages(source(here::here("load_packages.R")))
source(here::here("util", "train_test_functions.R"))

# set options
cancer_type <- opt$cancer_type
predictor <- opt$predictor
null_model <- opt$null_model
file_identifier <- ifelse(null_model,
                          str_c(cancer_type, predictor, "null", sep = "_"),
                          str_c(cancer_type, predictor, sep = "_"))

# define directories
mdl.dir <- here::here("models")
norm.dir <- here::here("normalized_data")
res.dir <- here::here("results")
rcn.dir <- file.path(norm.dir, "reconstructed_data")
rcn.res.dir <- file.path(res.dir, "reconstructed_data")

# define input files
supervised.model.files <- list.files(mdl.dir,
                                     pattern = paste0(file_identifier,
                                                      "_train_3_models"),
                                     full.names = TRUE)
recon.files <- list.files(rcn.dir,
                          pattern = paste0(file_identifier,
                                           "_reconstructed_data_"),
                          full.names = TRUE)

# get filename.seeds (identifiers for each replicate)
# from the reconstructed data files
filename.seeds <- unique(substr(recon.files,
                                (nchar(recon.files)-7),
                                (nchar(recon.files)-4)))

# define output files
cm.file.lead <- paste0(file_identifier,
                       "_prediction_reconstructed_data_confusionMatrices_")
kap.file.lead <- paste0(file_identifier,
                        "_prediction_reconstructed_data_kappa_")

#### main ----------------------------------------------------------------------

platforms <- c("array", "seq")
recon.methods <- c("PCA")

for (seed in filename.seeds) {

  # error-handling -- want to make sure there is a corresponding supervised
  # model file to current reconstructed data file (seed)
  check.model.file <- any(grepl(seed, supervised.model.files))
  if (!check.model.file) {
    stop(paste("There is no corresponding supervised model file for
          filename.seed:", seed))
  }

  rep.count <- grep(seed, filename.seeds)
  message(paste("\n\n#### CATEGORY PREDICTION",
                rep.count, "of", length(filename.seeds), "####\n\n"))

  # read in supervised models (LASSO, linear SVM, random forest)
  train.rds <- supervised.model.files[grep(seed, supervised.model.files)]
  train.list <- readRDS(train.rds)
  
  # remove Seurat data from this analysis because it is already in reduced space
  train.list$seurat <- NULL
  
  # need to read in corresponding sample.df
  sample.df.file <-
    file.path(res.dir,
              paste0(file_identifier,
                     "_matchedSamples_training_testing_split_labels_",
                     seed, ".tsv"))
  sample.df <- data.table::fread(sample.df.file, data.table = F)
  sample.df$category <- as.factor(sample.df$category)

  # initialize list to hold confusionMatrices & kappa statistics
  kappa.list <- list()
  cm.list <- list()

  for (plt in platforms) {
    plt.list <- list()
    plt.kap.list <- list()
    for (rcn in recon.methods) {

      # read in reconstructed data from current platform and reconstruction
      # method
      file.identifier <- paste(rcn, plt, seed, sep = "_")
      recon.rds <- recon.files[grep(file.identifier, recon.files)]
      recon.list <- readRDS(recon.rds)

      # return list of confusion matrix objects AND kappa statistics
      cm_kappa.list <- PredictWrapper(train.model.list = train.list,
                                      pred.list = recon.list,
                                      sample.df = sample.df,
                                      only.kap = FALSE,
                                      run.parallel = FALSE)

      # get confusionMatrix objects
      plt.list[[rcn]] <- cm_kappa.list$confusion_matrix_objects

      # get kappa statistics
      plt.kap.list[[rcn]] <- cm_kappa.list$kappa_statistics

      # remove reconstructed data and cm_kappa_list
      rm(recon.list, cm_kappa.list)
      gc()

    }

    cm.list[[plt]] <- plt.list
    kappa.list[[plt]] <- plt.kap.list

  }

  # save confusion matrices
  cm.file.name <- file.path(rcn.res.dir, paste0(cm.file.lead, seed, ".RDS"))
  saveRDS(cm.list, file = cm.file.name)

  # get kappa stats into data.frame from nested list and save as data.frame
  kappa.df <- reshape2::melt(kappa.list)
  colnames(kappa.df) <- c("Perc.seq", "Classifier", "Normalization",
                          "Measure", "Kappa", "Reconstruction", "Platform")
  kap.file.name <- file.path(rcn.res.dir, paste0(kap.file.lead, seed, ".tsv"))
  write.table(kappa.df, file = kap.file.name, row.names = F, quote = F,
              sep = "\t")

  rm(train.list, kappa.list, cm.list)
  gc()

}
