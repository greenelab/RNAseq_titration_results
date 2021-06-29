# J. Taroni Oct 2016
# The purpose of this script is to perform PAM50 subtype prediction
# (from 2-train_test_brca_subtype.R) on test/holdout data that has been
# reconstructed using the components from ICA and PCA on training data (the
# output of 4-ica_pca_feature_reconstruction.R). It outputs a list of
# confusionMatrix objects and a data.frame of Kappa statistics from these
# predictions.
# It should be run from the command line.
# USAGE: Rscript 5-predict_subtype_reconstructed_data.R

suppressMessages(source("load_packages.R"))
source(file.path("util", "train_test_functions.R"))

mdl.dir <- "models"
rcn.dir <- file.path("normalized_data", "reconstructed_data")
res.dir <- "results"
rcn.res.dir <- file.path(res.dir, "reconstructed_data")

model.files <- list.files(mdl.dir, full.names = TRUE)
supervised.model.files <- model.files[grep("3_models", model.files)]

recon.files <- list.files(rcn.dir, full.names = TRUE)

# get filename.seeds (identifiers for each of 10 replicates)
# from the reconstructed data files
filename.seeds <- unique(substr(recon.files,
                                (nchar(recon.files)-7),
                                (nchar(recon.files)-4)))

cm.file.lead <-
  "BRCA_subtype_prediction_reconstructed_data_confusionMatrices_"
kap.file.lead <- sub("confusionMatrices", "kappa", cm.file.lead)

#### main ----------------------------------------------------------------------

platforms <- c("array", "seq")
recon.methods <- c("ICA", "PCA")

for (seed in filename.seeds) {

  # error-handling -- want to make sure there is a corresponding supervised
  # model file to current reconstructed data file (seed)
  check.model.file <- any(grepl(seed, supervised.model.files))
  if (!check.model.file) {
    stop(paste("There is no corresponding supervised model file for
          filename.seed:", seed))
  }

  rep.count <- grep(seed, filename.seeds)
  message(paste("\n\n#### SUBTYPE PREDICTION",
                rep.count, "of", length(filename.seeds), "####\n\n"))

  # read in supervised models (LASSO, linear SVM, random forest)
  train.rds <-
    supervised.model.files[grep(seed, supervised.model.files)]
  train.list <- readRDS(train.rds)

  # need to read in corresponding sample.df
  sample.df.file <-
    file.path(res.dir,
           paste0(
             "BRCA_matchedSamples_PAM50Array_training_testing_split_labels_",
              seed, ".tsv"))
  sample.df <- data.table::fread(sample.df.file, data.table = F)
  sample.df$subtype <- as.factor(sample.df$subtype)

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

      # get confusionMatrix objects
      plt.list[[rcn]] <- PredictWrapper(train.model.list = train.list,
                                        pred.list = recon.list,
                                        sample.df = sample.df,
                                        return.kap = FALSE,
                                        run.parallel = FALSE)
      # get just Kappa statistic
      plt.kap.list[[rcn]] <- PredictWrapper(train.model.list = train.list,
                                            pred.list = recon.list,
                                            sample.df = sample.df,
                                            return.kap = TRUE,
                                            run.parallel = FALSE)


      # remove reconstructed datat
      rm(recon.list)
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
