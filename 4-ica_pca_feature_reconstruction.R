# J. Taroni Aug 2016
# The purpose of this script is to perform unsupervised learning on BRCA train-
# ing data (output of 1-normalize_titrated_data.R), specifically ICA and PCA,
# and to transform test data into the training data reduced dimensional space,
# and back out ('reconstruction') and to then calculate the 'reconstruction
# error' (MASE).
#
# It should be run from the command line.
# USAGE: Rscript 4-ica_pca_feature_reconstruction.R <n.comp> <initial.seed>
# n.comp refers to the number of components (PC/IC) that should be used
# for reconstruction.
#

suppressMessages(source("load_packages.R"))
source(file.path("util", "train_test_functions.R"))
source(file.path("util", "ICA_PCA_reconstruction_functions.R"))

args <- commandArgs(trailingOnly = TRUE)
n.comp <- as.integer(args[1])
initial.seed <- as.integer(args[2])
if (is.na(initial.seed)) {
  message("\nInitial seed set to default: 346")
  initial.seed <- 346
} else {
  message(paste("\nInitial seed set to:", initial.seed))
}

res.dir <- "results"
norm.dir <- "normalized_data"
mdl.dir <- "models"
rcn.dir <- file.path("normalized_data", "reconstructed_data")
rcn.res.dir <- file.path(res.dir, "reconstructed_data")
lf <- list.files(norm.dir, full.names = TRUE)
train.files <- lf[grepl("BRCA_array_seq_train_titrate_normalized_list_", lf)]
test.files <- lf[grepl("BRCA_array_seq_test_data_normalized_list_", lf)]
filename.seeds <- substr(train.files,
                         (nchar(train.files)-7),
                         (nchar(train.files)-4))
df.file.lead <- paste0("BRCA_reconstruction_error_", n.comp,
                       "_components_")
mdl.file.lead <- paste0("BRCA_array_seq_train_", n.comp,
                        "_components_object_")
rcn.file.lead <- paste0("BRCA_reconstructed_data_", n.comp,
                        "_components_")

#### main ----------------------------------------------------------------------
platforms <- c("array", "seq")
#recon.methods <- c("ICA", "PCA")
recon.methods <- c("PCA") # July 2021 update to no longer run ICA

for (seed in filename.seeds) {
  rep.count <- grep(seed, filename.seeds)
  message(paste("\n\n#### RECONSTRUCTION ROUND",
                rep.count, "of", length(filename.seeds), "####\n\n"))

  # set seed
  set.seed(initial.seed)

  #### read in data ####
  message("Reading in data...")
  train.rds <- train.files[grepl(seed, train.files)]
  test.rds <- test.files[grepl(seed, test.files)]
  train.data <- readRDS(train.rds)
  test.data <- readRDS(test.rds)
  train.data <- RestructureNormList(train.data)

  # get rid of TDM at 0 & 100% seq levels
  train.data$tdm$`0` <- NULL
  train.data$tdm$`100` <- NULL

  # for each method to be used for reconstruction
  for (rcn in recon.methods) {
    message(paste("  ", rcn, "on training set"))

    # perform reconstruction method on the training data
    train.comp.list <- TrainSetCompAnalysis(train.list = train.data,
                                            num.comp = n.comp,
                                            comp.method = rcn)

    # write the component objects to file in the models directory
    comp.rds.name <- paste0(mdl.file.lead, rcn, "_", seed, ".RDS")
    saveRDS(train.comp.list, file = file.path(mdl.dir, comp.rds.name))

    # reconstruction on the holdout data
    for(plt in platforms) { # for the two platforms -- microarray and RNA-seq
      message(paste("\t Performing", plt, "reconstruction"))
      if (plt == "seq" & is.null(train.comp.list$tdm$`0`)){
        # At the 0% RNA-seq level, TDM RNA-seq test data is transformed using the
        # log-transformed 100% array data on the reference. So, use
        # log-transformed 100% array data as the training set for evaluating the
        # TDM method at 0% RNA-seq level.
        train.comp.list$tdm$`0` <- train.comp.list$log$`0`
        train.comp.list$tdm <- train.comp.list$tdm[c(10, 1:9)]
      }

      # perform the reconstruction experiment, which will return reconstructed
      # holdout out data in data.table format suitable for subtype prediction
      # and calculate the reconstruction error (MASE) to be returned as a
      # data.frame
      results <- ReconstructionWrapper(train.list = train.comp.list,
                                       test.list = test.data[[plt]],
                                       num.comps = n.comp)

      # save recon objects
      message("\t   Saving reconstructed holdout data")
      recon.rds <- paste0(rcn.file.lead, rcn, "_", plt, "_", seed,".RDS")
      saveRDS(results$recon, file = file.path(rcn.dir, recon.rds))

      # write error data.frame to file
      error.df <- results$mase.df
      error.df <- cbind(error.df, rep(plt, nrow(error.df)))
      colnames(error.df)[ncol(error.df)] <- "platform"
      error.df.name <- paste0(df.file.lead, rcn, "_", plt, "_", seed, ".tsv")
      message("\t   Saving MASE data.frame")
      write.table(error.df, file = file.path(rcn.res.dir, error.df.name),
                  quote = FALSE, row.names = FALSE, sep = "\t")

      rm(results, error.df)
      gc()

    }

  }

  rm(train.data, test.data)
  gc()

}
