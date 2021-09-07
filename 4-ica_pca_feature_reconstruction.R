# J. Taroni Aug 2016
# The purpose of this script is to perform unsupervised learning on TCGA train-
# ing data (output of 1-normalize_titrated_data.R), PCA or ICA,
# and to transform test data into the training data reduced dimensional space,
# and back out ('reconstruction') and to then calculate the 'reconstruction
# error' (MASE).
#
# It should be run from the command line.
# USAGE: Rscript 4-ica_pca_feature_reconstruction.R --cancer_type --predictor --n_components --seed
# n_components refers to the number of components (PC/IC) that should be used
# for reconstruction.

option_list <- list(
  optparse::make_option("--cancer_type",
                        default = NULL,
                        help = "Cancer type"),
  optparse::make_option("--predictor",
                        default = NULL,
                        help = "Predictor used"),
  optparse::make_option("--n_components",
                        default = 50,
                        help = "Number of compenents [default: %default]"),
  optparse::make_option("--seed",
                        default = 346,
                        help = "Random seed [default: %default]")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
source(here::here("util/option_functions.R"))
check_options(opt)

# load libraries
suppressMessages(source(here::here("load_packages.R")))
source(here::here("util", "train_test_functions.R"))
source(here::here("util", "ICA_PCA_reconstruction_functions.R"))

# set options
cancer_type <- opt$cancer_type
predictor <- opt$predictor
file_identifier <- str_c(cancer_type, predictor, sep = "_")
n.comp <- as.integer(opt$n_components)

# set seed
initial.seed <- as.integer(opt$seed)
set.seed(initial.seed)
message(paste("\nInitial seed set to:", initial.seed))

# define directories
res.dir <- here::here("results")
norm.dir <- here::here("normalized_data")
mdl.dir <- here::here("models")
rcn.dir <- file.path(norm.dir, "reconstructed_data")
rcn.res.dir <- file.path(res.dir, "reconstructed_data")

# define input files
lf <- list.files(norm.dir, full.names = TRUE)
train.files <- lf[grepl(paste0(file_identifier,
                               "_array_seq_train_titrate_normalized_list_"), lf)]
test.files <- lf[grepl(paste0(file_identifier,
                              "_array_seq_test_data_normalized_list_"), lf)]

# parse filename seeds
filename.seeds <- substr(train.files,
                         (nchar(train.files)-7),
                         (nchar(train.files)-4))

# define output files
df.file.lead <- paste0(file_identifier,
                       "_reconstruction_error_", n.comp, "_components_")
mdl.file.lead <- paste0(file_identifier,
                        "_array_seq_train_", n.comp, "_components_object_")
rcn.file.lead <- paste0(file_identifier,
                        "_reconstructed_data_", n.comp, "_components_")

#### main ----------------------------------------------------------------------
platforms <- c("array", "seq")
#recon.methods <- c("ICA", "PCA")
recon.methods <- c("PCA") # July 2021 update to no longer run ICA

for (seed in filename.seeds) {
  rep.count <- grep(seed, filename.seeds)
  message(paste("\n\n#### RECONSTRUCTION ROUND",
                rep.count, "of", length(filename.seeds), "####\n\n"))

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
      # holdout out data in data.table format suitable for category prediction
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
