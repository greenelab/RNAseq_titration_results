# J. Taroni Aug 2016
# The purpose of this script is to perform unsupervised learning on BRCA train-
# ing data (output of 1-normalize_titrated_data.R), specifically ICA and PCA,
# and to transform test data into the training data reduced dimensional space,
# and back out ('reconstruction') and to then calculate the 'reconstruction 
# error' (MASE).
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

set.seed(initial.seed)

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
 recon.methods <- c("ICA", "PCA")

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
  
  # for array data (first to be analyzed)
  # get rid of TDM at 0 & 100% seq levels
  train.data$tdm$`0` <- NULL
  train.data$tdm$`100` <- NULL

  for (plt in platforms) {
    # deal with TDM normalization in training data
    if (plt == "seq" & is.null(train.data$tdm$`0`)){
      # At the 0% RNA-seq level, TDM RNA-seq test data is transformed using the
      # log-transformed 100% array data on the reference. So, use 
      # log-transformed 100% array data as the training set for evaluating the 
      # TDM method at 0% RNA-seq level. 
      train.data$tdm$`0` <- train.data$log$`0`
      train.data$tdm <- train.data$tdm[c(10, 1:9)]
    } 
    # evaluate platform test set PCA and ICA
    for (rcn in recon.methods) {
      message(paste("\t", plt, rcn, "reconstruction"))
      results <- 
        CompAnalysisEvalWrapper(train.list = train.data, 
                                test.list = test.data[[plt]],
                                no.comp = n.comp,
                                method = rcn,
                                platform = plt)
      # add error data.frame to list that holds all error data.frames
      # for this round of reconstruction
      message("\t\t writing error data.frame to file...")
      df.file.name <- paste0(df.file.lead, "_", plt, "_", rcn, "_", seed, 
      		   ".tsv")
      write.table(results$error.df, 
                  file = file.path(rcn.res.dir, df.file.name), quote = F,
      			      row.names = F, sep = "\t")
      # save prcomp or fastICA objects to model directory
      message("\t\t saving prcomp/fastICA objects RDS...")
      comp.rds.name <- paste0(mdl.file.lead, plt, "_", rcn, "_", 
                              seed, ".RDS")
      saveRDS(results$COMP, file = file.path(mdl.dir, comp.rds.name))
      # save reconstructed data to reconstructed data directory
      message("\t\t saving reconstructed data RDS...")
      recon.rds.name <- paste0(rcn.file.lead, plt, "_", rcn, "_", 
                               seed, ".RDS")
      saveRDS(results$RECON, file = file.path(rcn.dir, recon.rds.name))
    
      rm(results)	
      gc()
	   
    }  
    
  }
  
  rm(train.data, test.data)
  gc()

}
