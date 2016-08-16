# J. Taroni Aug 2016
# The purpose of this script is to perform unsupervised learning on BRCA train-
# ing data (output of 1-normalize_titrated_data.R), specifically ICA and PCA,
# and to transform test data into the training data reduced dimensional space,
# and back out ('reconstruction') and to then calculate the 'reconstruction 
# error' (MASE).
# It should be run from the command line.
# USAGE: Rscript 4-ica_pca_feature_reconstruction.R <n.comp>
# n.comp refers to the number of components (PC/IC) that should be used
# for reconstruction.
# 

suppressMessages(source("load_packages.R"))
source("util/train_test_functions.R")
source("util/ICA_PCA_reconstruction_functions.R")

args <- commandArgs(trailingOnly = TRUE)
n.comp <- as.integer(args[1])

kInitialSeed <- 346
set.seed(kInitialSeed)

res.dir <- "results/"
norm.dir <- "normalized_data/"
mdl.dir <- "models/"
lf <- list.files(norm.dir, full.names = TRUE)
train.files <- lf[grepl("BRCA_array_seq_train_titrate_normalized_list_", lf)]
test.files <- lf[grepl("BRCA_array_seq_test_data_normalized_list_", lf)]
filename.seeds <- substr(train.files,
                         (nchar(train.files)-7),
                         (nchar(train.files)-4))
df.file.lead <- paste0("BRCA_array_seq_ICA_PCA_reconstruction_error_", n.comp,
                       "_components_")
mdl.file.lead <- paste0("BRCA_array_seq_train_ICA_PCA_objects_", n.comp,
                        "_components_")

#### main ----------------------------------------------------------------------
rep.count <- 1
for (seed in filename.seeds) {
  message(paste("\n\n#### RECONSTRUCTION ROUND", 
                rep.count, "of", length(filename.seeds), "####\n\n"))
  
  #### read in data ####
  train.rds <- train.files[grepl(seed, train.files)]
  test.rds <- test.files[grepl(seed, test.files)]
  train.data <- readRDS(train.rds)
  test.data <- readRDS(test.rds)
  train.data <- RestructureNormList(train.data)
  
  # initialize list for all PCA & ICA objects
  comp.list <- list()
  
  # for array data, get rid of TDM at 0 & 100% seq levels
  train.data$tdm$`0` <- NULL
  train.data$tdm$`100` <- NULL
  
  #### PCA array test set ####
  message("\tArray PCA reconstruction")
  pca.array.recon <- CompAnalysisEvalWrapper(train.list = train.data, 
                                             test.list = test.data$array,
                                             no.comp = n.comp,
                                             method = "PCA",
                                             platform = "array")

  # save prcomp objects into a list
  comp.list$array$PCA <- pca.array.recon$COMP

  #### ICA array test set ####
  message("\tArray ICA reconstruction")
  ica.array.recon <- CompAnalysisEvalWrapper(train.list = train.data, 
                                             test.list = test.data$array,
                                             no.comp = n.comp,
                                             method = "ICA",
                                             platform = "array")

  # save fastICA objects into a list
  comp.list$array$ICA <- ica.array.recon$COMP

  # for seq reconstruction
  train.data$tdm$`0` <- train.data$log$`0`
  train.data$tdm <- train.data$tdm[c(10, 1:9)]
  #### PCA seq test set ####
  message("\tRNA-seq PCA reconstruction")
  pca.seq.recon <- CompAnalysisEvalWrapper(train.list = train.data, 
                                           test.list = test.data$seq,
                                           no.comp = n.comp,
                                           method = "PCA",
                                           platform = "seq")
  
  # save prcomp objects into a list
  comp.list$seq$PCA <- pca.seq.recon$COMP

  #### ICA seq test set ####
  message("\tRNA-seq ICA reconstruction")
  ica.seq.recon <- CompAnalysisEvalWrapper(train.list = train.data, 
                                           test.list = test.data$seq,
                                           no.comp = n.comp,
                                           method = "ICA",
                                           platform = "seq")
 
  # save fastICA objects into a list
  comp.list$seq$ICA <- ica.array.recon$COMP
  
  #### combine all results ####
  master.df <- rbind(pca.array.recon$MASE, 
                     ica.array.recon$MASE, 
                     pca.seq.recon$MASE, 
                     ica.seq.recon$MASE)
  
  mstr.filename <- paste0(res.dir, df.file.lead, seed, ".tsv")
  write.table(master.df, file = mstr.filename, sep="\t", row.names = F,
              quote = F)
  
  # save prcomp and fastICA objects to model directory
  comp.rds.name <- paste0(mdl.dir, mdl.file.lead, seed, ".RDS")
  saveRDS(comp.list, file = comp.rds.name)
  
  rep.count <- rep.count + 1  

}
