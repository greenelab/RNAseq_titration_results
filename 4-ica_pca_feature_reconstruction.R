# J. Taroni Aug 2016
# The purpose of this script is to perform unsupervised learning on BRCA train-
# ing data (output of 1-normalize_titrated_data.R), specifically ICA and PCA,
# and to transform test data into the training data reduced dimensional space,
# and back out ('reconstruction') and to then calculate the 'reconstruction 
# error' (RMSE).
# It should be run from the command line.
# USAGE: Rscript 4-ica_pca_feature_reconstruction.R <n.comp>
# n.comp refers to the number of components (PC/IC) that should be used
# for reconstruction.
# 

suppressMessages(source("load_packages.R"))
source("util/train_test_functions.R")
source("util/ICA_PCA_reconstruction_functions.R")

args <- commandArgs(trailingOnly = TRUE)
n.comp <- args[1]

kInitialSeed <- 346
set.seed(kInitialSeed)

res.dir <- "results/"
norm.dir <- "normalized_data/"
lf <- list.files(norm.dir, full.names = TRUE)
train.files <- lf[grepl("BRCA_array_seq_train_titrate_normalized_list_", lf)]
test.files <- lf[grepl("BRCA_array_seq_test_data_normalized_list_", lf)]
filename.seeds <- substr(train.files,
                         (nchar(train.files)-7),
                         (nchar(train.files)-4))
df.file.lead <- "BRCA_array_seq_ICA_PCA_reconstruction_error_"

#### main ----------------------------------------------------------------------
rep.count <- 1
for (seed in filename.seeds) {
  message(paste("\n\n#### RECONSTRUCTION ROUND", 
                rep.count, "of", length(filename.seeds), "####\n\n"))
  # read in data
  train.rds <- train.files[grepl(seed, train.files)]
  test.rds <- test.files[grepl(seed, test.files)]
  train.data <- readRDS(train.rds)
  test.data <- readRDS(test.rds)
  train.data <- RestructureNormList(train.data)
  train.data$tdm$`0` <- NULL
  train.data$tdm$`100` <- NULL
  #### PCA array test set ####
  message("\tArray PCA reconstruction")
  pca.array.recon <- CompAnalysisEvalWrapper(train.list = train.data, 
                                             test.list = test.data$array,
                                             no.comp = n.comp,
                                             method = "PCA",
                                             platform = "array")
  pca.array.recon <- cbind(pca.array.recon, rep("array", nrow(pca.array.recon)))
  colnames(pca.array.recon)[ncol(pca.array.recon)] <- "platform"
  #### ICA array test set ####
  message("\tArray ICA reconstruction")
  ica.array.recon <- CompAnalysisEvalWrapper(train.list = train.data, 
                                             test.list = test.data$array,
                                             no.comp = n.comp,
                                             method = "ICA",
                                             platform = "array")
  ica.array.recon <- cbind(ica.array.recon, rep("array", nrow(ica.array.recon)))
  colnames(ica.array.recon)[ncol(ica.array.recon)] <- "platform"
  
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
  pca.seq.recon <- cbind(pca.seq.recon, rep("seq", nrow(pca.seq.recon)))
  colnames(pca.seq.recon)[ncol(pca.seq.recon)] <- "platform"
  #### ICA seq test set ####
  message("\tRNA-seq ICA reconstruction")
  ica.seq.recon <- CompAnalysisEvalWrapper(train.list = train.data, 
                                           test.list = test.data$seq,
                                           no.comp = n.comp,
                                           method = "ICA",
                                           platform = "seq")
  ica.seq.recon <- cbind(ica.seq.recon, rep("seq", nrow(ica.seq.recon)))
  colnames(ica.seq.recon)[ncol(ica.seq.recon)] <- "platform"
  #### combine all results ####
  mstr.filename <- paste0(res.dir, df.file.lead, seed, ".tsv")
  master.df <- rbind(pca.array.recon, ica.array.recon, 
                     pca.seq.recon, ica.seq.recon)
  write.table(master.df, file = mstr.filename, sep="\t", row.names = F,
              quote = F)
  rep.count <- rep.count + 1  
}
