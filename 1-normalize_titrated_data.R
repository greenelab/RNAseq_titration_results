# J. Taroni Jun 2016
# The purpose of this script is to read in TGCA array and sequencing data,
# already pre-processed to only include test tumor samples,
# (output of 0-expression_data_overlap_and_split.R) and to normalize
# the data.
# It should be run from the command line through the run_experiments.R script

option_list <- list(
  optparse::make_option("--cancer_type",
                        default = NA_character_,
                        help = "Cancer type"),
  optparse::make_option("--predictor",
                        default = NA_character_,
                        help = "Predictor used"),
  optparse::make_option("--seed1",
                        default = NA_integer_,
                        help = "Random seed"),
  optparse::make_option("--seed2",
                        default = NA_integer_,
                        help = "Random seed"),
  optparse::make_option("--null_model",
                        action = "store_true",
                        default = FALSE,
                        help = "Refer to models with permuted dependent variable (within subtype if predictor is a gene)"),
  optparse::make_option("--ncores",
                        default = NA_integer_,
                        help = "Set the number of cores to use")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
source(here::here("util/option_functions.R"))
check_options(opt)

# load libraries
suppressMessages(source(here::here("load_packages.R")))
source(here::here("util", "normalization_functions.R"))

# set options
cancer_type <- opt$cancer_type
predictor <- opt$predictor
null_model <- opt$null_model
file_identifier <- ifelse(null_model,
                          str_c(cancer_type, predictor, "null", sep = "_"),
                          str_c(cancer_type, predictor, sep = "_"))
ncores <- min(parallel::detectCores() - 1,
              opt$ncores,
              na.rm = TRUE)

# set seed
filename.seed <- as.integer(opt$seed1)
initial.seed <- as.integer(opt$seed2)
set.seed(initial.seed)

# define directories
data.dir <- here::here("data")
norm.data.dir <- here::here("normalized_data")
res.dir <- here::here("results")

# name input files
seq.file <- paste0(cancer_type, "RNASeq_matchedOnly_ordered.pcl")
array.file <- paste0(cancer_type, "array_matchedOnly_ordered.pcl")
train.test.file <- paste0(file_identifier,
                          "_matchedSamples_training_testing_split_labels_",
                          filename.seed, ".tsv")

# name output files
norm.test.object <- paste0(file_identifier,
                           "_array_seq_test_data_normalized_list_",
                           filename.seed, ".RDS")
norm.train.object <- paste0(file_identifier,
                            "_array_seq_train_titrate_normalized_list_",
                            filename.seed, ".RDS")

#### read in data --------------------------------------------------------------

seq.data <- fread(file.path(data.dir, seq.file), data.table = FALSE)
array.data <- fread(file.path(data.dir, array.file), data.table = FALSE)
sample.train.test <- fread(file.path(res.dir, train.test.file), data.table = FALSE)

#### split samples, titrate ----------------------------------------------------

train.sample.names <- as.character(sample.train.test$sample[
  which(sample.train.test$split == "train")])
test.sample.names <- as.character(sample.train.test$sample[
  which(sample.train.test$split == "test")])

# get samples for 'titration'
titration.seed <- sample(1:10000, 1)
message(paste("Random seed for titration:",
              titration.seed), appendLF = TRUE)

set.seed(titration.seed)
titrate.sample.list <- lapply(seq(0, 1, by = 0.1),
                              function(x) GetTitratedSampleNames(train.sample.names,
                                                                 x))
names(titrate.sample.list) <- as.character(seq(0, 100, by = 10))

# these samples will be the RNA-seq samples in any given 'titration' experiment
# remove rows that are equal to all ones -- for any combination + test data
# z-score processing will not work on such rows
seq.dt.list <- lapply(titrate.sample.list,
                      function(x) seq.data[, c(1, which(colnames(seq.data) %in% x))])
seq.dt.list[["test"]] <-
  seq.data[, c(1, which(colnames(seq.data) %in% test.sample.names))]
all.same.list <- lapply(seq.dt.list[2:12],
                        function(x){
                          vals <- x[, 2:ncol(x)]
                          indx <- which(apply(vals, 1, check_all_same))
                          return(indx)
                        } )
all.same.indx <- unique(unlist(all.same.list))
# if no rows have all same value (in previous lapply), all.same.indx is integer(0)
# subsetting data frames by -integer(0) results in no rows
# so check that integer vector has length > 0 before subsetting
if (length(all.same.indx) > 0) {
  array.data <- array.data[-all.same.indx, ]
  seq.data <- seq.data[-all.same.indx, ]
}

#### get datatables to mix -----------------------------------------------------

# get a list that contains an
# array data.table and seq data.table for each level of 'titration'
array.train <-
  data.table(array.data[,
                        c(1, which(colnames(array.data) %in% train.sample.names))])

seq.train <-
  data.table(seq.data[,
                      c(1, which(colnames(seq.data) %in% train.sample.names))])

titrate.mix.dt.list <- lapply(titrate.sample.list,
                              function(x) GetDataTablesForMixing(array.train,
                                                                 seq.train, x))

#### normalize train data ------------------------------------------------------

# initialize in the list to hold normalized data
norm.titrate.list <- list()

# single platform array normalization
norm.titrate.list[["0"]] <-
  SinglePlatformNormalizationWrapper(titrate.mix.dt.list[[1]]$array,
                                     platform = "array",
                                     add.untransformed = TRUE,
                                     add.qn.z = TRUE)

# parallel backend
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

# 'mixed' both platform normalization
norm.titrate.list[2:10] <-
  foreach(n = 2:10, .packages = "tidyverse") %dopar% {
    NormalizationWrapper(titrate.mix.dt.list[[n]]$array,
                         titrate.mix.dt.list[[n]]$seq,
                         add.untransformed = TRUE,
                         add.qn.z = TRUE,
                         add.cn = TRUE,
                         add.seurat.training = TRUE)
  }

# stop parallel backend
parallel::stopCluster(cl)
# sort out names
names(norm.titrate.list)[2:10] <- names(titrate.mix.dt.list)[2:10]

# single platform seq normalization
norm.titrate.list[["100"]] <-
  SinglePlatformNormalizationWrapper(titrate.mix.dt.list[[11]]$seq,
                                     platform = "seq",
                                     add.untransformed = TRUE,
                                     add.qn.z = TRUE)

#### normalize test data -------------------------------------------------------
array.test <-
  data.table(array.data[,
                        c(1, which(colnames(array.data) %in% test.sample.names))])
seq.test <-
  data.table(seq.data[, c(1, which(colnames(seq.data) %in% test.sample.names))])

# array normalization
array.test.norm.list <-
  SinglePlatformNormalizationWrapper(array.test,
                                     platform = "array",
                                     add.untransformed = TRUE,
                                     add.qn.z = TRUE,
                                     add.cn.test = TRUE,
                                     add.seurat.test = TRUE,
                                     training.list = norm.titrate.list)

# seq normalization
# initialize list to hold normalized seq data
seq.test.norm.list <- list()

# LOG normalization
seq.test.norm.list[["log"]] <- LOGSeqOnly(seq.test)
# NPN
seq.test.norm.list[["npn"]] <- NPNSingleDT(seq.test)

# start parallel backend
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

# QN -- requires reference data
# initialize list to hold QN data
seq.qn.list <- list()

# for 0% seq - use 0% LOG array data
seq.qn.list[["0"]] <- QNSingleWithRef(ref.dt = norm.titrate.list$`0`$log,
                                      targ.dt = seq.test)

# for 10-90% seq - use the "raw array" training data at each level of sequencing
# data (this is LOG data, but only the array samples)
seq.qn.list[2:10] <-
  foreach(i = 2:10) %dopar% {
    QNSingleWithRef(ref.dt = norm.titrate.list[[i]]$raw.array,
                    targ.dt = seq.test)
  }
names(seq.qn.list)[2:10] <- names(norm.titrate.list)[2:10]

# stop parallel back end
parallel::stopCluster(cl)

# QN 100% seq by itself (preProcessCore::normalize.quantiles)
seq.qn.list[["100"]] <- QNSingleDT(seq.test)

# add QN seq data to list of normalized test data
seq.test.norm.list[["qn"]] <- seq.qn.list
rm(seq.qn.list)

# start parallel backend
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

# QN-Z -- requires reference data
# initialize list to hold QN data
seq.qnz.list <- list()

# for 0% seq - use 0% LOG array data
seq.qnz.list[["0"]] <- QNZSingleWithRef(ref.dt = norm.titrate.list$`0`$log,
                                        targ.dt = seq.test)

# for 10-90% seq - use the "raw array" training data at each level of sequencing
# data (this is LOG data, but only the array samples)
seq.qnz.list[2:10] <-
  foreach(i = 2:10) %dopar% {
    QNZSingleWithRef(ref.dt = norm.titrate.list[[i]]$raw.array,
                    targ.dt = seq.test)
  }
names(seq.qnz.list)[2:10] <- names(norm.titrate.list)[2:10]

# stop parallel back end
parallel::stopCluster(cl)

# QNZ 100% seq by itself (preProcessCore::normalize.quantiles)
seq.qnz.list[["100"]] <- QNZSingleDT(seq.test)

# add QNZ seq data to list of normalized test data
seq.test.norm.list[["qn-z"]] <- seq.qnz.list
rm(seq.qnz.list)

# start parallel back end
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

# TDM normalization -- requires references
# initialize list to hold TDM data
seq.tdm.list <- list()

# for 0% seq - use 0% LOG array data
seq.tdm.list[["0"]] <- TDMSingleWithRef(ref.dt = norm.titrate.list$`0`$log,
                                        targ.dt = seq.test)
# for 10-90% seq - use the "raw array" training data at each level of sequencing
# data (this is LOG data, but only the array samples)
seq.tdm.list[2:10] <-
  foreach(i = 2:10) %dopar% {
    TDMSingleWithRef(ref.dt = norm.titrate.list[[i]]$raw.array,
                     targ.dt = seq.test)
  }
names(seq.tdm.list)[2:10] <- names(norm.titrate.list)[2:10]

# stop parallel backend
parallel::stopCluster(cl)

# 100% is not applicable for TDM
seq.tdm.list["100"] <- list(NULL)

# add TDM seq data to list of normalized test data
seq.test.norm.list[["tdm"]] <- seq.tdm.list
rm(seq.tdm.list)

# z-score seq test data
seq.test.norm.list[["z"]] <- ZScoreSingleDT(seq.test)

# untransformed seq test data
seq.test.norm.list[["un"]] <- seq.test

# CrossNorm RNA-seq test
# Rescale each column, quantile normalize, then rescale each row
seq.test.norm.list[["qn (cn)"]] <- rescale_datatable(seq.test,
                                                by_column = TRUE) %>%
  QNSingleDT(zero.to.one = TRUE)

# Seurat RNA-seq test
# for 10-90% seq - use the integrated training data at each %RNA-seq

# parallel backend
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

seq.seurat.list <- foreach(i = 2:10, .packages = "tidyverse") %dopar% { # 2:10 corresponds to 10%-90%
  
  if (!is.null(norm.titrate.list[[i]][["seurat_model"]])) {
    
    SeuratProjectPCATestData(seq.test,
                             norm.titrate.list[[i]][["seurat_model"]],
                             vbose = TRUE)
    
  } else {
    NULL
  }
  
}

names(seq.seurat.list) <- names(norm.titrate.list)[2:10] # 2:10 corresponds to 10%-90%

# stop parallel backend
parallel::stopCluster(cl)

# add Seurat RNA-seq test data to list of normalized test data
seq.test.norm.list[["seurat"]] <- seq.seurat.list
rm(seq.seurat.list)

# combine array and seq test data into a list
test.norm.list <- list(array = array.test.norm.list,
                       seq = seq.test.norm.list)

# save test data
saveRDS(test.norm.list, file = file.path(norm.data.dir, norm.test.object))

# save train data after removing Seurat models (just keep Seurat-normed data)
for (n in names(norm.titrate.list)) {
  if ("seurat_model" %in% names(norm.titrate.list[[n]])) {
    norm.titrate.list[[n]][["seurat_model"]] <- NULL
  }
}

saveRDS(norm.titrate.list, file = file.path(norm.data.dir, norm.train.object))
