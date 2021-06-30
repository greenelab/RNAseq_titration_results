# J. Taroni Jun 2016
# The purpose of this script is to read in TGCA array and sequencing data,
# already pre-processed to only include test tumor samples,
# (output of 0-expression_data_overlap_and_split.R) and to normalize
# the data.
# It should be run from the command line through the run_experiments.R script

suppressMessages(source("load_packages.R"))
source(file.path("util", "normalization_functions.R"))

args <- commandArgs(trailingOnly = TRUE)
filename.seed <- as.integer(args[1])
initial.seed <- as.integer(args[2])
set.seed(initial.seed)

data.dir <- "data"
norm.data.dir <- "normalized_data"
seq.file <- "BRCARNASeq_matchedOnly_ordered.pcl"
array.file <- "BRCAarray_matchedOnly_ordered.pcl"
norm.test.object <-
  paste0("BRCA_array_seq_test_data_normalized_list_", filename.seed, ".RDS")
norm.train.object <-
  paste0("BRCA_array_seq_train_titrate_normalized_list_", filename.seed, ".RDS")

res.dir <- "results"
train.test.file <-
  paste0("BRCA_matchedSamples_PAM50Array_training_testing_split_labels_",
         filename.seed, ".tsv")
train.test.labels <- file.path(res.dir, train.test.file)

#### read in data --------------------------------------------------------------

seq.data <- fread(file.path(data.dir, seq.file), data.table = FALSE)
array.data <- fread(file.path(data.dir, array.file), data.table = FALSE)
sample.train.test <- read.delim(train.test.labels)

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
all1.list <- lapply(seq.dt.list[2:12],
                    function(x){
                      vals <- x[, 2:ncol(x)]
                      indx <- which(apply(vals, 1, function(x) all(x == 1)))
                      return(indx)
                    } )
all.1.indx <- unique(unlist(all1.list))
array.data <- array.data[-all.1.indx, ]
seq.data <- seq.data[-all.1.indx, ]

#### get datatables to mix -----------------------------------------------------

# get a list that contains an
# array data.table and seq data.table for each level of 'titration'
array.train <-
  data.table(array.data[,
                        c(1,
                          which(colnames(array.data) %in% train.sample.names))])
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
                                     platform = "array")

# parallel backend
cl <- parallel::makeCluster(detectCores() - 1)
doParallel::registerDoParallel(cl)

# 'mixed' both platform normalization
norm.titrate.list[2:10] <-
  foreach(n = 2:10) %dopar% {
    NormalizationWrapper(titrate.mix.dt.list[[n]]$array,
                         titrate.mix.dt.list[[n]]$seq)
  }

# stop parallel backend
parallel::stopCluster(cl)
# sort out names
names(norm.titrate.list)[2:10] <- names(titrate.mix.dt.list)[2:10]

# single platform seq normalization
norm.titrate.list[["100"]] <-
  SinglePlatformNormalizationWrapper(titrate.mix.dt.list[[11]]$seq,
                                     platform = "seq")

# save train data
saveRDS(norm.titrate.list, file = file.path(norm.data.dir, norm.train.object))

#### normalize test data -------------------------------------------------------
array.test <-
  data.table(array.data[,
                        c(1,
                          which(colnames(array.data) %in% test.sample.names))])
seq.test <-
  data.table(seq.data[, c(1, which(colnames(seq.data) %in% test.sample.names))])

# array normalization
array.test.norm.list <-
  SinglePlatformNormalizationWrapper(array.test, platform = "array")

# seq normalization
# initialize list to hold normalized seq data
seq.test.norm.list <- list()

# LOG normalization
seq.test.norm.list[["log"]] <- LOGSeqOnly(seq.test)
# NPN
seq.test.norm.list[["npn"]] <- NPNSingleDT(seq.test)

# start parallel backend
cl <- parallel::makeCluster(detectCores() - 1)
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

# start parallel back end
cl <- parallel::makeCluster(detectCores() - 1)
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

# combine array and seq test data into a list
test.norm.list <- list(array = array.test.norm.list,
                       seq = seq.test.norm.list)

# save test data
saveRDS(test.norm.list, file = file.path(norm.data.dir, norm.test.object))
