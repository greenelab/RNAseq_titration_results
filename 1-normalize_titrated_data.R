# J. Taroni Jun 2016
# The purpose of this script is to read in TGCA array and sequencing data,
# already pre-processed to only include test tumor samples, 
# (output of 0-expression_data_overlap_and_split.R) and to normalize
# the data.  
# It should be run from the command line through the run_experiments.R script

suppressMessages(source("load_packages.R"))
source("normalization_functions.R")

args <- commandArgs(trailingOnly = TRUE)
filename.seed <- as.integer(args[1])
kInitialSeed <- as.integer(args[2])
set.seed(kInitialSeed)

data.dir <- "data/"
norm.data.dir <- "normalized_data/"
seq.file <- "BRCARNASeq_matchedOnly_ordered.pcl"
array.file <- "BRCAarray_matchedOnly_ordered.pcl"
norm.test.object <- 
  paste0("BRCA_array_seq_test_data_normalized_list_", filename.seed, ".RDS")
norm.train.object <- 
  paste0("BRCA_array_seq_train_titrate_normalized_list_", filename.seed, ".RDS")

train.test.labels <- 
  paste0("BRCA_matchedSamples_PAM50Array_training_testing_split_labels_", 
         filename.seed, ".tsv")

#### functions -----------------------------------------------------------------
seq.data <- fread(paste0(data.dir, seq.file), data.table = F)
array.data <- fread(paste0(data.dir, array.file), data.table = F)
sample.train.test <- read.delim(train.test.labels)
#### split samples, titrate ----------------------------------------------------
train.sample.names <- as.character(sample.train.test$sample[
  which(sample.train.test$split == "train")])
test.sample.names <- as.character(sample.train.test$sample[
  which(sample.train.test$split == "test")]) 
# get samples for 'titration' 
p.list <- list()
i <- 1
for (p in seq(0, 1, by = 0.1)) {
  p.list[i] <- p
  i <- i + 1
}

titration.seed <- sample(1:10000, 1)
message(paste("Random seed for titration:", 
              titration.seed), appendLF=TRUE)

set.seed(titration.seed)
titrate.sample.list <- lapply(p.list, 
                          function(x) GetTitratedSampleNames(train.sample.names,
                                                             x))
names(titrate.sample.list) <- as.character(seq(0, 100, by=10))
# these samples will be the RNA-seq samples in any given 'titration' experiment
# remove rows that are equal to all ones -- for any combination + test data
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
# array datatable and seq datatable for each level of 'titration'
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
norm.titrate.list[["0"]] <- SinglePlatformNormalizationWrapper(
  titrate.mix.dt.list[[1]]$array, platform = "array") 

# parallel backend
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)
# 'mixed' both platform normalization
norm.titrate.list[2:10] <- 
  foreach(n = 2:10) %dopar% {
    NormalizationWrapper(titrate.mix.dt.list[[n]]$array,
                       titrate.mix.dt.list[[n]]$seq)
  }
stopCluster(cl)
names(norm.titrate.list)[2:10] <- names(titrate.mix.dt.list)[2:10]
# single platform seq normalization
norm.titrate.list[["100"]] <- SinglePlatformNormalizationWrapper(
  titrate.mix.dt.list[[11]]$seq, platform = "seq")

# save train data
saveRDS(norm.titrate.list, file=paste0(norm.data.dir, norm.train.object))

#### normalize test data -------------------------------------------------------
array.test <- 
  data.table(array.data[, 
                        c(1, 
                          which(colnames(array.data) %in% test.sample.names))])
seq.test <- 
  data.table(seq.data[, c(1, which(colnames(seq.data) %in% test.sample.names))])

# array normalization
array.test.norm.list <- SinglePlatformNormalizationWrapper(
  array.test, platform = "array") 

# seq normalization
seq.test.norm.list <- list()

seq.test.norm.list[["log"]] <- LOGSeqOnly(seq.test)
seq.test.norm.list[["npn"]] <- NPNSingleDT(seq.test)

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)
seq.qn.list <- list()
seq.qn.list[["0"]] <- QNSingleWithRef(ref.dt = norm.titrate.list$`0`$log,
                                      targ.dt = seq.test)
seq.qn.list[2:10] <- foreach(i=2:10) %dopar% {
  QNSingleWithRef(ref.dt=norm.titrate.list[[i]]$raw.array,
                  targ.dt = seq.test)
}  
names(seq.qn.list)[2:10] <- names(norm.titrate.list)[2:10]
seq.qn.list[["100"]] <- QNSingleDT(seq.test)
seq.test.norm.list[["qn"]] <- seq.qn.list
rm(seq.qn.list)
stopCluster(cl)

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)
seq.tdm.list <- list()
seq.tdm.list[["0"]] <- TDMSingleWithRef(ref.dt = norm.titrate.list$`0`$log,
                                     targ.dt = seq.test)
seq.tdm.list[2:10] <- foreach(i=2:10) %dopar% {
  TDMSingleWithRef(ref.dt=norm.titrate.list[[i]]$raw.array,
                   targ.dt = seq.test)
}  
names(seq.tdm.list)[2:10] <- names(norm.titrate.list)[2:10]
seq.tdm.list["100"] <- list(NULL)
seq.test.norm.list[["tdm"]] <- seq.tdm.list
rm(seq.tdm.list)
stopCluster(cl)

seq.test.norm.list[["z"]] <- ZScoreSingleDT(seq.test)

test.norm.list <- list(array=array.test.norm.list,
                       seq=seq.test.norm.list)

saveRDS(test.norm.list, file=paste0(norm.data.dir, norm.test.object))
