# J. Taroni Jul 2016
# The purpose of this script is to train LASSO, linear SVM, and 
# predictive models on normalized and mixed array and RNA-seq data 
# (output of 1-normalized_titrated_data.R) and then to perform predictions on
# normalized test data. 
# It should be run from the command line through the run_experiments.R script

suppressMessages(source("load_packages.R"))
source(file.path("util", "train_test_functions.R"))

args <- commandArgs(trailingOnly = TRUE)
filename.seed <- as.integer(args[1])
initial.seed <- as.integer(args[2])
set.seed(initial.seed)

norm.data.dir <- "normalized_data"
mdl.dir <- "models"
res.dir <- "results"
norm.test.object <-
   paste0("BRCA_array_seq_test_data_normalized_list_", filename.seed, ".RDS")
norm.train.object <- 
  paste0("BRCA_array_seq_train_titrate_normalized_list_", filename.seed, ".RDS")

trained.models.object <- 
  paste0("BRCA_train_3_models_", filename.seed, ".RDS")
train.kappa.file <- 
  file.path(res.dir, 
            paste0("BRCA_train_3_models_training_set_total_kappa_", 
                   filename.seed, ".tsv"))
array.kappa.file <-
  file.path(res.dir, paste0("BRCA_train_3_models_array_kappa_", filename.seed, 
                            ".tsv"))
seq.kappa.file <- file.path(res.dir, paste0("BRCA_train_3_models_seq_kappa_", 
                                            filename.seed, ".tsv"))

train.test.labels <- 
  file.path(res.dir, 
            paste0("BRCA_matchedSamples_PAM50Array_training_testing_split_labels_", 
                  filename.seed, ".tsv"))

#### load data -----------------------------------------------------------------

sample.train.test <- fread(train.test.labels, data.table = FALSE)
norm.titrate.list <- readRDS(file.path(norm.data.dir, norm.train.object))
norm.test.list <- readRDS(file.path(norm.data.dir, norm.test.object))

# subtype levels for each perc of seq data
subtype.norm.list <- 
  lapply(norm.titrate.list, 
         function(x) GetOrderedSubtypeLabels(x$z, sample.train.test))

# restructure normalized list so that it's organized by normalization method
restr.train.list <- RestructureNormList(norm.titrate.list)
rm(norm.titrate.list)

#### training ------------------------------------------------------------------

folds.seed <- sample(1:10000, 1)
message(paste("Random seed for createFolds:", folds.seed), appendLF = TRUE)
set.seed(folds.seed)
folds.list <- lapply(subtype.norm.list, function(x) createFolds(x, k = 5))

# parallel backend
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

resample.seed <- sample(1:10000, 1)
message(paste("Random seed for resampling:", resample.seed), appendLF=TRUE)

train.model.list <- 
  foreach(n = 1:length(restr.train.list)) %do% {  # foreach norm method
#    foreach(m = 1:length(subtype.norm.list)) %dopar% {  # foreach % seq level
    foreach(m = 1:length(subtype.norm.list)) %do% {  # foreach % seq level
      print(c(n, m))
      TrainThreeModels(dt = restr.train.list[[n]][[m]], 
                       subtype =  subtype.norm.list[[m]], 
                       seed = resample.seed,
                       folds.list = folds.list[[m]])
      
    }
  }

# stop parallel backend
stopCluster(cl)

print(str(train.model.list))
print("gets here 1")

# get names 
names(train.model.list) <- names(restr.train.list)
train.model.list <- mapply(function(x, y){
                              names(x) <- names(y)
                              return(x)
                            }, x = train.model.list,
                            y = restr.train.list,
                            SIMPLIFY = TRUE)

print("gets here 2")

# restructure trained model list so from top to bottom: norm method -> model
# type -> % seq level (0 - 100)
train.model.list <- RestructureTrainedList(train.model.list)

# save predictive models
saveRDS(train.model.list, file = file.path(mdl.dir, trained.models.object))

print("gets here 3")

#### training kappa ---------------------------------------------------------
# get rid of 0, 100 tdm list, they're NULL
restr.train.list$tdm$`0` <- NULL
restr.train.list$tdm$`100` <- NULL

# get training kappa stats and write to file
train.kappa.df <- PredictWrapper(train.model.list = train.model.list,
                                 pred.list = restr.train.list,
                                 sample.df = sample.train.test,
                                 return.kap = TRUE)

print(paste("writing to ", train.kappa.file))

write.table(train.kappa.df, file = train.kappa.file, sep = "\t", 
            row.names = FALSE, quote = FALSE)
#### predictions - test data ---------------------------------------------------

# get predictions on array test data as a data frame
array.kappa.df <- PredictWrapper(train.model.list = train.model.list,
                                 pred.list = norm.test.list$array,
                                 sample.df = sample.train.test,
                                 return.kap = TRUE)

print(paste("writing to ", array.kappa.file))

write.table(array.kappa.df, file = array.kappa.file, sep = "\t", 
            row.names = FALSE, quote = FALSE)

# for the 0 perc seq level of the titration, the model tested on log transformed
# array data (100% array data) should be tested on the TDM transformed seq data
for(i in 1:length(train.model.list[[5]])){
  train.model.list[[5]][[i]]$`0` <- train.model.list[["log"]][[i]]$`0`
  train.model.list[[5]][[i]] <- train.model.list[[5]][[i]][c(10, 1:9)]
}

# get rid of 100 tdm list, it's NULL
norm.test.list$seq$tdm$`100` <- NULL

# get predictions on RNA-seq test data as a data frame
seq.kappa.df <- PredictWrapper(train.model.list = train.model.list,
                               pred.list = norm.test.list$seq,
                               sample.df = sample.train.test,
                               return.kap = TRUE)

print(paste("writing to ", seq.kappa.file))

write.table(seq.kappa.df, file = seq.kappa.file, sep = "\t", row.names = FALSE, 
            quote = FALSE)
