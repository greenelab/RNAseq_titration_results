# J. Taroni Jul 2016
# The purpose of this script is to train LASSO, linear SVM, and
# predictive models on normalized and mixed array and RNA-seq data
# (output of 1-normalized_titrated_data.R) and then to perform predictions on
# normalized test data.
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
  optparse::make_option("--seed3",
                        default = NA_integer_,
                        help = "Random seed")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
source(here::here("util/option_functions.R"))
check_options(opt)

# load libraries
suppressMessages(source(here::here("load_packages.R")))
source(here::here("util", "train_test_functions.R"))

# set options
cancer_type <- opt$cancer_type
predictor <- opt$predictor
file_identifier <- str_c(cancer_type, predictor, sep = "_")

# set seed
filename.seed <- opt$seed1
initial.seed <- opt$seed3
set.seed(initial.seed)

# define directories
norm.data.dir <- here::here("normalized_data")
mdl.dir <- here::here("models")
res.dir <- here::here("results")

# define input files
norm.test.object <- paste0(file_identifier,
                           "_array_seq_test_data_normalized_list_",
                           filename.seed, ".RDS")
norm.train.object <- paste0(file_identifier,
                            "_array_seq_train_titrate_normalized_list_",
                            filename.seed, ".RDS")
train.test.labels <- paste0(file_identifier,
                            "_matchedSamples_training_testing_split_labels_",
                            filename.seed, ".tsv")

# define output files
trained.models.object <- paste0(file_identifier,
                                "_train_3_models_",
                                filename.seed, ".RDS")
train.kappa.file <- file.path(res.dir,
                              paste0(file_identifier,
                                     "_train_3_models_training_set_total_kappa_",
                                     filename.seed, ".tsv"))
array.kappa.file <- file.path(res.dir,
                              paste0(file_identifier,
                                     "_train_3_models_array_kappa_",
                                     filename.seed, ".tsv"))
seq.kappa.file <- file.path(res.dir,
                            paste0(file_identifier,
                                   "_train_3_models_seq_kappa_",
                                   filename.seed, ".tsv"))

#### load data -----------------------------------------------------------------

sample.train.test <- fread(file.path(res.dir, train.test.labels), data.table = FALSE)
norm.titrate.list <- readRDS(file.path(norm.data.dir, norm.train.object))
norm.test.list <- readRDS(file.path(norm.data.dir, norm.test.object))

# set each category as a factor
sample.train.test$category <- as.factor(sample.train.test$category)

# category levels for each perc of seq data
category.norm.list <- lapply(norm.titrate.list,
                             function(x) GetOrderedCategoryLabels(x$z,
                                                                  sample.train.test))

# restructure normalized list so that it's organized by normalization method
restr.train.list <- RestructureNormList(norm.titrate.list)
rm(norm.titrate.list)

#### training ------------------------------------------------------------------

folds.seed <- sample(1:10000, 1)
message(paste("Random seed for createFolds:", folds.seed), appendLF = TRUE)
set.seed(folds.seed)
folds.list <- lapply(category.norm.list, function(x) createFolds(x, k = 5))

# parallel backend
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

resample.seed <- sample(1:10000, 1)
message(paste("Random seed for resampling:", resample.seed), appendLF=TRUE)

train.model.list <-
  foreach(n = 1:length(restr.train.list)) %do% {  # foreach norm method
    foreach(m = 1:length(category.norm.list)) %dopar% {  # foreach % seq level
      TrainThreeModels(dt = restr.train.list[[n]][[m]],
                       category = category.norm.list[[m]],
                       seed = resample.seed,
                       folds.list = folds.list[[m]])

    }
  }

# stop parallel backend
stopCluster(cl)

# get names
names(train.model.list) <- names(restr.train.list)
train.model.list <- mapply(
  function(x, y){
    names(x) <- names(y)
    return(x)
  }, x = train.model.list,
  y = restr.train.list,
  SIMPLIFY = TRUE)

# restructure trained model list so from top to bottom: norm method -> model
# type -> % seq level (0 - 100)
train.model.list <- RestructureTrainedList(train.model.list)

# save predictive models
saveRDS(train.model.list, file = file.path(mdl.dir, trained.models.object))

#### training kappa ---------------------------------------------------------
# get rid of 0, 100 tdm list, they're NULL
restr.train.list$tdm$`0` <- NULL
restr.train.list$tdm$`100` <- NULL

# get training kappa stats and write to file
train.kappa.df <- PredictWrapper(train.model.list = train.model.list,
                                 pred.list = restr.train.list,
                                 sample.df = sample.train.test,
                                 only.kap = TRUE)

write.table(train.kappa.df, file = train.kappa.file, sep = "\t",
            row.names = FALSE, quote = FALSE)

#### predictions - test data ---------------------------------------------------

# get predictions on array test data as a data frame
array.kappa.df <- PredictWrapper(train.model.list = train.model.list,
                                 pred.list = norm.test.list$array,
                                 sample.df = sample.train.test,
                                 only.kap = TRUE)

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
                               only.kap = TRUE)

write.table(seq.kappa.df, file = seq.kappa.file, sep = "\t",
            row.names = FALSE, quote = FALSE)
