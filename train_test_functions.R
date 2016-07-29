GetKappa <- function(model, dt.mat, subtype, model.type = NULL){
  # This function takes a model, a normalized gene expression matrix,
  # and performs prediction and returns the Kappa statistic based on the 
  # observed subtype labels supplied
  #
  # Args:
  #   model: a predictive model of the class train (output of train() in the
  #   caret package)
  #   dt.mat: a matrix where genes are columns and rows are samples
  #   subtype: a vector of the true subtype labels
  #   model.type: what kind of predictive model will be used?
  #
  # Returns:
  #   Kappa: the Kappa statistic associated with the prediction
  #
  if (model.type == "glmnet") {
    prd <- predict(model, dt.mat, s=model$lambda.1se, type="class")
    cm <- confusionMatrix(as.factor(as.vector(prd)), subtype)
  } else {
    prd <- predict(model, dt.mat)
    tbl <- table(prd, subtype)
    cm <- confusionMatrix(tbl)
  }
  kap <- as.numeric(cm$overall["Kappa"])
  return(kap)
}

PredictKappa <- function(model, dt, sample.df, model.type=NULL){
  # This function takes a model, a normalized gene expression matrix,
  # and performs prediction and returns the total accuracy based on the 
  # observed subtype labels supplied
  #
  # Args:
  #   model: a predictive model of the class train OR class cv.glmnet
  #   dt: a data table where the columns are the samples and the rows are genes
  #   sample.df: the data frame that maps sample name/header to subtype and
  #              train/test set labels 
  #              output of 0-expression_data_overlap_and_split.R
  #   model.type: is the model from glmnet (class: cv.glmnet) or from caret 
  #               class 
  #
  # Returns:
  #   Kappa: the Kappa statistic of the prediction
  #
  require(caret)
  require(glmnet)
  require(ranger)
  require(kernlab)
  require(data.table)
   
  subtype <- GetOrderedSubtypeLabels(dt, sample.df)
  dt.mat <- t(dt[, 2:ncol(dt), with = F])
  pred.kappa <- GetKappa(model, dt.mat, subtype, model.type)
  return(pred.kappa)
}


GetOrderedSubtypeLabels <- function(exp.dt, sample.df){
  # This function takes a data.table of gene expression and a data.frame that
  # maps sample names to subtype labels 
  #
  # Args:
  #   exp.dt: a data table where the columns are the samples and the rows 
  #           are genes
  #   sample.df: the data frame that maps sample name/header to subtype and
  #              train/test set labels 
  #              output of 0-expression_data_overlap_and_split.R
  #
  # Returns:
  #   subtype.labels: corresponding subtype labels ordered to match the columns/
  #                   sample name order
  #
  exp.samples <- colnames(exp.dt)[2:ncol(exp.dt)]
  indx <- vector()
  for (i in 1:length(exp.samples)) {
    indx[i] <- which(sample.df$sample == exp.samples[i])
  }
  subtype.labels <- sample.df$subtype[indx]
  return(subtype.labels)
}

TrainThreeModels <- function(dt, subtype, seed, folds.list){
  # This function takes a data.table of gene expression and sample 
  # subtype labels in the same order as the columns. Seed and a list of folds
  # from createFolds are included for reproducibility purposes. Returns a list 
  # of 3 predictive models - LASSO (glmnet), random forest (caret/ranger),
  # and linear SVM (caret/kernlab)
  #
  # Args:
  #   dt: a data.table where the columns are the samples and the rows 
  #       are genes
  #   subtype: sample subtype labels in the same order as the columns
  #   seed: integer to be used to set seed and general seed list to supply to
  #         caret::train for parallel processing
  #   folds.list: a list with indices for each fold, from caret::createFolds
  #
  # Returns:
  #   train.list: a list of three
  #
  require(caret)
  require(glmnet)
  require(ranger)
  require(kernlab)
  require(data.table)
  
  if (!(is.null(dt))) {
    
    
    set.seed(seed)
    # need a seed list for parallel processing purposes?
    seed.list <- vector(mode = "list", length = 6)
    for (i in 1:5) seed.list[[i]]<- sample.int(n=1000, 3) # tuneLength is 3
    seed.list[[6]] <- sample.int(n=1000, 1) # last model
    
    fit.control <- trainControl(method = "cv",
                                number = 5, # 5-fold cross-validation
                                seeds = seed.list,
                                index = folds.list, # list of 5 sets of indices
                                allowParallel = T) # use parallel processing
    
    fold.vector <- vector()  # foldids for cv.glmnet so the folds are the same
    # as caret models
    for (i in 1:length(subtype)) {
      fold.vector[i] <- which(unlist(lapply(folds.list, function(x) i %in% x)))
    }
    
    # initialize list that will hold the predictive models
    train.list <- list()
    
    # LASSO
    train.list[["glmnet"]] <- cv.glmnet(t(dt[, 2:ncol(dt), with = F]), 
                                        subtype,
                                        family = "multinomial",
                                        foldid = fold.vector, # fold 'labels'
                                        parallel = T,
                                        type.measure="class") 
    # Random Forest    
    train.list[["rf"]] <- train(t(dt[, 2:ncol(dt), with = F]), 
                                subtype,
                                method = "ranger", 
                                trControl = fit.control,
                                tuneLength = 3)
    # Linear SVM
    train.list[["svm"]] <- train(t(dt[, 2:ncol(dt), with = F]), 
                                 subtype,
                                 method = "svmLinear", 
                                 trControl = fit.control,
                                 tuneLength = 3)
    # stopCluster(cl)
    train.list[["seeds"]] <- seed.list
    return(train.list)
  } else {
    return(list(NULL))
  }
}

RestructureNormList <- function(norm.list){
  # This function takes the normalized training data, the output of 
  # 1-normalized_titrate_data.R (SinglePlattformNormalizationWrapper/
  # NormalizationWrapper) which is organized by the 'titration level' or percent
  # RNA-seq and restructures the list such that it is organized by normalization
  # method for use to train predictive models
  # 
  # Args:
  #   norm.list: a list of normalized and concatenated expression data.tables
  # 
  # Returns:
  #   restr.norm.list: a list of normalized and concatenated expression 
  #                    data.tables, now reorganized by normalization method
  # 

  # identify the different methods used
  unique.methods <- unique(unlist(lapply(norm.list, names)))
  unique.methods <- unique.methods[-(which(unique.methods == "raw.array"))]
  restr.norm.list <- list()
  for (unq in unique.methods) {
    restr.norm.list[[unq]] <- lapply(norm.list, 
                                     function(x) if (unq %in% names(x)){
                                       return(x[[unq]])
                                     })
  }
  return(restr.norm.list)
}


RestructureTrainedList <- function(train.list){
  # This function takes a train.list (from applying TrainThreeModels), 
  # which is organized from top to bottom: normalization method -> percent 
  # sequencing data -> model type and restructures the list so that the 
  # structure from top to bottom: normalization method -> model type -> % seq.
  # This restructuring aids in prediction on test data and subsequent plotting.
  # 
  # Args:
  #   train.list: a list of predictive models, organized as stated above
  # 
  # Returns:
  #   norm.list: a list of predictive models that is reorganized as stated above
  # 
  norm.list <- list()
  for (i in 1:ncol(train.list)) { # for each normalization method
    new.list <- list() 
    for (meth in c("glmnet", "rf", "svm", "seeds")) { # for each model type
      # and also include seeds
      meth.list <- list()
      for (j in 1:length(train.list[, i])) { # get the models of model type
        # meth in order of % seq (0 - 100)
        meth.list[[names(train.list[, i][j])]] <- 
          train.list[, i][[j]][[meth]]
      }
      new.list[[meth]] <- meth.list
    }
    norm.list[[colnames(train.list)[i]]] <- as.matrix(new.list)
  }
  
  return(norm.list)
  
}

PredictArrayDataWrapper <- function(norm.array.list, train.list, sample.df){
  # This function takes a list of normalized array data (test set), a 
  # restructured train.list (from RestructuredTrainedList), and the data.frame
  # that maps sample names to subtype labels and train/test set labels. It
  # returns the Kappa statistic from the the various models
  # 
  # Args:
  #   norm.array.list: a list of normalized array data (test set)
  #   train.list: a list of predictive models, restructured using 
  #               RestructuredTrainedList
  #   sample.df: the data frame that maps sample name/header to subtype and
  #              train/test set labels 
  #              output of 0-expression_data_overlap_and_split.R
  # 
  # Returns:
  #   pred.list: a list of Kappa statistics from predictions on the test 
  #              array data
  # 
  
  # parallel backend 
  cl <- makeCluster(2)
  registerDoParallel(cl)
  
  pred.list <- foreach(n = 1:length(train.list)) %do% { # for each norm method
    
    mthd <- names(train.list)[n]  # what is the name of the norm method?
    if (mthd == "tdm") {  # when TDM is used to normalize the % seq data 
                          # included, we should test on log transformed array 
                          # data test set
      foreach(m = 1:3) %do% {  # exclude the seeds element of list (#4)
        foreach(o = 1:length(train.list[[n]][[m]]),
                .export=c("PredictKappa",
                          "GetKappa",
                          "GetKappaGlmnet",
                          "GetOrderedSubtypeLabels")) %dopar% {
                            PredictKappa(model = train.list[[n]][[m]][[o]],
                                         dt = norm.array.list[["log"]],
                                         sample.df = sample.df,
                                       model.type = rownames(train.list[[n]])[m]
                            )
                          }
      }
    } else {  # all other norm methods can use the corresponding array data 
      foreach(m = 1:3) %do% { # exclude the seeds element of list (#4)
        foreach(o = 1:length(train.list[[n]][[m]]), # for each %seq level
                .export=c("PredictKappa",
                          "GetKappa",
                          "GetKappaGlmnet",
                          "GetOrderedSubtypeLabels")) %dopar% {
                            PredictKappa(model = train.list[[n]][[m]][[o]],
                                         dt = norm.array.list[[mthd]],
                                         sample.df = sample.df,
                                       model.type = rownames(train.list[[n]])[m]
                            )
                          }
      }
    }
  }
  
  # sort out names, foreach does not retain element names like the apply family
  for (n in 1:length(train.list)) { # norm methods
    names(pred.list)[n] <- names(train.list)[n]
    for (m in 1:3) { # model types
      names(pred.list[[n]])[m] <- rownames(train.list[[n]])[m]
      for (o in 1:length(train.list[[n]][[m]])) { # %seq levels
        names(pred.list[[n]][[m]])[o] <- names(train.list[[n]][[m]])[o]
      }
    }
  }
  
  # stop parallel back end
  stopCluster(cl)
  
  return(pred.list)
  
}

PredictSeqDataWrapper <- function(norm.seq.list, train.list, sample.df){
  # This function takes a list of normalized RNA-seq data (test set), a 
  # restructured train.list (from RestructuredTrainedList), and the data.frame
  # that maps sample names to subtype labels and train/test set labels. It
  # returns the Kappa statistic from the the various models
  # 
  # Args:
  #   norm.seq.list: a list of normalized RNA-seq data (test set)
  #   train.list: a list of predictive models, restructured using 
  #               RestructuredTrainedList
  #   sample.df: the data frame that maps sample name/header to subtype and
  #              train/test set labels 
  #              output of 0-expression_data_overlap_and_split.R
  # 
  # Returns:
  #   pred.list: a list of Kappa statistics from predictions on the test 
  #              seq data
  # 
  
  # parallel backend 
  cl <- makeCluster(2)
  registerDoParallel(cl)
  
  pred.list <- foreach(n = 1:length(train.list)) %do% {
    
    mthd <- names(train.list)[n]
    # for the methods that require arrays as reference to normalize the 
    # sequencing data (TDM and QN), the seq test data has been normalized using
    # whatever arrays were used in the training set 
    if (mthd == "tdm" | mthd == "qn") {
      foreach(m = 1:3) %do% {  # exclude the seeds element of list (#4)
        foreach(o = 1:length(train.list[[n]][[m]]),
                .export=c("PredictKappa",
                          "GetKappa",
                          "GetKappaGlmnet",
                          "GetOrderedSubtypeLabels")) %dopar% {
                            PredictKappa(model = train.list[[n]][[m]][[o]],
                                         dt = norm.seq.list[[mthd]][[o]],
                                         sample.df = sample.df,
                                       model.type = rownames(train.list[[n]])[m]
                            )
                          }
      }
    } else { # the gene level normalization methods only have one seq test 
             # data.table
      foreach(m = 1:3) %do% {  # exclude the seeds element of list (#4)
        foreach(o = 1:length(train.list[[n]][[m]]),
                .export=c("PredictKappa",
                          "GetKappa",
                          "GetKappaGlmnet",
                          "GetOrderedSubtypeLabels")) %dopar% {
                            PredictKappa(model = train.list[[n]][[m]][[o]],
                                         dt = norm.seq.list[[mthd]],
                                         sample.df = sample.df,
                                       model.type = rownames(train.list[[n]])[m]
                            )
                          }
      }
    }
  }
  
  # sort out names
  for (n in 1:length(train.list)) {  # norm methods
    names(pred.list)[n] <- names(train.list)[n]
    for (m in 1:3) {  # model type
      names(pred.list[[n]])[m] <- rownames(train.list[[n]])[m]
      for (o in 1:length(train.list[[n]][[m]])) {  # %seq level
        names(pred.list[[n]][[m]])[o] <- names(train.list[[n]][[m]])[o]
      }
    }
  }
  
  stopCluster(cl)
  
  return(pred.list)
  
}

GetTrainingSetKappa <- function(model.list, train.data.list, subtype.list){
  # This function takes a list of normalized training data (from 
  # RestructuredNormList), a restructured train.list (from 
  # RestructuredTrainedList), and the list of subtype labels used for training
  # 
  # Args:
  #   train.data.list: a list of normalized training data, restructured using
  #                    RestructureNormList
  #   model.list: a list of predictive models, restructured using 
  #               RestructuredTrainedList
  #   subtype.list: list of subtype labels that were used for training
  #   
  # Returns:
  #   train.kappa.mlt: a data.frame of Kappa statistics 
  #                   from predictions on the training data
  
  # Error-handling 
  equal.length <- (length(train.data.list) != length(model.list))
  if (equal.length) {
    stop("train.data.list and model.list must be of equal length")
  }
  
  train.kappa.list <- list() # initialize list to hold all Kappa results
  for (i in 1:length(train.data.list)) { # norm methods
    norm.list <- list() # initialize lost to hold all results from the
                        # current normalization method 
    for (j in 1:3) {  # exclude the seeds element of list (#4)
      classif.list <- list() # initialize list to hold all results from the 
                             # current model type
      model.type <- rownames(model.list[[i]])[j]
      for (k in 1:length(train.data.list[[i]])) { # for each seq level
        dt <- train.data.list[[i]][[k]]
        dt.mat <- t(dt[, 2:ncol(dt), with=F])
        perc.seq <- names(train.data.list[[i]])[k]
        if (model.type == "glmnet") {
          classif.list[[perc.seq]] <-
            GetKappaGlmnet(model = model.list[[i]][[j]][[k]],
                              dt.mat = dt.mat,
                              subtype = subtype.list[[perc.seq]])
        } else {
          classif.list[[perc.seq]] <-
            GetKappa(model = model.list[[i]][[j]][[k]],
                        dt.mat = dt.mat,
                        subtype = subtype.list[[perc.seq]])
        } 
      }  
      norm.list[[model.type]] <- classif.list
    }
    train.kappa.list[[names(train.data.list)[i]]] <- norm.list
  }
  
  train.kappa.mlt <- melt(train.kappa.list)
  colnames(train.kappa.mlt) <- c("kappa", "perc.seq", "classifier", 
                               "norm.method")
  return(train.kappa.mlt)

}
