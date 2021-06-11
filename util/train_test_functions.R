GetCM <- function(model, dt.mat, subtype, 
                  model.type = NULL, return.kappa = TRUE){
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
  #   return.kappa: logical, should the Kappa statistic associated with
  #                 the prediction be returned (TRUE) or should the entire 
  #                 confusionMatrix be returned (FALSE)?
  #
  # Returns:
  #   if return.kappa = TRUE
  #     kap: the Kappa statistic associated with the prediction
  #   if return.kappa != TRUE
  #     cm: the confusionMatrix (caret) associated with the prediction
  #
  if (model.type == "glmnet") {
    prd <- predict(model, dt.mat, s=model$lambda.1se, type="class")
    cm <- confusionMatrix(as.factor(as.vector(prd)), subtype)
  } else {
    prd <- predict(model, dt.mat)
    tbl <- table(prd, subtype)
    cm <- confusionMatrix(tbl)
  }
  
  # should we return the Kappa statistic (return.kappa = TRUE), or the
  # full confusionMatrix for the prediction?
  if (return.kappa) {
    kap <- as.numeric(cm$overall["Kappa"])
    return(kap)
  } else {
    return(cm)
  }
}

PredictCM <- function(model, dt, sample.df, 
                      model.type = NULL, return.kappa = TRUE){
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
  #   return.kappa: logical, should the Kappa statistic associated with
  #                 the prediction be returned (TRUE) or should the entire 
  #                 confusionMatrix be returned (FALSE)?
  #                 
  # Returns:
  #   if return.kappa = TRUE
  #     kap: the Kappa statistic associated with the prediction
  #   if return.kappa != TRUE
  #     cm: the confusionMatrix (caret) associated with the prediction
  #
  require(caret)
  require(glmnet)
  require(ranger)
  require(kernlab)
  require(data.table)
  
  subtype <- GetOrderedSubtypeLabels(dt, sample.df)
  dt.mat <- t(dt[, 2:ncol(dt), with = F])
  
  # check if dt.mat is character -- if so, make numeric
  if (any(apply(dt.mat, 1, is.character))) {
    dt.mat <- apply(dt.mat, 2, as.numeric)
  }
  
  pred.cm <- GetCM(model, dt.mat, subtype, model.type, return.kappa)
  return(pred.cm)
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
    for (i in 1:5) seed.list[[i]]<- sample.int(n=1000, 6) #3) # tuneLength is 3
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
    
    t_dt <- t(dt[,-1, with = F])
    colnames(t_dt) <- paste0("c", seq(1:ncol(t_dt)))

    # LASSO
#    train.list[["glmnet"]] <- cv.glmnet(t(dt[, 2:ncol(dt), with = F]), 
#                                        subtype,
#                                        family = "multinomial",
#                                        foldid = fold.vector, # fold 'labels'
#                                        parallel = T,
#                                        type.measure="class") 
    # Random Forest
    #train.list[["rf"]] <- train(t(dt[, 2:ncol(dt), with = F]),

    print("about to train.list[[rf]]")

    train.list[["rf"]] <- train(t_dt, 
                                subtype,
                                method = "ranger", 
                                trControl = fit.control,
                                tuneLength = 3)
    # Linear SVM
    #train.list[["svm"]] <- train(t(dt[, 2:ncol(dt), with = F]), 

    print("about to train.list[[svm]]")

    train.list[["svm"]] <- train(t_dt,
                                 subtype,
                                 method = "svmLinear", 
                                 trControl = fit.control,
                                 tuneLength = 3)

    # LASSO
    #train.list[["glmnet"]] <- cv.glmnet(t(dt[, 2:ncol(dt), with = F]),

    print("about to train.list[[glmnet]]")

    train.list[["glmnet"]] <- cv.glmnet(t_dt,
                                        subtype,
                                        family = "multinomial",
                                        foldid = fold.vector, # fold 'labels'
                                        parallel = T,
                                        type.measure="class")
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
  for (norm in 1:ncol(train.list)) { # for each normalization method
    new.list <- list() 
    for (meth in c("glmnet", "rf", "svm", "seeds")) { # for each model type
      # and also include seeds
      meth.list <- list()
      for (seq.indx in 1:length(train.list[, norm])) { # get the models of model 
        # type
        # meth in order of % seq (0 - 100)
        meth.list[[names(train.list[, norm][seq.indx])]] <- 
          train.list[, norm][[seq.indx]][[meth]]
      }
      new.list[[meth]] <- meth.list
    }
    norm.list[[colnames(train.list)[norm]]] <- as.matrix(new.list)
  }
  
  return(norm.list)
  
}

PredictWrapper <- function(train.model.list, pred.list, sample.df, 
                           return.kap = TRUE, run.parallel = TRUE) {
  # This function is a wrapper for performing subtype prediction on 
  # normalized expression data.tables (training data, holdout data, or
  # reconstructed data) using the models trained on training data
  # in the supervised analysis (train.list from 2-train_test_brca_subtype.R).
  # 
  # Args:
  #   train.model.list: a list of predictive models 
  #                     (LASSO, linear SVM, random forest) 
  #   pred.list: list of normalized expression data.tables, labels for this
  #              data are going to be predicted
  #   sample.df: the data frame that maps sample name/header to subtype and
  #              train/test set labels 
  #              output of 0-expression_data_overlap_and_split.R
  #   return.kap: logical; should the entire confusionMatrix (FALSE) or just
  #               the Kappa statistic associated with the prediction be 
  #               returned?
  #   run.parallel: logical; should predictions be run in parallel?
  # 
  # Returns:
  #   if return.kap = TRUE
  #     kappa.df: a data.frame of Kappa statistics
  #   else  
  #     norm.list: a list of confusionMatrix objects (if return.kap = FALSE) 

  
  # since parallelizing -- other requirements are captured in PredictCM
  require(foreach)
  
  # function for parallel prediction of seq levels
  ParallelPredictFunction <- function(model.list,
                                      pred.data.list,
                                      sample.dataframe,
                                      mdl.type,
                                      rtrn.kappa = TRUE) {
    
    return.list <-
      foreach(seq.iter = seq_along(model.list)) %dopar% {
        seq.lvl <- names(model.list)[seq.iter]
        # get the model for % seq level
        trained.model <- model.list[[seq.lvl]]
        # if pred.data.list is a list (has different sequencing levels),
        # loop through the list -- this will be the case for RNA-seq hold out 
        # data for TDM and QN methods, as well as training data and 
        # reconstructed data
        if (class(pred.data.list) == "list") {
          # for some cases of TDM normalized data, not every seq level is 
          # present
          if (!is.null(pred.data.list[[seq.lvl]])) {
            pred.dt <- pred.data.list[[seq.lvl]]
          } else {
            pred.dt <- NULL
          }
        } else {  # pred.data.list is a single data.table otherwise
          pred.dt <- pred.data.list
        }
        
        # if pred.dt is not null, then perform prediction
        if (!is.null(pred.dt)) {
          PredictCM(model = trained.model,
                    dt = pred.dt,
                    sample.df = sample.dataframe,
                    model.type = mdl.type,
                    return.kappa = rtrn.kappa)
          
        }
        
      }
    
    names(return.list) <- names(model.list)
    
    return(return.list)
    
  }
  
  if (run.parallel) {  # if run.parallel is false, %dopar% in
    # ParallelPredictionFunction will run sequentially without parallel
    # backend
    
    # start parallel backend
    cl <- parallel::makeCluster(2)
    doParallel::registerDoParallel(cl)
    # these functions are required for prediction
    parallel::clusterExport(cl,
                            c("PredictCM",
                              "GetCM",
                              "GetOrderedSubtypeLabels"))
  }
  
  norm.methods <- names(train.model.list) 
  model.names <- c("glmnet", "rf", "svm")
  norm.list <-
    # for each normalization method
    foreach(mthd.iter = seq_along(train.model.list)) %do% {
      mthd <- norm.methods[mthd.iter]  # use method name rather than index
      # if the method is tdm and there's no tdm data to use for prediction
      # use log data -- this will be the case for array hold out sets
      if (mthd == "tdm" & !("tdm" %in% names(pred.list))) {
        input.list <- pred.list[["log"]]
      } else {  # otherwise, just use the corresponding data for prediction
        input.list <- pred.list[[mthd]]
      }
      # for each model (glmnet, rf, svm)
      foreach(mdl.iter = seq_along(model.names)) %do% {
        mdl <- model.names[mdl.iter]  # use model name rather than model index
        # do parallel prediction
        ParallelPredictFunction(model.list = 
                                  train.model.list[[mthd]][mdl, ][[mdl]],
                                pred.data.list = input.list,
                                sample.dataframe = sample.df,
                                mdl.type = mdl,
                                rtrn.kappa = return.kap)
        
      }
      
    }
  
  if (run.parallel) {
    # stop parallel backend
    parallel::stopCluster(cl)
  }
  
  # sort out names
  names(norm.list) <- names(train.model.list)
  for (nrm.it in seq_along(norm.list)) {
    names(norm.list[[nrm.it]]) <- model.names
  }
  
  # if returning Kappa -- melt the list into a data.frame and return
  # the data.frame
  if (return.kap) {
    kappa.df <- reshape2::melt(norm.list)
    colnames(kappa.df) <- c("kappa", "perc.seq", "classifier", "norm.method")
    return(kappa.df)
  } else {  # otherwise, return the list of confusionMatrix objects
    return(norm.list)
  }
  
}
