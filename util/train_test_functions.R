GetCM <- function(model, dt.mat, category,
                  model.type = NULL, return.only.kappa = TRUE){
  # This function takes a model, a normalized gene expression matrix,
  # and performs prediction and returns the Kappa statistic based on the
  # observed category labels supplied
  #
  # Args:
  #   model: a predictive model of the class train (output of train() in the
  #   caret package)
  #   dt.mat: a matrix where genes are columns and rows are samples
  #   category: a vector of the true category labels
  #   model.type: what kind of predictive model will be used?
  #   return.only.kappa: logical, should the Kappa statistic associated with
  #                      the prediction be returned (TRUE) or should the entire
  #                      confusionMatrix be returned (FALSE)?
  #
  # Returns:
  #   if return.only.kappa = TRUE
  #     kap: the Kappa statistic associated with the prediction
  #   if return.only.kappa != TRUE
  #     cm: the confusionMatrix (caret) associated with the prediction
  #
  if (model.type == "glmnet") {
    prd <- predict(model, dt.mat, s=model$lambda.1se, type="class")
    cm <- confusionMatrix(as.factor(as.vector(prd)), category)
  } else {
    prd <- predict(model, dt.mat)
    tbl <- table(prd, category)
    cm <- confusionMatrix(tbl)
  }

  # should we return only the Kappa statistic (return.only.kappa = TRUE), or the
  # full confusionMatrix for the prediction?
  if (return.only.kappa) {
    kap <- as.numeric(cm$overall["Kappa"])
    return(kap)
  } else {
    return(cm)
  }
}

PredictCM <- function(model, dt, sample.df,
                      model.type = NULL, return.only.kappa = TRUE){
  # This function takes a model, a normalized gene expression matrix,
  # and performs prediction and returns the total accuracy based on the
  # observed category labels supplied
  #
  # Args:
  #   model: a predictive model of the class train OR class cv.glmnet
  #   dt: a data table where the columns are the samples and the rows are genes
  #   sample.df: the data frame that maps sample name/header to category and
  #              train/test set labels
  #              output of 0-expression_data_overlap_and_split.R
  #   model.type: is the model from glmnet (class: cv.glmnet) or from caret
  #               class
  #   return.only.kappa: logical, should the Kappa statistic associated with
  #                      the prediction be returned (TRUE) or should the entire
  #                      confusionMatrix be returned (FALSE)?
  #
  # Returns:
  #   if return.only.kappa = TRUE
  #     kap: the Kappa statistic associated with the prediction
  #   if return.only.kappa != TRUE
  #     cm: the confusionMatrix (caret) associated with the prediction
  #
  require(caret)
  require(glmnet)
  require(ranger)
  require(kernlab)
  require(data.table)

  category <- GetOrderedCategoryLabels(dt, sample.df)
  #dt.mat <- t(dt[, 2:ncol(dt), with = F])
  dt.mat <- t(dt[, -1, with = F])
  colnames(dt.mat) <- paste0("c", seq(1:ncol(dt.mat)))

  # check if dt.mat is character -- if so, make numeric
  if (any(apply(dt.mat, 1, is.character))) {
    dt.mat <- apply(dt.mat, 2, as.numeric)
  }

  pred.cm <- GetCM(model, dt.mat, category, model.type, return.only.kappa)
  return(pred.cm)
}

GetOrderedCategoryLabels <- function(exp.dt, sample.df){
  # This function takes a data.table of gene expression and a data.frame that
  # maps sample names to category labels
  #
  # Args:
  #   exp.dt: a data table where the columns are the samples and the rows
  #           are genes
  #   sample.df: the data frame that maps sample name/header to category and
  #              train/test set labels
  #              output of 0-expression_data_overlap_and_split.R
  #
  # Returns:
  #   category.labels: corresponding category labels ordered to match the columns/
  #                   sample name order
  #
  exp.samples <- colnames(exp.dt)[2:ncol(exp.dt)]
  indx <- vector()
  for (i in 1:length(exp.samples)) {
    indx[i] <- which(sample.df$sample == exp.samples[i])
  }
  category.labels <- sample.df$category[indx]
  return(category.labels)
}

TrainThreeModels <- function(dt, category, seed, folds.list){
  # This function takes a data.table of gene expression and sample
  # category labels in the same order as the columns. Seed and a list of folds
  # from createFolds are included for reproducibility purposes. Returns a list
  # of 3 predictive models - LASSO (glmnet), random forest (caret/ranger),
  # and linear SVM (caret/kernlab)
  #
  # Args:
  #   dt: a data.table where the columns are the samples and the rows
  #       are genes
  #   category: sample category labels in the same order as the columns
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
    for (i in 1:5) seed.list[[i]]<- sample.int(n=1000, 6) # needs 6 even though tuneLength is 3?
    seed.list[[6]] <- sample.int(n=1000, 1) # last model

    fit.control <- trainControl(method = "cv",
                                number = 5, # 5-fold cross-validation
                                #savePredictions = TRUE,
                                classProbs = TRUE,
                                seeds = seed.list,
                                index = folds.list, # list of 5 sets of indices
                                allowParallel = T) # use parallel processing

    fold.vector <- vector()  # foldids for cv.glmnet so the folds are the same
    # as caret models
    for (i in 1:length(category)) {
      fold.vector[i] <- which(unlist(lapply(folds.list, function(x) i %in% x)))
    }

    # initialize list that will hold the predictive models
    train.list <- list()

    t_dt <- t(dt[,-1, with = F])
    colnames(t_dt) <- paste0("c", seq(1:ncol(t_dt)))

    # LASSO
    train.list[["glmnet"]] <- cv.glmnet(t_dt,
                                        category,
                                        family = "multinomial",
                                        foldid = fold.vector, # fold 'labels'
                                        parallel = T,
                                        type.measure = "class")
    # Random Forest
    train.list[["rf"]] <- train(t_dt,
                                category,
                                method = "ranger",
                                trControl = fit.control,
                                tuneLength = 3)
    # Linear SVM
    train.list[["svm"]] <- train(t_dt,
                                 category,
                                 method = "svmLinear",
                                 trControl = fit.control,
                                 tuneLength = 3)

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
                           only.kap = TRUE, run.parallel = TRUE) {
  # This function is a wrapper for performing category prediction on
  # normalized expression data.tables (training data, holdout data, or
  # reconstructed data) using the models trained on training data
  # in the supervised analysis (train.list from 2-train_test_category.R).
  #
  # Args:
  #   train.model.list: a list of predictive models
  #                     (LASSO, linear SVM, random forest)
  #   pred.list: list of normalized expression data.tables, labels for this
  #              data are going to be predicted
  #   sample.df: the data frame that maps sample name/header to category and
  #              train/test set labels
  #              output of 0-expression_data_overlap_and_split.R
  #   only.kap: logical; (TRUE) should only the Kappa statistics associated with
  #             the predictions be returned? Or a list of all confusionMatrix
  #             objects and the associated Kappa statistics (FALSE)
  #
  #   run.parallel: logical; should predictions be run in parallel?
  #
  # Returns:
  #   if only.kap = TRUE
  #     kappa.df: a data.frame of Kappa statistics
  #   else
  #     a list of all confusionMatrix objects ($confusion_matrix_objects) and
  #     and the associated Kappa statistics ($kappa_statistics)


  # since parallelizing -- other requirements are captured in PredictCM
  require(foreach)

  # function for parallel prediction of seq levels
  ParallelPredictFunction <- function(model.list,
                                      pred.data.list,
                                      sample.dataframe,
                                      mdl.type,
                                      only.kappa = TRUE) {

    return.list <-
      foreach(seq.iter = seq_along(model.list)) %dopar% {
        seq.lvl <- names(model.list)[seq.iter]
        # get the model for % seq level
        trained.model <- model.list[[seq.lvl]]
        # if pred.data.list is a list (has different sequencing levels),
        # loop through the list -- this will be the case for RNA-seq hold out
        # data for TDM, QN, and QN-Z methods, as well as training data and
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
          
          if (only.kappa) {
            kappa <- PredictCM(model = trained.model,
                               dt = pred.dt,
                               sample.df = sample.dataframe,
                               model.type = mdl.type,
                               return.only.kappa = only.kappa)
            
            auc <- PredictAUC(model = trained.model,
                              dt = pred.dt,
                              sample.df = sample.dataframe,
                              model.type = mdl.type)  
            
            return(data.frame("kappa" = kappa,
                              "auc" = auc))
            
          } else {
            
            CM <- PredictCM(model = trained.model,
                            dt = pred.dt,
                            sample.df = sample.dataframe,
                            model.type = mdl.type,
                            return.only.kappa = only.kappa)
            
            return(CM)
            
          }
          
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
                              "GetOrderedCategoryLabels",
                              "PredictAUC",
                              "mean_one_versus_all_AUC"))
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
                                only.kappa = only.kap)
        
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

  # if only returning Kappa -- melt the list into a data.frame and return
  # the data.frame
  if (only.kap) {
    kappa.df <- norm.list %>%
      # when there is null test data at a particular %RNA-seq, discard that null
      purrr::modify_depth(2, function(x) discard(x, is.null)) %>%
      reshape2::melt(id.vars = NULL) %>%
      tidyr::pivot_wider(names_from = "variable",
                         values_from = "value")
    colnames(kappa.df) <- c("perc.seq", "classifier", "norm.method", "kappa", "auc")
    return(kappa.df)
  } else {  # otherwise, return two objects:
    # 1. the list of confusionMatrix objects (norm.list)
    # 2. kappa data frame (kappa.df)

    # create kappa.df from norm.list
    # level 3 of norm.list is confusion matrix associated with %seq, classifier, and normalization method
    # purrr::modify_depth applies a function to confusion matrix to return kappa
    # list can then be melted and columns renamed
    kappa.df <- norm.list %>%
      # when there is null test data at a particular %RNA-seq, discard that null
      purrr::modify_depth(2, function(x) discard(x, is.null)) %>%
      purrr::modify_depth(3, function(x) x$overall["Kappa"]) %>%
      reshape2::melt() %>%
      dplyr::rename("kappa" = "value",
                    "perc.seq" = "L3",
                    "classifier" = "L2",
                    "norm.method" = "L1")

    # create list and return
    return(list("confusion_matrix_objects" = norm.list,
                "kappa_statistics" = kappa.df))
  }

}

mean_one_versus_all_AUC <- function(probabilities_matrix,
                                    true_subtypes) {
  
  true_subtypes <- factor(true_subtypes)
  
  if (sort(levels(true_subtypes)) != sort(colnames(probabilities_matrix))) {
    
    # if true subtypes do not match the categories with probabilities
    stop("True subtype label levels do not match predicted subtype levels in mean_one_versus_all_AUC().")
    
  } else {
    
    # calculate the AUC value for each subtype vs. all others
    auc_vector <- purrr::map_dbl(.x = levels(true_subtypes),
                                 .f = function(subtype) MLmetrics::AUC(y_pred = probabilities_matrix[,subtype],
                                                                       y_true = as.numeric(true_subtypes == subtype)))
    
  }
  
  # take the unweighted mean of the AUC values
  # (this matches the caret::multiClassSummary value of AUC)
  return(mean(auc_vector))
  
}

PredictAUC <- function(model, dt, sample.df,
                       model.type = NULL) {
  
  category <- GetOrderedCategoryLabels(dt, sample.df)
  dt.mat <- t(dt[, -1, with = F])
  colnames(dt.mat) <- paste0("c", seq(1:ncol(dt.mat)))
  
  # check if dt.mat is character -- if so, make numeric
  if (any(apply(dt.mat, 1, is.character))) {
    dt.mat <- apply(dt.mat, 2, as.numeric)
  }
  
  if (model.type == "glmnet") {
    
    predicted_probabilities <- predict(model,
                                       dt.mat,
                                       s = model$lambda.1se,
                                       type = "response")[,,1]
    
    mean_auc <- mean_one_versus_all_AUC(probabilities_matrix = predicted_probabilities,
                                        true_subtypes = category)
    
  } else { # if model.type == "rf" or "svm"
    
    predicted_probabilities <- predict(model,
                                       dt.mat,
                                       type = "prob")
    
    mean_auc <- mean_one_versus_all_AUC(probabilities_matrix = predicted_probabilities,
                                        true_subtypes = category)
  }
  
  return(mean_auc)
  
}
