
# to be shared between multiple methods
GetRMSE <- function(true.mat, recon.mat){
  # This function takes two gene expression matrices and compares them 
  # the reconstruction error between the two (RMSE) is calculated on a per
  # sample basis.
  # 
  # Args:
  #   true.mat: gene expression matrix before reconstruction, samples are rows
  #             genes are columns
  #   recon.mat: gene expression matrix after reconstruction, samples are rows
  #              genes are columns
  # 
  # Returns:
  #   rmse: root squared mean error between samples before and after (a vector)
  # 
  # Error-handing
  if (all.equal(dim(true.mat), dim(recon.mat)) != TRUE){
    stop("true.mat and recon.mat must be of equal dimensions")
  }
  diff.mat <- true.mat-recon.mat
  # calculate RMSE on a per sample basis
  rmse <- apply(diff.mat, 1, function(x) sqrt(mean(x^2)))
  return(rmse)
}

ExpressionDataTableToMatrixT <- function(dataTable){
  # This function transposes a data.table such that then samples are rows,
  # columns are genes
  # 
  # Args:
  #   dataTable: data.table where the first column contains gene identifiers,
  #              columns are samples
  #              
  # Returns:
  #   data.mat: transposed data matrix
  # 
  require(data.table)
  data.mat <- data.frame(dataTable[, 2:ncol(dataTable), with=F])
  data.mat <- apply(data.mat, 2, as.numeric)
  return(t(data.mat))  # transpose
}

EvaluatePCARecon <- function(train.dt, test.dt, no.comp = 50){
  # This function takes training and test set data.tables
  # performs PCA on training, projects test data into training PC space
  # and back out - quantifies 'reconstruction error' using RMSE
  # 
  # Args:
  #   train.dt: a data.table of training data
  #   test.dt: a data.table of test data
  #   no.comp: number of principal components to use for reconstruction
  # 
  # Returns:
  #   rmse: root mean squared error of reconstruction
  #   
  require(data.table)
  
  train.is.dt <- "data.table" %in% class(train.dt)
  test.is.dt <- "data.table" %in% class(test.dt)
  any.not.dt <- !(any(c(train.is.dt, test.is.dt)))
  if (any.not.dt) {
    stop("train.dt and test.dt must both be data.tables")
  }
  if (!(all(train.dt[[1]] %in% test.dt[[1]]))) {
    stop("Gene identifiers in data.tables must match")
  }
  
  PrcompDataTable <- function(dataTable){
    mat.t <- ExpressionDataTableToMatrixT(dataTable)
    pca <- prcomp(mat.t)
    return(pca)
  }
  
  PCNewDataTable <- function(pc, dataTable, no.comp){
    # Error-handling 
    if (class(pc) != "prcomp"){
      stop("pc must be prcomp object")
    }
    # get new data (dataTable) into correct form (matrix, transposed)
    mat.t <- ExpressionDataTableToMatrixT(dataTable)
    proj <- scale(mat.t, pc$center, pc$scale) %*% pc$rotation[, 1:no.comp] 
    recon <- proj %*% t(pc$rotation[, 1:no.comp])
    return(recon)
  }
  
  # PCA on training data
  pc <- PrcompDataTable(train.dt)
  # project test data into training PC space & back out
  recon <- PCNewDataTable(pc, test.dt, no.comp)
  # transpose mat for test data, to be used as "true.mat"
  test.dm.t <- ExpressionDataTableToMatrixT(test.dt)
  
  # quantification: how well did we do at reconstruction?
  rmse <- GetRMSE(true.mat = test.dm.t, recon.mat = recon)
  return(rmse)
  
}

EvaluateICARecon <- function(train.dt, test.dt, no.comp = 50){
  # This function takes training and test set data.tables
  # performs PCA on training, projects test data into PC space, then IC space
  # and back out - quantifies 'reconstruction error' using RMSE
  # 
  # Args:
  #   train.dt: a data.table of training data
  #   test.dt: a data.table of test data
  #   no.comp: number of independent components to use
  # 
  # Returns:
  #   rmse: root mean squared error of reconstruction
  #   
  require(data.table)
  # Error-handling
  train.is.dt <- "data.table" %in% class(train.dt)
  test.is.dt <- "data.table" %in% class(test.dt)
  any.not.dt <- !(any(c(train.is.dt, test.is.dt)))
  if (any.not.dt) {
    stop("train.dt and test.dt must both be data.tables")
  }
  if (!(all(train.dt[[1]] %in% test.dt[[1]]))) {
    stop("Gene identifiers in data.tables must match")
  }
  
  FastICADataTable <- function(dataTable, n.comp){
    mat.t <- ExpressionDataTableToMatrixT(dataTable)
    ica <- fastICA::fastICA(mat.t, n.comp = no.comp, alg.typ = "parallel",
                            fun = "exp", method = "C")
    return(ica)
  }
  
  ICNewDataTable <- function(ic, dataTable){
    # Error-handling
    # Are the correct matrices in the list ic
    correct.matrices <- all(names(ic) == c("X", "K", "W", "A", "S"))
    if (!(correct.matrices)) {
      stop("ic must be output of fastICA::fastICA")
    }
    # get new data (dataTable) into correct form (matrix, transposed)
    mat.t <- ExpressionDataTableToMatrixT(dataTable)
    recon <- (scale(mat.t, center=TRUE, scale=FALSE) %*% ic$K %*% ic$W %*% ic$A)
    return(recon)
  }
  
  # ICA on training data
  ic <- FastICADataTable(train.dt, no.comp)
  # project test data into training IC space & back out
  recon <- ICNewDataTable(ic, test.dt)
  # transpose mat for test data, to be used as "true.mat"
  test.dm.t <- ExpressionDataTableToMatrixT(test.dt)
  
  # quantification: how well did we do at reconstruction?
  rmse <- GetRMSE(true.mat = test.dm.t, recon.mat = recon)
  return(rmse)
  
}

CompAnalysisEvalWrapper <- function(train.list, 
                                    test.list,
                                    no.comp = 50,
                                    platform = "array",
                                    method = "PCA"){
  # This is a wrapper function for ICA/PCA reconstruction error
  # 
  # Args:
  #   train.list: list of training data.tables - output of 
  #               1-normalize_titrated_data.R
  #   test.list: list of test data.tables - output of 
  #              1-normalize_titrated_data.R
  #   no.comp: number of independent or principal components for reconstruction
  #   platform: is the test data 'array' or 'seq' data? character
  #   method: perform 'ICA' or 'PCA'? character
  # 
  # Return:
  #   mlt.c: a data.frame of the rmse results, with the seq level, 
  #          normalization method, component analysis type included
  #          
  check.method <- any(c(method == "PCA", method == "ICA"))
  if (!check.method){
    stop("method must be 'PCA' or 'ICA'")
  }
  
  suppressMessages(require(reshape))
  cl <- parallel::makeCluster(detectCores() - 1)
  doParallel::registerDoParallel(cl)
  
  if (platform == "seq"){
    comp.list <- foreach(i = 1:length(train.list)) %do% {
      mthd <- names(train.list)[i]
      if (mthd == "qn" | mthd == "tdm"){
        foreach(j = 1:length(train.list[[i]]),
                .export = c("EvaluatePCARecon",
                            "EvaluateICARecon",
                            "ExpressionDataTableToMatrixT",
                            "GetRMSE")) %dopar% {
                              if (method == "PCA") {
                                EvaluatePCARecon(train.list[[i]][[j]], 
                                                 test.list[[mthd]][[j]])
                              } else {
                                EvaluateICARecon(train.list[[i]][[j]], 
                                                 test.list[[mthd]][[j]])
                              }
                            }
      } else {
        foreach(j = 1:length(train.list[[i]]),
                .export = c("EvaluatePCARecon",
                            "EvaluateICARecon",
                            "ExpressionDataTableToMatrixT",
                            "GetRMSE")) %dopar% {
                              if (method == "PCA") {
                                EvaluatePCARecon(train.list[[i]][[j]], 
                                                 test.list[[mthd]])
                              } else {
                                EvaluateICARecon(train.list[[i]][[j]], 
                                                 test.list[[mthd]])
                              }
                            }
      }
    }
  } else if (platform == "array") {
    comp.list <- foreach(i = 1:length(train.list)) %do% {
      mthd <- names(train.list)[i]
      if (mthd == "tdm") {  
        foreach(j = 1:length(train.list[[i]]),
                .export = c("EvaluatePCARecon",
                            "EvaluateICARecon",
                            "ExpressionDataTableToMatrixT",
                            "GetRMSE")) %dopar% {
                              if (method == "PCA") {
                                EvaluatePCARecon(train.list[[i]][[j]], 
                                                 test.list[["log"]])
                              } else {
                                EvaluateICARecon(train.list[[i]][[j]], 
                                                 test.list[["log"]])
                              }
                            }
      } else {
        foreach(j = 1:length(train.list[[i]]),
                .export = c("EvaluatePCARecon",
                            "EvaluateICARecon",
                            "ExpressionDataTableToMatrixT",
                            "GetRMSE")) %dopar% {
                              if (method == "PCA") {
                                EvaluatePCARecon(train.list[[i]][[j]], 
                                                 test.list[[mthd]])
                              } else {
                                EvaluateICARecon(train.list[[i]][[j]], 
                                                 test.list[[mthd]])
                              }
                            }
      }
    }
  } else {
    stop("platform must be 'array' or 'seq'")
  }
  
  parallel::stopCluster(cl)
  
  # sort out names
  names(comp.list) <- names(train.list)
  for (i in 1:length(train.list)){
    names(comp.list[[i]]) <- names(train.list[[i]])
  }
  
  mlt.c <- reshape::melt(comp.list)
  mlt.c <- cbind(mlt.c, rep(method, nrow(mlt.c)))
  colnames(mlt.c) <- c("value", "seq.level", "norm.method", "recon.method")
  return(mlt.c)
  
}

