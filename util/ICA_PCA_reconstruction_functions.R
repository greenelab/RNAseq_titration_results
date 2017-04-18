GetMASE <- function(true.mat, recon.mat){
  # This function takes two gene expression matrices (before and after 
  # reconstruction) and compares them .
  # The reconstruction error between the two (MASE: mean absolute scaled error) 
  # is calculated on a per gene basis and then summarized using the mean.
  # 
  # Args:
  #   true.mat: gene expression matrix before reconstruction, samples are rows
  #             genes are columns
  #   recon.mat: gene expression matrix after reconstruction, samples are rows
  #              genes are columns
  # 
  # Returns:
  #   MASE: MASE before and after reconstruction, per gene basis,
  #         a vector
  # 
  
  # Error-handling  
  if (all.equal(dim(true.mat), dim(recon.mat)) != TRUE){
    stop("true.mat and recon.mat must be of equal dimensions")
  }
  
  CalculateMASE <- function(y, y.pred){
    # For a pair of vectors, a gene's measured expression value (y) and a gene's
    # reconstructed expression value (y.pred), calculate the genes mean absolute
    # scaled error.
    # number of observations
    n <- length(y)
    # difference in gene expression and reconstructed gene expression
    abs.errors <- abs(y - y.pred)
    # mean of true expression values
    y.bar <- mean(y)
    # calculate absolute scaled error
    scaled.errors <- abs.errors/(sum(abs(y - y.bar)) / n)
    # calculate mean absolute scaled error
    mase <- mean(scaled.errors)
    return(mase)
  }
  
  # looks like for loop is actually faster than apply family of functions here
  mase <- vector()
  for (col in 1:ncol(true.mat)){
    # for each gene (column), calculate the MASE between the true expression 
    # values and the expression values after reconstruction
    mase[col] <- CalculateMASE(y = true.mat[, col], 
                               y.pred = recon.mat[, col])
  }
  
  # return a vector of mean absolute scaled errors
  return(mase)
  
}

CalculateReconstructionError <- function(recon.test, test.dt) {
  # This function calculates the reconstruction error, here mean absolute 
  # scaled error (MASE) on a per-gene basis, between a matrix of reconstructed 
  # data and a data.table of test (holdout) data. 
  #
  # Args:
  #   recon.test: gene expression matrix after reconstruction, samples are rows
  #               genes are columns
  #   test.dt: data.table of gene expression data before reconstruction, 
  #            where the first column contains gene identifiers, columns are 
  #            samples
  # 
  # Returns:
  #   mase: MASE before and after reconstruction, per gene basis,
  #         a vector
  
  test.dm.t <- ExpressionDataTableToMatrixT(test.dt)
  # quantification: how well did we do at reconstruction?
  mase <- GetMASE(true.mat = test.dm.t, recon.mat = recon.test)
  return(mase)
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
  # include gene identifiers
  rownames(data.mat) <- dataTable[[1]]
  return(t(data.mat))  # transpose
}

GetReconstructedDataTable <- function(t.data.mat){
  # This function takes a matrix (rows - samples; columns - genes) and 
  # returns a data.table where the gene identifiers (data matrix colnames) are 
  # the first column and genes are rows, samples are columns. It is intended to
  # be used to get reconstructed test data into a form that can be used with
  # models trained for BRCA 
  # 
  # Args:
  #   t.data.mat: a matrix where genes are columns and rows are samples
  #              
  # Returns:
  #   dataTable: a data.table where rows are genes and columns are samples; 
  #              the first column contains gene identifiers
  # 
  dataTable <- 
    data.table::data.table(cbind(colnames(t.data.mat), t(t.data.mat)))
  colnames(dataTable) <- chartr(".", "-", colnames(dataTable))
  colnames(dataTable)[1] <- "gene"
  return(dataTable)
}

PerformCompAnalysis <- function(train.dt, no.components = 50, method = "PCA") {
  # This function performs components analysis, either Principal Components
  # Analysis (PCA) or Independent Components Analysis (ICA) (as determined
  # by the method arg), on a single data.table of gene expression using
  # a user-specified number of components.
  #
  # Args:
  #   train.dt: a data.table of gene expression data where the first column 
  #             contains gene identifiers, columns are samples
  #   no.components: an integer that determines the number of components to
  #                  be extracted with ICA (default is 50)
  #   method: which method should be used - either PCA (default) or ICA
  #
  # Returns:
  #   comp.obj: either a prcomp object (method = "PCA") or the output of
  #             fastICA::fastICA (a list of matrices; method = "ICA")    
  
  # internal functions   
  PrcompDataTable <- function(dataTable){
    # perform PCA on a data.table of gene expression values 
    # where genes are rows and columns are samples
    mat.t <- ExpressionDataTableToMatrixT(dataTable)
    pca <- prcomp(mat.t)
    return(pca)
  }
  
  FastICADataTable <- function(dataTable, no.comp){
    # perform ICA on a data.table of gene expression values 
    # where genes are rows and columns are samples
    mat.t <- ExpressionDataTableToMatrixT(dataTable)
    ica <- fastICA::fastICA(mat.t, n.comp = no.comp, alg.typ = "parallel",
                            fun = "exp", method = "C")
    return(ica)
  }

  # error-handling
  train.is.dt <- "data.table" %in% class(train.dt)
  if (!(train.is.dt)) {
    stop("train.dt must be a data.table")
  }
  
  if (!(any(c(method == "PCA", method == "ICA")))) {
    stop("method must be 'ICA' or 'PCA'")
  }
 
  if (method == "PCA") {
    comp.obj <- PrcompDataTable(train.dt)
  } else {
    comp.obj <- FastICADataTable(train.dt, no.comp = no.components)
  }
  
  return(comp.obj)
   
}


TrainSetCompAnalysis <- function(train.list, 
                                 num.comp = 50, 
                                 comp.method = "PCA") {
  # This function performs components analysis, either ICA or PCA, on the list
  # of training data that contains different normalization methods and 
  # percentages of RNA-seq samples. 
  #
  # Args:
  #   train.list: a nested list of normalized training data, organized first by
  #               normalization method and then by % sequencing data -- the
  #               output of 1-normalize_titrated_data.R
  #   num.comp: an integer that determines the number of components to
  #             be extracted with ICA (default is 50)
  #   comp.method: which method should be used - either PCA (default) or ICA
  #
  # Returns:
  #   comp.list: a list of component objects (either prcomp or list of matrices, 
  #              for PCA and ICA, respectively) -- one for each data.table of 
  #              training data
    
  # start parallel backend
  cl <- parallel::makeCluster(2)
  doParallel::registerDoParallel(cl)
  parallel::clusterExport(cl,
                          c("PerformCompAnalysis",
                            "ExpressionDataTableToMatrixT"))
  
  suppressMessages(require(foreach))
  
  comp.list <- 
    # for each normalization method
    foreach(norm.iter = seq_along(train.list)) %do% {
      # identify by method name, rather than index
      norm.mthd <- names(train.list)[norm.iter]
      # for each % seq data - 0 - 100 for the majority of methods
      foreach(seq.iter = seq_along(train.list[[norm.mthd]])) %dopar% {
        # identify by % seq name, rather than index
        seq.lvl <- names(train.list[[norm.mthd]])[seq.iter]
        # get the input data.table
        input.train.dt <- train.list[[norm.mthd]][[seq.lvl]]
        # perform ICA or PCA on the input training data
        PerformCompAnalysis(train.dt = input.train.dt,
                            no.components = num.comp,
                            method = comp.method)
      }  
    }
 
  # stop parallel backend
  parallel::stopCluster(cl)
  
  # sort out names
  names(comp.list) <- names(train.list)
  for (nrm.iter in 1:length(train.list)){
    names(comp.list[[nrm.iter]]) <- names(train.list[[nrm.iter]])
  }
  
  # return the list of components objects
  return(comp.list)
  
}

PCNewDataTable <- function(pc, dataTable, no.comp){
  # This function reconstructs test expression data using training PC space,
  # for internal use by PerformReconstructionExperiment
  #
  # Args:
  #   pc: prcomp object
  #   dataTable: a data.table of normalized test expression data, typically
  #              test or holdout data, to be projected onto the training space
  #              (pc) and back out, i.e., reconstructed
  #   no.comp: the number of PCs to use for the reconstruction
  #
  # Returns:
  #   recon.test: reconstructed test data
  
  # Error-handling 
  if (class(pc) != "prcomp"){
    stop("pc must be prcomp object")
  }
  # get new data (dataTable) into correct form (matrix, transposed)
  mat.t <- ExpressionDataTableToMatrixT(dataTable)
  proj.test <- scale(mat.t, pc$center, pc$scale) %*% pc$rotation[, 1:no.comp] 
  # recon.test <- proj.test %*% t(pc$rotation[, 1:no.comp])
  # use tcrossprod, has been optimized for matrix calculations
  recon.test <- tcrossprod(proj.test, pc$rotation[, 1:no.comp])
  return(recon.test)
}

ICNewDataTable <- function(ic, dataTable){
  # This function reconstructs test expression data using training IC space,
  # for internal use by PerformReconstructionExperiment
  #
  # Args:
  #   pc: list returned by fastICA::fastICA
  #   dataTable: a data.table of normalized test expression data, typically
  #              test or holdout data, to be projected onto the training space
  #              (ic) and back out, i.e., reconstructed
  #
  # Returns:
  #   recon.test: reconstructed test data
  
  # Error-handling
  # Are the correct matrices in the list ic
  correct.matrices <- all(names(ic) == c("X", "K", "W", "A", "S"))
  if (!(correct.matrices)) {
    stop("ic must be output of fastICA::fastICA")
  }
  # get new data (dataTable) into correct form (matrix, transposed)
  mat.t <- ExpressionDataTableToMatrixT(dataTable)
  whiten.mat <- (scale(mat.t, center=TRUE, scale=FALSE) %*% ic$K)
  recon.mat <- whiten.mat %*% (ic$W %*% ic$A)
  return(recon.mat)
}

PerformReconstructionExperiment <- function(train.comp, test.dt, 
                                            n.comps = 50,
                                            comp.method = "PCA") {
  # This function is the main lift of reconstruction experiment -- it will 
  # perform reconstruction and calculate reconstruction error (MASE). For
  # internal use by ParallelReconstructionFunction.
  #
  # Args:
  #   train.comp: training component object (output of prcomp or 
  #               fastICA::fastICA)
  #   n.comps: number of components to use for reconstruction (required for PCA 
  #            only)
  #   comp.method: method to be used, either PCA (default) or ICA
  #
  # Returns:
  #   a list with the following elements:
  #     recon: the reconstructed data in data.table form (columns are samples)
  #            this is suitable for use with the models trained in 
  #            2-train_test_brca_subtype.R
  #     mase: the reconstruction error (MASE) on a per-gene basis
  
  # project test data (test.dt) onto training space and reconstruct
  if (comp.method == "PCA") {
    recon.test <- PCNewDataTable(train.comp, test.dt, no.comp = n.comps)
  } else {
    recon.test <- ICNewDataTable(train.comp, test.dt)
  }
  # get reconstructed data in a form that will be amenable to 
  # prediction later
  recon.dt <- GetReconstructedDataTable(recon.test)
  # calculate MASE
  error.measure <- CalculateReconstructionError(recon.test, test.dt)
  # return a list that contains a data.table of reconstructed data
  # and the MASE
  return(list("recon" = recon.dt, "mase" = error.measure))
}

ParallelReconstructionFunction <- function(comp.list,
                                           test.data.list,
                                           number.components = 50,
                                           c.method = "PCA"){
  # This function performs the reconstruction experiment
  # (PerformReconstructionExperiment) on all the data from a particular 
  # normalization method (i.e., all % of sequencing data included) in parallel
  # for all sequencing levels. For internal use by ReconstructionWrapper.
  #
  # Args:
  #   comp.list: list of training component objects (output of prcomp or
  #              fastICA::fastICA) for the normalization method under 
  #              consideration (typically with seq levels 0-100)
  #   test.data.list: a list or data.table of test (holdout) data for the 
  #                   normalization method under consideration; can be either 
  #                   array or RNA-seq data 
  #   number.components: number of components to be used for reconstruction 
  #                      (default is 50)
  #   c.method: method to be used, PCA (default) or ICA
  #
  # Returns:
  #   return.list: list of reconstructed data and error measures (the output
  #                 of PerformReconstructionExperiment
    
  return.list <-
    # for each % sequencing data
    foreach(seq.iter = seq_along(comp.list)) %dopar% {  
      # identify by % seq name, rather than index
      seq.lvl <- names(comp.list)[seq.iter]
      # get corresponding ICA or PCA results from training
      comp.obj <- comp.list[[seq.lvl]]
      # if the test.data.list data has seq levels, i.e., is a 
      # list, loop through it. This will happen in the case 
      # of RNA-seq test data for TDM and QN norm methods.
      # Otherwise, the 'test.list' argument is a data.table, 
      # make it 'test.dt'
      if (class(test.data.list) == "list") {
        test.dt <- test.data.list[[seq.lvl]]
      } else {
        test.dt <- test.data.list
      }
      # perform the reconstruction experiment for each seq level
      PerformReconstructionExperiment(train.comp = comp.obj,
                                      test.dt = test.dt,
                                      n.comps = number.components,
                                      comp.method = c.method)
    }
  
  names(return.list) <- names(comp.list)
  return(return.list)
  
}

ReconstructionWrapper <- function(train.list, test.list, 
                                  num.comps = 50) {
  # This function performs reconstruction experiment on all test data -- i.e.,
  # all normalization methods which will be input to 
  # ParallelReconstructionFunction which will PerformReconstructionExperiment
  # for all test data.
  #
  # Args: 
  #   train.list: list of training component objects; the output of
  #               TrainSetCompAnalysis
  #   test.list: list of test/holdout data, the top-level should be 
  #              the normalization method (i.e., one platform [array or seq]
  #              should be input at at time) -- test list is typically
  #              output from 1-normalize_titrated_data.R
  #   num.comps: number of components to be used for reconstruction (for PCA
  #              reconstruction)
  # 
  # Returns: 
  #   a list with the following elements:
  #     recon: a list of the reconstructed test/holdout data as data.tables;
  #            this is appropriate for use with PredictWrapper
  #     mase.df: a data.frame of the reconstruction errors that were calculated
  #              for all norm methods and % seq (mean absolute scaled error)
  #   
  
  suppressMessages(require(foreach))
  
  # get method from class of objects in train.list
  # 'method detection'
  is.pca <- class(train.list[[1]][[1]]) == "prcomp"
  is.ica <- all(names(train.list[[1]][[1]]) == c("X", "K", "W", "A", "S"))
  # if the train.list isn't pca or ica output
  if (!(any(c(is.pca, is.ica)))) {
    stop("train.list must contain output of prcomp or fastICA::fastICA")
  }
  # if any individual element is not one of these objects, it will be 
  # caught by PCNewDataTable or ICNewDataTable
  
  # assign comp.mthd for use with ParallelReconstructionFunction
  if (is.pca) {
    comp.mthd <- "PCA" 
  } else {  # is.ica = TRUE, because of error handling above
    comp.mthd <- "ICA"
  }

  # start parallel backend
  cl <- parallel::makeCluster(2)
  doParallel::registerDoParallel(cl)
  # these functions are required for prediction
  parallel::clusterExport(cl,
                          c("GetMASE",
                            "GetReconstructedDataTable",
                            "ExpressionDataTableToMatrixT",
                            "PCNewDataTable", 
                            "ICNewDataTable",
                            "PerformReconstructionExperiment",
                            "CalculateReconstructionError"))
  
  norm.methods <- names(train.list)
  # for each normalization method
  recon.return.list <- 
    foreach(mthd.iter = seq_along(train.list)) %do% {
      mthd <- norm.methods[mthd.iter]  # use method name rather than index
      # if the method is tdm and there's no tdm data to use for prediction
      # use log data -- this will be the case for array hold out sets
      if (mthd == "tdm" & !("tdm" %in% names(test.list))) {
        input.list <- test.list[["log"]]
      } else {  # otherwise, just use the corresponding data for prediction
        input.list <- test.list[[mthd]]
      }
      # do parallel prediction for all seq levels
      ParallelReconstructionFunction(comp.list = train.list[[mthd]],
                                     test.data.list = input.list,
                                     number.components = num.comps,
                                     c.method = comp.mthd)
    }
  
  # names
  names(recon.return.list) <- norm.methods
  
  # save MASE into a data.frame
  mase.df <- 
    reshape2::melt(lapply(recon.return.list,
                          function(x) parallel::mclapply(x, 
                                                         function(y) y$mase)))
  # add comp method used, using cbind as rep produces a character vector
  mase.df <- cbind(mase.df, rep(comp.mthd, nrow(mase.df)))
  
  # include gene ids
  gene.ids <- test.list[[1]][[1]]
  n.reps <- nrow(mase.df) / length(gene.ids)
  mase.df <- dplyr::bind_cols(data.frame(rep(gene.ids, n.reps)),
                               mase.df)
  colnames(mase.df) <- c("gene", "MASE", "perc.seq", "norm.method", 
                         "comp.method")
  
  # restructure recon objects (without error measures) to be written
  # to RDS
  recon.list <- lapply(recon.return.list,
                       function(x) parallel::mclapply(x, 
                                                      function(y) y$recon))
  
  # stop parallel backend
  parallel::stopCluster(cl)
  
  return(list("recon" = recon.list, "mase.df" = mase.df))
  
}
