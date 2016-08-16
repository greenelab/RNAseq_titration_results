# to be shared between multiple methods
GetMASE <- function(true.mat, recon.mat){
  # This function takes two gene expression matrices (before and after 
  # reconstruction) and compares them .
  # The reconstruction error between the two (MASE: mean squared absolute error) 
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

EvaluatePCARecon <- function(train.dt, test.dt, no.comp = 50){
  # This function takes training and test set data.tables
  # performs PCA on training, projects test data into training PC space and 
  # back out - quantifies 'reconstruction error' using MASE
  # 
  # Args:
  #   train.dt: a data.table of training data
  #   test.dt: a data.table of test data
  #   no.comp: number of principal components to use for reconstruction
  # 
  # Returns:
  #   A list with the values
  #     MASE: mean(mean absolute scaled error between test values and 
  #                reconstructed test values, per gene basis)
  #     PC: the PrcompDataTable output (principal components analysis 
  #         on the training data)
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
    # perform PCA on a data.table of gene expression values 
    # where genes are rows and columns are samples
    mat.t <- ExpressionDataTableToMatrixT(dataTable)
    pca <- prcomp(mat.t)
    return(pca)
  }
  
  PCNewDataTable <- function(pc, dataTable, no.comp){
    # Reconstruct test expression data using training PC space
    # Error-handling 
    if (class(pc) != "prcomp"){
      stop("pc must be prcomp object")
    }
    # get new data (dataTable) into correct form (matrix, transposed)
    mat.t <- ExpressionDataTableToMatrixT(dataTable)
    proj.test <- scale(mat.t, pc$center, pc$scale) %*% pc$rotation[, 1:no.comp] 
    recon.test <- proj.test %*% t(pc$rotation[, 1:no.comp])
    return(recon.test)
  }
  
  # PCA on training data
  pc <- PrcompDataTable(train.dt)
  # project test data into training PC space & back out
  recon.test <- PCNewDataTable(pc, test.dt, no.comp)
  # transpose mat for test data, to be used as "true.mat"
  test.dm.t <- ExpressionDataTableToMatrixT(test.dt)
  
  # quantification: how well did we do at reconstruction?
  mase <- GetMASE(true.mat = test.dm.t, recon.mat = recon.test)
  names(mase) <- colnames(test.dm.t)
  
  # return reconstruction error and PCA output
  return(list("MASE" = mase, "COMP" = pc))
  
}

EvaluateICARecon <- function(train.dt, test.dt, no.comp = 50){
  # This function takes training and test set data.tables
  # performs ICA on training, projects test data into PC space, then IC space
  # and back out - quantifies 'reconstruction error' using MASE
  # 
  # Args:
  #   train.dt: a data.table of training data
  #   test.dt: a data.table of test data
  #   no.comp: number of independent components to use
  # 
  # Returns:
  #   A list with the values
  #     MASE: mean(mean absolute scaled error between test values and 
  #                reconstructed test values, per gene basis)
  #     IC: the FastICADataTable output (independent components analysis 
  #         on the training data)
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
    # perform ICA on a data.table of gene expression values 
    # where genes are rows and columns are samples
    mat.t <- ExpressionDataTableToMatrixT(dataTable)
    ica <- fastICA::fastICA(mat.t, n.comp = no.comp, alg.typ = "parallel",
                            fun = "exp", method = "C")
    return(ica)
  }
  
  ICNewDataTable <- function(ic, dataTable){
    # Reconstruct test expression data using training IC space
    # Error-handling
    # Are the correct matrices in the list ic
    correct.matrices <- all(names(ic) == c("X", "K", "W", "A", "S"))
    if (!(correct.matrices)) {
      stop("ic must be output of fastICA::fastICA")
    }
    # get new data (dataTable) into correct form (matrix, transposed)
    mat.t <- ExpressionDataTableToMatrixT(dataTable)
    whiten.mat <- (scale(mat.t, center=TRUE, scale=FALSE) %*% ic$K )
    recon.mat <- (whiten.mat %*% ic$W %*% ic$A)
    return(recon.mat)
  }
  
  # ICA on training data
  ic <- FastICADataTable(train.dt, no.comp)
  # project test data into training IC space & back out
  recon.test <- ICNewDataTable(ic, test.dt)
  # transpose mat for test data, to be used as "true.mat"
  test.dm.t <- ExpressionDataTableToMatrixT(test.dt)
  
  # quantification: how well did we do at reconstruction?
  mase <- GetMASE(true.mat = test.dm.t, recon.mat = recon.test)
  
  # return reconstruction error and ICA output
  return(list("MASE" = mase, "COMP" = ic))
  
}
