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
    errors <- y - y.pred
    # mean of true expression values
    y.bar <- mean(y)
    # calculate absolute scaled error
    scaled.errors <- errors/(sum(abs(y - y.bar)) / n)
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
