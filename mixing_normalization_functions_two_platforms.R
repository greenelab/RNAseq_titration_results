NAToZero <- function(dt, un=0) gdata::NAToUnknown(dt, un)

GetTitratedSampleNames <- function(sample.set, p){
  # This function takes a vector of sample names and selects some percentage p
  # of them -- to be used to concatenate array and RNA-seq data together
  #
  # Args:
  #   sample.set: a character vector of sample names
  #   p: the percentage of samples to select
  #
  # Returns:
  #   selected.samples: the selected subset of sample names, character
  #
  sample.set <- as.character(sample.set)
  no.select <- ceiling(length(sample.set) * p)  # number of samples to select 
  selected.samples <- sample(sample.set, no.select)
  return(selected.samples)
}

ZScoreProcessing <- function(array.dt, seq.dt){
  # This function takes array and RNA-seq data in the form of data.table 
  # to be 'mixed' (concatenated) and returns one z-scored, 
  # zero to one transformed data.table
  # 
  # Args:
  #   array.dt: data.table of array data where the first column contains 
  #             gene identifiers, the columns are samples, 
  #             rows are gene measurements
  #   seq.dt:   data.table of RNA-seq data where the first column contains 
  #             gene identifiers, the columns are samples, 
  #             rows are gene measurements
  # 
  # Returns:
  #   zto.z.dt: z-scored, zero to one transformed data.table that contains both
  #             array and RNA-seq samples
  require(TDM)
  require(data.table)
  # Error-handling
  array.is.dt <- "data.table" %in% class(array.dt)
  seq.is.dt <- "data.table" %in% class(seq.dt)
  any.not.dt <- !(any(c(array.is.dt, seq.is.dt)))
  if (any.not.dt) {
    stop("array.dt and seq.dt must both be data.tables")
  }
  if (!(all(array.dt[[1]] %in% seq.dt[[1]]))) {
    stop("Gene identifiers in data.tables must match")
  }
  cat("Z-score...\n")
  array.mat <- data.frame(array.dt[, 2:ncol(array.dt), with = F])
  z.array <- t(apply(array.mat, 1, function(x) scale(as.numeric(x))))
  seq.mat <- data.frame(seq.dt[, 2:ncol(seq.dt), with = F])
  z.seq <- t(apply(seq.mat, 1, function(x) scale(as.numeric(x))))
  cat("\tConcatenation...\n")
  z.dt <- data.table(cbind(array.dt[[1]], z.array, z.seq))
  colnames(z.dt) <- c("gene", 
                      chartr(".", "-", colnames(array.mat)),
                      chartr(".", "-", colnames(seq.mat)))
  # z.dt <- NAToZero(z.dt)
  cat("\tZero to one transformation...\n")
  zto.z.dt <- zero_to_one_transform(z.dt)
  return(zto.z.dt)
}

QNProcessing <- function(array.dt, seq.dt){
  # This function takes array and RNA-seq data in the form of data.table 
  # to be 'mixed' (concatenated) and returns one quantile normalized, 
  # zero to one transformed data.table. The array data is used as the target 
  # distribution.
  # 
  # Args:
  #   array.dt: data.table of array data where the first column contains 
  #             gene identifiers, the columns are samples, 
  #             rows are gene measurements
  #   seq.dt:   data.table of RNA-seq data where the first column contains 
  #             gene identifiers, the columns are samples, 
  #             rows are gene measurements
  # 
  # Returns:
  #   zto.qn.cat:  quantile normalized, zero to one transformed data.table 
  #                that contains both array and RNA-seq samples
  require(preprocessCore)
  require(TDM)
  require(data.table)
  # Error-handling
  array.is.dt <- "data.table" %in% class(array.dt)
  seq.is.dt <- "data.table" %in% class(seq.dt)
  any.not.dt <- !(any(c(array.is.dt, seq.is.dt)))
  if (any.not.dt) {
    stop("array.dt and seq.dt must both be data.tables")
  }
  if (!(all(array.dt[[1]] %in% seq.dt[[1]]))) {
    stop("Gene identifiers in data.tables must match")
  }
  ref.values <- data.frame(array.dt[, 2:ncol(array.dt), with = F])
  target.values <- data.frame(seq.dt[, 2:ncol(seq.dt), with = F])
  cat("Quantile normalization...\n")
  qn.targ <- normalize.quantiles.determine.target(
    data.matrix(ref.values), target.length = nrow(ref.values))
  qn.seq <- normalize.quantiles.use.target(data.matrix(target.values), qn.targ,
                                           copy = F)
  qn.seq <- data.table(cbind(seq.dt[[1]], qn.seq))
  colnames(qn.seq) <- chartr(".", "-", colnames(seq.dt))
  # array.dt <- NAToZero(array.dt)
  # qn.seq <- NAToZero(qn.seq)
  cat("\tZero to one transformation...\n")
  zto.array.dt <- zero_to_one_transform(array.dt)
  zto.qn.seq <- zero_to_one_transform(qn.seq)
  cat("\tConcatenation...\n")
  zto.qn.cat <- data.table(cbind(zto.array.dt, zto.qn.seq[, 2:ncol(zto.qn.seq),
                                                          with=F]))
  return(zto.qn.cat)
}

NPNProcessing <- function(array.dt, seq.dt){
  # This function takes array and RNA-seq data in the form of data.table 
  # to be 'mixed' (concatenated) and returns one nonparanormal normalized, 
  # zero to one transformed data.table
  # 
  # Args:
  #   array.dt: data.table of array data where the first column contains 
  #             gene identifiers, the columns are samples, 
  #             rows are gene measurements
  #   seq.dt:   data.table of RNA-seq data where the first column contains 
  #             gene identifiers, the columns are samples, 
  #             rows are gene measurements
  # 
  # Returns:
  #   npn.cat: NPN normalized, zero to one transformed data.table that contains both
  #            array and RNA-seq samples
  require(huge)
  require(TDM)
  require(data.table)
  # Error-handling
  array.is.dt <- "data.table" %in% class(array.dt)
  seq.is.dt <- "data.table" %in% class(seq.dt)
  any.not.dt <- !(any(c(array.is.dt, seq.is.dt)))
  if (any.not.dt) {
    stop("array.dt and seq.dt must both be data.tables")
  }
  if (!(all(array.dt[[1]] %in% seq.dt[[1]]))) {
    stop("Gene identifiers in data.tables must match")
  }
  ref.values <- data.frame(array.dt[, 2:ncol(array.dt), with = F])
  target.values <- data.frame(seq.dt[, 2:ncol(seq.dt), with = F])
  npn.ref <- data.matrix(ref.values)
  npn.array <- huge.npn(t(npn.ref), npn.func = "shrinkage", npn.thresh = NULL, 
                        verbose = TRUE)
  npn.targ <- data.matrix(target.values)
  npn.seq <- huge.npn(t(npn.targ), npn.func = "shrinkage", npn.thresh = NULL, 
                      verbose = TRUE)
  cat("\tConcatenation...\n")
  npn.cat <- data.table(cbind(array.dt[[1]], t(npn.array), t(npn.seq)))
  colnames(npn.cat) <- c("gene",
                         chartr(".", "-", colnames(ref.values)),
                         chartr(".", "-", colnames(target.values)))
  
  # npn.cat <- NAToZero(npn.cat)
  cat("\tZero to one transformation...\n")
  zto.npn.cat <- zero_to_one_transform(npn.cat)
  return(zto.npn.cat)
}

TDMProcessing <- function(array.dt, seq.dt){
  # This function takes array and RNA-seq data in the form of data.table 
  # to be 'mixed' (concatenated) and returns one TDM normalized, 
  # zero to one transformed data.table. The array data is used as the reference 
  # distribution.
  # 
  # Args:
  #   array.dt: data.table of array data where the first column contains 
  #             gene identifiers, the columns are samples, 
  #             rows are gene measurements
  #   seq.dt:   data.table of RNA-seq data where the first column contains 
  #             gene identifiers, the columns are samples, 
  #             rows are gene measurements
  # 
  # Returns:
  #   tdm.cat: TDM normalized, zero to one transformed data.table 
  #                that contains both array and RNA-seq samples
  require(TDM)
  require(data.table)
  # Error-handling
  array.is.dt <- "data.table" %in% class(array.dt)
  seq.is.dt <- "data.table" %in% class(seq.dt)
  any.not.dt <- !(any(c(array.is.dt, seq.is.dt)))
  if (any.not.dt) {
    stop("array.dt and seq.dt must both be data.tables")
  }
  if (!(all(array.dt[[1]] %in% seq.dt[[1]]))) {
    stop("Gene identifiers in data.tables must match")
  }
  cat("TDM transformation...\n")
  tdm.seq <- tdm_transform(target_data = seq.dt, 
                           ref_data = array.dt, 
                           negative = FALSE, 
                           filter_p = FALSE, 
                           inv_reference = TRUE, 
                           log_target=TRUE)
  # array.dt <- NAToZero(array.dt)
  # tdm.seq <- NAToZero(tdm.seq)
  cat("\tZero to one transformation...\n")
  zto.array <- zero_to_one_transform(array.dt)
  zto.tdm.seq <- zero_to_one_transform(tdm.seq)
  cat("\tConcatenation...\n")
  tdm.cat <- data.table(cbind(zto.array, zto.tdm.seq[, 2:ncol(zto.tdm.seq), 
                                                     with = F]))
  return(tdm.cat)
}

LOGProcessing <- function(array.dt, seq.dt){
  # This function takes array and RNA-seq data in the form of data.table 
  # to be 'mixed' (concatenated) and returns one log transformed, 
  # zero to one transformed data.table
  # 
  # Args:
  #   array.dt: data.table of array data where the first column contains 
  #             gene identifiers, the columns are samples, 
  #             rows are gene measurements
  #   seq.dt:   data.table of RNA-seq data where the first column contains 
  #             gene identifiers, the columns are samples, 
  #             rows are gene measurements
  # 
  # Returns:
  #   log.cat:  log transformed, zero to one transformed data.table 
  #             that contains both array and RNA-seq samples
  require(TDM)
  require(data.table)
  # Error-handling
  array.is.dt <- "data.table" %in% class(array.dt)
  seq.is.dt <- "data.table" %in% class(seq.dt)
  any.not.dt <- !(any(c(array.is.dt, seq.is.dt)))
  if (any.not.dt) {
    stop("array.dt and seq.dt must both be data.tables")
  }
  if (!(all(array.dt[[1]] %in% seq.dt[[1]]))) {
    stop("Gene identifiers in data.tables must match")
  }
  cat("Log transformation seq data...\n")
  log.seq <- log_transform_p1(seq.dt)
  # array.dt <- NAToZero(array.dt)
  # log.seq <- NAToZero(log.seq)
  cat("\tZero to one transformation...\n")
  zto.array <- zero_to_one_transform(array.dt)
  zto.log.seq <- zero_to_one_transform(log.seq)
  cat("\tConcatenation...\n")
  log.cat <- data.table(cbind(zto.array, zto.log.seq[, 2:ncol(zto.log.seq), 
                                                     with = F]))
  return(log.cat)
}

NormalizationWrapper <- function(array.dt, seq.dt){
  # This function takes array and RNA-seq data in the form of data.table 
  # to be 'mixed' (concatenated) and returns a list of normalized data.tables
  #
  # Args:
  #   array.dt: data.table of array data where the first column contains 
  #             gene identifiers, the columns are samples, 
  #             rows are gene measurements
  #   seq.dt:   data.table of RNA-seq data where the first column contains 
  #             gene identifiers, the columns are samples, 
  #             rows are gene measurements
  # 
  # Returns:
  #	  norm.list: a list of normalized, zero to one 'mixed' data.tables (log 
  #              transformation, nonparanormal normalized, quantile normalized, 
  #              TDM normalized, z-scored)
  require(data.table)
  require(TDM)
  # if any negative values are found, then inverse log transform and relog 
  # transform using x+1
  any.negative <- any(as.vector(as.matrix(array.dt[, 2:ncol(array.dt),
                                                   with=F])) < 0)
  if (any.negative) {
    cat("\nLog transformation array data...\n")
    array.dt <- inv_log_transform(array.dt)
    array.dt <- log_transform_p1(array.dt)
  }
  # convert NA to zero
  array.dt <- NAToZero(array.dt)  
  norm.list <- list()
  # save array data.table to be used as 'reference' for test data
  norm.list[["raw.array"]] <- array.dt
  # normalization methods, concatenation, zero to one transformation
  norm.list[["log"]] <- LOGProcessing(array.dt, seq.dt)
  norm.list[["npn"]] <- NPNProcessing(array.dt, seq.dt)
  norm.list[["qn"]] <- QNProcessing(array.dt, seq.dt)
  norm.list[["tdm"]] <- TDMProcessing(array.dt, seq.dt)
  norm.list[["z"]] <- ZScoreProcessing(array.dt, seq.dt)
  return(norm.list) 
}

GetDataTablesForMixing <- function(array.data, seq.data,
                                   titrate.sample.names){
  # This function takes two full data.tables that contain array and RNA-seq data,
  # matched samples and a vector of sample names. It returns a list that contains 
  # an array data.table (samples not in the sample names) and an RNA-seq data.table
  # (samples in the sample names).
  #
  # Args:
  #   array.data: data.table of array data where the first column contains 
  #               gene identifiers, the columns are samples, 
  #               rows are gene measurements
  #   seq.data:   data.table of RNA-seq data where the first column contains 
  #               gene identifiers, the columns are samples, 
  #               rows are gene measurements
  #   titrate.sample.names: samples to be taken from RNA-seq data.table
  #
  # Returns:
  #   mix.dt.list: list that contains the array data.table ("array") and 
  #                the RNA-seq data.table ("seq") to be normalized 
  #                and concatenated
  require(data.table)
  # Error-handling   
  array.is.dt <- "data.table" %in% class(array.data)
  seq.is.dt <- "data.table" %in% class(seq.data)
  any.not.df <- !(any(c(array.is.dt, seq.is.dt)))
  unequal.dim <- !(all(dim(array.data) == dim(seq.data)))
  if (any.not.df | unequal.dim) {
    stop("array.data and seq.data must both be data.tables of equal dimensions")
  }
  if (!(all(array.data[[1]] %in% seq.data[[1]]))) {
    stop("Gene identifiers in data.tables must match")
  }
  if ((length(titrate.sample.names) > 0) 
  	& (length(titrate.sample.names) != (ncol(seq.data) - 1))) {
    array.dt <- 
      array.data[, which(!(colnames(array.data) %in% titrate.sample.names)),
                 with = F]
    array.dt <- data.table(array.dt)
    seq.dt <- 
      seq.data[, c(1, which(colnames(seq.data) %in% titrate.sample.names)), 
               with = F]
    seq.dt <- data.table(seq.dt)
    mix.dt.list <- list("array" = array.dt, "seq" = seq.dt)
  } else if (length(titrate.sample.names) == (ncol(seq.data) - 1)) {
    mix.dt.list <- list("seq" = data.table(seq.data))
  } else {
    mix.dt.list <- list("array" = data.table(array.data))
  }
}
