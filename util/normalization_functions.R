NAToZero <- function(dt, un=0) suppressWarnings(gdata::NAToUnknown(dt, un))

#### single platform functions -------------------------------------------------
LOGArrayOnly <- function(array.dt, zero.to.one = TRUE){
  # This function takes array data in the form of a data.table and returns a
  # log transformed, zero to one transformed (if zero.to.one = TRUE) data.table
  # 
  # Args:
  #   array.dt: data.table where the first column contains gene identifiers,
  #             the columns are samples, rows are gene measurements
  #	  zero.to.one: logical - should the data be zero to one transformed?
  # 
  # Returns:
  #   array.dt: log, zero to one transformed if zero.to.one = TRUE, data.table
  require(data.table)
  # check to make sure object passed is a data.table
  array.is.dt <- "data.table" %in% class(array.dt)
  if (!(array.is.dt)) {
    stop("\nInput must be a data.table")
  }
  # if any negative values are found, then inverse log transform and relog 
  # transform using x+1
  any.negative <- any(as.vector(as.matrix(array.dt[, 2:ncol(array.dt),
                                                   with=F])) < 0)
  if (any.negative) {
    #  message("Log transformation array data...\n")
    array.dt <- TDM::inv_log_transform(array.dt)
    array.dt <- TDM::log_transform_p1(array.dt)
  }
  # convert NA to zero
  array.dt <- NAToZero(array.dt)
  #  message("\tZero to one transformation...\n")
  if (zero.to.one) {
  	array.dt <- TDM::zero_to_one_transform(array.dt)
  }
  return(array.dt)
}

LOGSeqOnly <- function(seq.dt, zero.to.one = TRUE){
  # This function takes RNA-seq data in the form of a data.table and returns a
  # log transformed, zero to one transformed (if zero.to.one = TRUE) data.table
  # 
  # Args:
  #   seq.dt: data.table where the first column contains gene identifiers,
  #             the columns are samples, rows are gene measurements
  #	  zero.to.one: logical - should the data be zero to one transformed?
  # 
  # Returns:
  #   log.dt: log, zero to one transformed if zero.to.one = TRUE, data.table
  require(data.table)
  # check to make sure object passed is a data.table
  seq.is.dt <- "data.table" %in% class(seq.dt)
  if (!(seq.is.dt)) {
    stop("\nInput must be a data.table")
  }
  #  message("Log transformation seq data...\n")
  log.dt <- TDM::log_transform_p1(seq.dt)
  # log.dt <- NAToZero(log.dt)
  #  message("\tZero to one transformation...\n")
  if (zero.to.one) {
  	log.dt <- TDM::zero_to_one_transform(log.dt)
  }
  return(log.dt)
}

QNSingleDT <- function(dt, zero.to.one = TRUE){
  # This function takes gene expression data in the form of a data.table 
  # and returns a quantile normalized, zero to one transformed (if zero.to.one
  # = TRUE) data.table
  # 
  # Args:
  #   dt: data.table where the first column contains gene identifiers,
  #             the columns are samples, rows are gene measurements
  #	  zero.to.one: logical - should the data be zero to one transformed?
  # 
  # Returns:
  #   qn.dt: quantile normalized, 
  #          zero to one transformed if zero.to.one = TRUE, data
  require(data.table)
  dt.is.dt <- "data.table" %in% class(dt)
  if (!(dt.is.dt)) {
    stop("\nInput must be a data.table")
  }
  val <- data.frame(dt[, 2:ncol(dt), with = F])
  #  message("Quantile normalization...\n")
  qn <- preprocessCore::normalize.quantiles(data.matrix(val), copy = F)
  qn.dt <- data.table(cbind(as.character(dt[[1]]), qn))
  colnames(qn.dt) <- chartr(".", "-", colnames(dt))
  # qn.dt <- NAToZero(qn.dt)
  #  message("\tZero to one transformation...\n")
  if (zero.to.one) {
  	qn.dt <- TDM::zero_to_one_transform(qn.dt)
  }
  return(qn.dt)
}

NPNSingleDT <- function(dt, zero.to.one = TRUE){
  # This function takes gene expression data in the form of a data.table 
  # and returns a nonparanormal normalized, zero to one transformed
  # (if zero.to.one = TRUE) data.table
  # 
  # Args:
  #   dt: data.table where the first column contains gene identifiers,
  #             the columns are samples, rows are gene measurements
  #	  zero.to.one: logical - should the data be zero to one transformed?
  #
  # Returns:
  #   npn.dt: nonparanormal normalized, 
  #           if zero.to.one = TRUE zero to one transformed, data.table
  require(data.table)
  dt.is.dt <- "data.table" %in% class(dt)
  if (!(dt.is.dt)) {
    stop("\nInput must be a data.table")
  }
  val <- data.frame(dt[, 2:ncol(dt), with = F])
  val.mat <- data.matrix(val)
  npn.mat <- huge.npn(t(val.mat), npn.func = "shrinkage", 
                             npn.thresh = NULL, 
                             verbose = FALSE)
  npn.dt <- data.table(cbind(as.character(dt[[1]]), t(npn.mat)))
  colnames(npn.dt) <- chartr(".", "-", colnames(dt))
  # npn.dt <- NAToZero(npn.dt)
  #  message("\tZero to one transformation...\n")
  if (zero.to.one) {
    npn.dt <- zero_to_one_transform(npn.dt)
  }
  return(npn.dt)
}

ZScoreSingleDT <- function(dt, zero.to.one = TRUE){
  # This function takes gene expression data in the form of a data.table 
  # and returns a z-scored, zero to one transformed data.table
  # 
  # Args:
  #   dt: data.table where the first column contains gene identifiers,
  #             the columns are samples, rows are gene measurements
  #	  zero.to.one: logical - should the data be zero to one transformed?
  #
  # Returns:
  #   z.dt: z-scored, if zero.to.one = TRUE zero to one transformed, data.table
  require(data.table)
  dt.is.dt <- "data.table" %in% class(dt)
  if (!(dt.is.dt)) {
    stop("\nInput must be a data.table")
  }
  val <- data.frame(dt[, 2:ncol(dt), with = F])
  #  message("Z-score...\n")
  z.val <- t(apply(val, 1, function(x) scale(as.numeric(x))))
  z.dt <- data.table(cbind(as.character(dt[[1]]), z.val))  
  colnames(z.dt) <- chartr(".", "-", colnames(dt))
  # z.dt <- NAToZero(z.dt)
  #  message("\tZero to one transformation...\n")
  if (zero.to.one) {
    z.dt <- zero_to_one_transform(z.dt)
  }
  return(z.dt)
}

QNSingleWithRef <- function(ref.dt, targ.dt, zero.to.one = TRUE){
  # This function takes array gene expression data.table as a reference and
  # an RNA-seq expression data.table ('target') and returns the quantile
  # normalized, zero to one transformed (if zero.to.one = TRUE) RNA-seq 
  # expression data.table
  # 
  # Args:
  #   ref.dt: array data.table where the first column contains gene identifiers,
  #           the columns are samples, rows are gene measurements
  #   targ.dt: RNA-seq data.table where the first column contains gene 
  #            identifiers, the columns are samples, rows are gene measurements
  #	  zero.to.one: logical - should the data be zero to one transformed?
  #
  # Returns:
  #   qn.targ: quantile normalized (quantiles from array data), zero to one 
  #            transformed if zero.to.one = TRUE, data.table
  # 
  require(data.table)  
  # Error-handling
  ref.is.dt <- "data.table" %in% class(ref.dt)
  targ.is.dt <- "data.table" %in% class(targ.dt)
  any.not.dt <- !(any(c(ref.is.dt, targ.is.dt)))
  if (any.not.dt) {
    stop("ref.dt and targ.dt must both be data.tables")
  }
  if (!(all(ref.dt[[1]] %in% targ.dt[[1]]))) {
    stop("Gene identifiers in data.tables must match")
  }
  ref.values <- data.frame(ref.dt[, 2:ncol(ref.dt), with = F])
  target.values <- data.frame(targ.dt[, 2:ncol(targ.dt), with = F])
  #  message("Quantile normalization...\n")
  # get target object "reference" for the quantile normalization
  qn.ref <- 
    preprocessCore::normalize.quantiles.determine.target(
                        data.matrix(ref.values), 
                        target.length = nrow(ref.values))
  
  # quantile normalize the data, against reference (array) distribution, using
  # replacement, not averaging
  qn.targ <- 
    preprocessCore::normalize.quantiles.use.target(data.matrix(target.values), 
                                                   qn.ref,
                                                   copy = F)
  qn.targ <- data.table(cbind(as.character(targ.dt[[1]]), qn.targ))
  colnames(qn.targ) <- chartr(".", "-", colnames(qn.targ))  
  #  message("\tZero to one transformation...\n")
  if (zero.to.one) {
    qn.targ <- TDM::zero_to_one_transform(qn.targ)
  }
  return(qn.targ) 
}

TDMSingleWithRef <- function(ref.dt, targ.dt, zero.to.one = TRUE){
  # This function takes array gene expression data.table as a reference and
  # an RNA-seq expression data.table ('target') and returns the TDM
  # normalized, zero to one transformed (if zero.to.one = TRUE) RNA-seq 
  # expression data.table
  # 
  # Args:
  #   ref.dt: array data.table where the first column contains gene identifiers,
  #           the columns are samples, rows are gene measurements
  #   targ.dt: RNA-seq data.table where the first column contains gene 
  #            identifiers, the columns are samples, rows are gene measurements
  #	  zero.to.one: logical - should the data be zero to one transformed?
  #
  # Returns:
  #   tdm.targ: TDM normalized (array data as reference), zero to one 
  #                transformed if zero.to.one = TRUE, data.table
  require(data.table)
  # Error-handling
  ref.is.dt <- "data.table" %in% class(ref.dt)
  targ.is.dt <- "data.table" %in% class(targ.dt)
  any.not.dt <- !(any(c(ref.is.dt, targ.is.dt)))
  if (any.not.dt) {
    stop("ref.dt and targ.dt must both be data.tables")
  }
  if (!(all(ref.dt[[1]] %in% targ.dt[[1]]))) {
    stop("Gene identifiers in data.tables must match")
  }
  #  message("TDM transformation...\n")
  tdm.targ <- TDM::tdm_transform(target_data = targ.dt, 
                                 ref_data = ref.dt, 
                                 negative = FALSE, 
                                 filter_p = FALSE, 
                                 inv_reference = TRUE, 
                                 log_target=TRUE)
  #  message("\tZero to one transformation...\n")
  if (zero.to.one) {
    tdm.targ <- TDM::zero_to_one_transform(tdm.targ)
  }
  return(tdm.targ)
}

SinglePlatformNormalizationWrapper <- function(dt, platform = "array", 
                                               zto = TRUE){
  # This function is a wrapper for processing expression data.tables that
  # contain only one RNA assay platform (array or seq). It returns a list of 
  # normalized data.tables.
  # 
  # Args:
  #   dt: data.table where the first column contains gene identifiers,
  #       the columns are samples, rows are gene measurements from a single
  #       platform
  #   platform: character, array or sequencing data?
  #	  zto: logical - should data be zero to one transformed?
  #   
  # Returns:
  #   norm.list: a list of normalized, 
  #              if zero.to.one = TRUE zero to one transformed, data.tables
  # 
  
  norm.list <- list()
  if (platform == "array") {
    norm.list[["log"]] <- LOGArrayOnly(dt, zto)
    norm.list[["npn"]] <- NPNSingleDT(norm.list$log, zto)
    norm.list[["qn"]] <- QNSingleDT(norm.list$log, zto)
    norm.list[["z"]] <- ZScoreSingleDT(norm.list$log, zto)
  } else if (platform == "seq") {
    norm.list[["log"]] <- LOGSeqOnly(dt, zto)
    norm.list[["npn"]] <- NPNSingleDT(dt, zto)
    norm.list[["qn"]] <- QNSingleDT(dt, zto)
    norm.list[["z"]] <- ZScoreSingleDT(dt, zto)
  } else {
    stop("platform parameter should be set to 'array' or 'seq'")
  }
  return(norm.list)
}


#### cross-platform functions --------------------------------------------------

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

ZScoreProcessing <- function(array.dt, seq.dt, zero.to.one = TRUE){
  # This function takes array and RNA-seq data in the form of data.table 
  # to be 'mixed' (concatenated) and returns one z-scored, 
  # zero to one transformed (if zero.to.one = TRUE) data.table
  # Args:
  #   array.dt: data.table of array data where the first column contains 
  #             gene identifiers, the columns are samples, 
  #             rows are gene measurements
  #   seq.dt:   data.table of RNA-seq data where the first column contains 
  #             gene identifiers, the columns are samples, 
  #             rows are gene measurements
  #   zero.to.one: logical - should data be zero to one transformed?
  # 
  # Returns:
  #   z.dt: z-scored, zero to one transformed if zero.to.one = TRUE, 
  #         data.table that contains both array and RNA-seq samples
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
  #  message("Z-score...\n")
  array.mat <- data.frame(array.dt[, 2:ncol(array.dt), with = F])
  z.array <- t(apply(array.mat, 1, function(x) scale(as.numeric(x))))
  seq.mat <- data.frame(seq.dt[, 2:ncol(seq.dt), with = F])
  z.seq <- t(apply(seq.mat, 1, function(x) scale(as.numeric(x))))
  #  message("\tConcatenation...\n")
  z.dt <- data.table(cbind(array.dt[[1]], z.array, z.seq))
  colnames(z.dt) <- c("gene", 
                      chartr(".", "-", colnames(array.mat)),
                      chartr(".", "-", colnames(seq.mat)))
  # z.dt <- NAToZero(z.dt)
  #  message("\tZero to one transformation...\n")
  if (zero.to.one) {
    z.dt <- TDM::zero_to_one_transform(z.dt)
  }
  return(z.dt)
}

QNProcessing <- function(array.dt, seq.dt, zero.to.one = TRUE){
  # This function takes array and RNA-seq data in the form of data.table 
  # to be 'mixed' (concatenated) and returns one quantile normalized, 
  # zero to one transformed (if zero.to.one = TRUE) data.table. The array data 
  # is used as the target distribution.
  # 
  # Args:
  #   array.dt: data.table of array data where the first column contains 
  #             gene identifiers, the columns are samples, 
  #             rows are gene measurements
  #   seq.dt:   data.table of RNA-seq data where the first column contains 
  #             gene identifiers, the columns are samples, 
  #             rows are gene measurements
  #   zero.to.one: logical - should data be zero to one transformed?
  # 
  # Returns:
  #   qn.cat:  quantile normalized, zero to one transformed if
  #	           zero.to one = TRUE, data.table that contains both array and 
  #            RNA-seq samples
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
  #  message("Quantile normalization...\n")
  qn.targ <- preprocessCore::normalize.quantiles.determine.target(
    data.matrix(ref.values), target.length = nrow(ref.values))
  qn.seq <- preprocessCore::normalize.quantiles.use.target(
    data.matrix(target.values), qn.targ, copy = F)
  qn.seq <- data.table(cbind(seq.dt[[1]], qn.seq))
  colnames(qn.seq) <- chartr(".", "-", colnames(seq.dt))
  # array.dt <- NAToZero(array.dt)
  # qn.seq <- NAToZero(qn.seq)
  #  message("\tZero to one transformation...\n")
  if (zero.to.one) {
    array.dt <- TDM::zero_to_one_transform(array.dt)
    qn.seq <- TDM::zero_to_one_transform(qn.seq)
  }
  #  message("\tConcatenation...\n")
  qn.cat <- data.table(cbind(array.dt, qn.seq[, 2:ncol(qn.seq), with = F]))
  return(qn.cat)
}

NPNProcessing <- function(array.dt, seq.dt, zero.to.one = TRUE){
  # This function takes array and RNA-seq data in the form of data.table 
  # to be 'mixed' (concatenated) and returns one nonparanormal normalized, 
  # zero to one transformed (if zero.to.one = TRUE) data.table
  # 
  # Args:
  #   array.dt: data.table of array data where the first column contains 
  #             gene identifiers, the columns are samples, 
  #             rows are gene measurements
  #   seq.dt:   data.table of RNA-seq data where the first column contains 
  #             gene identifiers, the columns are samples, 
  #             rows are gene measurements
  #   zero.to.one: logical - should data be zero to one transformed?
  #
  # Returns:
  #   npn.cat: NPN normalized, zero to one transformed if zero.to.one = TRUE, 
  #            data.table that contains both array and RNA-seq samples
  #             
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
  npn.array <- huge::huge.npn(t(npn.ref), npn.func = "shrinkage", 
                              npn.thresh = NULL, verbose = FALSE)
  npn.targ <- data.matrix(target.values)
  npn.seq <- huge::huge.npn(t(npn.targ), npn.func = "shrinkage", 
                            npn.thresh = NULL, verbose = FALSE)
  #  message("\tConcatenation...\n")
  npn.cat <- data.table(cbind(array.dt[[1]], t(npn.array), t(npn.seq)))
  colnames(npn.cat) <- c("gene",
                         chartr(".", "-", colnames(ref.values)),
                         chartr(".", "-", colnames(target.values)))
  
  # npn.cat <- NAToZero(npn.cat)
  #  message("\tZero to one transformation...\n")
  if (zero.to.one) {
    npn.cat <- TDM::zero_to_one_transform(npn.cat)
  }
  return(npn.cat)
}

TDMProcessing <- function(array.dt, seq.dt, zero.to.one = TRUE){
  # This function takes array and RNA-seq data in the form of data.table 
  # to be 'mixed' (concatenated) and returns one TDM normalized, 
  # zero to one transformed (if zero.to.one = TRUE) data.table. The array data 
  # is used as the reference distribution.
  # 
  # Args:
  #   array.dt: data.table of array data where the first column contains 
  #             gene identifiers, the columns are samples, 
  #             rows are gene measurements
  #   seq.dt:   data.table of RNA-seq data where the first column contains 
  #             gene identifiers, the columns are samples, 
  #             rows are gene measurements
  #   zero.to.one: logical - should data be zero to one transformed?
  # 
  # Returns:
  #   tdm.cat: TDM normalized, zero to one transformed if zero.to.one = TRUE, 
  #            data.table that contains both array and RNA-seq samples
  #
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
  #  message("TDM transformation...\n")
  tdm.seq <- TDM::tdm_transform(target_data = seq.dt, 
                                ref_data = array.dt, 
                                negative = FALSE, 
                                filter_p = FALSE, 
                                inv_reference = TRUE, 
                                log_target = TRUE)
  # array.dt <- NAToZero(array.dt)
  # tdm.seq <- NAToZero(tdm.seq)
  #  message("\tZero to one transformation...\n")
  if (zero.to.one) {
  	array.dt <- TDM::zero_to_one_transform(array.dt)
  	tdm.seq <- TDM::zero_to_one_transform(tdm.seq)
  }
  #  message("\tConcatenation...\n")
  tdm.cat <- data.table(cbind(array.dt, tdm.seq[, 2:ncol(tdm.seq), with = F]))
  
  return(tdm.cat)
}

LOGProcessing <- function(array.dt, seq.dt, zero.to.one = TRUE){
  # This function takes array and RNA-seq data in the form of data.table 
  # to be 'mixed' (concatenated) and returns one log transformed, 
  # zero to one transformed (if zero.to.one = TRUE) data.table
  # 
  # Args:
  #   array.dt: data.table of array data where the first column contains 
  #             gene identifiers, the columns are samples, 
  #             rows are gene measurements
  #   seq.dt:   data.table of RNA-seq data where the first column contains 
  #             gene identifiers, the columns are samples, 
  #             rows are gene measurements
  #   zero.to.one: logical - should data be zero to one transformed?
  #
  # Returns:
  #   log.cat:  log transformed, zero to one transformed if zero.to.one = TRUE,
  #             data.table that contains both array and RNA-seq samples
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
  #  message("Log transformation seq data...\n")
  log.seq <- TDM::log_transform_p1(seq.dt)
  # array.dt <- NAToZero(array.dt)
  # log.seq <- NAToZero(log.seq)
  #  message("\tZero to one transformation...\n")
  if (zero.to.one) {
    array.dt <- TDM::zero_to_one_transform(array.dt)
    log.seq <- TDM::zero_to_one_transform(log.seq)
  }
  #  message("\tConcatenation...\n")
  log.cat <- data.table(cbind(array.dt, log.seq[, 2:ncol(log.seq), with = F]))
  return(log.cat)
}

NormalizationWrapper <- function(array.dt, seq.dt, 
                                 zto = TRUE){
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
  #	  zto: logical - should data be zero to one transformed?
  #
  # Returns:
  #	  norm.list: a list of normalized, zero to one 'mixed' data.tables 
  #              if zto = TRUE (log transformation, nonparanormal normalized, 
  #              quantile normalized, TDM normalized, z-scored)
  #
  require(data.table)
  
  # if any negative values are found, then inverse log transform and relog 
  # transform using x+1
  any.negative <- any(as.vector(as.matrix(array.dt[, 2:ncol(array.dt),
                                                   with=F])) < 0)
  if (any.negative) {
    #  message("\nLog transformation array data...\n")
    array.dt <- TDM::inv_log_transform(array.dt)
    array.dt <- TDM::log_transform_p1(array.dt)
  }
  # convert NA to zero
  array.dt <- NAToZero(array.dt)  
  norm.list <- list()
  # save array data.table to be used as 'reference' for test data
  norm.list[["raw.array"]] <- array.dt
  # normalization methods, concatenation, zero to one transformation
  norm.list[["log"]] <- LOGProcessing(array.dt, seq.dt, zto)
  norm.list[["npn"]] <- NPNProcessing(array.dt, seq.dt, zto)
  norm.list[["qn"]] <- QNProcessing(array.dt, seq.dt, zto)
  norm.list[["tdm"]] <- TDMProcessing(array.dt, seq.dt, zto)
  norm.list[["z"]] <- ZScoreProcessing(array.dt, seq.dt, zto)
 
  return(norm.list) 
}

GetDataTablesForMixing <- function(array.data, seq.data,
                                   titrate.sample.names){
  # This function takes two full data.tables that contain array and RNA-seq 
  # data,matched samples and a vector of sample names. It returns a list that 
  # contains an array data.table (samples not in the sample names) and an 
  # RNA-seq data.table (samples in the sample names).
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

