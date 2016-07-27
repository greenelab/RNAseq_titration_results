NAToZero <- function(dt, un=0) gdata::NAToUnknown(dt, un)

LOGArrayOnly <- function(array.dt){
  # This function takes array data in the form of a data.table and returns a
  # log transformed, zero to one transformed data.table
  # 
  # Args:
  #   array.dt: data.table where the first column contains gene identifiers,
  #             the columns are samples, rows are gene measurements
  # 
  # Returns:
  #   zto.array.dt: log, zero to one transformed data.table
  require(TDM)
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
    cat("Log transformation array data...\n")
    array.dt <- inv_log_transform(array.dt)
    array.dt <- log_transform_p1(array.dt)
  }
  # convert NA to zero
  array.dt <- NAToZero(array.dt)
  cat("\tZero to one transformation...\n")
  zto.array.dt <- zero_to_one_transform(array.dt)
  return(zto.array.dt)
}

LOGSeqOnly <- function(seq.dt){
  # This function takes RNA-seq data in the form of a data.table and returns a
  # log transformed, zero to one transformed data.table
  # 
  # Args:
  #   seq.dt: data.table where the first column contains gene identifiers,
  #             the columns are samples, rows are gene measurements
  # 
  # Returns:
  #   zto.log.dt: log, zero to one transformed data.table
  require(TDM)
  require(data.table)
  # check to make sure object passed is a data.table
  seq.is.dt <- "data.table" %in% class(seq.dt)
  if (!(seq.is.dt)) {
    stop("\nInput must be a data.table")
  }
  cat("Log transformation seq data...\n")
  log.dt <- log_transform_p1(seq.dt)
  # log.dt <- NAToZero(log.dt)
  cat("\tZero to one transformation...\n")
  zto.log.dt <- zero_to_one_transform(log.dt)
  return(zto.log.dt)
}

QNSingleDT <- function(dt){
  # This function takes gene expression data in the form of a data.table 
  # and returns a quantile normalized, zero to one transformed data.table
  # 
  # Args:
  #   dt: data.table where the first column contains gene identifiers,
  #             the columns are samples, rows are gene measurements
  # 
  # Returns:
  #   zto.qn.dt: quantile normalized, zero to one transformed data.table
  require(TDM)
  require(preprocessCore)
  require(data.table)
  dt.is.dt <- "data.table" %in% class(dt)
  if (!(dt.is.dt)) {
    stop("\nInput must be a data.table")
  }
  val <- data.frame(dt[, 2:ncol(dt), with = F])
  cat("Quantile normalization...\n")
  qn <- normalize.quantiles(data.matrix(val), copy = F)
  qn.dt <- data.table(cbind(as.character(dt[[1]]), qn))
  colnames(qn.dt) <- chartr(".", "-", colnames(dt))
  # qn.dt <- NAToZero(qn.dt)
  cat("\tZero to one transformation...\n")
  zto.qn.dt <- zero_to_one_transform(qn.dt)
  return(zto.qn.dt)
}


NPNSingleDT <- function(dt){
  # This function takes gene expression data in the form of a data.table 
  # and returns a nonparanormal normalized, zero to one transformed data.table
  # 
  # Args:
  #   dt: data.table where the first column contains gene identifiers,
  #             the columns are samples, rows are gene measurements
  # 
  # Returns:
  #   zto.npn: nonparanormal normalized, zero to one transformed data.table
  require(huge)
  require(TDM)
  require(data.table)
  dt.is.dt <- "data.table" %in% class(dt)
  if (!(dt.is.dt)) {
    stop("\nInput must be a data.table")
  }
  val <- data.frame(dt[, 2:ncol(dt), with = F])
  val.mat <- data.matrix(val)
  npn.mat <- huge.npn(t(val.mat), npn.func = "shrinkage", 
                             npn.thresh = NULL, 
                             verbose = TRUE)
  npn.dt <- data.table(cbind(as.character(dt[[1]]), t(npn.mat)))
  colnames(npn.dt) <- chartr(".", "-", colnames(dt))
  # npn.dt <- NAToZero(npn.dt)
  cat("\tZero to one transformation...\n")
  zto.npn <- zero_to_one_transform(npn.dt)
  return(zto.npn)
}

ZScoreSingleDT <- function(dt){
  # This function takes gene expression data in the form of a data.table 
  # and returns a z-scored, zero to one transformed data.table
  # 
  # Args:
  #   dt: data.table where the first column contains gene identifiers,
  #             the columns are samples, rows are gene measurements
  # 
  # Returns:
  #   zto.z.dt: z-scored, zero to one transformed data.table
  require(TDM)
  require(data.table)
  dt.is.dt <- "data.table" %in% class(dt)
  if (!(dt.is.dt)) {
    stop("\nInput must be a data.table")
  }
  val <- data.frame(dt[, 2:ncol(dt), with = F])
  cat("Z-score...\n")
  z.val <- t(apply(val, 1, function(x) scale(as.numeric(x))))
  z.dt <- data.table(cbind(as.character(dt[[1]]), z.val))  
  colnames(z.dt) <- chartr(".", "-", colnames(dt))
  # z.dt <- NAToZero(z.dt)
  cat("\tZero to one transformation...\n")
  zto.z.dt <- zero_to_one_transform(z.dt)
  return(zto.z.dt)
}

QNSingleWithRef <- function(ref.dt, targ.dt){
  # This function takes array gene expression data.table as a reference and
  # an RNA-seq expression data.table ('target') and returns the quantile
  # normalized RNA-seq expression data.table
  # 
  # Args:
  #   ref.dt: array data.table where the first column contains gene identifiers,
  #           the columns are samples, rows are gene measurements
  #   targ.dt: RNA-seq data.table where the first column contains gene 
  #            identifiers, the columns are samples, rows are gene measurements
  # 
  # Returns:
  #   zto.qn.targ: quantile normalized (quantiles from array data), zero to one 
  #                transformed data.table
  # 
    
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
  require(preprocessCore)
  require(TDM)
  require(data.table)
  ref.values <- data.frame(ref.dt[, 2:ncol(ref.dt), with = F])
  target.values <- data.frame(targ.dt[, 2:ncol(targ.dt), with = F])
  cat("Quantile normalization...\n")
  # get target object "reference" for the quantile normalization
  qn.ref <- normalize.quantiles.determine.target(
    data.matrix(ref.values), 
    target.length = nrow(ref.values))
  
  # quantile normalize the data, against reference (array) distribution, using
  # replacement, not averaging
  qn.targ <- normalize.quantiles.use.target(data.matrix(target.values), qn.ref,
                                           copy = F)
  qn.targ <- data.table(cbind(as.character(targ.dt[[1]]), qn.targ))
  colnames(qn.targ) <- chartr(".", "-", colnames(qn.targ))  
  cat("\tZero to one transformation...\n")
  zto.qn.targ <- zero_to_one_transform(qn.targ)
  return(zto.qn.targ) 
}

TDMSingleWithRef <- function(ref.dt, targ.dt){
  # This function takes array gene expression data.table as a reference and
  # an RNA-seq expression data.table ('target') and returns the TDM
  # normalized RNA-seq expression data.table
  # 
  # Args:
  #   ref.dt: array data.table where the first column contains gene identifiers,
  #           the columns are samples, rows are gene measurements
  #   targ.dt: RNA-seq data.table where the first column contains gene 
  #            identifiers, the columns are samples, rows are gene measurements
  # 
  # Returns:
  #   zto.tdm.targ: TDM normalized (array data as reference), zero to one 
  #                transformed data.table

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
  
  require(TDM)
  require(data.table)
  cat("TDM transformation...\n")
  tdm.targ <- tdm_transform(target_data = targ.dt, 
                           ref_data = ref.dt, 
                           negative = FALSE, 
                           filter_p = FALSE, 
                           inv_reference = TRUE, 
                           log_target=TRUE)
  cat("\tZero to one transformation...\n")
  zto.tdm.targ <- zero_to_one_transform(tdm.targ)
  return(zto.tdm.targ)
}

SinglePlatformNormalizationWrapper <- function(dt, platform = "array"){
  # This function is a wrapper for processing expression data.tables that
  # contain only one RNA assay platform (array or seq). It returns a list of 
  # normalized data.tables.
  # 
  # Args:
  #   dt: data.table where the first column contains gene identifiers,
  #       the columns are samples, rows are gene measurements from a single
  #       platform
  #   platform: character, array or sequencing data?
  #   
  # Returns:
  #   norm.list: a list of normalized, zero to one transformed data.tables
  # 
  norm.list <- list()
  if (platform == "array") {
    norm.list[["log"]] <- LOGArrayOnly(dt)
    norm.list[["npn"]] <- NPNSingleDT(norm.list$log)
    norm.list[["qn"]] <- QNSingleDT(norm.list$log)
    norm.list[["z"]] <- ZScoreSingleDT(norm.list$log)
  } else if (platform == "seq") {
    norm.list[["log"]] <- LOGSeqOnly(dt)
    norm.list[["npn"]] <- NPNSingleDT(dt)
    norm.list[["qn"]] <- QNSingleDT(dt)
    norm.list[["z"]] <- ZScoreSingleDT(dt)
  } else {
    stop("platform parameter should be set to 'array' or 'seq'")
  }
  return(norm.list)
}
