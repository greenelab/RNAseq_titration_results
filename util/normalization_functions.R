source(here::here("util", "CrossNorm.R"))

NAToZero <- function(dt, un=0) suppressWarnings(gdata::NAToUnknown(dt, un, force = TRUE))

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
  
  array.dt <- ensure_numeric_gex(array.dt)
  
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
  
  array.dt <- ensure_numeric_gex(array.dt)
  
  #  message("\tZero to one transformation...\n")
  if (zero.to.one) {
    array.dt <- rescale_datatable(array.dt)
  }
  return(data.table(array.dt))
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
  
  seq.dt <- ensure_numeric_gex(seq.dt)
  
  #  message("Log transformation seq data...\n")
  log.dt <- TDM::log_transform_p1(seq.dt)
  
  log.dt <- ensure_numeric_gex(log.dt)
  
  # log.dt <- NAToZero(log.dt)
  #  message("\tZero to one transformation...\n")
  if (zero.to.one) {
    log.dt <- rescale_datatable(log.dt)
  }
  return(data.table(log.dt))
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
  
  dt <- ensure_numeric_gex(dt)
  
  val <- data.frame(dt[, 2:ncol(dt), with = F])
  #  message("Quantile normalization...\n")
  qn <- preprocessCore::normalize.quantiles(data.matrix(val), copy = F)
  qn.dt <- data.table(cbind(as.character(dt[[1]]), qn))
  colnames(qn.dt) <- chartr(".", "-", colnames(dt))
  
  qn.dt <- ensure_numeric_gex(qn.dt)
  
  # qn.dt <- NAToZero(qn.dt)
  #  message("\tZero to one transformation...\n")
  if (zero.to.one) {
    qn.dt <- rescale_datatable(qn.dt)
  }
  return(data.table(qn.dt))
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
  
  dt <- ensure_numeric_gex(dt)
  
  val <- data.frame(dt[, 2:ncol(dt), with = F])
  val.mat <- data.matrix(val)
  npn.mat <- huge.npn(t(val.mat), npn.func = "shrinkage",
                      npn.thresh = NULL,
                      verbose = FALSE)
  npn.dt <- data.table(cbind(as.character(dt[[1]]), t(npn.mat)))
  colnames(npn.dt) <- chartr(".", "-", colnames(dt))
  
  npn.dt <- ensure_numeric_gex(npn.dt)
  
  # npn.dt <- NAToZero(npn.dt)
  #  message("\tZero to one transformation...\n")
  if (zero.to.one) {
    npn.dt <- rescale_datatable(npn.dt)
  }
  return(data.table(npn.dt))
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
  
  dt <- ensure_numeric_gex(dt)
  
  val <- data.frame(dt[, 2:ncol(dt), with = F])
  #  message("Z-score...\n")
  z.val <- t(apply(val, 1, function(x) scale(as.numeric(x))))
  z.dt <- data.table(cbind(as.character(dt[[1]]), z.val))
  colnames(z.dt) <- chartr(".", "-", colnames(dt))
  # z.dt <- NAToZero(z.dt)
  #  message("\tZero to one transformation...\n")
  if (zero.to.one) {
    z.dt <- rescale_datatable(z.dt)
  }
  return(data.table(z.dt))
}

QNZSingleDT <- function(dt, zero.to.one = TRUE){
  # This function takes gene expression data in the form of a data.table
  # and returns a quantile normalized and then z-scored,  zero to one
  # transformed (if zero.to.one = TRUE) data.table
  #
  # Args:
  #   dt: data.table where the first column contains gene identifiers,
  #             the columns are samples, rows are gene measurements
  #	  zero.to.one: logical - should the data be zero to one transformed?
  #
  # Returns:
  #   qn.dt: quantile normalized, z-scored
  #          zero to one transformed if zero.to.one = TRUE, data
  require(data.table)
  dt.is.dt <- "data.table" %in% class(dt)
  if (!(dt.is.dt)) {
    stop("\nInput must be a data.table")
  }
  
  dt <- ensure_numeric_gex(dt)
  
  val <- data.frame(dt[, 2:ncol(dt), with = F])
  #  message("Quantile normalization...\n")
  qn <- preprocessCore::normalize.quantiles(data.matrix(val), copy = F)
  z.qn <- t(apply(qn, 1, function(x) scale(as.numeric(x))))
  z.dt <- data.table(cbind(as.character(dt[[1]]), z.qn))
  colnames(z.dt) <- chartr(".", "-", colnames(dt))
  
  z.dt <- ensure_numeric_gex(z.dt)
  
  # z.dt <- NAToZero(z.dt)
  #  message("\tZero to one transformation...\n")
  if (zero.to.one) {
    z.dt <- rescale_datatable(z.dt)
  }
  return(data.table(z.dt))
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
  
  ref.dt <- ensure_numeric_gex(ref.dt)
  targ.dt <- ensure_numeric_gex(targ.dt)
  
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
  
  qn.targ <- ensure_numeric_gex(qn.targ)
  
  #  message("\tZero to one transformation...\n")
  if (zero.to.one) {
    qn.targ <- rescale_datatable(qn.targ)
  }
  return(data.table(qn.targ))
}

QNZSingleWithRef <- function(ref.dt, targ.dt, zero.to.one = TRUE){
  # This function takes array gene expression data.table as a reference and
  # an RNA-seq expression data.table ('target') and returns the quantile
  # normalized, z-scored, zero to one transformed (if zero.to.one = TRUE) RNA-seq
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
  #   qnz.dt: quantile normalized (quantiles from array data), z-scored,
  #           zero to one transformed if zero.to.one = TRUE, data.table
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
  
  ref.dt <- ensure_numeric_gex(ref.dt)
  targ.dt <- ensure_numeric_gex(targ.dt)
  
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
  
  # z-score the quantile normalized seq values
  qnz.targ <- t(apply(qn.targ, 1, function(x) scale(as.numeric(x))))
  
  #  message("\tConcatenation...\n")
  qnz.dt <- data.table(cbind(targ.dt[[1]], qnz.targ))
  
  colnames(qnz.dt) <- chartr(".", "-", colnames(targ.dt))
  
  qnz.dt <- ensure_numeric_gex(qnz.dt)
  
  #  message("\tZero to one transformation...\n")
  if (zero.to.one) {
    qnz.dt <- rescale_datatable(qnz.dt)
  }
  return(data.table(qnz.dt))
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
  
  ref.dt <- ensure_numeric_gex(ref.dt)
  targ.dt <- ensure_numeric_gex(targ.dt)
  
  #  message("TDM transformation...\n")
  tdm.targ <- TDM::tdm_transform(target_data = targ.dt,
                                 ref_data = ref.dt,
                                 negative = FALSE,
                                 filter_p = FALSE,
                                 inv_reference = TRUE,
                                 log_target=TRUE)
  
  tdm.targ <- ensure_numeric_gex(tdm.targ)
  
  #  message("\tZero to one transformation...\n")
  if (zero.to.one) {
    tdm.targ <- rescale_datatable(tdm.targ)
  }
  return(data.table(tdm.targ))
}

SinglePlatformNormalizationWrapper <- function(dt, platform = "array",
                                               zto = TRUE,
                                               add.untransformed = FALSE,
                                               add.qn.z = FALSE,
                                               add.cn.test = FALSE,
                                               add.seurat.test = FALSE,
                                               norm.list = NULL){
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
  #	  add.untransformed: logical - should untransformed (counts/RSEM) RNA-seq
  #	                     be added to the list?
  #	  add.qn.z: logical - should quantile normalized data that is then z-scored
  #	            be added to the list?
  #   add.cn.test =  logical - should CrossNorm data be added to the list?
  #   add.seurat.test =  logical - should Seurat test data be added to the list? Does NOT get zero-to-one transformed
  #   training.list = list containing normalized training data, needed when add.seurat.test = TRUE
  #
  # Returns:
  #   norm.list: a list of normalized,
  #              if zero.to.one = TRUE zero to one transformed, data.tables
  #

  ### Commented out to allow for log2-scaled array data to be added as
  ### untransformed negative control at 0% RNA-seq
  # error-handling
  #if (platform == "array" & add.untransformed) {
  #  warning("If add.transformed = TRUE, must be RNA-seq data (platform = seq).\n
  #           Setting add.untransformed to FALSE...")
  #  add.untransformed <- FALSE
  #}

  dt <- ensure_numeric_gex(dt)
  
  norm.list <- list()
  if (platform == "array") {
    norm.list[["log"]] <- LOGArrayOnly(dt, zto)
    norm.list[["npn"]] <- NPNSingleDT(norm.list$log, zto)
    norm.list[["qn"]] <- QNSingleDT(norm.list$log, zto)
    norm.list[["z"]] <- ZScoreSingleDT(norm.list$log, zto)
    # should quantile normalized data followed by z-transformation be added?
    if (add.qn.z) {
      norm.list[["qn-z"]] <- QNZSingleDT(norm.list$log, zto)
    }
    # should untransformed (log2 scale, not zero_to_one) array data be added?
    if (add.untransformed){
      norm.list[["un"]] <- UnNoZTOProcessing(array.dt = dt)
    }
  } else if (platform == "seq") {
    norm.list[["log"]] <- LOGSeqOnly(dt, zto)
    norm.list[["npn"]] <- NPNSingleDT(dt, zto)
    norm.list[["qn"]] <- QNSingleDT(dt, zto)
    norm.list[["z"]] <- ZScoreSingleDT(dt, zto)
    # should quantile normalized data followed by z-transformation be added?
    if (add.qn.z) {
      norm.list[["qn-z"]] <- QNZSingleDT(dt, zto)
    }
    # should untransformed RNA-seq data be added?
    if (add.untransformed){
      # by design, untransformed data should not be zero to one transformed,
      # so just add the data.table (dt) that contains RNA-seq data (RSEM)
      # to the list
      norm.list[["un"]] <- dt
    }
  } else {
    
    stop("platform parameter should be set to 'array' or 'seq'")
    
  }
  
  # Whether platform == array or seq, do we add CrossNorm and/or Seurat?
  
  # should CrossNorm be added?
  # Rescale each column, quantile normalize, then rescale each row
  if (add.cn.test){
    norm.list[["cn"]] <- rescale_datatable(dt,
                                           by_column = TRUE) %>%
      QNSingleDT(zero.to.one = zto)
  }
  
  # should Seurat test data be added?
  # do this in parallel for 10-90% RNA-seq using integration from each %
  if (add.seurat.test) {
    
    # parallel backend
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    
    seurat_projection_list <- foreach(i = 2:10) %dopar% { #2:10 corresponds to 10%-90%
      SeuratProjectPCATestData(dt,
                               training.list[[i]][["seurat_model"]],
                               vbose = TRUE)
    }
    
    names(seurat_projection_list) <- names(training.list)[2:10] #10%-90%
    
    # stop parallel backend
    parallel::stopCluster(cl)
    
    # add Seurat RNA-seq test data to list of normalized test data
    norm.list[["seurat"]] <- seurat_projection_list
    rm(seurat_projection_list)
    
  }
  
  return(norm.list)
}

check_all_same <- function(x, my_tolerance = 1e-9){
  # This function returns TRUE if all the elements of the vector are the same
  # within a numerical tolerance levels
  # Thank you: https://stackoverflow.com/a/4752834
  #
  # Args:
  #   x: a numeric vector
  #   my_tolerance: how close must two numbers be for them to be considered equal?
  #
  # Returns:
  #   TRUE or FALSE
  if (is.numeric(x) & is.numeric(my_tolerance)) {
    return(all(abs(max(x) - min(x)) < my_tolerance))  
  } else {
    stop("Vector and tolerance given to check_all_same() must be numeric.")
  }
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
  
  array.dt <- ensure_numeric_gex(array.dt)
  seq.dt <- ensure_numeric_gex(seq.dt)
  
  #  message("Z-score...\n")
  array.mat <- data.frame(array.dt[, 2:ncol(array.dt), with = F])
  z.array <- t(apply(array.mat, 1, function(x) scale(as.numeric(x))))
  seq.mat <- data.frame(seq.dt[, 2:ncol(seq.dt), with = F])
  z.seq <- t(apply(seq.mat, 1, function(x) scale(as.numeric(x))))
  #  message("\tConcatenation...\n")
  z.dt <- data.table(cbind(as.character(array.dt[[1]]), z.array, z.seq))
  colnames(z.dt) <- c("gene",
                      chartr(".", "-", colnames(array.mat)),
                      chartr(".", "-", colnames(seq.mat)))
  
  z.dt <- ensure_numeric_gex(z.dt)
  
  # z.dt <- NAToZero(z.dt)
  #  message("\tZero to one transformation...\n")
  if (zero.to.one) {
    z.dt <- rescale_datatable(z.dt)
  }
  return(data.table(z.dt))
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
  
  array.dt <- ensure_numeric_gex(array.dt)
  seq.dt <- ensure_numeric_gex(seq.dt)
  
  ref.values <- data.frame(array.dt[, 2:ncol(array.dt), with = F])
  target.values <- data.frame(seq.dt[, 2:ncol(seq.dt), with = F])
  #  message("Quantile normalization...\n")
  qn.targ <- preprocessCore::normalize.quantiles.determine.target(
    data.matrix(ref.values), target.length = nrow(ref.values))
  qn.seq <- preprocessCore::normalize.quantiles.use.target(
    data.matrix(target.values), qn.targ, copy = F)
  qn.seq <- data.table(cbind(seq.dt[[1]], qn.seq))
  colnames(qn.seq) <- chartr(".", "-", colnames(seq.dt))

  qn.seq <- ensure_numeric_gex(qn.seq)
  
  # array.dt <- NAToZero(array.dt)
  # qn.seq <- NAToZero(qn.seq)
  #  message("\tZero to one transformation...\n")
  if (zero.to.one) {
    array.dt <- rescale_datatable(array.dt)
    qn.seq <- rescale_datatable(qn.seq)
  }
  #  message("\tConcatenation...\n")
  qn.cat <- data.table(cbind(array.dt, qn.seq[, 2:ncol(qn.seq), with = F]))
  return(data.table(qn.cat))
}

QNZProcessing <- function(array.dt, seq.dt, zero.to.one = TRUE){
  # This function takes array and RNA-seq data in the form of data.table
  # to be 'mixed' (concatenated) and returns one quantile normalized, z-scored
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
  #   z.cat:  quantile normalized, z-scored, zero to one transformed if
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
  
  array.dt <- ensure_numeric_gex(array.dt)
  seq.dt <- ensure_numeric_gex(seq.dt)
  
  ref.values <- data.frame(array.dt[, 2:ncol(array.dt), with = F])
  target.values <- data.frame(seq.dt[, 2:ncol(seq.dt), with = F])
  #  message("Quantile normalization...\n")
  qn.targ <- preprocessCore::normalize.quantiles.determine.target(
    data.matrix(ref.values), target.length = nrow(ref.values))
  qn.seq <- preprocessCore::normalize.quantiles.use.target(
    data.matrix(target.values), qn.targ, copy = F)
  # z-score the array values (ref.values) and the quantile normalized seq values
  z.array <- t(apply(ref.values, 1, function(x) scale(as.numeric(x))))
  z.seq <- t(apply(qn.seq, 1, function(x) scale(as.numeric(x))))
  #  message("\tConcatenation...\n")
  z.dt <- data.table(cbind(as.character(array.dt[[1]]), z.array, z.seq))
  colnames(z.dt) <- c("gene",
                      chartr(".", "-", colnames(ref.values)),
                      chartr(".", "-", colnames(target.values)))
  
  z.dt <- ensure_numeric_gex(z.dt)

  # z.dt <- NAToZero(z.dt)
  #  message("\tZero to one transformation...\n")
  if (zero.to.one) {
    z.dt <- rescale_datatable(z.dt)
  }
  return(data.table(z.dt))
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
  
  array.dt <- ensure_numeric_gex(array.dt)
  seq.dt <- ensure_numeric_gex(seq.dt)
  
  ref.values <- data.frame(array.dt[, 2:ncol(array.dt), with = F])
  target.values <- data.frame(seq.dt[, 2:ncol(seq.dt), with = F])
  npn.ref <- data.matrix(ref.values)
  npn.array <- huge::huge.npn(t(npn.ref), npn.func = "shrinkage",
                              npn.thresh = NULL, verbose = FALSE)
  npn.targ <- data.matrix(target.values)
  npn.seq <- huge::huge.npn(t(npn.targ), npn.func = "shrinkage",
                            npn.thresh = NULL, verbose = FALSE)
  #  message("\tConcatenation...\n")
  npn.cat <- data.table(cbind(as.character(array.dt[[1]]), t(npn.array), t(npn.seq)))
  colnames(npn.cat) <- c("gene",
                         chartr(".", "-", colnames(ref.values)),
                         chartr(".", "-", colnames(target.values)))

  npn.cat <- ensure_numeric_gex(npn.cat)
  
  # npn.cat <- NAToZero(npn.cat)
  #  message("\tZero to one transformation...\n")
  if (zero.to.one) {
    npn.cat <- rescale_datatable(npn.cat)
  }
  return(data.table(npn.cat))
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
  
  array.dt <- ensure_numeric_gex(array.dt)
  seq.dt <- ensure_numeric_gex(seq.dt)
  
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
  
  tdm.seq <- ensure_numeric_gex(tdm.seq)
  
  if (zero.to.one) {
    array.dt <- rescale_datatable(array.dt)
    tdm.seq <- rescale_datatable(tdm.seq)
  }
  #  message("\tConcatenation...\n")
  tdm.cat <- data.table(cbind(array.dt, tdm.seq[, 2:ncol(tdm.seq), with = F]))

  return(data.table(tdm.cat))
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
  
  array.dt <- ensure_numeric_gex(array.dt)
  seq.dt <- ensure_numeric_gex(seq.dt)
  
  #  message("Log transformation seq data...\n")
  log.seq <- TDM::log_transform_p1(seq.dt)
  # array.dt <- NAToZero(array.dt)
  # log.seq <- NAToZero(log.seq)
  
  log.seq <- ensure_numeric_gex(log.seq)
  
  #  message("\tZero to one transformation...\n")
  if (zero.to.one) {
    array.dt <- rescale_datatable(array.dt)
    log.seq <- rescale_datatable(log.seq)
  }
  #  message("\tConcatenation...\n")
  log.cat <- data.table(cbind(array.dt, log.seq[, 2:ncol(log.seq), with = F]))
  return(data.table(log.cat))
}

UnNoZTOProcessing <- function(array.dt = NULL, seq.dt = NULL) {
  # This function takes array data and RNA-seq count data and combines them
  # with no transformation to the RNA-seq data ("untransformed") and no
  # zero to one transformation. It should be regarded as a negative control.
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
  #   dt.cat: data.table that contains concatenated array data and untransformed
  #           RNA-seq data, zero to one transformation is not applied
  #
  
  array.dt.null <- is.null(array.dt)
  seq.dt.null <- is.null(seq.dt)
  if (all(array.dt.null, seq.dt.null)) {
    stop("Cannot have array.dt and seq.dt both NULL in UnNoZTOProcessing()")
  }

  # If the only input is seq data, there is nothing to be done -- just return it
  if (array.dt.null & !seq.dt.null) {
    
    return(data.table(seq.dt)) # don't need to do anything to to seq.dt
    
  } else { # if there is array data, we need to do something to it

    array.dt <- ensure_numeric_gex(array.dt)
        
    # extract the gene vector and array column names
    gene_vector <- array.dt[,1]  
    array_column_names <- colnames(array.dt)
    
    array.dt <- LOGArrayOnly(array.dt, # does nothing if already LOG data
                             zero.to.one = FALSE)
    
    array_matrix <- data.matrix(array.dt[, -1, with = F])
    
    # if there is no seq data, set up the returned object with just array
    if (seq.dt.null) {
    
      un_datatable <- data.table(data.frame(gene_vector, array_matrix))
      colnames(un_datatable) <- array_column_names
      
    } else { # if there is both seq and array data to combine
      
      seq.dt <- ensure_numeric_gex(seq.dt)
      
      # extract seq column names, without the gene column
      seq_column_names <- colnames(seq.dt)[-1]
  
      # combine gene, array, and seq data    
      un_datatable <- data.table(data.frame(gene_vector,
                                            array_matrix,
                                            seq.dt[ , -1, with = F]))
      colnames(un_datatable) <- c(array_column_names, seq_column_names)
      
    }
    
    return(data.table(un_datatable))
    
  }
}

CNProcessing <-  function(array.dt, seq.dt,
                          paired = FALSE,
                          scale_inputs = TRUE,
                          zero.to.one = TRUE) {
  
  # This function takes array and RNA-seq data in the form of data.table
  # to be 'mixed' (concatenated) and applies the CrossNorm normalization method.
  # There are two flavors of CrossNorm: General (default) and Paired
  # If the columns of the two inputs are paired (e.g. from same patient) use paired = TRUE.
  # The two input data sets should be scaled 0-1 on the sample level (column).
  # To be consistent with other methods, we also scale the final matrix gene-wise 0-1.
  # The normalization method used is "quantile normalization".
  #
  # Args:
  #   array.dt: data.table of array data where the first column contains
  #             gene identifiers, the columns are samples,
  #             rows are gene measurements
  #   seq.dt:   data.table of RNA-seq data where the first column contains
  #             gene identifiers, the columns are samples,
  #             rows are gene measurements
  #   paired:   are the columns of array and seq paired (like before/after from same patient?)
  #   scale_inputs: logical - scale columns if inputs are on different scales to start with
  #   zero.to.one: logical - should data be gene-wise zero to one transformed?
  #
  # Returns:
  #   cn.mat:   CrossNorm'd data.table, zero to one transformed if zero.to.one = TRUE,
  #             that contains both array and RNA-seq samples
  
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
  
  gene_names <- array.dt[[1]]
  n_genes <- length(gene_names)
  
  array.dt <- ensure_numeric_gex(array.dt)
  seq.dt <- ensure_numeric_gex(seq.dt)
  
  array.colnames <- names(array.dt)
  seq.colnames <- names(seq.dt)
  
  # if the inputs are on different scales (array and RNA-seq usually are), scale
  # each sample (column) to be on the same 0-1 scale
  if (scale_inputs) {
    array.dt <- rescale_datatable(array.dt, by_column = TRUE)
    seq.dt <- rescale_datatable(seq.dt, by_column = TRUE)
  }

  # The number of array and RNA-seq samples
  n_array <- ncol(array.dt) - 1
  n_seq <- ncol(seq.dt) - 1
  
  combined.dt <- cbind(array.dt[,-1], seq.dt[,-1])
  array_seq_labels <- c(rep(0, n_array), rep(1, n_seq))
  
  # Create the combined data matrix that will be normalized
  if (paired) { # if the inputs are paired, simply stack the two inputs
    
    if (n_array == n_seq) {
      
      cn.dt <- PairedCrossNorm(combined.dt, array_seq_labels)
      
    } else {
      
      stop("Number of array and seq columns must be same for paired CrossNorm")
      
    }
    
  } else { # if the inputs are not paired, then each array sample is matched
    # with each seq sample to create a matrix with n_array*n_seq columns
    
    cn.dt <- GeneralCrossNorm(combined.dt, array_seq_labels)
    
  }
  
  cn.dt <- data.table(cbind(gene_names, cn.dt))
  names(cn.dt) <- c(array.colnames, seq.colnames[-1])
  
  cn.dt <- ensure_numeric_gex(cn.dt)
  
  # rescale each gene to be 0-1
  if (zero.to.one) {
    cn.dt <- rescale_datatable(cn.dt)
  }
  
  return(cn.dt)
  
}


SeuratIntegration <-  function(array.dt, seq.dt,
                               n_dims = 50,
                               vbose = TRUE) {
  
  # This function takes array and RNA-seq data in the form of data.tables
  # to be normalized and concatenated through Seurat integration.
  #
  # Args:
  #   array.dt: data.table of array data where the first column contains
  #             gene identifiers, the columns are samples,
  #             rows are gene measurements
  #   seq.dt:   data.table of RNA-seq data where the first column contains
  #             gene identifiers, the columns are samples,
  #             rows are gene measurements
  #   n_dims:   number of dimensions for reduction (default: 50)
  #   vbose:    use verbose screen output (default: TRUE, may give error if FALSE)
  #
  # Returns:
  #   array_seq.integrated: an integrated Seurat object including array and seq
  
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
  
  array.dt <- ensure_numeric_gex(array.dt)
  seq.dt <- ensure_numeric_gex(seq.dt)
  
  # List of seurat objects to be integrated
  array_seq_list <- list(array = array.dt[,-1],
                         seq = seq.dt[,-1])
  
  array_seq_list <- lapply(X = array_seq_list,
                           FUN = Seurat::CreateSeuratObject)
  
  # SCTransform
  array_seq_list <- lapply(X = array_seq_list,
                           FUN = Seurat::SCTransform,
                           verbose = vbose)
  
  features <- Seurat::SelectIntegrationFeatures(object.list = array_seq_list,
                                                verbose = vbose)
  
  array_seq_list <- Seurat::PrepSCTIntegration(object.list = array_seq_list,
                                               anchor.features = features,
                                               verbose = vbose)
    
  # integration
  array_seq.anchors <- Seurat::FindIntegrationAnchors(object.list = array_seq_list,
                                                      normalization.method = "SCT",
                                                      anchor.features = features,
                                                      dims = 1:n_dims,
                                                      verbose = vbose)
  
  array_seq.integrated <- Seurat::IntegrateData(anchorset = array_seq.anchors,
                                                normalization.method = "SCT",
                                                dims = 1:n_dims,
                                                verbose = vbose)
  
  # PCA
  array_seq.integrated <- Seurat::RunPCA(array_seq.integrated,
                                         npcs = n_dims,
                                         nfeatures.print = n_dims,
                                         verbose = vbose)
  
  return(array_seq.integrated)
  
}

SeuratPCATrainingData <- function(integrated_seurat_object) {
  
  # This function returns the PCA reduction of training data from the integrated Seurat object
  #
  # Args:
  #   integrated_seurat_object:  output from SeuratIntegration, an integrated
  #             Seurat object with both array and RNA-seq data combined
  #
  # Returns:
  #   properly formatted data.table (with dimension names as column 1, sample names as column names)
  
  reduced.dt <- integrated_seurat_object@reductions$pca@cell.embeddings %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("dimension") %>%
    data.table::as.data.table() %>%
    ensure_numeric_gex()
  
  return(reduced.dt)
  
}

SeuratProjectPCATestData <- function(test_data.dt,
                                     integrated_seurat_object,
                                     vbose = TRUE) {
  
  # This function transforms test data (SCTransform), then projects it to the
  # reduced dimension PCA space of the integrated_seurat_object.
  #
  # Args:
  #   test_data.dt: a data.table of test data (column 1 = gene)
  #   integrated_seurat_object:  output from SeuratIntegration, an integrated
  #             Seurat object with both array and RNA-seq data combined
  #   vbose:    use verbose screen output (default: TRUE, may give error if FALSE)
  #
  # Returns:
  #   properly formatted data.table (with dimension names as column 1, sample names as column names)
  
  require(data.table)
  # Error-handling
  test_data.is.dt <- "data.table" %in% class(test_data.dt)
  if (!test_data.is.dt) {
    stop("test_data must be data.table")
  }
  
  test_data.dt <- ensure_numeric_gex(test_data.dt)

  # Create query single cell object based on test data
  query.sct <- Seurat::CreateSeuratObject(test_data.dt[,-1]) %>%
    Seurat::SCTransform(verbose = vbose)
  
  # Detect transfer anchors between training and test
  ref_dims <- dim(integrated_seurat_object@reductions$pca@cell.embeddings)[2]
  anchors <- Seurat::FindTransferAnchors(reference = integrated_seurat_object,
                                         query = query.sct,
                                         dims = 1:ref_dims,
                                         reference.reduction = "pca",
                                         verbose = vbose)
  
  # Add PCA mapping to query 
  query.sct <- Seurat::MapQuery(anchorset = anchors,
                                reference = integrated_seurat_object,
                                query = query.sct,
                                reference.reduction = "pca",
                                verbose = vbose)
  
  # Convert samples x PCs to PCs by sample and name rows same as training
  test_reduced.dt <- query.sct@reductions$ref.pca@cell.embeddings %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("dimension") %>%
    mutate(dimension = stringr::str_replace(dimension, "refpca", "PC")) %>%
    data.table::as.data.table() %>%
    ensure_numeric_gex()
  
  return(test_reduced.dt)
  
}

NormalizationWrapper <- function(array.dt, seq.dt,
                                 zto = TRUE,
                                 add.untransformed = FALSE,
                                 add.qn.z = FALSE,
                                 add.cn = FALSE,
                                 add.seurat.training = FALSE){
  
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
  #	  add.untransformed: logical - should untransformed (counts/RSEM) RNA-seq
  #	                     be concatenated to array data and added to the list?
  #	  add.qn.z: logical - should quantile normalized data that is then z-scored
  #	            be added to the list?
  #   add.cn =  logical - should CrossNorm data be added to the list?
  #   add.seurat.training =  logical - should Seurat training data be added to the list? Does NOT get zero-to-one transformed
  #
  # Returns:
  #	  norm.list: a list of normalized, zero to one 'mixed' data.tables
  #              if zto = TRUE (log transformation, nonparanormal normalized,
  #              quantile normalized, TDM normalized, z-scored, quantile
  #              normalized + z-transformed [if add.qn.z = TRUE], untransformed
  #              [if add.untransformed = TRUE])
  #
  require(data.table)

  # if zero to one transformation is to be performed AND untransformed RNA-seq
  # data is to be added, warn the user that there is no zero to one
  # transformation in the untransformed step
  if (zto & add.untransformed) {
    warning("Zero to one transformation will be performed for most normalization
            methods. Untransformed RNA-seq data step does not include zero to
            one transformation.")
  }

  array.dt <- ensure_numeric_gex(array.dt)
  seq.dt <- ensure_numeric_gex(seq.dt)
  
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
  
  array.dt <- ensure_numeric_gex(array.dt)
  
  norm.list <- list()
  # save array data.table to be used as 'reference' for test data
  norm.list[["raw.array"]] <- array.dt
  # normalization methods, concatenation, zero to one transformation
  norm.list[["log"]] <- LOGProcessing(array.dt, seq.dt, zto)
  norm.list[["npn"]] <- NPNProcessing(array.dt, seq.dt, zto)
  norm.list[["qn"]] <- QNProcessing(array.dt, seq.dt, zto)
  norm.list[["tdm"]] <- TDMProcessing(array.dt, seq.dt, zto)
  norm.list[["z"]] <- ZScoreProcessing(array.dt, seq.dt, zto)
  # should quantile normalized data that then z-transformed be added?
  if (add.qn.z) {
    norm.list[["qn-z"]] <- QNZProcessing(array.dt, seq.dt, zto)
  }
  # should untransformed data be added?
  if (add.untransformed) {
    norm.list[["un"]] <- UnNoZTOProcessing(array.dt, seq.dt)
  }
  # should CrossNorm be added?
  if (add.cn) {
    norm.list[["cn"]] <- CNProcessing(array.dt = array.dt,
                                      seq.dt = seq.dt)
  }

  # should Seurat training data be added?
  if (add.seurat.training) {
    
    norm.list[["seurat_model"]] <- SeuratIntegration(array.dt,
                                                     seq.dt,
                                                     n_dims = 50,
                                                     vbose = TRUE)
    
    
    norm.list[["seurat"]] <- SeuratPCATrainingData(norm.list[["seurat_model"]])
    
  }
  
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
  
  array.data <- ensure_numeric_gex(array.data)
  seq.data <- ensure_numeric_gex(seq.data)
  
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
    
    array.dt <- ensure_numeric_gex(array.dt)
    seq.dt <- ensure_numeric_gex(seq.dt)
    
    mix.dt.list <- list("array" = array.dt, "seq" = seq.dt)
  } else if (length(titrate.sample.names) == (ncol(seq.data) - 1)) {
    mix.dt.list <- list("seq" = data.table(seq.data))
  } else {
    mix.dt.list <- list("array" = data.table(array.data))
  }
}

rescale_01 <- function(data_vector){
  # rescale values in a vector to [0,1]
  # Inputs: vector of numeric values
  # Returns: rescaled vector

  # if all the values are the same, return 0 vector
  if (check_all_same(data_vector)) {
    return(rep(0, length(data_vector)))
  } else {
    min_value <- min(data_vector)
    max_value <- max(data_vector)
    rescaled_values <- (data_vector - min_value)/(max_value - min_value)
    return(rescaled_values)
  }
}

rescale_datatable <- function(data_table, by_column = FALSE){
  # rescale each row (or column) of a data table to [0,1]
  # applies rescale_01() to each row (or column)
  # Inputs: gene expression data table
  #   first column of input is genes
  #   remaining columns are expression values
  #   by_column = FALSE rescale each row, if TRUE rescale each column
  # Returns: scaled gene expression data table
  
  data_table <- ensure_numeric_gex(data_table)
  
  data_matrix = data.matrix(data_table[, -1, with = F])
  
  # Rescale each row or column [0,1]
  if (by_column) {
    rescaled_data_matrix = apply(data_matrix, 2, rescale_01)
  } else {
    rescaled_data_matrix = t(apply(data_matrix, 1, rescale_01))
  }
  
  # Include gene symbols in result
  result = data.table(data.frame(data_table[,1], rescaled_data_matrix))
  colnames(result) <- colnames(data_table)
  
  result <- ensure_numeric_gex(result)
  
  return(data.table(result))
  
}

ensure_numeric_gex <- function(input_data) {
  # Numeric version of data.table
  #
  # Ensure gene expression values are numeric in a given data.table
  #
  # input_data: a data.table with gene in the first column and gene
  # expression in the remaining columns
  #
  # returns a data.table with numeric values for the gene expression columns

  if ("data.table" %in% class(input_data)) {

    # save gene names as character vector
    gene_names <- as.character(input_data[[1]])

    # force gene expression values to be numeric
    gex_values <- t(apply(input_data[,-1], 1, as.numeric))

    # create data frame of gene names and gene expression values
    # set column names to be same as input data
    return_df <- data.frame(gene_names, gex_values)
    names(return_df) <- names(input_data)

    # return as data.table
    return(data.table(return_df))

  } else {

    stop("\nInput must be a data.table")

  }
}
