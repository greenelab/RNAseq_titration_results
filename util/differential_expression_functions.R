source(file.path("util", "train_test_functions.R"))

GetDesignMatrices <- function(norm.list, sample.df, subtype = "Basal") {
  # This function takes a list of normalized expression data and returns
  # a list of design matrices to be used to detect differentially expressed 
  # genes between BRCA subtype (by default: Basal) and all other subtypes with
  # the limma package 
  # 
  # Args:
  #   norm.list: a nested list of normalized expression data; the top level of 
  #              the list corresponds to % seq included in the expression data
  #   sample.df: a data.frame that maps the sample names to the subtype labels;
  #              'sample' and 'subtype' columns respectively
  #             
  #   
  # Returns:
  #   design.list: a list of design matrices corresponding to each % seq in 
  #                norm.list
  #

  # error-handling
  if (!("sample" %in% colnames(sample.df) 
        & "subtype" %in% colnames(sample.df))) {
    stop("sample.df must have columns named sample and subtype")
  }
  
  # initialize list to hold design matrices
  design.list <- list()
  
  # for each 'titration level' in the normalized data list (norm.list)
  seq.lvls <- names(norm.list)
  for (amt.seq in seq.lvls) {
    
    # get a factor vector of labels, ordered to match the columns of
    # the amt.seq under consideration
    group.factor <- 
      as.character(GetOrderedSubtypeLabels(norm.list[[amt.seq]]$log,
                                           sample.df))
    group.factor[which(group.factor != subtype)] <- "Other"
    group.factor <- as.factor(group.factor)
    
    # build design matrix
    design.mat <- model.matrix(~0 + group.factor)
    colnames(design.mat) <- levels(group.factor)
    rownames(design.mat) <- 
      colnames(norm.list[[amt.seq]]$log)[2:ncol(norm.list[[amt.seq]]$log)]
  
    # add to list
    design.list[[amt.seq]] <- design.mat
    
  }
  
  return(design.list)
  
}

GetFiteBayes <- function(exprs, design.mat) {
  # function for performing differential expression analysis to be passed
  # takes an expression matrix and a design matrix & returns a limma 
  # fit
  # 
  # Args:
  #   exprs: either a matrix, data.table, or EList of expression data
  #   design.mat: a design matrix to be used for differential expression
  #               analysis
  # 
  # Returns:
  #   fit2: differential expression fit results, ready for use with 
  #         limma::topTable
  #         
  # get data.table into correct format for use with limma 
  if ("data.table" %in% class(exprs) & colnames(exprs)[1] == "gene") {
    # error-handling 
    if (!(all.equal(colnames(exprs)[2:ncol(exprs)], rownames(design.mat)))) {
      stop("exprs colnames must match design.mat rownames")
    } 
    gene.ids <- exprs[["gene"]]
    tmp.mat <- as.data.frame(exprs[, 2:ncol(exprs), with = FALSE])   
    rownames(tmp.mat) <- gene.ids
    tmp.mat <- t(apply(tmp.mat , 1, as.numeric))
    colnames(tmp.mat) <- colnames(exprs)[2:ncol(exprs)]
    exprs <- tmp.mat
    rm(tmp.mat)
  } else if (class(exprs) == "EList") {
    # error-handling
    if (!(all.equal(colnames(exprs$E), rownames(design.mat)))) {
      stop("exprs colnames must match design.mat rownames")
    } 
  } else {
    # error-handling 
    if (!(all.equal(colnames(exprs), rownames(design.mat)))) {
      stop("exprs colnames must match design.mat rownames")
    } 
  }
  
  #### differential expression ####
  fit <- limma::lmFit(exprs, design.mat)
  
  # detect subtype to be used for differential expression
  subtype <- colnames(design.mat)[which(colnames(design.mat) != "Other")]
  if (subtype == "Her2") {
    cont.eval <- "Her2vsOther= Her2-Other"
  } else if (subtype == "LumA") {
    cont.eval <- "LumAvsOther= LumA-Other"
  } else if (subtype == "LumB") {
    cont.eval <- "LumBvsOther= LumB-Other"
  } else if (subtype == "Normal") {
    cont.eval <- "NormalvsOther= Normal-Other"
  } else {
    cont.eval <- "BasalvsOther= Basal-Other"
  }
  
  cont.matrix <-
     eval(parse(
       text = paste0(
         'limma::makeContrasts(', cont.eval, ', levels = design.mat)')))
  fit2 <- limma::contrasts.fit(fit, cont.matrix)
  fit2 <- limma::eBayes(fit2)
  return(fit2)
  
}  

GetFiteBayesList <- function(norm.list, design.list) {
  # This function takes a nested list of normalized expression data and 
  # the corresponding list design matrices (output from GetDesignMatrices) and
  # will return differential expression results (eBayes fit from limma). 
  # 
  # We treat the 100% untransformed seq data (RSEM) -- one of the 'silver 
  # standards' -- with limma::voom, which is "intended to process RNA-seq ...
  # data prior to linear modelling in limma" (limma Reference Manual,
  # 28 Jan 2017; see also Law, et al. Genome Biology. 2014.)
  # 
  # Args:
  #   norm.list: a nested list of normalized expression data, the top level 
  #              should be %seq
  #   design.list: a list of design matrices that correspond to each %seq;
  #                output of GetDesignMatrices
  #   
  # Returns:
  #   eBayes.fit.list: nested list of eBayes fits (differential expression 
  #                    results) for each %seq and each normalization method
  #   
  
  # want to assess the following normalization methods 
  norm.methods <- c("log", "npn", "qn", "tdm", "z", "qn-z", "un")
  
  # initialize list for eBayes fits
  eBayes.fit.list <- list()
  # for each %seq data
  for (seq.iter in seq_along(norm.list)) {
    seq.lvl <- names(norm.list)[seq.iter]
    # 100% RNA-seq data has special considerations wrt processing untransformed
    # count data (RSEM)
    if (seq.lvl == "100") {
      nrm.mthds <- head(norm.methods, -1)
      eBayes.fit.list[[seq.lvl]] <-
        lapply(norm.list[[seq.lvl]][names(norm.list[[seq.lvl]]) %in% 
                                      nrm.mthds],
               function(x) GetFiteBayes(exprs = x, 
                                        design.mat = design.list[[seq.lvl]]))
      # special case -- untransformed count data (RSEM)
      exprs <- norm.list[[seq.lvl]]$un
      exprs <- as.data.frame(exprs[, 2:ncol(exprs), with = FALSE])
      rownames(exprs) <- norm.list[[seq.lvl]]$un$gene
      # voom
      voom.counts <- limma::voom(counts = exprs)
      eBayes.fit.list[[seq.lvl]][["un"]] <- GetFiteBayes(voom.counts, 
                                                         design.list[[seq.lvl]])
      
    } else {
      eBayes.fit.list[[seq.lvl]] <-  # for every norm method of interest
        lapply(norm.list[[seq.lvl]][names(norm.list[[seq.lvl]]) %in% 
                                      norm.methods],
               function(x) GetFiteBayes(exprs = x, 
                                        design.mat = design.list[[seq.lvl]]))
    }
  }
  
  # return list of eBayes fits
  return(eBayes.fit.list)
  
}
