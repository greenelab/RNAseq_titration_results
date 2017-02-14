source(file.path("util", "train_test_functions.R"))

GetDesignMat <- function(norm.dt, sample.df, subtype) {
  # This function takes a data.table of expression data (1st column is gene ids; 
  # rows are genes, columns are samples), a data.frame mapping sample names to
  # subtype labels, and returns a design matrix specifying for use with limma
  # comparing the user-specified subtype to all other samples
  # 
  # Args:
  #   norm.dt: a data.table of normalized expression data, 1st column is
  #            gene identifiers. rows are genes; samples are columns
  #   sample.df: a data.frame that maps the sample names to the subtype labels;
  #              'sample' and 'subtype' columns respectively
  #   subtype: which subtype should be compared to all others?
  #   
  # Returns:
  #   design.mat: a design matrix (2-group comparison) for use with limma
  #
  
  # error-handling
  if (!("sample" %in% colnames(sample.df) 
        & "subtype" %in% colnames(sample.df))) {
    stop("sample.df must have columns named sample and subtype")
  }
  
  # get a factor vector of labels, ordered to match the columns of
  # norm.dt
  group.factor <- 
    as.character(GetOrderedSubtypeLabels(norm.dt,
                                         sample.df))
  group.factor[which(group.factor != subtype)] <- "Other"
  group.factor <- as.factor(group.factor)
  
  # build design matrix
  design.mat <- model.matrix(~0 + group.factor)
  colnames(design.mat) <- levels(group.factor)
  rownames(design.mat) <- colnames(norm.dt)[2:ncol(norm.dt)]
  
  return(design.mat)
}

GetDesignMatrixList <- function(norm.list, sample.df, subtype = "Basal") {
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
  #   subtype: which subtype should be compared to all others? default is Basal
  #             
  #   
  # Returns:
  #   design.list: a list of design matrices corresponding to each % seq in 
  #                norm.list
  #
  
  # initialize list to hold design matrices
  design.list <- list()
  
  # for each 'titration level' in the normalized data list (norm.list)
  seq.lvls <- names(norm.list)
  for (amt.seq in seq.lvls) {
    design.list[[amt.seq]] <- GetDesignMat(norm.dt = norm.list[[amt.seq]]$log,
                                           sample.df = sample.df,
                                           subtype = subtype)
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

#### plotting functions --------------------------------------------------------

GetAllGenesTopTable <- function(fit.object, adjust = "BH") {
  # Get all genes into a top table -- reports the t-statistic, p-value, and 
  # FDR (adj.P.Val) if adjust = "BH", see limma::topTable adjust.methods for 
  # more information
  # 
  # Args:
  #   fit.object: a limma fit object (post eBayes) to extract the topTable from
  #   adjust: character; how should correction for multiple hypothesis testing 
  #           be performed? FDR ("BH") by default
  #
  # Returns:
  #   top.table: the corresponding topTable that includes all genes tested
  #   
  top.table <- limma::topTable(fit.object, number = dim(fit.object$t)[1],
                               adjust.method = adjust)
  
  return(top.table)
}

GetGeneSetJaccard <- function(silver.set, top.table, 
                              cutoff = 0.05) {
  # Given a "silver standard" of differentially expressed genes (a character
  # vector), find the Jaccard similarity between the silver standard and another
  # set of differentially expressed genes (DEGs; from limma::topTable output) 
  # using a user-specified adjusted p-value threshold
  # 
  # Args:
  #   silver.set: a character vector of the "silver standard" genes 
  #   top.table: output of limma::topTable (typically from GetAllGenesTopTable)
  #              for the other differential expression experiment being
  #              considered
  #   cutoff: the adjusted p-value threshold - all genes with an adjusted p
  #           below this threshold will be considered differentially expressed
  #   
  # Returns:
  #   jacc: Jaccard similarity of silver standard and other set of DEGs
  
  # get differentially expressed gene set from experiment topTable 
  tt.set <- rownames(top.table)[which(top.table$adj.P.Val < cutoff)]
  # find Jaccard similarity between experiment and silver standard
  jacc <-
    length(intersect(silver.set, tt.set)) / length(union(silver.set, tt.set))
  # return Jaccard similarity
  return(jacc)
}


PlotProportionDE <- function(fit.list, adjust.method = "BH", cutoff = 0.05) {
  # from a list of limma fits, plot the proportion of genes that are 
  # considered differentially expressed under a user-specified adjusted p 
  # cutoff (FDR < 5% by default); return topTable(s) as well
  #
  # Args:
  #   fit.list: a nested list of limma eBayes fits
  #   adjust.method: how should the p-values be adjusted for multiple 
  #                  hypothesis testing?
  #
  # Returns
  #   A list containing:
  #     top.table.list: the limma::topTable output (includes all genes)
  #                     for fits
  #     p: the bar plot of the proportion of genes that are differentially 
  #        expressed
  #        
  
  # colorblind friendly palette
  cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                 "#0072B2", "#D55E00", "#CC79A7")
  
  # get topTables
  top.table.list <- 
    lapply(fit.list,  # for each level of % seq
           function(x)   
             lapply(x, # for each normalization method
                    function(y) GetAllGenesTopTable(y, adjust = adjust.method)))
  
  # how many genes are differentially expressed @ cutoff
  deg.count.count.list <- 
    lapply(top.table.list,  # for each level of % seq
           function(x) lapply(x,  # for each normalization method  
                              function(y) sum(y$adj.P.Val < cutoff)))
  # get DEG counts as data.frame
  deg.count.df <- reshape2::melt(deg.count.count.list)
  
  # get proportion differentially expressed genes
  deg.count.df$value <- deg.count.df$value / dim(top.table.list[[1]][[1]])[1]
  colnames(deg.count.df) <- c("perc.ltdeg.count", "normalization", "perc.seq")
  
  # order % seq so bar plot displays 0-100
  deg.count.df$perc.seq <- factor(deg.count.df$perc.seq, 
                                  levels = seq(0, 100, 10))
  
  # capitalize normalization methods for display 
  deg.count.df$normalization <- as.factor(toupper(deg.count.df$normalization))
  
  # bar plots of proportion of genes that are differentially expressed
  # in each experiment
  p <- ggplot(deg.count.df, aes(perc.seq, perc.ltdeg.count, 
                                fill = normalization)) + 
    geom_bar(stat = "identity", position = "dodge", colour = "black") +
    theme_classic() +
    xlab("% RNA-seq samples") +
    ylab("proportion of genes FDR < 5%") +
    scale_fill_manual(values = cbPalette)
  
  return(list("top.table.list" = top.table.list, "plot" = p))
}

PlotSilverStandardJaccard <- function(top.table.list, title,
                                      cutoff = 0.05){
  # Given a list of top tables, plot the Jaccard similarity between the 
  # microarray and RNA-seq "silver standards" and all other experiments
  #
  # Args:
  #   top.table.list: list of limma::topTable objects from PlotProportionDE
  #   title: a character vector to be used for the title of the plot
  #   cutoff: adjusted p-value threshold to be used
  #   
  # Returns:
  #   p: Jaccard similarity line plot
  # 
  
  # colorblind friendly palette
  cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                 "#0072B2", "#D55E00", "#CC79A7")
  
  ### "silver standards" ###
  
  # 100% RNA-seq data RSEM using limma::voom processing step
  seq.truth.set <- 
    rownames(top.table.list$`100`$un)[which(top.table.list$`100`$un$adj.P.Val < 
                                              cutoff)]
  # LOG 100% array data
  array.truth.set <- 
    rownames(top.table.list$`0`$log)[which(top.table.list$`0`$log$adj.P.Val < 
                                             cutoff)]
  
  # how similiar are DEG results to the RNA-seq silver standard? as data.frame
  seq.jacc.list <- 
    lapply(top.table.list, 
           function(x) lapply(x,
                              function(y) GetGeneSetJaccard(seq.truth.set, y)))
  seq.jacc.df <- reshape2::melt(seq.jacc.list)
  
  # how similiar are DEG results to the microarray silver standard? as data.frame
  array.jacc.list <- 
    lapply(top.table.list, 
           function(x) lapply(x,
                              function(y) GetGeneSetJaccard(array.truth.set, 
                                                            y)))
  array.jacc.df <- reshape2::melt(array.jacc.list)
  
  # combine seq and array similarity results 
  array.jacc.df <- cbind(array.jacc.df, rep("Microarray", nrow(array.jacc.df)))
  seq.jacc.df <- cbind(seq.jacc.df, rep("RNA-seq", nrow(seq.jacc.df)))
  colnames(seq.jacc.df) <- colnames(array.jacc.df) <- c("jaccard", 
                                                        "normalization",
                                                        "perc.seq",
                                                        "platform")
  mstr.df <- rbind(array.jacc.df, seq.jacc.df)
  
  # order % seq so plot displays 0-100
  mstr.df$perc.seq <- factor(mstr.df$perc.seq, levels = seq(0, 100, 10))
  
  # capitalize normalization methods for display 
  mstr.df$normalization <- as.factor(toupper(mstr.df$normalization))
  
  # line plot
  p <- ggplot(mstr.df, aes(perc.seq, jaccard, 
                           color = platform, fill = platform)) +
    facet_wrap(~normalization, ncol = 4) +
    geom_line(aes(group = platform), position = position_dodge(0.3)) +
    geom_point(aes(group = platform), position = position_dodge(0.3)) +
    theme_bw() +
    scale_colour_manual(values = cbPalette[c(2, 3)]) +
    ggtitle(title) +
    xlab("% RNA-seq") +
    ylab("Jaccard similarity") + 
    theme(axis.text.x=element_text(angle = 45, vjust = 0.5))
  return(p)
  
}

#### small n functions ---------------------------------------------------------

GetSmallNSilverStandardJaccard <- function(top.table.list, cutoff = 0.05){
  # This function takes a list of limma::topTable output, derives "silver
  # standards" from 100% array and 100% seq data, and finds the Jaccard 
  # similarity between the standards and 50-50 experiment differential 
  # expression results (quantile norm and z-score)
  # 
  # Args:
  #   top.table.list: a list of limma::topTable output with all genes -
  #                   output of GetAllGenesTopTable
  #   cutoff: FDR cutoff, defaults to FDR < 5%
  #   
  # Returns:
  #   jacc.df: a data.frame of the Jaccard results with columns corresponding to
  #            the Jaccard value, the silver standard used for the
  #            comparison ("platform"), the normalization method, and the
  #            n used (number of samples; "no.samples")
  #
  
  # initialize list to hold Jaccard similarity
  jacc.list <- list()
  
  # for each n (number of samples)
  for (smpl.no.iter in seq_along(top.table.list)) {
    
    current.smpl.tt <- top.table.list[[smpl.no.iter]]
    current.n <- names(top.table.list)[smpl.no.iter]
    
    array.top <- current.smpl.tt$log
    array.standard <- rownames(array.top)[which(array.top$adj.P.Val < cutoff)]
    seq.top <- current.smpl.tt$un
    seq.standard <- rownames(seq.top)[which(seq.top$adj.P.Val < cutoff)]
    
    jacc.list[[current.n]]$qn$array <-
      GetGeneSetJaccard(array.standard, 
                        top.table = current.smpl.tt$qn,
                        cutoff)
    jacc.list[[current.n]]$qn$seq <-
      GetGeneSetJaccard(seq.standard, 
                        top.table = current.smpl.tt$qn,
                        cutoff)
    jacc.list[[current.n]]$z$array <-
      GetGeneSetJaccard(array.standard, 
                        top.table = current.smpl.tt$z,
                        cutoff)
    jacc.list[[current.n]]$z$seq <-
      GetGeneSetJaccard(seq.standard, 
                        top.table = current.smpl.tt$z,
                        cutoff)
    
  }
  
  jacc.df <- reshape2::melt(jacc.list)
  colnames(jacc.df) <-c("jaccard", "platform", "normalization", "no.samples")
  
  # rename platforms
  plt.recode.str <- 
    "'array' = 'Microarray'; 'seq' = 'RNA-seq'"
  jacc.df$platform <- car::recode(jacc.df$platform,
                                  recodes = plt.recode.str)
  
  # order no.samples so plot displays from smallest n to largest
  jacc.df$no.samples <- 
    factor(jacc.df$no.samples, 
           levels = sort(unique(as.numeric(as.character(jacc.df$no.samples)))))
  
  # capitalize normalization methods for display
  jacc.df$normalization <- as.factor(toupper(jacc.df$normalization))
  
  return(jacc.df)
  
}

SmallNDEGWrapper <- function(norm.list, sample.df, subtype = "Basal") {
  # Perform differential expression analysis for the small n experiments
  # from a list of normalized data from SmallNNormWrapper
  # 
  # Args:
  #   norm.list: a list of normalized expression data from SmallNNormWrapper
  #   sample.df: a data.frame mapping sample names to subtype labels
  #   subtype: which subtype should be compared to all others? default is Basal
  #   
  # Returns:
  #    fit.list: a list of differential expression results for all 4 elements
  #              of norm.list (fits from limma)
  # 
  # error-handling
  if (!("sample" %in% colnames(sample.df) 
        & "subtype" %in% colnames(sample.df))) {
    stop("sample.df must have columns named sample and subtype")
  }
  
  full.design.mat <- GetDesignMat(norm.list$log, sample.df, subtype)
  half.design.mat <- GetDesignMat(norm.list$qn, sample.df, subtype)
  
  # initialize list to hold fits
  fit.list <- list()
  
  # 100% array data control
  fit.list[["log"]] <- GetFiteBayes(norm.list$log, full.design.mat)
  
  # 100% seq untransformed (RSEM) data (control) requires limma::voom 
  exprs <- norm.list$un
  exprs <- as.data.frame(exprs[, 2:ncol(exprs), with = FALSE])
  rownames(exprs) <- norm.list$un$gene
  voom.counts <- limma::voom(counts = exprs)
  fit.list[["un"]] <- GetFiteBayes(voom.counts, full.design.mat)
  
  # 50-50 split data experiment - quantile normalization and z-transformation
  fit.list[["qn"]] <- GetFiteBayes(norm.list$qn, half.design.mat)
  fit.list[["z"]] <- GetFiteBayes(norm.list$z, half.design.mat)
  
  return(fit.list)
  
}

GetSamplesforMixingSmallN <- function(n, sample.df, subtype = "Basal") {
  # This function is designed to identify sample names to be used in the "small
  # n" differential expression experiment 
  # 
  # Arg:
  #   n: number of samples (for each subtype - no. of replicates)
  #   sample.df: a data.frame mapping sample names to subtype labels
  #   subtype: which subtype should be compared to all others? defaults to Basal
  # 
  # Returns:
  #   A list comprised of the following:
  #     all: names of all samples to be used for the experiment (100% of
  #          a single platform)
  #     array: names of samples to be used for array half of experiment
  #     seq: names of samples to be used for sequencing half of experiment
  #   
  
  # get samples from the subtype of interest
  subtype.samples <-  
    as.character(
      sample.df$sample[sample(which(sample.df$subtype == subtype), n)]
    )
  
  # get comparator group samples 
  other.samples <- 
    as.character(
      sample.df$sample[sample(which(sample.df$subtype != subtype), n)]
    )  
  
  # all samples
  all.samples <- c(subtype.samples, other.samples)
  # array half
  array.samples <- c(sample(subtype.samples, floor(n * .5)),
                     sample(other.samples, ceiling(n * .5)))
  # seq half
  seq.samples <- all.samples[!(all.samples %in% array.samples)]
  
  return(list("all" = all.samples, "array" = array.samples, 
              "seq" = seq.samples))
  
}  

SmallNNormWrapper <- function(array.dt, seq.dt, mix.list, zto = FALSE) {
  # This function is a wrapper for normalizing data for the "small n" 
  # differential expression experiment -- it will return a list of normalized
  # data that includes two controls (100% log-transformed array data and 
  # and 100% untransformed RNA-seq count data (RSEM)) and two experimental 
  # datasets (50-50 array and seq data either quantile normalized or 
  # z-transformed)
  # 
  # Args:
  #   array.dt: a data.table of all microarray data, first column contains gene
  #             identifiers, rows are genes, samples are columns
  #   seq.dt: a data.table of all seq data, first column contains gene
  #             identifiers, rows are genes, samples are columns  
  #   mix.list: a list with elements "all", "array", and "seq" - the output 
  #             of GetSamplesforMixingSmallN
  #   zto: logical - should data be zero to one transformed? default is FALSE
  #   
  # Returns:
  #   norm.list: a normalized list of data with log (100% array data),
  #              un (100% untransformed RNA-seq data), qn (quantile normalized 
  #              data that is 50-50 array/seq data [the array data is used as 
  #              the reference]), and z (z-transformed data that is 50-50 
  #              array/seq data [z-scores within platform + concatenation])
  #              
  
  # error-handling 
  if (!(all(names(mix.list) == c("all", "array", "seq")))) {
    stop("mix.list should have all, array, and seq elements (output of 
         GetSamplesForMixingSmallN)")
  }
  
  # array and seq data.tables to be used as controls
  array.full.dt <- array.dt[, c(1, which(colnames(array.dt) %in% 
                                           mix.list$all)), with = FALSE] 
  seq.full.dt  <- seq.dt[, c(1, which(colnames(seq.dt) %in% 
                                        mix.list$all)), with = FALSE] 
  
  # array and seq data.tables to be used in the 50-50 experiment     
  array.half.dt <- array.dt[, c(1, which(colnames(array.dt) %in% 
                                           mix.list$array)), with = FALSE]
  seq.half.dt <- seq.dt[, c(1, which(colnames(seq.dt) %in% 
                                       mix.list$seq)), with = FALSE]
  
  # remove any rows (genes) that all have count = 1, will cause issues with
  # z-score
  GetAllOnesRowIndex <- function(x){
    vals <- x[, 2:ncol(x), with = FALSE]
    indx <- which(apply(vals, 1, function(x) all(x == 1)))
    return(indx)
  }
  
  all1.indx <- unique(c(GetAllOnesRowIndex(seq.full.dt), 
                        GetAllOnesRowIndex(seq.half.dt)))
  
  array.full.dt <- array.full.dt[-all1.indx, ]
  seq.full.dt <- seq.full.dt[-all1.indx, ]
  array.half.dt <- array.half.dt[-all1.indx, ]
  seq.half.dt <- seq.half.dt[-all1.indx, ]
  
  # initialize list to hold normalized data
  norm.list <- list()
  # control conditions 100% log-transformed array data & 100% RSEM seq data
  norm.list[["log"]] <- LOGArrayOnly(array.full.dt, zto)
  norm.list[["un"]] <- seq.full.dt
  # 50-50 experiments
  norm.list[["qn"]] <- QNProcessing(array.dt = array.half.dt,
                                    seq.dt = seq.half.dt,
                                    zero.to.one = zto)
  norm.list[["z"]] <- ZScoreProcessing(array.dt = array.half.dt,
                                       seq.dt = seq.half.dt, 
                                       zero.to.one = zto)
  return(norm.list)
  
}
