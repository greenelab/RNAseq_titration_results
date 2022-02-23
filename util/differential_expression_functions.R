source(here::here("util", "normalization_functions.R"))
source(here::here("util", "train_test_functions.R"))
source(here::here("util", "color_blind_friendly_palette.R"))

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
        & "category" %in% colnames(sample.df))) {
    stop("sample.df must have columns named sample and category")
  }

  # get a factor vector of labels, ordered to match the columns of
  # norm.dt
  group.factor <-
    as.character(GetOrderedCategoryLabels(norm.dt,
                                          sample.df))
  group.factor[which(group.factor != subtype)] <- "Other"
  group.factor <- as.factor(group.factor)

  # build design matrix
  design.mat <- model.matrix(~0 + group.factor)
  colnames(design.mat) <- levels(group.factor)
  rownames(design.mat) <- colnames(norm.dt)[2:ncol(norm.dt)]

  return(design.mat)
}

GetDesignMatrixList <- function(norm.list, sample.df, subtype) {
  # This function takes a list of normalized expression data and returns
  # a list of design matrices to be used to detect differentially expressed
  # genes between one subtype and all other subtypes with
  # the limma package
  #
  # Args:
  #   norm.list: a nested list of normalized expression data; the top level of
  #              the list corresponds to % seq included in the expression data
  #   sample.df: a data.frame that maps the sample names to the subtype labels;
  #              'sample' and 'subtype' columns respectively
  #   subtype: which subtype should be compared to all others?
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
  cont.eval <- paste0(subtype, "vsOther= ", subtype, "-Other")

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

GetGeneSetStats <- function(silver.set,
                           top.table,
                           cutoff = 0.05) {
  # Given a top table list of from RNA-seq or microarray-only data
  # (our "silver standard" of differentially expressed genes), compare those 
  # true positives and true negatives differentially expressed genes (DEGs)
  # against DEGs meeting out cutoff criteria in other experimental settings
  # using Jaccard similarity, Rand index, and Spearman rank correlation
  #
  # Args:
  #   silver.set: output of limma::topTable (typically from GetAllGenesTopTable)
  #               for the silver standard RNA-seq or array data being compared against
  #   top.table: output of limma::topTable (typically from GetAllGenesTopTable)
  #              for the other differential expression experiment being
  #              considered
  #   cutoff: the adjusted p-value threshold - all genes with an adjusted p
  #           below this threshold will be considered differentially expressed
  #
  # Returns:
  #   jacc: Jaccard similarity
  #   rand: Rand index
  #   spearman: Spearman correlation rho
  
  # start with silver standard
  silver_df <- silver.set %>%
    rownames_to_column(var = "gene") %>%
    select(gene, adj.P.Val) %>%
    rename("silver.adj.P.Val" = "adj.P.Val") %>%
    mutate(silver_group = as.integer(silver.adj.P.Val < cutoff))
  
  # do the same with our experimental setting
  experimental_df <- top.table %>%
    rownames_to_column(var = "gene") %>%
    select(gene, adj.P.Val) %>%
    rename("experimental.adj.P.Val" = "adj.P.Val") %>%
    mutate(experimental_group = as.integer(experimental.adj.P.Val < cutoff))
    
  # join those data sets together to match up gene names
  combined_df <- silver_df %>%
    left_join(experimental_df,
              by = "gene") %>%
    mutate(silver_group = factor(silver_group, levels = c(0, 1)),
           experimental_group = factor(experimental_group, levels = c(0, 1)))
  
  # calculate agreement between two results
  contingency_table <- combined_df %>%
    select(silver_group,
           experimental_group) %>%
    table()
  
  TN <- contingency_table[1,1]
  TP <- contingency_table[2,2]
  total <- sum(contingency_table)
  
  # calculate and return jaccard, rand index, and spearman
  if (total == TN) { # denominator of 0 leads to NaN
    jacc <- 0
  } else {
    jacc <- TP/(total - TN)  
  }
  rand <- (TP + TN)/total
  spearman <- as.vector(cor.test(combined_df$silver.adj.P.Val,
                                 combined_df$experimental.adj.P.Val,
                                 method = "spearman",
                                 exact = FALSE,
                                 continuity = TRUE)$estimate)
  return(data.frame("jaccard" = jacc,
                    "rand" = rand,
                    "spearman" = spearman))
}

GetDataProportionDE <- function(top.table.list,
                                adjust.method = "BH", cutoff = 0.05) {
  # from a top.table.list, return the proportion of genes that are
  # considered differentially expressed under a user-specified adjusted p
  # cutoff (FDR < 5% by default)
  #
  # Args:
  #   top.table.list: a nested list from GetAllGenesTopTable()
  #   adjust.method: how should the p-values be adjusted for multiple
  #                  hypothesis testing?
  #   cutoff = what adjusted p-value cutoff to use
  #
  # Returns
  #   A data frame with the proportion of genes differentially expressed
  #
  
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
  
  return(deg.count.df)
}

PlotProportionDE <- function(propDE_df,
                             subtypes,
                             cancer_type,
                             cutoff = 0.05) {
  # plot the proportion of DEGs by normalization method
  #
  # Args:
  #   propDE_df: output from GetDataProportionDE, a data frame with the
  #     proportion of genes differentially expressed
  #   subtypes: character string of subtypes being compared (for plot title)
  #   cancer_type: cancer type of data being used (for plot title)
  #   cutoff = what adjusted p-value cutoff was used
  #
  # Returns
  #   A plot object
  #

  # bar plots of proportion of genes that are differentially expressed
  # in each experiment
  p <- ggplot(propDE_df,
              aes(x = perc.seq,
                  y = perc.ltdeg.count,
                  #fill = normalization
                  )) +
    facet_wrap(~ normalization,
               ncol = 4) +
    geom_bar(stat = "identity",
             position = "dodge") +
    labs(x = "% RNA-seq Samples",
         y = str_c("Proportion of Genes FDR < ", cutoff),
         title = str_c(cancer_type, " ", subtypes, " DEGs")) +
    scale_x_discrete(labels = c("0", "", "", "", "",
                                "50", "", "", "", "",
                                "100"),
                     drop = FALSE) +
    expand_limits(y = 1) +
    theme_bw()
  
  return(p)
}

GetDataSilverStandardStats <- function(top.table.list,
                                       cutoff = 0.05){
  # Given a list of top tables, get data to plot the Jaccard similarity, Rand index, and
  # Spearman rho between the "silver standards" and all other experiments
  #
  # Args:
  #   top.table.list: list of limma::topTable objects from GetAllGenesTopTable
  #   cutoff: adjusted p-value threshold to be used
  #
  # Returns:
  #   Data frame to plot Jaccard similarity, Rand index, Spearman rank correlation
  #
  
  ### "silver standards" ###
  
  # 100% RNA-seq data RSEM using limma::voom processing step
  #top.table.list$`100`$un
  
  # LOG 100% array data
  #top.table.list$`0`$log
  
  # how similiar are DEG results to the RNA-seq silver standard?
  seq.stats.list <-
    lapply(top.table.list,
           function(x) lapply(x,
                              function(y) GetGeneSetStats(top.table.list$`100`$un,
                                                          y, cutoff = cutoff)))
  seq.stats.df <- reshape2::melt(seq.stats.list,
                                 id.vars = c("jaccard", "rand", "spearman"))
  
  # how similiar are DEG results to the microarray silver standard?
  array.stats.list <-
    lapply(top.table.list,
           function(x) lapply(x,
                              function(y) GetGeneSetStats(top.table.list$`0`$log,
                                                          y, cutoff = cutoff)))
  array.stats.df <- reshape2::melt(array.stats.list,
                                   id.vars = c("jaccard", "rand", "spearman"))
  
  # combine seq and array similarity results
  array.stats.df <- cbind(array.stats.df, rep("Microarray", nrow(array.stats.df)))
  seq.stats.df <- cbind(seq.stats.df, rep("RNA-seq", nrow(seq.stats.df)))
  colnames(seq.stats.df) <- colnames(array.stats.df) <- c("Jaccard",
                                                          "Rand",
                                                          "Spearman",
                                                          "Normalization",
                                                          "Perc.Seq",
                                                          "Platform")
  mstr.df <- rbind(array.stats.df, seq.stats.df)
  
  # order % seq so plot displays 0-100
  mstr.df$perc.seq <- factor(mstr.df$perc.seq, levels = seq(0, 100, 10))
  
  # capitalize normalization methods for display
  mstr.df$normalization <- as.factor(toupper(mstr.df$normalization))
  
  # gather jaccard rand and spearman
  mstr.df <- mstr.df %>%
    gather("Jaccard", "Rand", "Spearman",
           key = "measure", value = "value")
  
  return(mstr.df)
  
}

PlotSilverStandardStats <- function(measure_df,
                                    title,
                                    single_measure = FALSE){
  # Plot the Jaccard similarity, Rand index, or Spearman rho between the
  # "silver standards" and all other experiments
  #
  # Args:
  #   propDE_df: output from GetDataSilverStandardStats, data frame with each set of stats
  #   title: main title of plot
  #   single_measure: specify TRUE if only one measure (Jaccard, Rand, or Spearman) is being used as input
  #
  # Returns:
  #   Jaccard similarity line plot, and
  #   Rand index line plot, and
  #   Spearman rank correlation line plot
    
  # line plots
  plot_obj <- ggplot(measure_df,
                     aes(x = Perc.Seq,
                         y = value,
                         color = Platform,
                         fill = Platform)) +
    geom_line(aes(group = Platform),
              position = position_dodge(0.7)) +
    geom_point(aes(group = Platform),
               position = position_dodge(0.7)) +
    scale_x_discrete(labels = c("0", "", "", "", "",
                                "50", "", "", "", "",
                                "100")) +
    expand_limits(y = c(0,1)) +
    scale_y_continuous(breaks = seq(0, 1, 0.25)) +
    scale_colour_manual(values = cbPalette[c(2, 3)]) +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(x = "% RNA-seq in Experimental Dataset",
         y = ifelse(single_measure,
                    unique(measure_df$measure),
                    "Measure of similarity"),
         color = "Silver Standard Comparison Platform",
         fill = "Silver Standard Comparison Platform",
         title = title)
  
  if (single_measure) {
    plot_obj <- plot_obj +
      facet_wrap(~ Normalization,
                 ncol = 7)
    } else {
    plot_obj <- plot_obj +
      facet_grid(rows = vars(measure),
                 cols = vars(Normalization))
  }
  
  return(plot_obj)
  
}

#### small n functions ---------------------------------------------------------

GetSmallNSilverStandardStats <- function(top.table.list, cutoff = 0.05){
  # This function takes a list of limma::topTable output, derives "silver
  # standards" from 100% array and 100% seq data, and finds three stats:
  # Jaccard similarity, Rand index, and Spearman correlation
  # between the standards and any experiment differential expression results
  # (here, quantile norm and z-score are used for illustration)
  #
  # Args:
  #   top.table.list: a list of limma::topTable output with all genes -
  #                   output of GetAllGenesTopTable
  #   cutoff: FDR cutoff, defaults to FDR < 5%
  #
  # Returns:
  #   stats.df: a data.frame of three stats results with columns corresponding to
  #             the Jaccard value, Rand index, Spearman correlation,
  #             the silver standard used for the comparison ("platform"),
  #             the normalization method, and the n used (number of samples; "no.samples")
  #
  
  # initialize list to hold stats
  stats.list <- list()
  
  # for each n (number of samples)
  for (smpl.no.iter in seq_along(top.table.list)) {
    
    current.smpl.tt <- top.table.list[[smpl.no.iter]]
    current.n <- names(top.table.list)[smpl.no.iter]
    
    array.top <- current.smpl.tt$log
    seq.top <- current.smpl.tt$un
    
    stats.list[[current.n]]$qn$array <-
      GetGeneSetStats(silver.set = array.top,
                      top.table = current.smpl.tt$qn,
                      cutoff = cutoff)
    stats.list[[current.n]]$qn$seq <-
      GetGeneSetStats(silver.set = seq.top,
                      top.table = current.smpl.tt$qn,
                      cutoff = cutoff)
    stats.list[[current.n]]$z$array <-
      GetGeneSetStats(silver.set = array.top,
                      top.table = current.smpl.tt$z,
                      cutoff = cutoff)
    stats.list[[current.n]]$z$seq <-
      GetGeneSetStats(silver.set = seq.top,
                      top.table = current.smpl.tt$z,
                      cutoff = cutoff)
  }
  
  stats.df <- reshape2::melt(stats.list,
                             id.vars = c("jaccard", "rand", "spearman"))
  colnames(stats.df) <-c("Jaccard", "Rand", "Spearman",
                         "platform", "normalization", "no.samples")
  
  # rename platforms
  plt.recode.str <-
    "'array' = 'Microarray'; 'seq' = 'RNA-seq'"
  stats.df$platform <- car::recode(stats.df$platform,
                                   recodes = plt.recode.str)
  
  # order no.samples so plot displays from smallest n to largest
  stats.df$no.samples <-
    factor(stats.df$no.samples,
           levels = sort(unique(as.numeric(as.character(stats.df$no.samples)))))
  
  # capitalize normalization methods for display
  stats.df$normalization <- as.factor(toupper(stats.df$normalization))
  
  return(stats.df)
  
}

SmallNDEGWrapper <- function(norm.list, sample.df, subtype) {
  # Perform differential expression analysis for the small n experiments
  # from a list of normalized data from SmallNNormWrapper
  #
  # Args:
  #   norm.list: a list of normalized expression data from SmallNNormWrapper
  #   sample.df: a data.frame mapping sample names to subtype labels
  #   subtype: which subtype should be compared to all others?
  #
  # Returns:
  #    fit.list: a list of differential expression results for all 4 elements
  #              of norm.list (fits from limma)
  #
  # error-handling
  if (!("sample" %in% colnames(sample.df)
        & "category" %in% colnames(sample.df))) {
    stop("sample.df must have columns named sample and category")
  }

  full.design.mat <- GetDesignMat(norm.list$log, sample.df, subtype)
  part.design.mat <- GetDesignMat(norm.list$qn, sample.df, subtype)

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

  # (100-X)%/X% split data experiment - quantile normalization and z-transformation
  fit.list[["qn"]] <- GetFiteBayes(norm.list$qn, part.design.mat)
  fit.list[["z"]] <- GetFiteBayes(norm.list$z, part.design.mat)

  return(fit.list)

}

GetSamplesforMixingSmallN <- function(n, sample.df, subtype, seq_proportion) {
  # This function is designed to identify sample names to be used in the "small
  # n" differential expression experiment
  #
  # Arg:
  #   n: number of samples (for each subtype - no. of replicates)
  #   sample.df: a data.frame mapping sample names to subtype labels
  #   subtype: which subtype should be compared to all others
  #   seq_proportion: percentage of RNA-seq samples to include in mix
  #
  # Returns:
  #   A list comprised of the following:
  #     all: names of all samples to be used for the experiment (100% of
  #          a single platform)
  #     array: names of samples to be used for array part of experiment
  #     seq: names of samples to be used for sequencing part of experiment
  #

  # get samples from the subtype of interest
  subtype.samples <-
    as.character(
      sample.df$sample[sample(which(sample.df$category == subtype), n)]
    )

  # get comparator group samples
  other.samples <-
    as.character(
      sample.df$sample[sample(which(sample.df$category != subtype), n)]
    )

  # all samples
  all.samples <- c(subtype.samples, other.samples)
  # array part
  array.samples <- c(sample(subtype.samples, floor(n * (1 - seq_proportion))),
                     sample(other.samples, ceiling(n * (1 - seq_proportion))))
  # seq part
  seq.samples <- all.samples[!(all.samples %in% array.samples)]

  return(list("all" = all.samples, "array" = array.samples,
              "seq" = seq.samples))

}

SmallNNormWrapper <- function(array.dt, seq.dt, mix.list, zto = FALSE) {
  # This function is a wrapper for normalizing data for the "small n"
  # differential expression experiment -- it will return a list of normalized
  # data that includes two controls (100% log-transformed array data and
  # and 100% untransformed RNA-seq count data (RSEM)) and two experimental
  # datasets (titrated array and seq data either quantile normalized or
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
  #              data that is titrated array/seq data [the array data is used as
  #              the reference]), and z (z-transformed data that is titrated
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
  
  # array and seq data.tables to be used in the titrated experiment
  array.part.dt <- array.dt[, c(1, which(colnames(array.dt) %in%
                                           mix.list$array)), with = FALSE]
  seq.part.dt <- seq.dt[, c(1, which(colnames(seq.dt) %in%
                                       mix.list$seq)), with = FALSE]

  # remove any rows (genes) that all have same values, will cause issues with
  # z-score
  GetAllSameRowIndex <- function(x){
    if (ncol(x) > 1) { # check that there is at least one data column
      # there will be no data column when no seq or no array samples are included
      vals <- x[, 2:ncol(x), with = FALSE]
      indx <- which(apply(vals, 1, check_all_same))
      return(indx)
    } else {
      return(integer(0))
    }
  }

  all.same.indx <- unique(c(GetAllSameRowIndex(seq.full.dt),
                            GetAllSameRowIndex(seq.part.dt)))
  # if no rows are all same (in previous GetAllSameRowIndex), all.same.indx is integer(0)
  # subsetting data frames by -integer(0) results in no rows
  # so check that integer vector has length > 0 before subsetting
  if (length(all.same.indx) > 0) {
    array.full.dt <- array.full.dt[-all.same.indx, ]
    seq.full.dt <- seq.full.dt[-all.same.indx, ]
    array.part.dt <- array.part.dt[-all.same.indx, ]
    seq.part.dt <- seq.part.dt[-all.same.indx, ]
  }

  # initialize list to hold normalized data
  norm.list <- list()
  # control conditions 100% log-transformed array data & 100% RSEM seq data
  norm.list[["log"]] <- LOGArrayOnly(array.full.dt, zto)
  norm.list[["un"]] <- seq.full.dt
  # titrated experiments
  if (ncol(seq.part.dt) == 1) { # no seq samples added, array only
    norm.list[["qn"]] <- QNSingleDT(dt = array.part.dt,
                                    zero.to.one = zto)
    norm.list[["z"]] <- ZScoreSingleDT(dt = array.part.dt,
                                       zero.to.one = zto)
  } else if (ncol(array.part.dt) == 1) { # no array samnples added, seq only
    norm.list[["qn"]] <- QNSingleDT(dt = seq.part.dt,
                                    zero.to.one = zto)
    norm.list[["z"]] <- ZScoreSingleDT(dt = seq.part.dt,
                                       zero.to.one = zto)
    
  } else {
    norm.list[["qn"]] <- QNProcessing(array.dt = array.part.dt,
                                      seq.dt = seq.part.dt,
                                      zero.to.one = zto)
    norm.list[["z"]] <- ZScoreProcessing(array.dt = array.part.dt,
                                         seq.dt = seq.part.dt,
                                         zero.to.one = zto)
  }
  
  return(norm.list)
  
}
