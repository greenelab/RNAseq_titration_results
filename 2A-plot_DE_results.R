# J. Taroni Feb 2016
# This script compares differential expression "silver standards" which are
# differential expression analysis results from standard pipelines (i.e., 
# log transformed 100% array data and RSEM counts 100% RNA-seq [processed with
# limma::voom]) with differential expression results from RNA-seq titrated data 
# (0-100%) normalized various ways.
# 
# USAGE: Rscript 2A-plot_DE_results.R
 
library(ggplot2)

#### functions -----------------------------------------------------------------

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
    lapply(fit.list, 
           function(x) 
             lapply(x, 
                    function(y) GetAllGenesTopTable(y, adjust = adjust.method)))
  
  # how many genes are differentially expressed @ cutoff
  deg.count.count.list <- 
    lapply(top.table.list, 
           function(x) lapply(x, 
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

#### plot Basal v. Other results -----------------------------------------------

# basal v. other fit objects
basal.fit.list <- 
  readRDS(file.path(
    "results", "differential_expression",
    "BRCA_titration_differential_exp_eBayes_fits_BasalvOther.RDS"))

# plot proportion of genes that are diff expressed and get topTable(s)
basal.results <- PlotProportionDE(basal.fit.list)
ggsave(file.path("plots", 
                 "differential_expr_proportion_ltFDR5perc_BasalvOther.pdf"), 
       plot = basal.results$plot, width = 11, height = 4.25)

# plot jaccard similarity with silver standards
basal.jacc.plot <- PlotSilverStandardJaccard(basal.results$top.table.list,
                                             title = "Basal v. Other FDR < 5%")
ggsave(file.path("plots", "basal_v_other_jaccard_lineplot.pdf"),
       plot = basal.jacc.plot, width = 8.5, height = 4)


#### plot Her2 v. LumA results -------------------------------------------------

her2.fit.list <- 
  readRDS(file.path(
    "results", "differential_expression",
    "BRCA_titration_differential_exp_eBayes_fits_Her2vLumA.RDS"))

# plot proportion of genes that are diff expressed and get topTable(s)
her2.results <- PlotProportionDE(her2.fit.list)
ggsave(file.path("plots", 
                 "differential_expr_proportion_ltFDR5perc_Her2vLumA.pdf"), 
       plot = her2.results$plot, width = 11, height = 4.25)

# plot jaccard similarity with silver standards
her2.jacc.plot <- PlotSilverStandardJaccard(her2.results$top.table.list,
                                            title = "Her2 v. LumA FDR < 5%")
ggsave(file.path("plots", "her2_v_lumA_jaccard_lineplot.pdf"),
       plot = her2.jacc.plot, width = 8.5, height = 4)
