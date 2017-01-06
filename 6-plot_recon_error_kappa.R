# J. Taroni Oct 2016
# This script plots reconstruction errors (MASE and RMSE) from 
# 4-ica_pca_feature_reconstruction.R and the Kappa statistics associated with
# predictions on reconstructed data from 5-predict_subtype_reconstructed_data.R
# as violin and line plots, respectively.
# USAGE: Rscript 6-plot_recon_error_kappa.R

library(ggplot2)
library(dplyr)

plot.dir <- "plots"
rcn.res.dir <- file.path("results", "reconstructed_data")

kap.plot.file.lead <- file.path(plot.dir, "BRCA_kappa_reconstructed_data_")
err.plot.file.lead <- file.path(plot.dir, "BRCA_reconstruction_error_")

kappa.df.files <- list.files(rcn.res.dir, pattern = "kappa", full.names = TRUE)
error.files <- list.files(rcn.res.dir, pattern = "BRCA_reconstruction_error",
                              full.names = TRUE)

#### functions -----------------------------------------------------------------

DataSummary <- function(x) {
  # This function is supplied to ggplot2::stat_summary in order to plot the 
  # median value of a vector as a point and the "confidence interval on the 
  # median" used in notched boxplots as a vertical line. See boxplot.stats for
  # more information.
  m <- median(x)
  conf <- boxplot.stats(x)$conf
  ymin <- min(conf)
  ymax <- max(conf)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# color-blind friendly palette
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
               "#0072B2", "#D55E00", "#CC79A7")
#### plot kappa stats ----------------------------------------------------------

# read in kappa data.frames from each replicate and bind -- line plot with 
# boxplot "confidence intervals"
kappa.df.list <- list()
fl.iter <- 1
for (fl in kappa.df.files) {
  kappa.df.list[[fl.iter]] <- data.table::fread(fl, data.table = FALSE)
  fl.iter <- fl.iter + 1
}
kappa.master.df <- as.data.frame(data.table::rbindlist(kappa.df.list))
rm(kappa.df.list)

# order Perc.seq so line plot displays 0-100
kappa.master.df$Perc.seq <- factor(kappa.master.df$Perc.seq, 
                                   levels = seq(0, 100, 10))

# rename classifiers
kappa.master.df$Classifier[which(kappa.master.df$Classifier == "glmnet")] <- 
  "LASSO"
kappa.master.df$Classifier[which(kappa.master.df$Classifier == "rf")] <- 
  "Random Forest"
kappa.master.df$Classifier[which(kappa.master.df$Classifier == "svm")] <- 
  "Linear SVM"
kappa.master.df$Classifier <- as.factor(kappa.master.df$Classifier)

# get norm and reconstruction methods as factors 
kappa.master.df$Normalization <- as.factor(kappa.master.df$Normalization)
kappa.master.df$Reconstruction <- as.factor(kappa.master.df$Reconstruction)

# rename platforms
kappa.master.df$Platform[which(kappa.master.df$Platform == "array")] <-
  "Microarray"
kappa.master.df$Platform[which(kappa.master.df$Platform == "seq")] <- "RNA-seq"
kappa.master.df$Platform <- as.factor(kappa.master.df$Platform)

# for each normalization method, plot kappa stats 
norm.methods <- levels(kappa.master.df$Normalization)
for (norm in norm.methods) {
  plot.nm <- paste0(kap.plot.file.lead, norm, ".pdf")  # each norm method
  # line plot is saved as a PDF
  ggplot(kappa.master.df[which(kappa.master.df$Normalization == norm), ], 
         aes(x = Perc.seq, y = Kappa, color = Platform,
             shape = Platform)) + 
    facet_wrap(Reconstruction ~ Classifier, ncol=3) +
    stat_summary(fun.y = median, geom = "line", aes(group = Platform),
                 position = position_dodge(0.2)) +
    stat_summary(fun.data = DataSummary, aes(group = Platform),
                 position = position_dodge(0.2), size = 0.4) +
    ggtitle(toupper(norm)) + 
    xlab("% RNA-seq samples") +
    theme_bw() +
    scale_colour_manual(values = cbPalette[c(2, 3)])
  ggsave(plot.nm, plot = last_plot(), height = 8.5, width = 11)
}
rm(kappa.master.df)

#### plot error measures -------------------------------------------------------

# read in error measure data.frames from each replicate and bind -- violin plot
error.df.list <- list()
for(fl.iter in seq_along(error.files)){
  error.df.list[[fl.iter]] <- data.table::fread(error.files[fl.iter], 
                                             data.table = FALSE)
}
error.master.df <- as.data.frame(data.table::rbindlist(error.df.list))
rm(error.df.list)

# order perc.seq so plot displays 0-100
error.master.df$perc.seq <- factor(error.master.df$perc.seq, 
                                   levels = seq(0, 100, by = 10))

# get norm and reconstruction methods as factors 
error.master.df$normalization <- as.factor(error.master.df$normalization)
error.master.df$error.measure <- as.factor(error.master.df$error.measure)

# rename platforms
error.master.df$platform[which(error.master.df$platform == "array")] <-
  "Microarray"
error.master.df$platform[which(error.master.df$platform == "seq")] <- "RNA-seq"
error.master.df$platform <- as.factor(error.master.df$platform)

# reconstruction method as factor
error.master.df$method <- as.factor(error.master.df$method)

# take the average of each genes error across replicates 
error.mean.df <- error.master.df %>% 
                    group_by(gene, perc.seq, normalization, error.measure, 
                             platform, method) %>% 
                    summarise_each(funs(mean))
rm(error.master.df)
colnames(error.mean.df) <- c("Gene", "Perc_seq", "Normalization", 
                             "Error_Measure", "Platform", "Method", 
                             "Mean_Value")

# for each normalization method, plot reconstruction error
norm.methods <- levels(error.mean.df$Normalization)
for (norm in norm.methods) {
  plot.nm <- paste0(err.plot.file.lead, norm, ".pdf")  # each norm method
  # violin plot is saved as a PDF
  ggplot(error.mean.df[which(error.mean.df$Normalization == norm), ], 
         aes(Perc_seq, y = Mean_Value, color = Platform, fill = Platform)) +
    facet_wrap(Method ~ Error_Measure, ncol = 2) +
    theme_bw() +
    scale_colour_manual(values = cbPalette[c(2, 3)]) +
    geom_violin(colour = "black", position = position_dodge(0.8),
                alpha = 0.2) +
    stat_summary(fun.y = median, geom = "line", aes(group = Platform),
                 position=position_dodge(0.8)) +
    stat_summary(fun.y = median, geom = "point", aes(group = Platform),
                 position = position_dodge(0.8)) +
    xlab("% RNA-seq") +
    ylab("Mean Value (per gene)")
    ggtitle(toupper(norm))
  ggsave(plot.nm, plot = last_plot(), height = 8.5, width = 8.5)
}
