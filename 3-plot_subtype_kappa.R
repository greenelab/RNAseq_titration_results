# J. Taroni Jul 2016
# The purpose of this script is to plot Kappa statistics from predictions on 
# test data.It should be run from the command line through the 
# classifier_repeat_wrapper.R script

suppressMessages(source("load_packages.R"))
library(ggplot2)
library(data.table)

plot.dir <- "plots/"
res.dir <- "results/"
lf <- list.files(res.dir, full.names = TRUE)
array.files <- lf[grepl("BRCA_train_3_models_array_kappa_", lf)]
seq.files <- lf[grepl("BRCA_train_3_models_seq_kappa_", lf)]

plot.file.lead <- "BRCA_train_3_models_kappa_"

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
#### read in data --------------------------------------------------------------
# read in the tables that contain the kappa statistics for predictions on test
# data
array.list <- list()  # initialize list that will hold all array tables
seq.list <- list()  # initialize list that will hold all the RNA-seq tables 
for (i in 1:length(array.files)) {
  array.list[[i]] <- fread(paste0(res.dir, array.files[i]), 
                           data.table = F)
  seq.list[[i]] <- fread(paste0(res.dir, seq.files[i]), 
                         data.table = F)
}
# combine all tables from each platform into a data.frame
array.df <- rbind.fill(array.list)
seq.df <- rbind.fill(seq.list)
rm(list=c("array.list", "seq.list"))

#### plot test set results -----------------------------------------------------
# color-blind friendly palette
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
               "#0072B2", "#D55E00", "#CC79A7")
# bind all kappa stats together
test.df <- cbind(rbind(array.df, seq.df),
                 c(rep("Microarray", nrow(array.df)), 
                   rep("RNA-seq", nrow(seq.df))))
colnames(test.df) <- c("Kappa", "Perc.Seq", "Classifier", 
                       "Normalization", "Platform")
test.df$Perc.Seq <- factor(test.df$Perc.Seq, levels = seq(0, 100, 10))
test.df$Classifier[which(test.df$Classifier == "glmnet")] <- "LASSO"
test.df$Classifier[which(test.df$Classifier == "rf")] <- "Random Forest"
test.df$Classifier[which(test.df$Classifier == "svm")] <- "Linear SVM"
test.df$Normalization <- as.factor(toupper(test.df$Normalization))

# plot performance of three models for a particular normalization method in a
# single plot 
norm.methods <- levels(test.df$Normalization)
for (norm in norm.methods) {
  plot.nm <- paste0(plot.dir, plot.file.lead, norm, "_test.pdf")
  ggplot(test.df[which(test.df$Normalization == norm), ], 
         aes(x=Perc.Seq, y=Kappa, color=Platform,
             shape=Platform)) + 
    facet_wrap(~Classifier, ncol=1) +
    stat_summary(fun.y = median, geom="line", aes(group=Platform),
                 position=position_dodge(0.2)) +
    stat_summary(fun.data=DataSummary, aes(group=Platform),
                 position=position_dodge(0.2), size=0.4) +
    ggtitle(norm) + 
    xlab("% RNA-seq samples") +
    theme_bw() +
    scale_colour_manual(values=cbPalette[c(2, 3)])
  ggsave(plot.nm, plot=last_plot(), height = 14, width = 8)
}

# plot performance of a classifier/model type on all normalization method in a
# single plot 
test.df$Classifier <- as.factor(test.df$Classifier)
cls.methods <- unique(test.df$Classifier)
for (cls in cls.methods) {
  plot.nm <- paste0(plot.dir, plot.file.lead, cls, "_test.pdf")
  ggplot(test.df[which(test.df$Classifier == cls), ], 
         aes(x=Perc.Seq, y=Kappa, color=Platform,
             shape=Platform)) + 
    facet_wrap(~Normalization, ncol=5) +
    stat_summary(fun.y = median, geom="line", aes(group=Platform),
                 position=position_dodge(0.2)) +
    stat_summary(fun.data=DataSummary, aes(group=Platform),
                 position=position_dodge(0.2), size=0.4) +
    ggtitle(cls) + 
    xlab("% RNA-seq samples") +
    theme_bw() +
    scale_colour_manual(values=cbPalette[c(2, 3)])
  ggsave(plot.nm, plot=last_plot(), height = 3.5, width = 15)
}
