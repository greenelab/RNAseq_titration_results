# J. Taroni Feb 2016
# The purpose of this analysis is to examine how normalization methods
# (quantile normalization or z-transformation) perform wrt differential 
# expression when there are a small number of samples on each platform 
# (50-50 split microarray and RNA-seq).
# 
# USAGE: Rscript 3A-small_n_differential_expression.R

suppressMessages(source("load_packages.R"))
source(file.path("util", "normalization_functions.R"))
source(file.path("util", "differential_expression_functions.R"))
source(file.path("util", "color_blind_friendly_palette.R"))

args <- commandArgs(trailingOnly = TRUE)
initial.seed <- as.integer(args[1])
if (is.na(initial.seed)) {
  message("\nInitial seed set to default: 3255")
  initial.seed <- 3255
} else {
  message(paste("\nInitial seed set to:", initial.seed))
}
set.seed(initial.seed)

data.dir <- "data"
deg.dir <- file.path("results", "differential_expression")

seq.file <- file.path(data.dir, "BRCARNASeq_matchedOnly_ordered.pcl")
array.file <- file.path(data.dir, "BRCAarray_matchedOnly_ordered.pcl")

smpl.file <- 
  file.path("results", 
            "BRCA_matchedSamples_PAM50Array_training_testing_split_labels_3061.tsv")

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
  return(c(y = m, ymin = ymin, ymax = ymax))
}

#### read in data --------------------------------------------------------------

seq.data <- data.table::fread(seq.file, data.table = F)
array.data <- data.table::fread(array.file, data.table = F)
sample.df <- read.delim(smpl.file)

sample.names <- sample.df$sample

#### main ----------------------------------------------------------------------

# leave only Her2 and LumA samples to choose from & make data.table
# remove all samples that are not Her2 or LumA 
samples.to.keep <- 
  sample.df$sample[which(sample.df$subtype %in% c("LumA", "Her2"))]
  
array.dt <- data.table(array.data[, 
                                  c(1, which(colnames(array.data) %in% 
                                      samples.to.keep))])
seq.dt <- data.table(seq.data[, 
                                c(1, which(colnames(seq.data) %in% 
                                             samples.to.keep))])
sample.df <- sample.df[which(sample.df$sample %in% samples.to.keep), ]

# different sizes of n to test
no.samples <- c(3, 4, 5, 6, 8, 10, 15, 25, 50)

# initialize list to hold jaccard index data.frames from the 10 trials
jacc.df.list <- list()

# we're going to repeat the small n experiment 10 times
for (trial.iter in 1:10) {
  
  # for each n (3...50), get the sample names that will be included in the
  # experiment and on each platform
  sample.list <- 
    lapply(no.samples,  # for each n (3...50)
           function(x) GetSamplesforMixingSmallN(x, sample.df, 
                                                 subtype = "Her2"))
  
  # initialize list to hold differential expression results (eBayes output)
  master.deg.list <- list()
  
  for (smpl.no.iter in seq_along(sample.list)) {  # for each n (3...50)
    # normalize data
    norm.list <- SmallNNormWrapper(array.dt = array.dt,
                                   seq.dt = seq.dt,
                                   mix.list = sample.list[[smpl.no.iter]],
                                   zto = FALSE)
    # perform differential expression analysis
    master.deg.list[[as.character(no.samples[smpl.no.iter])]] <- 
      SmallNDEGWrapper(norm.list = norm.list, sample.df = sample.df, 
                       subtype = "Her2")
  }
  
  top.table.list <- 
    lapply(master.deg.list,  # for each n (3...5)
           function(x)  # for each normalization method
             lapply(x, function(y) GetAllGenesTopTable(y)))  # extract DEGs
  
  # how do the 50/50 array/seq differentially expressed genes compared to
  # the platform-specific standards?
  jacc.df.list[[trial.iter]] <- GetSmallNSilverStandardJaccard(top.table.list,
                                                               cutoff = 0.1) 
  
}

# combine jaccard similarity data.frames into one data.frame


jacc.df <- data.table::rbindlist(jacc.df.list)

write.table(jacc.df, 
            file = file.path("results", "differential_expression",
                             "small_n_Her2vLumA_50-50_jaccard_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# line plot is saved as a PDF
ggplot(jacc.df, aes(x = no.samples, y = jaccard, color = platform)) + 
  facet_wrap(~ normalization, ncol = 1) +
  stat_summary(fun.y = median, geom = "line", aes(group = platform),
               position = position_dodge(0.2)) +
  stat_summary(fun.data = DataSummary, aes(group = platform),
               position = position_dodge(0.2), size = 0.2) +
  theme_bw() +
  ggtitle("Her2 vs. LumA FDR < 10%") +
  ylab("Jaccard similarity to standard") +
  xlab("Number of samples (n)") +
  scale_colour_manual(values = cbPalette[c(2, 3)])
ggsave(filename = file.path("plots", 
                            "small_n_Her2vLumA_50-50_jaccard_lineplots.pdf"),
       plot = last_plot())
