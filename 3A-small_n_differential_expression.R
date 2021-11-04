# J. Taroni Feb 2016
# The purpose of this analysis is to examine how normalization methods
# (quantile normalization or z-transformation) perform wrt differential
# expression when there are a small number of samples on each platform
#
# USAGE: Rscript 3A-small_n_differential_expression.R --cancer_type --subtype_vs_subtype

option_list <- list(
  optparse::make_option("--cancer_type",
                        default = NA_character_,
                        help = "Cancer type"),
  optparse::make_option("--subtype_vs_subtype",
                        default = NA_character_,
                        help = "Subtypes used in head-to-head comparison (comma-separated without space e.g. Type1,Type2)"),
  optparse::make_option("--seed",
                        default = 3255,
                        help = "Random seed")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
source(here::here("util/option_functions.R"))
check_options(opt)

# load libraries
suppressMessages(source(here::here("load_packages.R")))
source(here::here("util", "normalization_functions.R"))
source(here::here("util", "differential_expression_functions.R"))
source(here::here("util", "color_blind_friendly_palette.R"))

# set options
cancer_type <- opt$cancer_type
subtype_vs_subtype <- opt$subtype_vs_subtype
two_subtypes <- as.vector(stringr::str_split(subtype_vs_subtype, pattern = ",", simplify = TRUE))
file_identifier <- str_c(cancer_type, "subtype", sep = "_") # we are only working with subtype models here

# set seed
initial.seed <- opt$seed
set.seed(initial.seed)
message(paste("\nInitial seed set to:", initial.seed))

# define directories
data.dir <- here::here("data")
res.dir <- here::here("results")
deg.dir <- file.path(res.dir, "differential_expression")


# define input files
seq.file <- file.path(data.dir,
                      paste0(cancer_type, "RNASeq_matchedOnly_ordered.pcl"))
array.file <- file.path(data.dir,
                        paste0(cancer_type, "array_matchedOnly_ordered.pcl"))
smpl.file <- file.path(res.dir,
                       list.files(res.dir, # this finds the first example of a subtypes file from cancer_type
                                  pattern = paste0(file_identifier, # and does not rely on knowing a seed
                                                   "_matchedSamples_training_testing_split_labels_"))[1])

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

# check that subtypes are in sample.df
for(subtype in two_subtypes) {
  if (!(subtype %in% sample.df$category)) {
    stop(paste("Subtype", subtype, "not found in sample file",
               smpl.file, "in 3A-small_n_differential_expression.R."))
  }
}

sample.names <- sample.df$sample

#### main ----------------------------------------------------------------------

# leave only subtypes of interest to choose from & make data.table
# remove all samples that are not subtypes of interest
samples.to.keep <-
  sample.df$sample[which(sample.df$category %in% two_subtypes)]

array.dt <- data.table(array.data[,
                                  c(1, which(colnames(array.data) %in%
                                               samples.to.keep))])
seq.dt <- data.table(seq.data[,
                              c(1, which(colnames(seq.data) %in%
                                           samples.to.keep))])
sample.df <- sample.df[which(sample.df$sample %in% samples.to.keep), ]

smaller_subtype_size <- min(table(droplevels(sample.df$category)))

# different sizes of n to test
no.samples <- c(3, 4, 5, 6, 8, 10, 15, 25, 50)
no.samples <- no.samples[which(no.samples <= smaller_subtype_size)]

message(paste("Smaller subtype has", smaller_subtype_size, "samples,",
              "so using up to", max(no.samples), "samples in 3A-small_n_differential_expression.R"))

# initialize list to hold Jaccard, Rand, Spearman data from the 10 trials
stats.df.list <- list()

# Do this at 0-100% RNA-seq titration levels
# parallel backend
cl <- parallel::makeCluster(detectCores() - 1)
doParallel::registerDoParallel(cl)

# at each titration level (0-100% RNA-seq)
stats.df.list[1:9] <- foreach(seq_prop = seq(0.1, .9, 0.1), .packages = c("tidyverse")) %dopar% {
  
  # we're going to repeat the small n experiment 10 times
  stats.df.iter_list <- list() # this is returned to stats.df.list each iteration
  for (trial.iter in 1:10) {
    
    # for each n (3...50), get the sample names that will be included in the
    # experiment and on each platform
    sample.list <-
      lapply(no.samples,  # for each n (3...50)
             function(x) GetSamplesforMixingSmallN(x, sample.df,
                                                   subtype = data.table::last(two_subtypes),
                                                   seq_proportion = seq_prop))
    
    # initialize list to hold differential expression results (eBayes output)
    master.deg.list <- list()
    
    for (smpl.no.iter in seq_along(sample.list)) {  # for each n (3...50)
      # normalize data
      n_array <- length(sample.list[[smpl.no.iter]]$array)
      n_seq <- length(sample.list[[smpl.no.iter]]$seq)
      
      if (n_array >= 3 & n_seq >= 3) { # require at least three array and seq samples
        norm.list <- SmallNNormWrapper(array.dt = array.dt,
                                       seq.dt = seq.dt,
                                       mix.list = sample.list[[smpl.no.iter]],
                                       zto = FALSE)
        # perform differential expression analysis
        master.deg.list[[as.character(no.samples[smpl.no.iter])]] <-
          SmallNDEGWrapper(norm.list = norm.list, sample.df = sample.df,
                           subtype = data.table::last(two_subtypes)) 
      }
    }
    
    top.table.list <-
      lapply(master.deg.list,  # for each n (3...50)
             function(x)  # for each normalization method
               lapply(x, function(y) GetAllGenesTopTable(y)))  # extract DEGs
    
    # how do the 50/50 array/seq differentially expressed genes compared to
    # the platform-specific standards?
    if (length(top.table.list) > 0) {
      stats.df.iter_list[[trial.iter]] <- GetSmallNSilverStandardStats(top.table.list,
                                                                       cutoff = 0.1)  
    }
  }
  stats.df.iter_list # return stats.df.iter_list to stats.df.list
}

# stop parallel backend
parallel::stopCluster(cl)

# renames list levels
names(stats.df.list)[1:9] <- as.character(seq(10, 90, 10))

# combine jaccard similarity data.frames into one data.frame
subtypes_combination <- stringr::str_c(two_subtypes, collapse = "v")
subtypes_combination_nice <- stringr::str_c(two_subtypes, collapse = " vs. ")

stats.df <- reshape2::melt(stats.df.list,
                           id.vars = c("platform", "normalization", "no.samples"))
names(stats.df) <- c("platform", "normalization", "no.samples", "metric", "value",
                     "iteration", "seq_prop")
stats.df <- pivot_wider(stats.df,
                        names_from = "metric",
                        values_from = "value")

write.table(stats.df,
            file = file.path(deg.dir,
                             paste0(file_identifier,
                                    "_small_n_",
                                    subtypes_combination,
                                    "_results.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

# line plot is saved as a PDF
# TODO future problem: make sure the x-axis values are consistent across plots
for (percent_rna_seq in as.integer(names(stats.df.list))) {
  stats.df.pct <- stats.df %>%
    filter(seq_prop == percent_rna_seq)
  
  ggplot(stats.df.pct, aes(x = no.samples,
                           y = jaccard,
                           color = platform)) +
    facet_wrap(~ normalization, ncol = 1) +
    stat_summary(fun = median, geom = "line", aes(group = platform),
                 position = position_dodge(0.2)) +
    stat_summary(fun.data = DataSummary, aes(group = platform),
                 position = position_dodge(0.2), size = 0.2) +
    theme_bw() +
    ggtitle(paste(subtypes_combination_nice, "FDR < 10%")) +
    ylab("Jaccard similarity") +
    xlab("Number of samples (n)") +
    scale_colour_manual(values = cbPalette[c(2, 3)]) +
    theme(text = element_text(size = 18))
  ggsave(filename = here::here("plots",
                               str_c(file_identifier, "small_n", subtypes_combination,
                                     percent_rna_seq, "pct_rna_seq_jaccard_lineplots.pdf", sep = "_")),
         plot = last_plot(), width = 5, height = 7)
  
  ggplot(stats.df.pct, aes(x = no.samples,
                           y = rand,
                           color = platform)) +
    facet_wrap(~ normalization, ncol = 1) +
    stat_summary(fun = median, geom = "line", aes(group = platform),
                 position = position_dodge(0.2)) +
    stat_summary(fun.data = DataSummary, aes(group = platform),
                 position = position_dodge(0.2), size = 0.2) +
    theme_bw() +
    ggtitle(paste(subtypes_combination_nice, "FDR < 10%")) +
    ylab("Rand index") +
    xlab("Number of samples (n)") +
    scale_colour_manual(values = cbPalette[c(2, 3)]) +
    theme(text = element_text(size = 18))
  ggsave(filename = here::here("plots",
                               str_c(file_identifier, "small_n", subtypes_combination,
                                     percent_rna_seq, "pct_rna_seq_rand_lineplots.pdf", sep = "_")),
         plot = last_plot(), width = 5, height = 7)
  
  ggplot(stats.df.pct, aes(x = no.samples,
                           y = spearman,
                           color = platform)) +
    facet_wrap(~ normalization, ncol = 1) +
    stat_summary(fun = median, geom = "line", aes(group = platform),
                 position = position_dodge(0.2)) +
    stat_summary(fun.data = DataSummary, aes(group = platform),
                 position = position_dodge(0.2), size = 0.2) +
    theme_bw() +
    ggtitle(paste(subtypes_combination_nice, "FDR < 10%")) +
    ylab("Spearman correlation") +
    xlab("Number of samples (n)") +
    scale_colour_manual(values = cbPalette[c(2, 3)]) +
    theme(text = element_text(size = 18))
  ggsave(filename = here::here("plots",
                               str_c(file_identifier, "small_n", subtypes_combination,
                                     percent_rna_seq, "pct_rna_seq_spearman_lineplots.pdf", sep = "_")),
         plot = last_plot(), width = 5, height = 7)
}
