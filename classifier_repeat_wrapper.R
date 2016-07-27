# J. Taroni Jul 2016
# The purpose of this script is to read in TGCA array and sequencing data,
# to preprocess leaving only overlapping genes and samples with complete 
# subtype information, and to split the data intotraining and testing sets
# It should be run from the command line.
# USAGE: Rscript classifier_repeat_wrapper.R

kInitialSeed <- 12
set.seed(kInitialSeed)

seeds <- sample(1:10000, 10)

for(seed in seeds){
  system(paste("Rscript run_experiments.R", seed))
}

system("Rscript 3-plot_subtype_kappa.R")
