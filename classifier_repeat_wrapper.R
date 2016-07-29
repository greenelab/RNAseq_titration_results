# J. Taroni Jul 2016
# This script is a wrapper for running the BRCA subtype pipeline repeatedly with
# different random seeds.
# It should be run from the command line.
# USAGE: Rscript classifier_repeat_wrapper.R <n.repeats>

args <- commandArgs(trailingOnly = TRUE)
n.repeats <- as.integer(args[1])
if (is.na(n.repeats)) {
  message("Number of repeats set to default: 10")
  n.repeats <- 10
} else {
  message(paste("Number of repeats set to", n.repeats))
}

kInitialSeed <- 12
set.seed(kInitialSeed)

seeds <- sample(1:10000, n.repeats)

for(seed in seeds){
  system(paste("Rscript run_experiments.R", seed))
}

system("Rscript 3-plot_subtype_kappa.R")
