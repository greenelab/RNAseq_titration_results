# J. Taroni Jul 2016
# This script is a wrapper for running the BRCA subtype pipeline repeatedly with
# different random seeds.
# It should be run from the command line.
# USAGE: Rscript classifier_repeat_wrapper.R <n.repeats>

args <- commandArgs(trailingOnly = TRUE)
n.repeats <- as.integer(args[1])
if (is.na(n.repeats)) {
  message("\nNumber of repeats set to default: 10")
  n.repeats <- 10
} else {
  message(paste("\nNumber of repeats set to", n.repeats))
}

initial.seed <- 12
set.seed(initial.seed)

seeds <- sample(1:10000, n.repeats)

rep.count <- 1
for(seed in seeds){
  message(paste("\n\n#### REPEAT NUMBER", rep.count, "####\n\n"))
  system(paste("Rscript run_experiments.R", seed))
  rep.count <- rep.count + 1
}

system("Rscript 3-plot_subtype_kappa.R")
