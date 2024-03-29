# J. Taroni Jul 2016
# This script is a wrapper for running the BRCA subtype pipeline repeatedly with
# different random seeds.
# It should be run from the command line.
# USAGE: Rscript classifier_repeat_wrapper.R --cancer_type [BRCA|GBM] --predictor [subtype|TP53|PIK3CA] --n_repeats (default: 10) --null_model --ncores

option_list <- list(
  optparse::make_option("--cancer_type",
                        default = NA_character_,
                        help = "Cancer type"),
  optparse::make_option("--predictor",
                        default = NA_character_,
                        help = "Predictor used"),
  optparse::make_option("--n_repeats",
                        default = 10,
                        help = "Number of times experiment is repeated [default: %default]"),
  optparse::make_option("--null_model",
                        action = "store_true",
                        default = FALSE,
                        help = "Permute dependent variable (within subtype if predictor is a gene)"),
  optparse::make_option("--ncores",
                        default = NA_integer_,
                        help = "Set the number of cores to use")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
source(here::here("util/option_functions.R"))
check_options(opt)

# set options
cancer_type <- opt$cancer_type
predictor <- opt$predictor
n.repeats <- opt$n_repeats
null_model <- opt$null_model
ncores <- min(parallel::detectCores() - 1,
              opt$ncores,
              na.rm = TRUE)

message(paste("\nPredicting", predictor,
              "in", cancer_type,
              ifelse(null_model,
                     "(null model) ...",
                     "...")))
message(paste("\nNumber of repeats set to", n.repeats))
message(paste("\nUsing", ncores, "out of", parallel::detectCores(), "cores"))

initial.seed <- 12
set.seed(initial.seed)

seeds <- sample(1:10000, n.repeats)

rep.count <- 1
for(seed in seeds){
  message(paste("\n\n#### REPEAT NUMBER", rep.count, "####\n\n"))
  system(paste("Rscript run_experiments.R",
               "--cancer_type", cancer_type,
               "--predictor", predictor,
               "--seed", seed,
               ifelse(null_model,
                      "--null_model",
                      ""),
               "--ncores", ncores))
  rep.count <- rep.count + 1
}

system(paste("Rscript 3-combine_category_kappa.R",
             "--cancer_type", cancer_type,
             "--predictor", predictor,
             ifelse(null_model,
                    "--null_model",
                    "")))
