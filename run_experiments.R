# J. Taroni Jul 2016
# The purpose of this script is to run the BRCA subtype classifier pipeline
# for RNA-seq 'titration.'
# It should be run from the command line.
# USAGE: Rscript run_experiments.R --cancer_type [BRCA|GBM] --seed integer
# It also may be run through the classifier_repeat_wrapper.R

option_list <- list(
  optparse::make_option("--cancer_type",
                        default = NULL,
                        help = "Cancer type"),
  optparse::make_option("--predictor",
                        default = NULL,
                        help = "Predictor used"),
  optparse::make_option("--seed",
                        default = NULL,
                        help = "Random seed")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
source(here::here("util/option_functions.R"))
check_options(opt)

# set options
cancer_type <- opt$cancer_type
predictor <- opt$predictor

# set seed
initial.seed <- opt$seed
set.seed(initial.seed)

# these seeds should be between 1000 and 9999 (be 4 digits) to match later file name parsing
seeds <- sample(1000:9999, 3)

message(paste("Initial seed:", initial.seed))
message(paste("Secondary seeds:", stringr::str_c(seeds, collapse = ", ")))

message("Getting overlap and splitting into training and testing sets...")
system(paste("Rscript 0-expression_data_overlap_and_split.R",
             "--cancer_type", cancer_type,
             "--predictor", predictor,
             "--seed1", seeds[1]))

message("\nNormalizing data...")
system(paste("Rscript 1-normalize_titrated_data.R",
             "--cancer_type", cancer_type,
             "--predictor", predictor,
             "--seed1", seeds[1],
             "--seed2", seeds[2]))

message("\nTraining and testing models...")
system(paste("Rscript 2-train_test_brca_subtype.R",
             "--cancer_type", cancer_type,
             "--predictor", predictor,
             "--seed1", seeds[1],
             "--seed3", seeds[3]))
