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
  optparse::make_option("--seed",
                        default = NULL,
                        help = "Random seed used to initiate this repeat")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
source("util/option_functions.R")
check_options(opt)

initial.seed <- opt$seed
set.seed(initial.seed)

# these seeds should be between 1000 and 9999 (be 4 digits) to match later file name parsing
seeds <- sample(1000:9999, 3)

message(paste("Initial seed:", initial.seed))
message(paste("Secondary seeds:", stringr::str_c(seeds, collapse = ", ")))

message("Getting overlap and splitting into training and testing sets...")
system(paste("Rscript 0-expression_data_overlap_and_split.R", seeds[1]))
message("\nNormalizing data...")
system(paste("Rscript 1-normalize_titrated_data.R", seeds[1], seeds[2]))
message("\nTraining and testing models...")
system(paste("Rscript 2-train_test_brca_subtype.R", seeds[1], seeds[3]))
