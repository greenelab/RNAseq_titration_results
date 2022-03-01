# J. Taroni Jul 2016
# The purpose of this script is to run the BRCA subtype classifier pipeline
# for RNA-seq 'titration.'
# It should be run from the command line.
# USAGE: Rscript run_experiments.R --cancer_type [BRCA|GBM] --predictor [subtype|TP53|PIK3CA] --seed integer --null_model --ncores
# It also may be run through the classifier_repeat_wrapper.R

option_list <- list(
  optparse::make_option("--cancer_type",
                        default = NA_character_,
                        help = "Cancer type"),
  optparse::make_option("--predictor",
                        default = NA_character_,
                        help = "Predictor used"),
  optparse::make_option("--seed",
                        default = NA_integer_,
                        help = "Random seed"),
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
null_model <- opt$null_model
ncores <- min(parallel::detectCores() - 1,
              opt$ncores,
              na.rm = TRUE)

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
             "--seed1", seeds[1],
             ifelse(null_model,
                    "--null_model",
                    "")))

message("\nNormalizing data...")
system(paste("Rscript 1-normalize_titrated_data.R",
             "--cancer_type", cancer_type,
             "--predictor", predictor,
             "--seed1", seeds[1],
             "--seed2", seeds[2],
             ifelse(null_model,
                    "--null_model",
                    ""),
             "--ncores", ncores))

message("\nTraining and testing models...")
system(paste("Rscript 2-train_test_category.R",
             "--cancer_type", cancer_type,
             "--predictor", predictor,
             "--seed1", seeds[1],
             "--seed3", seeds[3],
             ifelse(null_model,
                    "--null_model",
                    ""),
             "--ncores", ncores))
