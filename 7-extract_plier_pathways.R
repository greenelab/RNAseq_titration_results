# S. Foltz Nov 2021
# The purpose of this analysis is to use PLIER to identify expression pathways
# in our data, using pure microarray and RNA-seq data as comparison standards
# for data coming from different normalization methods and titration levels.
#
# USAGE: Rscript 7-extract_plier_pathways.R --cancer_type --seed

option_list <- list(
  optparse::make_option("--cancer_type",
                        default = NA_character_,
                        help = "Cancer type"),
  optparse::make_option("--seed",
                        default = 8934,
                        help = "Random seed")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
source(here::here("util/option_functions.R"))
check_options(opt)

# load libraries
suppressMessages(source(here::here("load_packages.R")))
source(here::here("util", "color_blind_friendly_palette.R"))

# set options
cancer_type <- opt$cancer_type
file_identifier <- str_c(cancer_type, "subtype", sep = "_") # assuming subtype

# set seed
initial.seed <- opt$seed
set.seed(initial.seed)
message(paste("\nInitial seed set to:", initial.seed))

# define directories
data.dir <- here::here("data")
norm.data.dir <- here::here("normalized_data")
res.dir <- here::here("results")

# define input files
# finds first example of a subtypes file from cancer_type, does not rely on seed
#norm.test.object <- file.path(norm.data.dir,
#                              list.files(norm.data.dir,
#                                         pattern = paste0(file_identifier,
#                                                          "_array_seq_test_data_normalized_list_"))[1])
norm.train.object <- file.path(norm.data.dir,
                               list.files(norm.data.dir,
                                          pattern = paste0(file_identifier,
                                                           "_array_seq_train_titrate_normalized_list_"))[1])
sample.file <- file.path(res.dir,
                         list.files(res.dir,
                                    pattern = paste0(file_identifier,
                                                     "_matchedSamples_training_testing_split_labels_"))[1])

#### read in data --------------------------------------------------------------

#norm.test.list <- read_rds(norm.test.object)
norm.train.list <- read_rds(norm.train.object)
sample.df <- read.delim(sample.file)

#### set up PLIER data ---------------------------------------------------------

data(bloodCellMarkersIRISDMAP)
data(canonicalPathways)
data(oncogenicPathways)
data(svmMarkers)

all.paths <- PLIER::combinePaths(bloodCellMarkersIRISDMAP,
                                 canonicalPathways,
                                 oncogenicPathways,
                                 svmMarkers)

common.genes <- PLIER::commonRows(all.paths,
                                  norm.train.list#TODO UPDATE)

#### main ----------------------------------------------------------------------

# repeatable PLIER function
run_plier <- function(expr_data, paths, genes){

  # minimum k for PLIER = 2*num.pc
  set.k <- 2*PLIER::num.pc(expr_data[genes, ])

  # PLIER main function
  PLIER::PLIER(expr_data[genes, ],
               paths[genes, ],
               k = set.k,
               trace = TRUE,
               scale = F)
}

# create an output list
plier_results_list <- list()

# PLIER at pure array
plier_results_list["array"] <- run_plier(DATA,
                                         all.paths,
                                         common.genes)

# PLIER at pure RNA-seq
plier_results_list["seq"] <- run_plier(DATA,
                                       all.paths,
                                       common.genes)

# Do this at 0-100% RNA-seq titration levels
# across each normalization method
# parallel backend
cl <- parallel::makeCluster(detectCores() - 1)
doParallel::registerDoParallel(cl)

# at each titration level (0-100% RNA-seq)
plier_results_list["titrated"][1:9] <- foreach(seq_prop = seq(0.1, .9, 0.1), .packages = c("tidyverse")) %dopar% {

  run_plier(DATA,
            all.paths,
            common.genes)

}

# stop parallel backend
parallel::stopCluster(cl)

# renames list levels
names(plier_results_list["titrated"])[1:9] <- as.character(seq(10, 90, 10))
