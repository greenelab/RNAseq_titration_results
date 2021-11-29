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
source(here::here("util", "normalization_functions.R"))
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
#norm.test.files <- file.path(norm.data.dir,
#                            list.files(norm.data.dir,
#                                       pattern = paste0(file_identifier,
#                                                        "_array_seq_test_data_normalized_list_")))
norm.train.files <- file.path(norm.data.dir,
                             list.files(norm.data.dir,
                                        pattern = paste0(file_identifier,
                                                         "_array_seq_train_titrate_normalized_list_")))
sample.files <- file.path(res.dir,
                         list.files(res.dir,
                                    pattern = paste0(file_identifier,
                                                     "_matchedSamples_training_testing_split_labels_")))

#### set up PLIER data ---------------------------------------------------------

data(bloodCellMarkersIRISDMAP)
data(canonicalPathways)
data(oncogenicPathways)
data(svmMarkers)

all.paths <- PLIER::combinePaths(bloodCellMarkersIRISDMAP,
                                 canonicalPathways,
                                 oncogenicPathways,
                                 svmMarkers)

#### loop over data for each seed and get PLIER results ------------------------

#for(seed_index in 1:length(norm.train.files)) {
for(seed_index in 1:2) {
  
  print(str_c("PLIER with data seed", seed_index,
              "out of", length(norm.train.files), "..."))
  
  #### read in data ------------------------------------------------------------
  
  #norm.test.list <- read_rds(norm.test.files[seed_index])
  norm.train.list <- read_rds(norm.train.files[seed_index])
  sample.df <- read.delim(sample.files[seed_index])
  
  # convert gene names column to row names
  # convert GBM gene names to SYMBOL
  convert_row_names <- function(expr, cancer_type){
    if (cancer_type == "GBM") {
      expr <- expr %>%
        mutate(gene = ensembldb::select(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
                                        keys = as.character(gene),
                                        keytype = "GENEID",
                                        columns = "SYMBOL")$SYMBOL)
    }
    column_to_rownames(expr,
                       var = "gene")
  }
  
  norm.train.list <- purrr::modify_depth(norm.train.list, 2,
                                         function(x) convert_row_names(expr = x,
                                                                       cancer_type = cancer_type))
  
  #### get common gene set -----------------------------------------------------
  
  common.genes <- PLIER::commonRows(all.paths,
                                    norm.train.list[["0"]][["z"]])
  
  #### main --------------------------------------------------------------------
  
  # create an output list
  plier_results_list <- list()
  
  # parallel backend
  cl <- parallel::makeCluster(detectCores() - 1)
  doParallel::registerDoParallel(cl)
  
  # at each titration level (0-100% RNA-seq)
  #perc_seq <- as.character(seq(0, 100, 10))
  #norm_methods <- c("log", "npn", "qn", "tdm", "z")
  perc_seq <- as.character(seq(0, 100, 50))
  norm_methods <- c("log", "z")
  plier_results_list <- foreach(ps = perc_seq, .packages = c("PLIER", "doParallel"), .export = c("check_all_same")) %dopar% {
    foreach(nm = norm_methods, .packages = c("PLIER", "doParallel"), .export = c("check_all_same")) %dopar% {
      
      if (nm %in% names(norm.train.list[[ps]])) {
        if(any(apply(norm.train.list[[ps]][[nm]], 1, check_all_same))) {
          c("Some rows all same value...")
        } else {
          # minimum k for PLIER = 2*num.pc
          set.k <- 2*PLIER::num.pc(norm.train.list[[ps]][[nm]][common.genes, ])
          
          # PLIER main function
          PLIER::PLIER(as.matrix(norm.train.list[[ps]][[nm]][common.genes, ]),
                       all.paths[common.genes, ],
                       k = set.k,
                       trace = FALSE,
                       scale = TRUE)
        }
      } else {
        c("No data for this combination...")
      }
    }
  }
  
  # stop parallel backend
  parallel::stopCluster(cl)
  
  # renames list levels
  names(plier_results_list) <- perc_seq
  for (i in perc_seq) {
    names(plier_results_list[[i]]) <- norm_methods
  }
  
  # TODO NOW DO COMPARISON METRIC  

}

# TODO PLOT THAT
