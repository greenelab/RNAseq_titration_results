# Steven Foltz Nov 2021
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

#### Function for converting column to row names -------------------------------

convert_row_names <- function(expr, cancer_type){
  # If the cancer type is GBM, convert ENSG to gene symbols
  # then convert the gene column to rownames
  #
  # Inputs: gene expression matrix with genes in first column, and cancer type
  # Returns: modified gene expression matrix with gene names as row names
  
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

#### Function to get jaccard values from a list of PLIER results ---------------

return_plier_jaccard <- function(test_PLIER, array_silver, seq_silver){
  # Given a set of PLIER results (which is a list), compare significant pathways
  # to two silver sets of pathways defined by array and RNA-seq data only
  # Jaccard similaritiy is defined as O(intersect)/O(union).
  #
  # Inputs: PLIER result, pathway set 1, pathway set 2
  # Returns: data frame with two rows (array, seq) with stats for each overlap
  
  if (is.list(test_PLIER)) {
    
    test_genes <- test_PLIER[["summary"]] %>%
      filter(FDR < 0.05) %>%
      pull(pathway) %>%
      unique()
    
    array_jaccard <- length(intersect(array_silver, test_genes))/length(union(array_silver, test_genes))
    seq_jaccard <- length(intersect(seq_silver, test_genes))/length(union(seq_silver, test_genes))  
    
    data.frame(silver = c("array", "seq"),
               n_silver = c(length(array_silver),
                            length(seq_silver)),
               n_test = length(test_genes),
               n_intersect = c(length(intersect(array_silver, test_genes)),
                               length(intersect(seq_silver, test_genes))),
               n_union = c(length(union(array_silver, test_genes)),
                           length(union(seq_silver, test_genes))),
               n_common_genes = nrow(test_PLIER[["Z"]]),
               k = ncol(test_PLIER[["Z"]]),
               jaccard = c(array_jaccard,
                           seq_jaccard))
    
    
  }
  
}

#### loop over data for each seed and get PLIER results ------------------------

jaccard_list <- list()

for(seed_index in 1:length(norm.train.files)) {
  
  message(str_c("PLIER with data seed", seed_index,
                "out of", length(norm.train.files), "...",
                sep = " "))
  
  #### read in data ------------------------------------------------------------
  
  #norm.test.list <- read_rds(norm.test.files[seed_index])
  norm.train.list <- read_rds(norm.train.files[seed_index])
  sample.df <- read.delim(sample.files[seed_index])
  
  # convert gene names column to row names
  # if GBM, also convert from GENEID to SYMBOL
  norm.train.list <- purrr::modify_depth(norm.train.list, 2,
                                         function(x) convert_row_names(expr = x,
                                                                       cancer_type = cancer_type))
  
  #### main --------------------------------------------------------------------
  
  # create an output list
  plier_results_list <- list()
  
  # parallel backend
  cl <- parallel::makeCluster(detectCores() - 1)
  doParallel::registerDoParallel(cl)
  
  # at each titration level (0-100% RNA-seq)
  # TODO remove test conditions
  #perc_seq <- as.character(seq(0, 100, 10))
  #norm_methods <- c("log", "npn", "qn", "tdm", "z")
  perc_seq <- as.character(seq(0, 100, 50))
  norm_methods <- c("log", "z")
  plier_results_list <- foreach(ps = perc_seq,
                                .packages = c("PLIER", "doParallel")) %do% { #par% {
    foreach(nm = norm_methods) %do% { #par% {
      
      if (nm %in% names(norm.train.list[[ps]])) {
        
        message(str_c(seed_index, ps, nm, sep = " "))
        
        # remove any rows with all the same value
        all.same.indx <- which(apply(norm.train.list[[ps]][[nm]], 1,
                                     check_all_same))
        message(length(all.same.indx))
        message(nrow(norm.train.list[[ps]][[nm]]))
        if (length(all.same.indx) > 0) {
          norm.train.list[[ps]][[nm]] <- norm.train.list[[ps]][[nm]][-all.same.indx, ]
        }
        message(nrow(norm.train.list[[ps]][[nm]]))
        
        # get common genes
        common.genes <- PLIER::commonRows(all.paths,
                                          norm.train.list[[ps]][[nm]])      
        message(length(common.genes))
        
        # minimum k for PLIER = 2*num.pc
        set.k <- 2*PLIER::num.pc(norm.train.list[[ps]][[nm]][common.genes, ])
        # TODO alternatively, should we just set one k for all data sets?
        #set.k <- 50 # set k the be the same arbitrary value for all runs
        message(set.k)
        message(nrow(as.matrix(norm.train.list[[ps]][[nm]][common.genes, ])))
        message(min(apply(norm.train.list[[ps]][[nm]], 1, sd)))
        message(nrow(all.paths[common.genes, ]))
        # PLIER main function
        
        tryCatch(expr = PLIER::PLIER(as.matrix(norm.train.list[[ps]][[nm]][common.genes, ]),
                                     all.paths[common.genes, ],
                                     k = set.k,
                                     scale = FALSE),
                 error = function(err) NA)
        
      } else {
        
        NA # return NA for easy check later 
        
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
  
  
  print(plier_results_list)
  
  # TODO remove this test: write out plier results
  readr::write_rds(x = plier_results_list,
                   path = str_c("plier_results_list.", seed_index, ".RDS"))
  
  # Jaccard comparison metric to array and seq silver standards
  # TODO what are best settings for silver standard? PLIER expects z-scored
  array_silver <- plier_results_list[["0"]][["z"]][["summary"]] %>%
    filter(FDR < 0.05) %>%
    pull(pathway) %>%
    unique()
  seq_silver <- plier_results_list[["100"]][["z"]][["summary"]] %>%
    filter(FDR < 0.05) %>%
    pull(pathway) %>%
    unique()
  
  # Check that silver standard pathways have non-zero length
  if (length(array_silver) > 0 & length(seq_silver) > 0) {
  
    jaccard_list[[seed_index]] <- purrr::modify_depth(plier_results_list, 2,
                                                      function(x) return_plier_jaccard(x, array_silver, seq_silver))
    
  } else {
    
    message(str_c("Silver standard array or seq significant pathways has non-zero length",
                  seed_index, percent_seq, normalization_method,
                  length(array_silver), length(seq_silver),
                  sep = " "))
    
  }
}

jaccard_df <- reshape2::melt(data = jaccard_list,
                             id.vars = c("silver", "n_silver", "n_test", "n_intersect", "n_union", "n_common_genes", "k"),
                             value.name = "jaccard") %>%
  rename("nmeth" = "L3",
         "pseq" = "L2",
         "seed_index" = "L1")

readr::write_tsv(x = jaccard_df,
                 path = here::here("test.tsv"))

# TODO PLOT THAT
