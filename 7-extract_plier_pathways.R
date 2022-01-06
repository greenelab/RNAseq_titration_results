# Steven Foltz Nov 2021
# The purpose of this analysis is to use PLIER to identify expression pathways
# in our data, and quantify the rate of return for significant PLIER pathways
# for data coming from different normalization methods and titration levels.
#
# USAGE: Rscript 7-extract_plier_pathways.R --cancer_type --seed

option_list <- list(
  optparse::make_option("--cancer_type",
                        default = NA_character_,
                        help = "Cancer type"
  ),
  optparse::make_option("--seed",
                        default = 8934,
                        help = "Random seed"
  )
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
message(paste("\nPLIER initial seed set to:", initial.seed))

# define directories
data.dir <- here::here("data")
norm.data.dir <- here::here("normalized_data")
res.dir <- here::here("results")

# define input files
# norm.test.files <- file.path(norm.data.dir,
#                            list.files(norm.data.dir,
#                                       pattern = paste0(file_identifier,
#                                                        "_array_seq_test_data_normalized_list_")))
norm.train.files <- file.path(
  norm.data.dir,
  list.files(norm.data.dir,
             pattern = paste0(
               file_identifier,
               "_array_seq_train_titrate_normalized_list_"
             )
  )
)
# sample.files <- file.path(res.dir,
#                         list.files(res.dir,
#                                    pattern = paste0(file_identifier,
#                                                     "_matchedSamples_training_testing_split_labels_")))

#### set up PLIER data ---------------------------------------------------------

data(bloodCellMarkersIRISDMAP)
data(canonicalPathways)
data(oncogenicPathways)
data(svmMarkers)

all.paths <- PLIER::combinePaths(
  bloodCellMarkersIRISDMAP,
  canonicalPathways,
  oncogenicPathways,
  svmMarkers
)
PLIER_pathways <- colnames(all.paths)

#### Function for converting column to row names -------------------------------

convert_row_names <- function(expr, cancer_type) {
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
                                      columns = "SYMBOL"
      )$SYMBOL)
  }
  
  column_to_rownames(expr,
                     var = "gene"
  )
}

#### Function to check if a PLIER model converged or not

check_failure_to_converge <- function(plier_result) {

  # if the plier result contained an error message
  if ("message" %in% names(plier_result)) {
    # and if that error massage refers to convergence failure
    if (str_detect(plier_result$message, "system is computationally singular")) {
      NA # return NA
    } else { # PLIER failed for another reason and we need to know about that
      stop("PLIER run failed for reason other than system is computationally singular")
    }
  } else {
    plier_result # return the plier result as is
  }
}

#### Functions to get jaccard values from a list of PLIER results ---------------

return_plier_jaccard_silver <- function(test_PLIER, array_silver, seq_silver) {
  # Given a set of PLIER results (which is a list), compare significant pathways
  # to two silver sets of pathways defined by array and RNA-seq data only
  # Jaccard similaritiy is defined as O(intersect)/O(union). If the input is not
  # a properly completed PLIER result, then return all NAs.
  # 
  # Inputs: PLIER result, pathway set 1, pathway set 2
  # Returns: data frame with two rows (array, seq) with stats for each overlap
  
  if ("summary" %in% names(test_PLIER)) {
    
    test_pathways <- test_PLIER[["summary"]] %>%
      filter(FDR < 0.05) %>%
      pull(pathway) %>%
      unique()
    
    array_jaccard <- length(intersect(array_silver, test_pathways)) / length(union(array_silver, test_pathways))
    seq_jaccard <- length(intersect(seq_silver, test_pathways)) / length(union(seq_silver, test_pathways))
    
    data.frame(
      silver = c("array", "seq"),
      n_silver = c(
        length(array_silver),
        length(seq_silver)
      ),
      n_test = length(test_pathways),
      n_intersect = c(
        length(intersect(array_silver, test_pathways)),
        length(intersect(seq_silver, test_pathways))
      ),
      n_union = c(
        length(union(array_silver, test_pathways)),
        length(union(seq_silver, test_pathways))
      ),
      n_common_genes = nrow(test_PLIER[["Z"]]),
      k = ncol(test_PLIER[["Z"]]),
      jaccard = c(
        array_jaccard,
        seq_jaccard
      )
    )
  } else {
    
    data.frame(
      silver = NA,
      n_silver = NA,
      n_test = NA,
      n_intersect = NA,
      n_union = NA,
      n_common_genes = NA,
      k = NA,
      jaccard = NA
    )  
    
  }
}

return_plier_jaccard_global <- function(test_PLIER, global_pathways) {
  # Given a set of PLIER results (which is a list), compare significant pathways
  # to a global set of pathways defined by the existing PLIER pathways
  # Jaccard similaritiy is defined as O(intersect)/O(union). If the input is not
  # a properly completed PLIER result, then return all NAs.
  #
  # Inputs: PLIER result, global pathways
  # Returns: data frame with one row with stats for each overlap
  
  if ("summary" %in% names(test_PLIER)) {
    
    test_pathways <- test_PLIER[["summary"]] %>%
      filter(FDR < 0.05) %>%
      pull(pathway) %>%
      unique()
    
    global_jaccard <- length(intersect(global_pathways, test_pathways)) / length(global_pathways)
    
    data.frame(
      n_global = length(global_pathways),
      n_test = length(test_pathways),
      n_intersect = length(intersect(global_pathways, test_pathways)),
      n_common_genes = nrow(test_PLIER[["Z"]]),
      k = ncol(test_PLIER[["Z"]]),
      jaccard = global_jaccard
    ) 
    
  } else {
    
    data.frame(
      n_global = NA,
      n_test = NA,
      n_intersect = NA,
      n_common_genes = NA,
      k = NA,
      jaccard = NA
    )
    
  }
}

#### loop over data for each seed and get PLIER results ------------------------

jaccard_list <- list()

for (seed_index in 1:length(norm.train.files)) {
  message(str_c("PLIER with data seed", seed_index,
                "out of", length(norm.train.files), "...",
                sep = " "
  ))
  
  #### read in data ------------------------------------------------------------
  
  # norm.test.list <- read_rds(norm.test.files[seed_index])
  norm.train.list <- read_rds(norm.train.files[seed_index])
  # sample.df <- read.delim(sample.files[seed_index])
  
  # convert gene names column to row names
  # if GBM, also convert from GENEID to SYMBOL
  norm.train.list <- purrr::modify_depth(
    norm.train.list, 2,
    function(x) {
      convert_row_names(
        expr = x,
        cancer_type = cancer_type
      )
    }
  )
  
  #### main --------------------------------------------------------------------
  
  # create an output list
  plier_results_list <- list()
  
  # parallel backend
  cl <- parallel::makeCluster(detectCores() - 1)
  doParallel::registerDoParallel(cl)
  
  # at each titration level (0-100% RNA-seq)
  perc_seq <- as.character(seq(0, 100, 10))
  norm_methods <- c("log", "npn", "qn", "qn-z", "tdm", "un", "z")
  #perc_seq <- as.character(seq(40, 50, 10))
  #norm_methods <- c("un", "z")
  plier_results_list <- foreach(
    ps = perc_seq,
    .packages = c("PLIER", "doParallel")
  ) %dopar% {
    foreach(
      nm = norm_methods,
      .errorhandling = "pass" # let pass on inside loop
    ) %dopar% {
      if (nm %in% names(norm.train.list[[ps]])) {
        
        # remove any rows with all the same value
        all.same.indx <- which(apply(
          norm.train.list[[ps]][[nm]], 1,
          check_all_same
        ))
        if (length(all.same.indx) > 0) {
          norm.train.list[[ps]][[nm]] <- norm.train.list[[ps]][[nm]][-all.same.indx, ]
        }
        
        # get common genes
        common.genes <- PLIER::commonRows(
          all.paths,
          norm.train.list[[ps]][[nm]]
        )
        
        # minimum k for PLIER = 2*num.pc
        set.k <- 2 * PLIER::num.pc(PLIER::rowNorm(norm.train.list[[ps]][[nm]][common.genes, ]))
        
        # PLIER main function
        PLIER::PLIER(as.matrix(norm.train.list[[ps]][[nm]][common.genes, ]),
                     all.paths[common.genes, ],
                     k = set.k,
                     scale = TRUE # PLIER z-scores input values by row
        )
      } else {
        NULL # return NULL for empty result; purrr will ignore this list element
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
  
  # write test file
  
  write_rds(x = plier_results_list,
            path = str_c("plier.", seed_index, ".rds"))
  
  # Check for failure to converge, and set to NA
  
  plier_results_list <- purrr::modify_depth(plier_results_list, 2,
                                            check_failure_to_converge
  )
  
  # Return pathway comparison for appropriate level of PLIER results list
  jaccard_list[[seed_index]] <- purrr::modify_depth(
    plier_results_list, 2,
    function(x) return_plier_jaccard_global(x, PLIER_pathways)
  )
}

if (length(jaccard_list) > 0) {
  
  # melt jaccard list elements into one data frame
  jaccard_df <- reshape2::melt(
    data = jaccard_list,
    id.vars = c(
      "n_global", "n_test", "n_intersect",
      "n_common_genes", "k"
    ),
    value.name = "jaccard"
  ) %>%
    rename(
      "nmeth" = "L3", # normalization method
      "pseq" = "L2", # percentage RNA-seq
      "seed_index" = "L1"
    )
  
  readr::write_tsv(
    x = jaccard_df,
    path = file.path(
      res.dir,
      str_c(file_identifier, "_PLIER_jaccard.tsv")
    )
  )
  
  # TODO PLOT THAT
}
