# J. Taroni Jan 2017
# The purpose of this analysis is to identify differentially expressed genes
# between one subtype, specified by the user, and all
# other subtypes using the limma package for varying amounts of RNA-seq data
# (0-100%, 10% added at a time; termed 'RNA-seq titration') and normalization
# methods. It takes RNA-seq and microarray data from matched samples as input,
# and performs RNA-seq titration and differential expression analysis.
#
# USAGE: Rscript 1A-detect_differentially_expressed_genes.R --cancer_type --subtype_vs_others --subtype_vs_subtype --seed

option_list <- list(
  optparse::make_option("--cancer_type",
                        default = NULL,
                        help = "Cancer type"),
  optparse::make_option("--subtype_vs_others",
                        default = NULL,
                        help = "Subtype used for comparison against all others"),
  optparse::make_option("--subtype_vs_subtype",
                        default = NULL,
                        help = "Subtypes used in head-to-head comparison (comma-separated without space e.g. Type1,Type2)"),
  optparse::make_option("--seed",
                        default = 98,
                        help = "Random seed [default: %default]")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
source("util/option_functions.R")
check_options(opt)

# load libraries
suppressMessages(source("load_packages.R"))
source(file.path("util", "normalization_functions.R"))
source(file.path("util", "differential_expression_functions.R"))

# set options
cancer_type <- opt$cancer_type
subtype_vs_others <- opt$subtype_vs_others
subtype_vs_subtype <- opt$subtype_vs_subtype
# really this could be any number of subtypes
two_subtypes <- as.vector(stringr::str_split(subtype_vs_subtype, pattern = ",", simplify = TRUE))

# set seed
initial.seed <- opt$seed
set.seed(initial.seed)
message(paste("\nInitial seed set to:", initial.seed))

# define directories
data.dir <- "data"
deg.dir <- file.path("results", "differential_expression")

# define input files
seq.file <- file.path(data.dir, paste0(cancer_type, "RNASeq_matchedOnly_ordered.pcl"))
array.file <- file.path(data.dir, paste0(cancer_type, "array_matchedOnly_ordered.pcl"))
smpl.file <- file.path("results",
                       list.files("results", # this finds the first example of a subtypes file from cancer_type
                                  pattern = paste0(cancer_type, # and does not rely on knowing a seed
                                                   "_matchedSamples_subtypes_training_testing_split_labels_"))[1])

# define output files
subtype_vs_others.rds <- file.path(deg.dir,
                                   paste0(cancer_type,
                                          "_titration_differential_exp_eBayes_fits_",
                                          subtype_vs_others, "vOther.RDS"))
two_subtypes.rds <- file.path(deg.dir,
                              paste0(cancer_type,
                                     "_titration_differential_exp_eBayes_fits_",
                                     stringr::str_c(two_subtypes, collapse = "v"), ".RDS"))
norm.rds <- file.path("normalized_data",
                      paste0(cancer_type, "_titration_no_ZTO_transform_with_UN.RDS"))

#### read in data --------------------------------------------------------------

seq.data <- data.table::fread(seq.file, data.table = F)
array.data <- data.table::fread(array.file, data.table = F)
sample.df <- read.delim(smpl.file)

# check that subtypes are in sample.df
for(subtype in c(subtype_vs_others, two_subtypes)) {
  if (!(subtype %in% sample.df$subtype)) {
    stop(paste("Subtype", subtype, "not found in sample file",
               smpl.file, "in 1A-detect_differentially_expressed_genes.R."))
  }
}

sample.names <- sample.df$sample

#### RNA-seq 'titration' -------------------------------------------------------

titration.seed <- sample(1:10000, 1)
message(paste("Random seed for titration:",
              titration.seed), appendLF=TRUE)

set.seed(titration.seed)
# these samples will be the RNA-seq samples in any given 'titration' experiment
titrate.sample.list <-
  lapply(seq(0, 1, by = 0.1),
         function(x) GetTitratedSampleNames(sample.names, x))

# remove rows that are equal to all ones in sequencing data -- these are
# essentially missing values and cause issues with z-transformation
seq.dt.list <-
  lapply(titrate.sample.list,
         function(x) seq.data[, c(1, which(colnames(seq.data) %in% x))])
all.same.list <- lapply(seq.dt.list[2:11],
                    function(x){
                      vals <- x[, 2:ncol(x)]
                      indx <- which(apply(vals, 1, all_same))
                      return(indx)
                    } )
all.same.indx <- unique(unlist(all.same.list))
# if no rows are all same (in previous lapply), all.same.indx is integer(0)
# subsetting data frames by -integer(0) results in no rows
# so check that integer vector has length > 0 before subsetting
if (length(all.same.indx) > 0) {
  array.data <- array.data[-all.same.indx, ]
  seq.data <- seq.data[-all.same.indx, ]  
}

# get a list that contains an array data.table and seq data.table for each
# each level of 'titration'
titrate.mix.dt.list <-
  lapply(titrate.sample.list,
         function(x) GetDataTablesForMixing(data.table(array.data),
                                            data.table(seq.data),
                                            x))
names(titrate.mix.dt.list) <- as.character(seq(0, 100, by=10))

#### normalize data ------------------------------------------------------------

# initialize in the list to hold normalized data
norm.titrate.list <- list()

# single platform array normalization
norm.titrate.list[["0"]] <-
  SinglePlatformNormalizationWrapper(titrate.mix.dt.list[[1]]$array,
                                     platform = "array",
                                     zto = FALSE,
                                     add.qn.z = TRUE)

# parallel backend
cl <- parallel::makeCluster(detectCores() - 1)
doParallel::registerDoParallel(cl)

# 'mixed' both platform normalization
norm.titrate.list[2:10] <-
  foreach(n = 2:10) %dopar% {
    NormalizationWrapper(titrate.mix.dt.list[[n]]$array,
                         titrate.mix.dt.list[[n]]$seq,
                         zto = FALSE,
                         add.untransformed = TRUE,
                         add.qn.z = TRUE)
  }
names(norm.titrate.list)[2:10] <- names(titrate.mix.dt.list)[2:10]

# stop parallel backend
parallel::stopCluster(cl)

# single platform seq normalization
norm.titrate.list[["100"]] <-
  SinglePlatformNormalizationWrapper(titrate.mix.dt.list[[11]]$seq,
                                     platform = "seq",
                                     zto = FALSE,
                                     add.untransformed = TRUE,
                                     add.qn.z = TRUE)

# save normalized data
saveRDS(norm.titrate.list, file = norm.rds)

#### Subtype v. Others  --------------------------------------------------------
# design matrices
design.mat.list <- GetDesignMatrixList(norm.titrate.list, sample.df,
                                       subtype = subtype_vs_others)
# differential expression
fit.results.list <- GetFiteBayesList(norm.list = norm.titrate.list,
                                     design.list = design.mat.list)
# save fit results to RDS
saveRDS(fit.results.list, file = subtype_vs_others.rds)

#### Subtype v. Subtype --------------------------------------------------------
# remove all samples that are not in these subtypes
samples.to.keep <-
  sample.df$sample[which(sample.df$subtype %in% two_subtypes)]

pruned.norm.list <-
  lapply(norm.titrate.list,
         function(x) lapply(x,
                            function(y) y[,
                                          c(1, which(colnames(y) %in%
                                                       samples.to.keep)),
                                          with = FALSE]))

# get design matrices
last_subtype.design.list <- GetDesignMatrixList(pruned.norm.list,
                                                sample.df,
                                                subtype = last(two_subtypes))
# differential expression
last_subtype.fit.results.list <- GetFiteBayesList(norm.list = pruned.norm.list,
                                                  design.list = last_subtype.design.list)

# save fit results to file
saveRDS(last_subtype.fit.results.list, file = two_subtypes.rds)
