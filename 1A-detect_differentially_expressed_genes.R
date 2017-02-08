# J. Taroni Jan 2017
# The purpose of this analysis is to identify differentially expressed genes
# between one BRCA subtype, specified by the user (default: Basal), and all
# other subtypes using the limma package for varying amounts of RNA-seq data
# (0-100%, 10% added at a time; termed 'RNA-seq titration') and normalization 
# methods. It takes RNA-seq and microarray data from matched samples as input, 
# and performs RNA-seq titration and differential expression analysis.
# 
# USAGE: Rscript 1A-detect_differentially_expressed_genes.R

suppressMessages(source("load_packages.R"))
source(file.path("util", "normalization_functions.R"))
source(file.path("util", "differential_expression_functions.R"))

args <- commandArgs(trailingOnly = TRUE)
initial.seed <- as.integer(args[1])
if (is.na(initial.seed)) {
  message("\nInitial seed set to default: 98")
  initial.seed <- 98
} else {
  message(paste("\nInitial seed set to:", initial.seed))
}
set.seed(initial.seed)

data.dir <- "data"
deg.dir <- file.path("results", "differential_expression")

seq.file <- file.path(data.dir, "BRCARNASeq_matchedOnly_ordered.pcl")
array.file <- file.path(data.dir, "BRCAarray_matchedOnly_ordered.pcl")
smpl.file <- 
  file.path("results", 
            "BRCA_matchedSamples_PAM50Array_training_testing_split_labels_3061.tsv")
basal.rds <- 
  file.path(deg.dir, 
            "BRCA_titration_differential_exp_eBayes_fits_BasalvOther.RDS")
her2.rds <- 
  file.path(deg.dir, 
            "BRCA_titration_differential_exp_eBayes_fits_Her2vLumA.RDS")
norm.rds <- file.path("normalized_data",
                      "BRCA_titration_no_ZTO_transform_with_UN.RDS")

#### read in data --------------------------------------------------------------

seq.data <- data.table::fread(seq.file, data.table = F)
array.data <- data.table::fread(array.file, data.table = F)
sample.df <- read.delim(smpl.file)

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

# remove rows that are equal to all ones in sequencing data
seq.dt.list <- 
  lapply(titrate.sample.list,
         function(x) seq.data[, c(1, which(colnames(seq.data) %in% x))])
all1.list <- lapply(seq.dt.list[2:11],
                    function(x){
                      vals <- x[, 2:ncol(x)]
                      indx <- which(apply(vals, 1, function(x) all(x == 1)))
                      return(indx)
                    } )
all.1.indx <- unique(unlist(all1.list))
array.data <- array.data[-all.1.indx, ]
seq.data <- seq.data[-all.1.indx, ]

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

#### Basal v. Other  -----------------------------------------------------------
# design matrices
design.mat.list <- GetDesignMatrices(norm.titrate.list, sample.df, 
                                     subtype = "Basal")
# differential expression
fit.results.list <- GetFiteBayesList(norm.list = norm.titrate.list,
                                     design.list = design.mat.list)
# save fit results to RDS
saveRDS(fit.results.list, file = basal.rds)

#### Her2 v. LumA --------------------------------------------------------------
# remove all samples that are not Her2 or LumA 
samples.to.keep <- sample.df$sample[which(sample.df$subtype == "LumA" |
                                            sample.df$subtype == "Her2")]
pruned.norm.list <- 
  lapply(norm.titrate.list, 
         function(x) lapply(x, 
                            function(y) y[, 
                                          c(1, which(colnames(y) %in% 
                                                          samples.to.keep)),
                                         with = FALSE]))

# get design matrices
her2.design.list <- GetDesignMatrices(pruned.norm.list,
                                          sample.df, subtype = "Her2")
# differential expression
her2.fit.results.list <- GetFiteBayesList(norm.list = pruned.norm.list,
                                          design.list = her2.design.list)

# save fit results to file
saveRDS(her2.fit.results.list, file = her2.rds)
