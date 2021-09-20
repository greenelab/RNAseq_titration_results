# J. Taroni Jun 2016
# The purpose of this script is to read in TGCA array and sequencing data,
# to preprocess leaving only overlapping genes and samples with complete
# category information, and to split the data into training and testing sets
# It should be run from the command line through the run_experiments.R script

option_list <- list(
  optparse::make_option("--cancer_type",
                        default = NA_character_,
                        help = "Cancer type"),
  optparse::make_option("--predictor",
                        default = NA_character_,
                        help = "Predictor used"),
  optparse::make_option("--seed1",
                        default = NA_integer_,
                        help = "Random seed"),
  optparse::make_option("--null_model",
                        action = "store_true",
                        default = FALSE,
                        help = "Permute dependent variable (within subtype if predictor is a gene)")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
source(here::here("util/option_functions.R"))
check_options(opt)

# load libraries
suppressMessages(source(here::here("load_packages.R")))

# set options
cancer_type <- opt$cancer_type
predictor <- opt$predictor
null_model <- opt$null_model
file_identifier <- ifelse(null_model,
                          str_c(cancer_type, predictor, "null", sep = "_"),
                          str_c(cancer_type, predictor, sep = "_"))

# set seed
initial.seed <- as.integer(opt$seed1)
set.seed(initial.seed)
# set seed for spliting into train/test here, before null_model scramble
split.seed <- sample(1:10000, 1)

# define directories
data.dir <- here::here("data")
plot.dir <- here::here("plots")
res.dir <- here::here("results")

# name input files
seq.exprs.filename <- paste0(cancer_type, "RNASeq.pcl")
array.exprs.filename <- paste0(cancer_type, "array.pcl")
clin.filename <- paste0("combined_clinical_data.", cancer_type, ".tsv")

# name output files
category.distribtion.plot <- paste0(file_identifier,
                                    "_dist_split_stacked_bar_",
                                    initial.seed, ".pdf")
train.test.labels <- paste0(file_identifier,
                            "_matchedSamples_training_testing_split_labels_",
                            initial.seed, ".tsv")

#### read in expression and clinical data --------------------------------------

# read in expression data as data.frame
seq.data <- fread(file.path(data.dir, seq.exprs.filename),
                  data.table = FALSE)
array.data <- fread(file.path(data.dir, array.exprs.filename),
                    data.table = FALSE)
clinical <- fread(file.path(data.dir, clin.filename),
                  data.table = FALSE)

# filter clinical data to keep tumor samples with complete data
# if the predictor is subtype, we only select subtype (twice, but once)
# if the predictor is a gene, we select subtype and the gene
# this ensures downstream mutation predictions will have subtype available as covariate
clinical <- clinical %>%
  select(Sample, Type, "subtype", all_of(predictor)) %>%
  rename("category" = all_of(predictor)) %>%
  filter(Type == "tumor") %>%
  tidyr::drop_na()

# change first column name to "gene"
colnames(array.data)[1] <- colnames(seq.data)[1] <- "gene"

# remove tumor-adjacent samples from the array data set
array.tumor.smpls <- clinical$Sample
array.tumor.smpls <- substr(array.tumor.smpls, 1, 15)

array.category <- clinical$category

# filter array data only to include tumor samples
array.data <- array.data[, c(1, which(colnames(array.data) %in%
                                        array.tumor.smpls))]

# what are the overlapping sample names -- "matched" samples?
# includes "gene" column
sample.overlap <- intersect(colnames(array.data), colnames(seq.data))

# what are the overlapping genes between the two platforms?
gene.overlap <- intersect(array.data$gene, seq.data$gene)

# filter the expression data for matched samples and overlapping genes
array.matched <- array.data[which(array.data$gene %in% gene.overlap),
                            sample.overlap]
seq.matched <- seq.data[which(seq.data$gene %in% gene.overlap),
                        sample.overlap]

# reorder genes on both platforms
array.matched <- array.matched[order(array.matched$gene), ]
seq.matched <- seq.matched[order(seq.matched$gene), ]

# reorder samples on both platforms
array.matched <- array.matched[, c(1, (order(colnames(array.matched)[-1]) + 1))]
seq.matched <- seq.matched[, c(1, (order(colnames(seq.matched)[-1]) + 1))]

# check reording sample names worked as expected
if (any(colnames(array.matched) != colnames(seq.matched))) {
  stop("Column name reordering did not work as expected in 0-expression_data_overlap_and_split.R")
}

# keep category labels for samples with expression data
array.category <- as.factor(array.category[which(array.tumor.smpls %in%
                                                   colnames(array.matched))])

array.tumor.smpls <- array.tumor.smpls[which(array.tumor.smpls %in%
                                               colnames(array.matched))]

# remove "unmatched" / "raw" expression data
rm(array.data, seq.data)

# write matched only samples to pcl files
array.output.nm <- sub(".pcl", "_matchedOnly_ordered.pcl", array.exprs.filename)
array.output.nm <- file.path(data.dir, array.output.nm)
write.table(array.matched, file = array.output.nm,
            row.names = FALSE, quote = FALSE, sep = "\t")

seq.output.nm <- sub(".pcl", "_matchedOnly_ordered.pcl", seq.exprs.filename)
seq.output.nm <- file.path(data.dir, seq.output.nm)
write.table(seq.matched, file = seq.output.nm,
            row.names = FALSE, quote = FALSE, sep = "\t")

#### split data into balanced training and testing sets ------------------------

# order array category to match the expression data order
array.category <- array.category[order(array.tumor.smpls)]

message(paste("\nRandom seed for splitting into testing and training:",
              split.seed), appendLF = TRUE)

set.seed(split.seed)
train.index <- unlist(createDataPartition(array.category, times = 1, p = (2/3)))

#### write training/test labels to file ----------------------------------------

lbl <- rep("test", length(array.tumor.smpls))
lbl[train.index] <- "train"
lbl.df <- tibble(sample = colnames(array.matched)[2:ncol(array.matched)],
                 split = lbl,
                 category = as.character(array.category))

# add back subtype if predicting gene
if (predictor != "subtype") {
  lbl.df <- lbl.df %>% 
    left_join(clinical %>%
                select(Sample, subtype),
              by = c("sample" = "Sample"))
}

#### permute category labels for null model ------------------------------------
# this comes after createDataPartition() to ensure same samples go to train/test
# grouping by split ensure labels remain balanced within train and test
# if null_model is specified and predicting subtype, permute subtype labels
# if null_model is specified and predicting mutation status,
#   permute mutation labels WITHIN subtype

if (null_model) {
  if (predictor == "subtype") { # here, subtype = category
    lbl.df <- lbl.df %>%
      group_by(split) %>%
      mutate(category = sample(category)) %>%
      ungroup()
  } else { # if predictor not subtype, then must be mutation
    lbl.df <- lbl.df %>% # subtype = subtype, category = TP53 or PIK3CA 0/1
      group_by(split, subtype) %>% # sample within subtype
      mutate(category = sample(category)) %>%
      ungroup()
  }
}

write.table(lbl.df,
            file = file.path(res.dir, train.test.labels),
            quote = FALSE, sep = "\t", row.names = FALSE)

#### plot category distributions ------------------------------------------------
cbPalette <- c("#000000", "#E69F00", "#56B4E9",
               "#009E73", "#F0E442","#0072B2", "#D55E00", "#CC79A7")

plot.df <- lbl.df %>%
  mutate(split = case_when(split == "train" ~ "train (2/3)",
                           split == "test" ~ "test (1/3)")) %>%
  bind_rows(lbl.df %>% mutate(split = "whole"))

plot.nm <- file.path(plot.dir, category.distribtion.plot)
ggplot(plot.df, aes(x = split, fill = category)) +
  geom_bar() +
  theme_classic() +
  scale_fill_manual(values = cbPalette) +
  ggsave(plot.nm,
         height = 6,
         width = 6)
