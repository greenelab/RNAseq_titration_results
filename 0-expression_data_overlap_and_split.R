# J. Taroni Jun 2016
# The purpose of this script is to read in TGCA array and sequencing data,
# to preprocess leaving only overlapping genes and samples with complete
# category information, and to split the data into training and testing sets
# It should be run from the command line through the run_experiments.R script

option_list <- list(
  optparse::make_option("--cancer_type",
                        default = NULL,
                        help = "Cancer type"),
  optparse::make_option("--predictor",
                        default = NULL,
                        help = "Predictor used"),
  optparse::make_option("--seed1",
                        default = NULL,
                        help = "Random seed")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
source(here::here("util/option_functions.R"))
check_options(opt)

# load libraries
suppressMessages(source(here::here("load_packages.R")))

# set options
cancer_type <- opt$cancer_type
predictor <- opt$predictor
file_identifier <- str_c(cancer_type, predictor, sep = "_")

# set seed
initial.seed <- as.integer(opt$seed1)
set.seed(initial.seed)

# define directories
data.dir <- here::here("data")
plot.dir <- here::here("plots")
res.dir <- here::here("results")

# name input files
seq.exprs.filename <- paste0(cancer_type, "RNASeq.pcl")
array.exprs.filename <- paste0(cancer_type, "array.pcl")
clin.filename <- paste0(cancer_type, "Clin.tsv")

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

if (cancer_type == "BRCA") { # rename from PAM50
  colnames(clinical)[4] <- "subtype"
}

# filter clinical data to keep tumor samples with complete data
# if the predictor is subtype, we only select subtype (twice, but once)
# if the predictor is a gene, we select subtype and the gene
# this ensures downstream mutation predictions will have subtype available as covariate
clinical <- clinical %>%
  select(Sample, Type, "subtype", all_of(predictor)) %>%
  rename("category" = predictor) %>%
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

split.seed <- sample(1:10000, 1)
message(paste("\nRandom seed for splitting into testing and training:",
              split.seed), appendLF = TRUE)

set.seed(split.seed)
train.index <- unlist(createDataPartition(array.category, times = 1, p = (2/3)))

#### plot category distributions ------------------------------------------------
whole.df <- cbind(as.character(array.category),
                  rep("whole", length(array.category)))
train.df <- cbind(as.character(array.category[train.index]),
                  rep("train (2/3)", length(train.index)))
test.df <- cbind(as.character(array.category[-train.index]),
                 rep("test (1/3)",
                     length(array.category)-length(train.index)))
mstr.df <- rbind(whole.df, train.df, test.df)

colnames(mstr.df) <- c("category", "split")
cbPalette <- c("#000000", "#E69F00", "#56B4E9",
               "#009E73", "#F0E442","#0072B2", "#D55E00", "#CC79A7")

plot.nm <- file.path(plot.dir, category.distribtion.plot)
ggplot(as.data.frame(mstr.df), aes(x = split, fill = category)) +
  geom_bar() +
  theme_classic() +
  scale_fill_manual(values = cbPalette) +
  ggsave(plot.nm,
         height = 6,
         width = 6)

#### write training/test labels to file ----------------------------------------

lbl <- rep("test", length(array.tumor.smpls))
lbl[train.index] <- "train"
lbl.df <- cbind(colnames(array.matched)[2:ncol(array.matched)],
                lbl,
                as.character(array.category))
colnames(lbl.df) <- c("sample", "split", "category")

write.table(lbl.df,
            file = file.path(res.dir, train.test.labels),
            quote = FALSE, sep = "\t", row.names = FALSE)
