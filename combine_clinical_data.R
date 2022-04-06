# Script combines clinical data from all cancer types to one data frame
# Clinical data includes: subtype and TP53/PIK3CA mutation status
# Steven Foltz August 2021

option_list <- list(
  optparse::make_option("--cancer_type",
                        default = NA_character_,
                        help = "Cancer type"),
  optparse::make_option("--clinical_input",
                        default = NA_character_,
                        help = "Clinical information input file path (.tsv)"),
  optparse::make_option("--mutation_input",
                        default = NA_character_,
                        help = "Mutation input file path (.tsv)"),
  optparse::make_option("--combined_output",
                        default = NA_character_,
                        help = "Combined subtype and mutation output file path (.tsv)"),
  optparse::make_option("--overwrite",
                        action = "store_true",
                        default = FALSE,
                        help = "Overwrite existing output files [default: %default]")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
source(here::here("util/option_functions.R"))
check_options(opt)

# load libraries
suppressMessages(library(tidyverse))

# set options
cancer_type <- opt$cancer_type
clinical_input_filepath <- opt$clinical_input
mutation_input_filepath <- opt$mutation_input
combined_output_filepath <- opt$combined_output

################################################################################
# Read in clinical and mutation data
################################################################################
clinical_df <- read_tsv(clinical_input_filepath, # treat all columns equally
                        col_types = cols(.default = col_character())) %>%
  mutate(Sample = substr(Sample, 1, 15)) # remove extra parts of TCGA ID

mutation_df <- read_tsv(mutation_input_filepath, # treat all columns equally
                        col_types = cols(.default = col_character())) %>%
  mutate(tcga_id = substr(tcga_id, 1, 15)) # remove extra parts of TCGA ID

################################################################################
# Combine clinical and mutation data
################################################################################

# combine data frames with left_join() to get the left side of venn diagram
# start join with clinical_df because later scripts expect column name = Sample
combined_df <- clinical_df %>%
  left_join(mutation_df,
            by = c("Sample" = "tcga_id"))

################################################################################
# Save output file
################################################################################

write_tsv(combined_df,
          combined_output_filepath)
