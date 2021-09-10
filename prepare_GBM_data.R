# Script prepares GBM expression data for use in pipelines
# For GBM array data: convert sample names and remove duplicate individuals
# For GBM seq data: filter all seq data for GBM samples only, convert genes IDs

# Steven Foltz July 2021

option_list <- list(
  optparse::make_option("--seq_input",
                        default = NA_character_,
                        help = "TCGA sequencing expression input file path"),
  optparse::make_option("--array_input",
                        default = NA_character_,
                        help = "refine.bio microarray expression input file path"),
  optparse::make_option("--metadata_input",
                        default = NA_character_,
                        help = "refine.bio aggregated metadata JSON file path"),
  optparse::make_option("--array_output",
                        default = NA_character_,
                        help = "Processed microarray data output file path"),
  optparse::make_option("--seq_output",
                        default = NA_character_,
                        help = "Processed sequencing data output file path"),
  optparse::make_option("--clinical_input",
                        default = NA_character_,
                        help = "Clinical information input file path (Excel file)"),
  optparse::make_option("--clinical_output",
                        default = NA_character_,
                        help = "Clinical information output file path (.tsv)"),
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
tcga_seq_expression_input_filepath <- opt$seq_input
gbm_array_expression_input_filepath <- opt$array_input
metadata_json_input_filepath <- opt$metadata_input
gbm_array_output_filepath <- opt$array_output
gbm_seq_output_filepath <- opt$seq_output
clinical_xlxs_input_filepath <- opt$clinical_input
clinical_tsv_output_filepath <- opt$clinical_output

################################################################################
# Array data
################################################################################

# read in refine.bio GBM array expression data
gbm_array_expression <- read_tsv(gbm_array_expression_input_filepath,
                                 col_types = cols(
                                   .default = col_double(),
                                   Gene = col_character()
                                 ))

# load up aggregated metadata json file
metadata_json <- jsonlite::fromJSON(metadata_json_input_filepath,
                                    simplifyVector = FALSE)

# accession IDs present in expression data
available_array_accession_ids <- colnames(gbm_array_expression)[-1]

# All about TCGA barcodes
# TCGA barcodes are defined in this format: TCGA-XX-YYYY-ZZ*
# XX is the two character tissue source site (TSS) (defines a combination of hospital system and cancer type -- so multiple TSSs map to single cancer type)
# YYYY is the four digit participant ID specific to a TSS (so TCGA-XX-YYYY defines an individual patient -- patients from different TSSs may have same participant ID)
# ZZ is the sample type (starts with 0 for tumor samples; specifically 01 for primary solid tumors)
# After ZZ (*) is more specific information not relevant in this context
# Links to more info:
# https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
# https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tissue-source-site-codes
# https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes

# starting with flattened metadata, parse out raw TCGA IDs and filter for tumors
all_array_tumor_samples <- tibble(accession = names(metadata_json$samples)) %>%
  filter(accession %in% available_array_accession_ids) %>% # check these are ==
  rowwise() %>%
  mutate(tcga_id_raw = metadata_json$samples[[accession]]$refinebio_annotations[[1]]$characteristics_ch1[[9]] %>%
           str_remove("sample: ")) %>%
  # tcga_id_raw is the entire TCGA barcode (includes information after ZZ which we need to drop)
  mutate(tcga_id = str_sub(tcga_id_raw, 1, 15), # so tcga_id includes TCGA-XX-YYYY-ZZ
         sample = str_sub(tcga_id_raw, 14, 15)) %>% # and sample refers to ZZ
  filter(sample == "01") %>% # filter for primary solid tumors only
  ungroup()

# keep one (first) accession per TCGA ID
array_accession_tcga_id_keep <- all_array_tumor_samples %>%
  group_by(tcga_id) %>%
  summarize(accession = sort(accession)[1])
accession_colnames_keep <- colnames(gbm_array_expression)[-1][colnames(gbm_array_expression)[-1] %in% array_accession_tcga_id_keep$accession]

# select columns to keep and rename with TCGA IDs
gbm_array_expression_renamed <- gbm_array_expression %>%
  select(c("Gene",
           array_accession_tcga_id_keep$accession))
colnames(gbm_array_expression_renamed) <- c("sample",
                                            array_accession_tcga_id_keep$tcga_id)

################################################################################
# Sequencing data
################################################################################

# read in column names of entire TCGA seq expression file
tcga_seq_expression_column_names <- read_tsv(tcga_seq_expression_input_filepath,
                                             col_types = cols(
                                               .default = col_double(),
                                               gene_id = col_character()),
                                             n_max = 0) %>%
  names()

# identify sequencing TCGA IDs of samples present in array data
# (a more inclusive approach would be selecting GBM samples based on TSS codes)
gbm_seq_tumor_samples <- tibble(tcga_id_raw = tcga_seq_expression_column_names[-1]) %>%
  mutate(tcga_patient = str_sub(tcga_id_raw, 1, 12), # TCGA-XX-YYYY to match with clinical
         tcga_id = str_sub(tcga_id_raw, 1, 15), # as with array, tcga_id refers to TCGA-XX-YYYY-ZZ
         sample = str_sub(tcga_id_raw, 14, 15)) %>% # and sample is ZZ
  filter(sample == "01") %>% # require sample to be primary solid tumor
  filter(tcga_id %in% array_accession_tcga_id_keep$tcga_id) %>% # keep array GBMs
  group_by(tcga_patient, tcga_id) %>%
  summarize(tcga_id_raw = sort(tcga_id_raw)[1]) # keep one raw ID per person

# now read in GBM subset of entire TCGA seq expression file
# this is faster and uses less memory than reading in entire file and then subsetting
# read these GBM columns only
tcga_seq_gbm_tf <- tcga_seq_expression_column_names[-1] %in% gbm_seq_tumor_samples$tcga_id_raw
# use these column types
# first column is 'c' for gene_id, then '-' for non-GBM samples, then 'd' for GBM
# setting column type to '-' skips over that column when reading file
tcga_seq_gbm_col_types <- str_c(c("c", c("-", "d")[tcga_seq_gbm_tf + 1]), collapse = "")
# read in my defined subset of columns with column types
gbm_seq_expression <- read_tsv(tcga_seq_expression_input_filepath,
                               col_types = tcga_seq_gbm_col_types)
colnames(gbm_seq_expression) <- c("gene_id",
                                  gbm_seq_tumor_samples$tcga_id)

# Detour to make gene ids consistent between array and seq files
# will convert seq format (SYMBOL|ENTREZ) to array format (ENSG)

# separate gene symbols from entrez ids (delimiter = "|")
symbol_entrez_ids <- gbm_seq_expression %>%
  select(gene_id) %>%
  separate(gene_id,
           into = c("SYMBOL", "ENTREZID"),
           sep = "\\|",
           remove = FALSE)

# map entrez ids to ensembl ids (GENEID)
entrez_ensembl_ids <- ensembldb::select(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
                                        keys= symbol_entrez_ids$ENTREZID,
                                        keytype = "ENTREZID",
                                        columns = "GENEID") %>%
  as_tibble() %>%
  mutate(ENTREZID = as.character(ENTREZID))

# collate gene name schemes
# filter for those that mapped and exist in array
# filter for ENSGs with a one-to-one mapping with entrez
gene_id_mapping_in_array <- symbol_entrez_ids %>%
  left_join(entrez_ensembl_ids,
            by = "ENTREZID") %>%
  filter(!is.na(GENEID)) %>%
  add_count(GENEID) %>%
  filter(n == 1) %>%
  filter(GENEID %in% gbm_array_expression$Gene) %>%
  select(gene_id, GENEID)

# starting with acceptable genes, left join with seq expression and select cols
gbm_seq_expression_renamed <- gene_id_mapping_in_array %>%
  left_join(gbm_seq_expression,
            by = "gene_id") %>%
  select(-gene_id) %>%
  rename("sample" = "GENEID")

################################################################################
# subtype information
################################################################################

# read in Table S7 from flagship GBM landscape paper (Brennan et al., Cell 2013)
# select and rename interesting columns
gbm_subtypes <- readxl::read_xlsx(path = clinical_xlxs_input_filepath,
                                  sheet = "Clinical Data",
                                  skip = 1) %>%
  right_join(gbm_seq_tumor_samples,
             by = c("Case ID" = "tcga_patient")) %>%
  select("tcga_id",
         "MGMT Status",
         "G-CIMP\r\n methylation",
         "IDH1\r\n status",
         "Expression\r\nSubclass") %>%
  rename("Sample" = "tcga_id",
         "MGMT_methylation_status" = "MGMT Status",
         "G-CIMP_methylation" = "G-CIMP\r\n methylation",
         "IDH1_mutation_status" = "IDH1\r\n status",
         "subtype" = "Expression\r\nSubclass") %>%
  mutate(subtype = na_if(subtype, "NA")) %>%
  mutate(Type = "tumor")

missing_clinical <- gbm_subtypes %>%
  filter(is.na(subtype)) %>%
  pull(Sample)

################################################################################
# Write to file, excluding samples without clinical info
################################################################################

write_tsv(gbm_array_expression_renamed %>%
            select(-all_of(missing_clinical)),
          path = gbm_array_output_filepath)

write_tsv(gbm_seq_expression_renamed %>%
            select(-all_of(missing_clinical)),
          path = gbm_seq_output_filepath)

write_tsv(gbm_subtypes %>%
            filter(!is.na(subtype)),
          path = clinical_tsv_output_filepath)
