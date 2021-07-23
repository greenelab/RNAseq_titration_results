# Script prepares GBM expression data for use in pipelines
# For GBM array data: convert sample names and remove duplicate individuals
# For GBM seq data: filter all seq data for GBM samples only, convert genes IDs

# Steven Foltz July 2021

library(tidyverse)

# use optparse here
tcga_seq_expression_filepath <- "data/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv"
gbm_array_expression_filepath
metadata_json_filepath
gbm_array_output_filename
gbm_seq_output_filename
clinical_xlxs_filepath
clinical_xlxs_output_filepath

################################################################################
# Array data
################################################################################

gbm_array_expression <- read_tsv(gbm_array_expression_filepath,
                                 col_types = cols(
                                   .default = col_double(),
                                   Gene = col_character()
                                 ))

metadata_json <- jsonlite::fromJSON(metadata_json_filepath,
                                    simplifyVector = FALSE)

available_array_accession_ids <- colnames(gbm_array_expression)[-1]
all_array_tumor_samples <- tibble(accession = names(metadata_json$samples)) %>%
  filter(accession %in% available_array_accession_ids) %>% # check these are ==
  rowwise() %>%
  mutate(tcga_id_raw = metadata_json$samples[[accession]]$refinebio_annotations[[1]]$characteristics_ch1[[9]] %>%
           str_remove("sample: ")) %>%
  mutate(tcga_id = str_sub(tcga_id_raw, 1, 12),
         sample = str_sub(tcga_id_raw, 14, 15)) %>%
  filter(str_starts(sample, "0")) %>% # remove non-tumor samples
  ungroup()

array_accession_tcga_id_keep <- all_array_tumor_samples %>%
  group_by(tcga_id) %>%
  summarize(accession = sort(accession)[1])

accession_colnames_keep <- colnames(gbm_array_expression)[-1][colnames(gbm_array_expression)[-1] %in% array_accession_tcga_id_keep$accession]

map_accession_colnames_to_tcga_id <- tibble(accession = accession_colnames_keep) %>%
  left_join(array_accession_tcga_id_keep,
            by = "accession")

gbm_array_expression_renamed <- gbm_array_expression %>%
  select(c("Gene", map_accession_colnames_to_tcga_id$accession))
colnames(gbm_array_expression_renamed) <- c("sample",
                                            map_accession_colnames_to_tcga_id$tcga_id)

# write to file
write_tsv(gbm_array_expression_renamed,
          file = array_output_filename)

################################################################################
# Seq data
################################################################################

# read in column names of TCGA seq expression file
tcga_seq_expression_column_names <- read_tsv(tcga_seq_expression_filepath,
                                                   col_types = cols(
                                                     .default = col_double(),
                                                     gene_id = col_character()),
                                                   n_max = 0) %>%
  names()

gbm_seq_tumor_samples <- tibble(tcga_id_raw = tcga_seq_expression_column_names[-1]) %>%
  mutate(tcga_id = str_sub(tcga_id_raw, 1, 12),
         sample = str_sub(tcga_id_raw, 14, 15)) %>%
  filter(str_starts(sample, "0")) %>% # remove non-tumor samples
  filter(tcga_id %in% array_accession_tcga_id_keep$tcga_id) %>% # keep array GBMs
  group_by(tcga_id) %>%
  summarize(tcga_id_raw = sort(tcga_id_raw)[1]) # keep one ID per person

gbm_seq_expression <- read_tsv(tcga_seq_expression_filepath,
                               col_types = cols(
                                 .default = col_double(),
                                 gene_id = col_character()),
                               col_select = c("gene_id",
                                              gbm_seq_tumor_samples$tcga_id_raw))

colnames(gbm_seq_expression) <- c("gene_id",
                                  gbm_seq_tumor_samples$tcga_id)

symbol_entrez_ids <- gbm_seq_expression %>%
  select(gene_id) %>%
  separate(gene_id,
           into = c("SYMBOL", "ENTREZID"),
           sep = "\\|",
           remove = FALSE)

entrez_ensembl_ids <- ensembldb::select(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
                                        keys= symbol_entrez_ids$ENTREZID,
                                        keytype = "ENTREZID",
                                        columns = "GENEID") %>%
  mutate(ENTREZID = as.character(ENTREZID))

gene_id_mapping_in_array <- symbol_entrez_ids %>%
  left_join(entrez_ensembl_ids,
            by = "ENTREZID") %>%
  filter(!is.na(GENEID)) %>%
  rowwise() %>%
  filter(GENEID %in% gbm_array_expression$Gene)

gene_id_single_mapping_in_array <- gene_id_mapping_in_array %>%
  count(gene_id) %>%
  filter(n == 1) %>%
  select(gene_id) %>%
  left_join(gene_id_mapping_in_array,
            by = "gene_id") %>%
  select(gene_id, GENEID)

gbm_seq_expression_renamed <- gene_id_single_mapping_in_array %>%
  left_join(gbm_seq_expression,
            by = "gene_id") %>%
  select(-gene_id) %>%
  rename("sample" = "GENEID")

# write to file
write_tsv(gbm_seq_expression_renamed,
          file = gbm_seq_output_filename)

################################################################################
# subtype information
################################################################################

gbm_subtypes <- readxl::read_xlsx(path = clinical_xlxs_filepath,
                                  sheet = "Clinical Data",
                                  skip = 1) %>%
  select("Case ID",
         "MGMT Status",
         "G-CIMP\r\n methylation",
         "IDH1\r\n status",
         "Expression\r\nSubclass") %>%
  rename("sample" = "Case ID",
         "MGMT_methylation_status" = "MGMT Status",
         "G-CIMP_methylation" = "G-CIMP\r\n methylation",
         "IDH1_mutation_status" = "IDH1\r\n status",
         "subtype" = "Expression\r\nSubclass") %>%
  write_tsv(file = clinical_xlxs_output_filepath)
