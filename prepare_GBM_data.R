#install.packages("jsonlite")
#BiocManager::install("EnsDb.Hsapiens.v86")

# Script prepares GBM expression data for use in pipelines
# For GBM array data: convert sample names and remove duplicate individuals
# For GBM seq data: filter all seq data for GBM samples only, convert genes IDs

# Steven Foltz 2021-07

library(tidyverse)
library(EnsDb.Hsapiens.v86)

gbm_array_expression <- read_tsv("data/GSE83130/GSE83130/GSE83130.tsv",
                                 col_types = cols(
                                   .default = col_double(),
                                   Gene = col_character()
                                 ))
metadata_json <- jsonlite::fromJSON("data/GSE83130/aggregated_metadata.json", simplifyVector = FALSE)

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
  dplyr::select(c("Gene", map_accession_colnames_to_tcga_id$accession))
colnames(gbm_array_expression_renamed) <- c("sample",
                                            map_accession_colnames_to_tcga_id$tcga_id)

# write to file




tcga_seq_expression <- read_tsv("data/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv")

gbm_seq_tumor_samples <- tibble(tcga_id_raw = colnames(tcga_seq_expression)[-1]) %>%
  mutate(tcga_id = str_sub(tcga_id_raw, 1, 12),
         sample = str_sub(tcga_id_raw, 14, 15)) %>%
  filter(str_starts(sample, "0")) %>% # remove non-tumor samples
  filter(tcga_id %in% array_accession_tcga_id_keep$tcga_id) # keep array GBMs

seq_tcga_id_raw_keep <- gbm_seq_tumor_samples %>%
  group_by(tcga_id) %>%
  summarize(tcga_id_raw = sort(tcga_id_raw)[1])

gbm_seq_expression_renamed <- tcga_seq_expression %>%
  select(c("gene_id", seq_tcga_id_raw_keep$tcga_id_raw))
colnames(gbm_seq_expression_renamed) <- c("gene_id",
                                          seq_tcga_id_raw_keep$tcga_id)


symbol_entrez_ids <- tcga_seq_expression %>%
  dplyr::select(gene_id) %>%
  tidyr::separate(gene_id,
                  into = c("SYMBOL", "ENTREZID"),
                  sep = "\\|",
                  remove = FALSE)

entrez_ensembl_ids <- ensembldb::select(EnsDb.Hsapiens.v86,
                                        keys= symbol_entrez_ids$ENTREZID,
                                        keytype = "ENTREZID",
                                        columns = "GENEID") %>%
  mutate(ENTREZID = as.character(ENTREZID))

gene_id_mapping_in_array <- symbol_entrez_ids %>%
  left_join(entrez_ensembl_ids,
            by = "ENTREZID") %>%
  dplyr::filter(!is.na(GENEID)) %>%
  rowwise() %>%
  dplyr::filter(GENEID %in% gbm_array_expression_renamed$sample)

gene_id_single_mapping_in_array <- gene_id_mapping_in_array %>%
  count(gene_id) %>%
  dplyr::filter(n == 1) %>%
  dplyr::select(gene_id) %>%
  left_join(gene_id_mapping_in_array,
            by = "gene_id") %>%
  dplyr::select(gene_id, GENEID)

gbm_seq_expression_renamed <- gene_id_single_mapping_in_array %>%
  left_join(gbm_seq_expression_renamed,
            by = "gene_id") %>%
  dplyr::select(-gene_id) %>%
  dplyr::rename("sample" = "GENEID")


# write to file
