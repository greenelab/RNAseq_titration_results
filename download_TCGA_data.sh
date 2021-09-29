#!/bin/bash
set -euo pipefail

# change to the directory of this script
cd "$(dirname "${BASH_SOURCE[0]}")"

# set data directory
data="data"
mkdir -p $data

# downlaod BRCA array and seq data from URLs
wget -nc -i brca_data_urls.txt '--directory-prefix='$data

# Obtain TCGA data freeze manifest file
# See here for more info: https://gdc.cancer.gov/about-data/publications/pancanatlas
manifest_url="https://gdc.cancer.gov/files/public/file/PanCan-General_Open_GDC-Manifest_2.txt"
manifest_basename=$(basename $manifest_url)
if [ -f $data/$manifest_basename ]; then
  echo TCGA file $data/$manifest_basename already exists and was not overwritten.
else
  echo Downloading $manifest_basename
  curl -o $data/$manifest_basename --silent $manifest_url
fi

# download specific files from TCGA manifest
mutations="mc3.v0.2.8.PUBLIC.maf.gz"
copy_number="broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.seg"
clinical="TCGA-CDR-SupplementalTableS1.xlsx"
rnaseq="EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv"

filename_array=($mutations \
                $copy_number \
                $clinical \
                $rnaseq)

for filename in ${filename_array[@]}; do
  if [ -f $data/$filename ]; then
    echo TCGA file $data/$filename already exists and was not overwritten.
  else
    echo Downloading $filename
    id=$(grep -w $filename $data/$manifest_basename | cut -f1)
    curl -o $data/$filename https://api.gdc.cancer.gov/data/$id
  fi
done

# get TCGA array expression using refine.bio client
# GSE83130 GBM
for accession in GSE83130; do
  if [ -d $data/$accession ]; then
    echo refine.bio download for $accession already exists and was not overwritten.
  else
    echo Downloading $accession
    refinebio create-token -s
    refinebio download-dataset \
      --email-address steven.foltz@ccdatalab.org \
      --path $data/$accession\.zip \
      --experiments $accession \
      --aggregation EXPERIMENT \
      --transformation NONE \
      --skip-quantile-normalization True
    unzip -d $data/$accession $data/$accession\.zip && rm -f $data/$accession\.zip
  fi
done

# download TCGA GBM clinical data including subtypes
# Publication: Brennan, C. W. et al. The somatic genomic landscape of glioblastoma. Cell 155, 462–477 (2013)
# Link to paper: https://doi.org/10.1016/j.cell.2013.09.034
gbm_clinical_link="https://www.cell.com/cms/10.1016/j.cell.2013.09.034/attachment/9cefc2e8-caac-4225-bcdd-70f105ccf568/mmc7.xlsx"
if [ -f $data/gbm_clinical_table_S7.xlsx ]; then
  echo GBM clinical spreadsheet $data/gbm_clinical_table_S7.xlsx already exists and was not overwritten.
else
  wget -O $data/gbm_clinical_table_S7.xlsx $gbm_clinical_link
fi

# check md5 sums of downloaded files
echo Checking md5 sums of downloaded files ...
md5sum --check check_sums.tsv
echo All files downloaded match expected md5 sums!

# modify BRCA clinical file column PAM50 to be subtype
sed -i 's/PAM50/subtype/' $data/BRCAClin.tsv

# process GBM data via script
echo Processing GBM data ...
Rscript prepare_GBM_data.R \
  --seq_input $data/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv \
  --array_input $data/GSE83130/GSE83130/GSE83130.tsv \
  --metadata_input $data/GSE83130/aggregated_metadata.json \
  --array_output $data/GBMarray.pcl \
  --seq_output $data/GBMRNASeq.pcl \
  --clinical_input $data/gbm_clinical_table_S7.xlsx \
  --clinical_output $data/GBMClin.tsv

# retrieve BRCA and GBM mutations in PIK3CA and TP53 from TCGA MC3
# output is stored in data/mutations.* TSV and MAF files
echo Retrieving mutation data from MC3 for BRCA and GBM ...
python3 retrieve_MC3_mutations.py

# combine clinical and mutation data into one data frame
echo Combining clinical and mutation data for BRCA ...
Rscript combine_clinical_data.R \
  --cancer_type BRCA \
  --clinical_input $data/BRCAClin.tsv \
  --mutation_input $data/mutations.BRCA.tsv \
  --combined_output $data/combined_clinical_data.BRCA.tsv
echo Combining clinical and mutation data for GBM ...
Rscript combine_clinical_data.R \
  --cancer_type GBM \
  --clinical_input $data/GBMClin.tsv \
  --mutation_input $data/mutations.GBM.tsv \
  --combined_output $data/combined_clinical_data.GBM.tsv

# get BRCA array expression data from TCGA Legacy Archive
# data/gdc_legacy_archive_brca_manifest.txt obtained from https://portal.gdc.cancer.gov/legacy-archive
# with search parameters
#   Cases
#     Cancer Program = TCGA
#     Project = TCGA-BRCA
#   Files
#     Data Category = Raw microarray data
#     Data Type = Raw intensities
#     Experimental Strategy = Gene expression array
################################################################################
# UNCOMMENT TO DOWNLOAD TCGA LEGACY ARCHIVE BRCA EXPRESSION ARRAY DATA
# Need to rebuild docker image with gdc-client uncommented
#brca_array_dir=$data/BRCA_array
#if [ -d $brca_array_dir ]; then
#  echo TCGA Legacy Archive data for BRCA already exists and was not overwritten.
#else
#  mkdir -p $brca_array_dir
#  gdc-client download --manifest gdc_legacy_archive_brca_manifest.txt --dir $brca_array_dir
#fi
################################################################################
