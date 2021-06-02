#!/bin/bash
set -euo pipefail

# change to the directory of this script
cd "$(dirname "${BASH_SOURCE[0]}")"

# create raw and processed data directories
data="data"
mkdir -p $data

manifest_url="https://gdc.cancer.gov/files/public/file/PanCan-General_Open_GDC-Manifest_2.txt"

# Obtain TCGA manifest file
manifest_basename=$(basename $manifest_url)
if [ -f $data/$manifest_basename ]; then
  echo TCGA file $data/$manifest_basename already exists and was not overwritten.
else
  echo Downloading $manifest_basename
  curl -o $data/$manifest_basename --silent $manifest_url
fi

# download specific files from TCGA manifest
copy_number="broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.seg"
mutations="mc3.v0.2.8.PUBLIC.maf.gz"
clinical="TCGA-CDR-SupplementalTableS1.xlsx"
rnaseq="EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv"

for filename in $copy_number $mutations $clinical $rnaseq; do
  if [ -f $data/$filename ]; then
    echo TCGA file $data/$filename already exists and was not overwritten.
  else
    echo Downloading $filename
    id=$(grep -w $filename $data/$manifest_basename | cut -f1)
    curl -o $data/$filename --silent https://api.gdc.cancer.gov/data/$id
  fi
done

# get data from refinebio
for accession in GSE83130 GSE68833; do

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
      --transformation NONE
    unzip -d $data/$accession data/$accession\.zip && rm -f data/$accession\.zip
  fi
done

# get BRCA data from TCGA Legacy Archive
# gdc_manifest.2021-06-01.txt was obtained from https://portal.gdc.cancer.gov/legacy-archive
# with search parameters
#   Cases
#     Cancer Program = TCGA
#     Project = TCGA-BRCA
#   Files
#     Data Category = Raw microarray data
#     Data Type = Raw intensities
#     Experimental Strategy = Gene expression array
brca_array_dir=$data/BRCA_array
if [ -d $brca_array_dir ]; then
  echo TCGA Legacy Archive data for BRCA already exists and was not overwritten.
else
  mkdir -p $brca_array_dir
  gdc-client download --manifest gdc_manifest.2021-06-01.txt --dir $brca_array_dir
fi

