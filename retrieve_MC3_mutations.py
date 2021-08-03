import gzip
import sys

# define output filepaths
output_tsv_filepath = "data/mutations.tsv"
output_maf_filepath = "data/mutations.maf"

# TCGA MC3 MAF from https://gdc.cancer.gov/about-data/publications/pancanatlas
mc3_filename = "data/mc3.v0.2.8.PUBLIC.maf.gz"

# cancer types and genes of interest
cancer_type_abbrevs = {"Breast invasive carcinoma": "BRCA",
                       "Glioblastoma multiforme": "GBM"}
cancer_types_of_interest = cancer_type_abbrevs.keys()
genes_of_interest = ["PIK3CA", "PTEN", "TP53"]

############################################################
# Tissue source sites define the cancer type of the sample #
############################################################

# TSS codes from https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tissue-source-site-codes
tcga_tss_codes = open("tcga_tss_codes.csv", "r")
tcga_tss_codes.readline()
tcga_tss_codes_dict = {}

# set up TSS dictionary with {tss: cancer_type}
for line in tcga_tss_codes:
  k,v = line.strip().split(",")
  tcga_tss_codes_dict[k] = v

tcga_tss_codes.close()

###############################
# Retrieve mutations from MC3 #
###############################

# simple mutation dictionary {cancer_type: {tcga_id: ["PTEN", "TP53"]}}
# this will be used at end to create simple 0/1 mutation status data frame
mutation_dict = {x:{} for x in cancer_types_of_interest}

# open up MC3 and define header lines
mc3 = gzip.open(mc3_filename, "rb")
maf_header = mc3.readline().decode('UTF-8').strip()
tsv_header = "\t".join(["tcga_id", "cancer_type", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "Hugo_Symbol", "Chromosome", "Start_Position", "HGVSc", "VARIANT_CLASS"])

output_tsv = open(output_tsv_filepath, "w")
output_maf = open(output_maf_filepath, "w")

output_tsv.write(tsv_header + "\n")
output_maf.write(maf_header + "\n")

# progress through each line of MC3
# check if sample is primary solid tumor ("01") and
# from out genes and cancer types of interest
for line in mc3:
  record = line.decode('UTF-8').strip().split("\t")
  hugo_symbol = record[0] # gene name
  tcga_id_raw = record[15] # TCGA-XX-YYYY-ZZ...
  tcga_id_raw_normal = record[16]
  is_tumor = tcga_id_raw.split("-")[3].startswith("01")
  tss_code = tcga_id_raw.split("-")[1]
  cancer_type = tcga_tss_codes_dict[tss_code]

  if is_tumor and cancer_type in cancer_types_of_interest:
    tcga_id = tcga_id_raw[0:15]
    
    # add TCGA ID to mutation dict
    if tcga_id not in mutation_dict[cancer_type]:
      mutation_dict[cancer_type][tcga_id] = []

    # if gene of interest, add to mutation dict list for that ID and outputs
    if hugo_symbol in genes_of_interest:
      chromosome = record[4]
      start_position = record[5]
      hgvsc = record[34]
      variant_class = record[94]
    
      mutation_dict[cancer_type][tcga_id].append(hugo_symbol)

      output_tsv.write("\t".join([tcga_id, cancer_type_abbrevs[cancer_type], tcga_id_raw, tcga_id_raw_normal, hugo_symbol, chromosome, start_position, hgvsc, variant_class]) + "\n")
      output_maf.write("\t".join(record) + "\n")
  
output_tsv.close()
output_maf.close()

# write cancer-type-specific simple output data frames (0/1 for gene mutation status)
for cancer_type in cancer_types_of_interest:

  simple_output_filename = ".".join(["data/mutations", cancer_type_abbrevs[cancer_type], "tsv"])
  simple_output = open(simple_output_filename, "w")
  simple_output_header = "\t".join(["tcga_id", "PIK3CA_mutation", "PTEN_mutation", "TP53_mutation"]) + "\n"

  simple_output.write(simple_output_header)

  for tcga_id, mutation_list in mutation_dict[cancer_type].items():
    mutation_status_list = [str(int(x in mutation_list)) for x in genes_of_interest]
    simple_output.write(tcga_id + "\t" + "\t".join(mutation_status_list) + "\n")

  simple_output.close()
