import gzip
import sys

output_df_filepath = "data/detailed_mutations.tsv"
output_maf_filepath = "data/mutations.maf"
simple_output_filepaths = {"Breast invasive carcinoma": "data/BRCA_mutations.tsv",
                           "Glioblastoma multiforme": "data/GBM_mutations.tsv"}

genes_of_interest = ["PIK3CA", "PTEN", "TP53"]
cancer_types_of_interest = simple_output_filepaths.keys()

# TSS codes from https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tissue-source-site-codes
tcga_tss_codes = open("tcga_tss_codes.csv", "r")
tcga_tss_codes.readline()
tcga_tss_codes_dict = {}

# set up TSS dictionary with {tss: cancer_type}
for line in tcga_tss_codes:
  k,v = line.strip().split(",")
  tcga_tss_codes_dict[k] = v

tcga_tss_codes.close()

# simple mutation dictionary {cancer_type: {tcga_id: ["PTEN", "TP53"]}}
mutation_dict = {x:{} for x in cancer_types_of_interest}

# MC3 from https://gdc.cancer.gov/about-data/publications/pancanatlas
mc3 = gzip.open("data/mc3.v0.2.8.PUBLIC.maf.gz", "rb")
maf_header = mc3.readline().decode('UTF-8').strip()
df_header = "\t".join(["tcga_id", "cancer_type", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "Hugo_Symbol", "Chromosome", "Start_Position", "HGVSc", "VARIANT_CLASS"])

output_df = open(output_df_filepath, "w")
output_maf = open(output_maf_filepath, "w")

output_df.write(df_header + "\n")
output_maf.write(maf_header + "\n")

for line in mc3:
  record = line.decode('UTF-8').strip().split("\t")
  hugo_symbol = record[0]
  tcga_id_raw = record[15]
  tcga_id_raw_normal = record[16]
  is_tumor = tcga_id_raw.split("-")[3].startswith("01")
  tss_code = tcga_id_raw.split("-")[1]
  cancer_type = tcga_tss_codes_dict[tss_code]

  if is_tumor and cancer_type in cancer_types_of_interest:
    tcga_id = tcga_id_raw[0:15]
    if tcga_id not in mutation_dict[cancer_type]:
      mutation_dict[cancer_type][tcga_id] = []

    if hugo_symbol in genes_of_interest:
      chromosome = record[4]
      start_position = record[5]
      hgvsc = record[34]
      variant_class = record[94]
    
      mutation_dict[cancer_type][tcga_id].append(hugo_symbol)

      output_df.write("\t".join([tcga_id, cancer_type, tcga_id_raw, tcga_id_raw_normal, hugo_symbol, chromosome, start_position, hgvsc, variant_class]) + "\n")
      output_maf.write("\t".join(record) + "\n")
  
output_df.close()
output_maf.close()

for cancer_type in cancer_types_of_interest:

  simple_output = open(simple_output_filepaths[cancer_type], "w")
  simple_output_header = "\t".join(["tcga_id", "PIK3CA_mutation", "PTEN_mutation", "TP53_mutation"]) + "\n"

  simple_output.write(simple_output_header)

  for tcga_id, mutation_list in mutation_dict[cancer_type].items():
    mutation_status_list = [str(int(x in mutation_list)) for x in genes_of_interest]
    simple_output.write(tcga_id + "\t" + "\t".join(mutation_status_list) + "\n")

  simple_output.close()
