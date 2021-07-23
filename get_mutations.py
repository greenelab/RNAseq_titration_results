import gzip

tcga_tss_codes = open("tcga_tss_codes.csv", "r")
tcga_tss_codes.readline()
tcga_tss_codes_dict = {}

for line in tcga_tss_codes:
  k,v = line.strip().split(",")
  tcga_tss_codes_dict[k] = v

mc3 = gzip.open("mc3.v0.2.8.PUBLIC.maf.gz", "rb")
mc3.readline()

total_n_dict = {}

for line in mc3:
  record = line.decode('ascii').strip().split("\t")
  gene = record[0]
  id = record[15]
  person = id[0:12]
  tumor_or_normal = id.split("-")[3]
  cancer_type = tcga_tss_codes_dict[id.split("-")[1]]
  if cancer_type in total_n_dict:
    total_n_dict[cancer_type].append(person)
  else:
    total_n_dict[cancer_type] = [person]

print(len(set(total_n_dict["Breast invasive carcinoma"])))
print(len(set(total_n_dict["Glioblastoma multiforme"])))

  #if cancer_type in ["Breast invasive carcinoma", "Glioblastoma multiforme"] and gene in ["PIK3CA", "PTEN", "TP53"] and tumor_or_normal.startswith("0"):
  #  print("\t".join([cancer_type, person, gene]))

