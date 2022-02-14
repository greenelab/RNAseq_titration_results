from datetime import datetime
import requests
import xmltodict

#### GEO -----------------------------------------------------------------------

# search for results matching term
search_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
geo_array_search_term = "homo+sapiens[Organism]+AND+expression+profiling+by+array[DataSet+Type]"
geo_rnaseq_search_term = "homo+sapiens[Organism]+AND+expression+profiling+by+high+throughput+sequencing[DataSet+Type]"

geo_array_initial_url = search_base + "&".join(["db=gds", "term="+geo_array_search_term])
geo_rnaseq_initial_url = search_base + "&".join(["db=gds", "term="+geo_rnaseq_search_term])

geo_dict = {"array": [geo_array_initial_url, 0],
            "rnaseq": [geo_rnaseq_initial_url, 0]}

for platform in geo_dict:
  
  initial_url = geo_dict[platform][0]
  initial_xml = requests.get(initial_url)
  initial_dict = xmltodict.parse(initial_xml.content)

  n_results = initial_dict['eSearchResult']['Count']

  # search again using n_results, save a query_key and WebEnv
  second_url = initial_url + "&RetMax=" + n_results + "&usehistory=y"
  second_xml = requests.get(second_url)
  second_dict = xmltodict.parse(second_xml.content)

  query_key = second_dict['eSearchResult']['QueryKey']
  webenv = second_dict['eSearchResult']['WebEnv']

  # now fetch the records of the search results
  fetch_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
  skip_these = ["(Submitter supplied) This SuperSeries is composed of the SubSeries listed below."]
  retstart = 0

  while retstart < int(n_results):
    fetch_parameters = "&".join(["db=gds",
      "query_key="+query_key,
      "WebEnv="+webenv,
      "retmax=10000",
      "retstart="+str(retstart)])
    fetch_url = fetch_base + fetch_parameters
    fetch_text = requests.get(fetch_url).text
    for result in fetch_text.split("\n\n"):
      record = result.split("\n")
      for entry in record:
        if entry in skip_these:
          continue
        else:  
          if entry.startswith("Platform"):
            geo_dict[platform][1] += int(entry.split(" ")[-2])
    retstart += 10000

#### Array Express -------------------------------------------------------------

# Do not include GEO results in AE (directsub=on)
ae_array_url =  "https://www.ebi.ac.uk/arrayexpress/ArrayExpress-Experiments.txt?keywords=&organism=Homo+sapiens&exptype%5B%5D=%22rna+assay%22&exptype%5B%5D=%22array+assay%22&array=&directsub=on"
ae_rnaseq_url = "https://www.ebi.ac.uk/arrayexpress/ArrayExpress-Experiments.txt?keywords=&organism=Homo+sapiens&exptype%5B%5D=%22rna+assay%22&exptype%5B%5D=%22sequencing+assay%22&array=&directsub=on"

ae_dict = {"array": [ae_array_url, 0],
           "rnaseq": [ae_rnaseq_url, 0]}

for platform in ae_dict:
  url = ae_dict[platform][0]
  results = requests.get(url).text.split("\n")[1:-1]
  for entry in results:
    n_assays = int(entry.split('\t')[4])
    ae_dict[platform][1] += n_assays

#### Print results -------------------------------------------------------------
total_array = geo_dict["array"][1] + ae_dict["array"][1]
total_rnaseq = geo_dict["rnaseq"][1] + ae_dict["rnaseq"][1]
ratio = total_array/total_rnaseq

print('\t'.join(["Platform", "GEO", "AE", "Total"]))
print('\t'.join([str(x) for x in ["Array",
                                  geo_dict["array"][1],
                                  ae_dict["array"][1],
                                  total_array]]))
print('\t'.join([str(x) for x in ["RNA-seq",
                                  geo_dict["rnaseq"][1],
                                  ae_dict["rnaseq"][1],
                                  total_rnaseq]]))
print("Ratio: "+str(ratio))
print("Date: "+datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC"))
                                  
