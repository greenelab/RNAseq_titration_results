from datetime import datetime
import requests
import xmltodict

###############################################################################
# This script queries two databases (GEO, ArrayExpress) to find human samples
# analyzed on array or RNA-seq platforms. It parses information from each data
# set and adds up the total number of samples from each platform. The output is
# table showing the number of samples from each database and platform, plus a
# ratio of the number of array samples to the number of RNA-seq samples. Also,
# since this information will change over time, a timestamp representing the
# moment these queries took place is included in the output.
#
# Usage: python3 util/search_geo_arrayexpress.py > output_file
#
# S. Foltz February 2022
###############################################################################

###############################################################################
# GEO - Gene Expression Omnibus
###############################################################################

# set up search terms and dictionary to track n_samples
search_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
geo_array_search_term = "homo+sapiens[Organism]+AND+expression+profiling+by+array[DataSet+Type]"
geo_rnaseq_search_term = "homo+sapiens[Organism]+AND+expression+profiling+by+high+throughput+sequencing[DataSet+Type]"

geo_array_initial_url = search_base +
"&".join(["db=gds", "term=" + geo_array_search_term])
geo_rnaseq_initial_url = search_base +
"&".join(["db=gds", "term=" + geo_rnaseq_search_term])

geo_dict = {"array": [geo_array_initial_url, 0],
            "rnaseq": [geo_rnaseq_initial_url, 0]}

# for each platform in the dictionary, search twice and then fetch samples
for platform in geo_dict:

    # first search to retrieve the total number of results
    initial_url = geo_dict[platform][0]
    initial_xml = requests.get(initial_url)
    initial_dict = xmltodict.parse(initial_xml.content)
    n_results = initial_dict['eSearchResult']['Count']

    # search again using n_results, save query_key and WebEnv for fetch
    second_url = initial_url + "&RetMax=" + n_results + "&usehistory=y"
    second_xml = requests.get(second_url)
    second_dict = xmltodict.parse(second_xml.content)
    query_key = second_dict['eSearchResult']['QueryKey']
    webenv = second_dict['eSearchResult']['WebEnv']

    # now fetch the records of the search results
    fetch_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
    # for any data sets with an entry in this list, we will skip that data set
    skip_these = ["(Submitter supplied) This SuperSeries is composed of the SubSeries listed below."]
    # fetch returns up to 10K results, so we need to define the start position
    # so we can increase the start position on subsequent fetches
    retstart = 0

    while retstart < int(n_results):
        fetch_parameters = "&".join(
          ["db=gds",
           "query_key=" + query_key,
           "WebEnv=" + webenv,
           "retmax=10000",
           "retstart=" + str(retstart)])
        fetch_url = fetch_base + fetch_parameters
        fetch_text = requests.get(fetch_url).text
        for result in fetch_text.split("\n\n"):  # split by \n
            record = result.split("\n")  # split by \n
            for entry in record:
                if entry in skip_these:  # data set should be skipped
                    continue
                else:  # otherwise, look for the line starting with "Platform"
                    if entry.startswith("Platform"):  # parse 2nd last element
                        n_samples = int(entry.split(" ")[-2])  # (n_samples)
                        geo_dict[platform][1] += n_samples  # increment count
        retstart += 10000  # increment the start position

###############################################################################
# ArrayExpress
###############################################################################

# set up search terms and dictionary to track n_samples
# do not include GEO results in AE (directsub=on)
ae_base_url = "https://www.ebi.ac.uk/arrayexpress/ArrayExpress-Experiments.txt?keywords="
ae_array_url = ae_base_url + "&organism=Homo+sapiens&exptype%5B%5D=%22rna+assay%22&exptype%5B%5D=%22array+assay%22&array=&directsub=on"
ae_rnaseq_url = ae_base_url + "&organism=Homo+sapiens&exptype%5B%5D=%22rna+assay%22&exptype%5B%5D=%22sequencing+assay%22&array=&directsub=on"

ae_dict = {"array": [ae_array_url, 0],
           "rnaseq": [ae_rnaseq_url, 0]}

# for each platform in the dictionary, get results from url
for platform in ae_dict:
    url = ae_dict[platform][0]
    results = requests.get(url).text.split("\n")[1:-1]  # skip first, last
    for entry in results:
        n_assays = int(entry.split('\t')[4])  # fifth column is n_assays
        ae_dict[platform][1] += n_assays  # increment the count

###############################################################################
# Print results
###############################################################################

total_array = geo_dict["array"][1] + ae_dict["array"][1]  # total number array
total_rnaseq = geo_dict["rnaseq"][1] + ae_dict["rnaseq"][1]  # total n RNA-seq
ratio = total_array/total_rnaseq  # array:RNA-seq

print('\t'.join(["Platform", "GEO", "AE", "Total"]))
print('\t'.join([str(x) for x in ["Array",
                                  geo_dict["array"][1],
                                  ae_dict["array"][1],
                                  total_array]]))
print('\t'.join([str(x) for x in ["RNA-seq",
                                  geo_dict["rnaseq"][1],
                                  ae_dict["rnaseq"][1],
                                  total_rnaseq]]))
print("Ratio: " + str(ratio))
print("Date: " + datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC"))
