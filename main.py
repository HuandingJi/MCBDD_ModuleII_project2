from chembl_webresource_client.new_client import new_client
import pandas as pd
import statistics
import uniprot as uni
import requests
import xml.etree.ElementTree as ET
from collections import Counter


def main():
    pd.set_option('display.max_columns', None)
    # retrieve drugs dataset
    drug = new_client.drug

    # 1. rank the drugs by year and then by name
    drug_approved = drug.filter(max_phase=4)
    drug_approved_rank_year_name = drug_approved.order_by("first_approval__gte").order_by("usan_stem")
    print("There are a total of " + str(len(drug_approved_rank_year_name)) + " drugs")
    # print("Drugs ranked by their approval year and name:")
    # i = 1
    # for item in drug_rank_year_name:
    #     print(str(i), item)
    #     i += 1
    # print("============= End of ranking =============")

    # 2. get approved drugs since 2013, calculate median number of protein targets associated wth each drug
    drug_after2013 = drug.filter(first_approval__gte=2013)
    print("There are a total of " + str(len(drug_after2013)) + " drugs approved after 2013")
    print(drug_after2013[0].keys())

    # Create a dictionary to store the UniProt accession numbers for each compound
    accession_numbers = {}
    # Loop through each recent approved compound
    count = 0
    for compound in drug_after2013:
        # Get the molecule ID for the compound
        molregno = compound['molecule_chembl_id']
        print("got molregno")
        # Get the protein targets associated with the compound
        targets = new_client.target.filter(molecules__chembl_id=molregno)
        print("got targets")
        # Get the UniProt accession numbers for each target
        uniprot_list = []
        for target in targets:
            try:
                uniprot_list.append(target['target_components'][0]['accession'])
            except:
                pass
        uniprot_accessions = set(uniprot_list)

        # uniprot_accessions = set([target['target_components'][0]['accession'] for target in targets])
        print("got uniprot accession")
        # Store the UniProt accession numbers in the dictionary
        accession_numbers[molregno] = uniprot_accessions
        print("number of uniprots: " + str(len(uniprot_accessions)))
        print(count)
        count += 1

    # Get the number of UniProt accession numbers for each compound
    num_targets = [len(accession_numbers[molregno]) for molregno in accession_numbers]

    # Calculate the median number of UniProt accession numbers
    median_num_targets = statistics.median(num_targets)

    # Print the median number of UniProt accession numbers
    print("The median number of UniProt accession numbers associated with each drug developed after 2013 is:",
          str(median_num_targets))


    # 3. Retrieve UniProt keywords associated with it, and identify most frequent key word

    # Create a list to store the keywords for each UniProt accession number
    keywords = []
    # Loop through each UniProt accession number
    for accession in set.union(*accession_numbers.values()):
        # Retrieve the keywords for the protein from uniprot web service
        url = f"https://www.uniprot.org/uniprot/{accession}.xml"
        response = requests.get(url)
        # Parse the XML response using ElementTree
        root = ET.fromstring(response.content)
        # Extract the keywords from the XML response and append them to the list
        keyword_elements = root.findall(".//{http://uniprot.org/uniprot}keyword")
        keyword_names = [element.attrib['id'] for element in keyword_elements]
        keywords.extend(keyword_names)

    # Count the frequency of each keyword in the list
    keyword_counts = Counter(keywords)

    # Get the most common keyword(s)
    most_common_keywords = [keyword for keyword, count in keyword_counts.most_common() if
                            count == keyword_counts.most_common(1)[0][1]]

    # Print the most common keyword(s)
    if len(most_common_keywords) == 1:
        print(f"The keyword associated with the most drugs approved since 2013 is: {most_common_keywords[0]}")
    else:
        print("The following keywords are associated with the most drugs approved since 2013:")
        for keyword in most_common_keywords:
            print(keyword)

    
# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()


