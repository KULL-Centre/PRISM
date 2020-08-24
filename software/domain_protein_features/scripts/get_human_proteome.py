"""get_human_proteome.py obtaines ids of all human proteins deposited in uniprot

Author: Johanna K.S. Tiemann

Date of last major changes: 2020-08-21

"""

# Standard library imports
import urllib.request
import sys


def extract_uniprot_human_proteome_ids(reviewed="yes"):
    """Extract humand proteome from uniprot"""
    url = (
        "https://www.uniprot.org/uniprot/"
        f"?query=reviewed:{reviewed}"
        "+AND+proteome:up000005640"
        "+&format=tab&columns=id"
    )
    data_array = []
    for line in urllib.request.urlopen(url):
        line = line.decode('utf-8')
        unprocessed = line[:-1].split('\t')
        data_array.append(unprocessed)
    
    flat_list = [item[0] for item in data_array]
    return flat_list[1:]


if __name__ == '__main__':
    if len(sys.argv) > 1:
        reviewed = sys.argv[4]
    else:
        reviewed = 'yes'
    extract_uniprot_human_proteome_ids(reviewed=reviewed)