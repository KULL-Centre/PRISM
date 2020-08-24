"""get_db_features obtains information from different databases

Author: Johanna K.S. Tiemann

Date of last major changes: 2020-08-24

"""

# Standard library imports
import urllib.request
import json


def get_disprot_disordered(uniprot_id):
    """
    Extracts disorder content (%) from disprot (https://disprot.org/) 
    or present uniprot ID
    """
    disprot_url = f'https://disprot.org/api/{uniprot_id}'
    request_header = {"format": "json"} 
    request = urllib.request.Request(disprot_url, headers=request_header)
    contents = urllib.request.urlopen(request).read()
    disordered_region = []
    if contents:
        all_json = json.loads(contents.decode("utf-8"))
        for elem in all_json['disprot_consensus']['full']:
            disordered_region.append([elem['start'],elem['end'],elem['type']])
    return disordered_region


def get_mobidb_disordered(uniprot_id):
    """
    Extracts disorder content (%) from mobidb (https://mobidb.bio.unipd.it/) 
    for present uniprot_ID
    """
    mobidb_url = f'https://mobidb.bio.unipd.it/ws/{uniprot_id}/consensus'
    request_header = {"format": "json"} 
    request = urllib.request.Request(mobidb_url, headers=request_header)
    contents = urllib.request.urlopen(request).read()
    all_json = json.loads(contents.decode("utf-8"))
    disordered_region = []
    if 'db' in all_json['mobidb_consensus']['disorder'].keys():
        for elem in all_json['mobidb_consensus']['disorder']['db'][0]['regions']:
            if elem[2]!='S':
                disordered_region.append(elem)
    return disordered_region