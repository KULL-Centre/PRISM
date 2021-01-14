"""get_domains.py script to extract domains and 
sequence/structural features from uniprot ids

Author: Johanna K.S. Tiemann

Date of last major changes: 2020-10-19

"""

# Standard library imports
import gzip
import json
import logging as log
import os
import urllib.parse
import urllib.request
import xmltodict

# Third party imports
import pandas as pd


def extract_single_protein_pfam( uniprot_id, verbose=False ):
    """
    Extracts domain information from pfam for present uniprot ID. Return type changed from dictionary to list of dictionaries.
    """

    disprot_url = (
        f'http://pfam.xfam.org/protein/{uniprot_id}?output=xml'
    )
    request_header = {"format": "xml"} 
    request = urllib.request.Request(disprot_url, headers=request_header)
    contents = urllib.request.urlopen(request).read()
    
    dictString = xmltodict.parse(contents.decode("utf-8"), attr_prefix='', cdata_key='#text')

    result_ls=[]
    #debug
    #print(dictString.keys())
    #print(dictString['error'])
    #debug
    if 'error' in dictString:
        print('pfam error:', dictString['error'])
        return()
    
    if 'matches' in dictString['pfam']['entry'].keys():
        subdic = dictString['pfam']['entry']['matches']['match']
        for keys in subdic:
            if isinstance(keys,str):
                result_ls.append({'acc': subdic['accession'], 'id':subdic['id'], 'start':int(subdic['location']['start']), 'end':int(subdic['location']['end'])})
            else:
                result_ls.append({'acc': keys['accession'], 'id':keys['id'], 'start':int(keys['location']['start']), 'end':int(keys['location']['end'])})

    if verbose:
        jsonString = json.dumps(xmltodict.parse(contents.decode("utf-8"), attr_prefix='', cdata_key='#text'), indent=4)
        return result_ls, jsonString
    else:
        return result_ls


def extract_pfam_nested(write_dir=''):
    """
    Extracts pfam regions which are nested/overlapping
    """
    url = 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.1/database_files/nested_locations.txt.gz'
    
    data_array = [['pfamA_acc','nested_pfamA_acc','pfamseq_acc','seq_version','seq_start','seq_end']]
    req = urllib.request.Request(url)
    req.add_header('Accept-Encoding', 'gzip')
    response = urllib.request.urlopen(req)
    content = gzip.decompress(response.read())
    
    decomp_req = content.splitlines()
    for line in decomp_req:
        data_array.append(line.decode('utf-8').split('\t'))
    all_pfam_df = pd.DataFrame(data=data_array[1:], columns=data_array[0])
    if write_dir:
        all_pfam_df.to_csv(os.path.join(write_dir,'nested_domains.csv'), index=False)
    return all_pfam_df


def extract_pfam_pdb_mapping():
    # extract domains with available pdbs from pfam
    url = 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/mappings/pdb_pfam_mapping.txt'

    data_array = []
    for line in urllib.request.urlopen(url):
        line = line.decode('utf-8')
        unprocessed = line[:-1].split('\t')
        data_array.append(unprocessed)
    all_pfam_df = pd.DataFrame(data=data_array[1:], columns=data_array[0])

    return all_pfam_df



def extract_pfam_all_release():
    """
    Extracts pdb regions for specific domains from pfam
    """
    url = 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.1/database_files/pdb_pfamA_reg.txt.gz'
    
    data_array = []
    req = urllib.request.Request(url)
    req.add_header('Accept-Encoding', 'gzip')
    response = urllib.request.urlopen(req)
    content = gzip.decompress(response.read())
    
    decomp_req = content.splitlines()
    for line in decomp_req:
        data_array.append(line.decode('utf-8').split('\t'))
    all_pfam_df = pd.DataFrame(data=data_array)

    return all_pfam_df



def extract_single_infos(uniprot_id, human_proteome_info_df):
    """
    Extract additional information like TM region, topology domains, 
    disulfid bonds, disorder from uniprot and other databases
    """
    headers = ['Entry', 'Entry name', 'Gene names  (primary )', 'Organism',
           'Cross-reference (PDB)', 'Cross-reference (DisProt)',
           'Cross-reference (MobiDB)', 'Length', 'Transmembrane', 'Intramembrane',
           'Topological domain', 'Sequence similarities', 'Domain [CC]',
           'Domain [FT]', 'Protein families', 'Coiled coil', 'Motif', 'Region',
           'Repeat', 'Zinc finger', 'Disulfide bond', 'Active site',
           'Binding site', 'Cross-reference (Pfam)',
           'Cross-reference (InterPro)', 'Cross-reference (SUPFAM)',
           'Cross-reference (PROSITE)']
    r = human_proteome_info_df.loc[human_proteome_info_df[headers[0]] == uniprot_id]
    get_transmembrane = list(human_proteome_info_df.loc[human_proteome_info_df[headers[0]] == uniprot_id]['Transmembrane'])[0].strip().split('; ')
    tm_region = []
    for index, elem in enumerate(get_transmembrane):
        if elem.startswith('TRANSMEM'):
            tm_str = get_transmembrane[index].split()[1].split('..')
            tm_str.append(get_transmembrane[index+1].split('"')[1])
            tm_region.append(tm_str)
    get_topdomain = list(human_proteome_info_df.loc[human_proteome_info_df[headers[0]] == uniprot_id]['Topological domain'])[0].strip().split('; ')
    top_domain = []
    for index, elem in enumerate(get_topdomain):
        if elem.startswith('TOPO_DOM'):
            top_str = get_topdomain[index].split()[1].split('..')
            top_str.append(get_topdomain[index+1].split('"')[1])
            top_domain.append(top_str)
    get_ss = list(human_proteome_info_df.loc[human_proteome_info_df[headers[0]] == uniprot_id]['Disulfide bond'])[0].strip().split('; ')
    ss_list = []
    for index, elem in enumerate(get_ss):
        if elem.startswith('DISULFID'):
            ss_str = get_ss[index].split()[1].split('..')
            if len(ss_str) == 1:
                ss_str.append(ss_str[0])
                ss_str.append('disulfid:special')
            else:
                ss_str.append('disulfid')
            ss_list.append(ss_str)
    get_activesite = list(human_proteome_info_df.loc[human_proteome_info_df[headers[0]] == uniprot_id]['Active site'])[0].strip().split('; ')
    activesite_list = []
    for index, elem in enumerate(get_activesite):
        if elem.startswith('ACT_SITE'):
            activesite_str = get_activesite[index].split()[1].split('..')
            activesite_str.append('func_site:active')
            activesite_list.append(activesite_str)
    get_bindingsite = list(human_proteome_info_df.loc[human_proteome_info_df[headers[0]] == uniprot_id]['Binding site'])[0].strip().split('; ')
    bindingsite_list = []
    for index, elem in enumerate(get_bindingsite):
        if elem.startswith('BINDING'):
            bindingsite_str = get_bindingsite[index].split()[1].split('..')
            bindingsite_str.append('func_site:binding')
            bindingsite_list.append(bindingsite_str)

    mobidb_dis_content = get_mobidb_disordered(uniprot_id)
    disprot_dis_content = get_disprot_disordered(uniprot_id)
    if len(mobidb_dis_content) < len(disprot_dis_content):
        disorder = disprot_dis_content
    else:
        disorder = mobidb_dis_content
    if tm_region:
        tm_region = [[f'tm:{elem[2]}',elem[0],elem[1]] for elem in tm_region]
    if top_domain:
        tmp_top_domain = []
        for elem in top_domain:
            if len(elem) == 3:
                tmp_top_domain.append([elem[2],elem[0],elem[1]])
            else:
                tmp_top_domain.append([elem[1],elem[0],elem[0]])
        top_domain = tmp_top_domain
    if disorder:
        disorder = [[f'disorder:{elem[2]}',elem[0],elem[1]] for elem in disorder]
    if ss_list:
        ss_list = [[elem[2],elem[0],elem[1]] for elem in ss_list]
    func_sites = bindingsite_list + activesite_list
    if func_sites:
        func_sites = [[elem[1],elem[0]] for elem in func_sites]
    result_dict = {
        'tm_region': tm_region,
        'topol_domain': top_domain,
        'disorder': disorder,
        'disulfid': ss_list,
        'func_sites': func_sites,
    }
    return result_dict

