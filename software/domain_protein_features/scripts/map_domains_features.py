"""map_domains_features.py maps the domains, the features and the pdbs together

Author: Johanna K.S. Tiemann

Date of last major changes: 2020-08-24

"""

# Standard library imports
import logging as log
import os
import sys
import requests
import time
import urllib.request
import urllib.parse
import json
import xmltodict
import gzip

# Third party imports
from bs4 import BeautifulSoup as bs
import numpy as np
import pandas as pd

# Local application imports
from get_db_features import get_disprot_disordered, get_mobidb_disordered
from get_domains import extract_single_protein_pfam


log_message = "verbose"

if log_message == "verbose":
    log.basicConfig(
        format='%(levelname)s:%(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=log.INFO
    )
elif log_message == "debug":
    log.basicConfig(
        format='%(levelname)s:%(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=log.WARNING
    )
else:
    log.basicConfig(
        format='%(levelname)s:%(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=log.ERROR
    )

logger = log.getLogger(__name__)


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


def map_pfam_pdb(uniprot_id, all_pfam_df, pdb_pfam_df, human_proteome_info_df, write_dir=''):
    """
    map all together into output format
    """
    result_dic = extract_single_protein_pfam(uniprot_id)
    tp_results = extract_single_infos(uniprot_id, human_proteome_info_df)
    tp_results_to_list = []
    for elem in tp_results:
        tp_results_to_list += tp_results[elem] 
    combined_dict = {}
    error_list = []
    for key in result_dic:
        result_arr = []
        result_arr.append(result_dic[key])
        results_info = []
        for key2 in tp_results:
            for entry in tp_results[key2]:
                if len(entry) == 2:
                    if int(result_dic[key][1]) <= int(entry[1]) <= int(result_dic[key][2]):
                        results_info.append(entry)
                        error_list.append(entry)
                else:
                    range1 = range(int(entry[1]),int(entry[2]))
                    range2 = range(int(result_dic[key][1]), int(result_dic[key][2]))
                    if set(range1).intersection(range2):
                        results_info.append(entry)
                        error_list.append(entry)
        result_arr.append(results_info)
        results_pdb = all_pfam_df.loc[all_pfam_df[4] == uniprot_id]
        results_pdb = results_pdb.loc[results_pdb[3] == key]
        arr = []
        for line in results_pdb.to_numpy():
            arr.append([line[2], line[5], line[6], line[8], line[10], line[11]])
        domain_pdb_array = []
        if isinstance(pdb_pfam_df,(pd.core.frame.DataFrame)):
            new = pdb_pfam_df["PFAM_ACC"].str.split(".", n = 1, expand = True) 
            pdb_pfam_df["PFAM_ACC_short"] = new[0] 
            domain_pdb = pdb_pfam_df.loc[pdb_pfam_df['PFAM_ACC_short'] == key]
            for line in domain_pdb.to_numpy():
                domain_pdb_array.append([line[0], line[1], line[2], line[3], '', ''])
            #add_domain_pdb = []
            #for elem in domain_pdb_array:
            #    if not any(elem[0] in x[0] for x in arr):
            #        add_domain_pdb.append(elem)
        result_arr.append([arr,domain_pdb_array])
        combined_dict[key] = result_arr
    if write_dir:
        with open(os.path.join(write_dir,f'domains_features_{uniprot_id}.json'), 'w') as json_file:
            json.dump(human_genome, json_file, indent=4)
    
    #information not covered by pfam regions
    uncovered_by_pfam = (set(tuple(i) for i in error_list) | set(tuple(i) for i in tp_results_to_list)) - (set(tuple(i) for i in error_list) & set(tuple(i) for i in tp_results_to_list))
    if len(uncovered_by_pfam) == 0:
        return combined_dict, False
    else:
        if write_dir:
            with open(os.path.join(write_dir,f'domains_features_uncovered_{uniprot_id}.json'), 'w') as json_file:
                json.dump(human_genome, json_file, indent=4)
        return combined_dict, uncovered_by_pfam

def get_domain_features_human_proteome(human_proteome_ids, all_pfam_df, pdb_pfam_df, human_proteome_info_df, write_dir=''):

    human_genome = {}
    uncovered_info = {}
    error_list = []
    for uniprot_id in human_proteome_ids:
        try:
            result, uncovered = map_pfam_pdb(uniprot_id, all_pfam_df, pdb_pfam_df, human_proteome_info_df, write_dir='')
            human_genome[uniprot_id] = result
            if uncovered:
                uncovered_info[uniprot_id] = uncovered
        except:
            error_list.append(uniprot_id)
    jsonString = json.dumps(human_genome, indent=4) 

    if write_dir:
        with open(os.path.join(write_dir,'domains_features_human_genome.json'), 'w') as json_file:
            json.dump(human_genome, json_file, indent=4)
        if uncovered_info:
            with open(os.path.join(write_dir,'domains_features_human_genome_uncovered.json'), 'w') as json_file:
                json.dump(uncovered_info, json_file, indent=4)
        if error_list:
            with open(os.path.join(write_dir,'domains_features_human_genome_error_ids.txt'), 'w') as json_file:
                json_file.write(f'{error_list}')
    if uncovered_info:
        logger.info(f'uncovered info by pfam: {uncovered_info}')
    if error_list:
        logger.info(f'ERROR: {error_list}')
    return human_genome
