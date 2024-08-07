"""prism_rosetta_parser.py contains functions to convert prism to mut & ddG to prism

Author: Johanna K.S. Tiemann
Date of last major changes: 2020-05-01

"""

# Standard library imports
from datetime import datetime
import json
import logging as logger
import re
import subprocess
import sys
import time


# Third party imports
import pandas as pd


# Local application imports
import rosetta_paths
sys.path.insert(1, rosetta_paths.prism_parser)
from PrismData import PrismParser, VariantData


def span_multi(x, region):
    out_tms = []
    for resi in x:
        if int(resi) in region:
            out_tms.append('True')
        else:
            out_tms.append('False')
    return ":".join(out_tms)

def deepTMHMM_multi(x, deepTM_start, deepTM_end, info_name):
    print(x)
    if not x in ['outside', 'inside', 'periplasm', 'Beta sheet', 'Beta_sheet', 'TMhelix', 'signal']:
        resi = x[0]
        if (int(resi) >= int(deepTM_start)) and (int(resi) <= int(deepTM_end)):
            return "_".join(info_name.split(' '))
        else:
            return x
    else:
        return "_".join(x.split(' '))


def rosetta_to_prism(ddg_file, prism_file, sequence, rosetta_info=None, version=1, sys_name='', 
    first_residue_number=1, sha_tag='', MP=False, span_file='', lipid_file='', deepTMHMM_file='', scale=2.9):
    
    sequence = sequence.replace('-', 'X')
    # create prism file with rosetta values
    logger.info('Create prism file with rosetta ddG values')
    variant = []
    norm_ddG_value = []
    std_ddG_value = []
    ddG_value = []
    with open(ddg_file, 'r') as fp:
        for line in fp:
            split_line = line.split(',')
            variant.append(split_line[0])
            norm_ddG_value.append(split_line[1])
            std_ddG_value.append(split_line[2].strip())
     #       ddG_value.append(split_line[1])

    def split_multiple(variant):
        aa_ref = []
        resi = []
        aa_var = []
        for variants in variant:
            tmp_aa_ref = []
            tmp_resi = []
            tmp_aa_var = []
            for seg in variants.split(':'):
                tmp_aa_ref.append(seg[0])
                tmp_resi.append(seg[1:-1])
                tmp_aa_var.append(seg[-1])
            aa_ref.append(tmp_aa_ref)
            resi.append(tmp_resi)
            aa_var.append(tmp_aa_var)
        return aa_ref, resi, aa_var

    aa_ref, resi, aa_var = split_multiple(variant)
    data = {
        'variant': pd.Series(variant),
        'mean_ddG': pd.Series(norm_ddG_value),
        'std_ddG': pd.Series(std_ddG_value),
        #      'ddG': pd.Series(ddG_value),
        'n_mut': 1,  # pd.Series([1 for x in range(len(variant))]),
        'aa_ref': aa_ref,
        "resi": resi,
        "aa_var": aa_var,
    }
    dataframeset = pd.DataFrame(data)

    if span_file!='':
        TM_regions = []
        with open(span_file) as fp:
            next(fp)
            next(fp)
            next(fp)
            next(fp)
            for line in fp:
                line = line.strip().split()
                TM_regions += range(int(line[0]), int(line[1])+1)
        dataframeset['TMspan'] = dataframeset['resi'].apply(lambda x: span_multi(x, TM_regions))

    if lipid_file!='':
        with open(lipid_file, 'r') as fp:
            data = json.load(fp)
        lipid_df = pd.DataFrame.from_dict( data, orient='index', columns=['LAR']).reset_index(drop=False)#.T.set_index('index')
        lipid_df = lipid_df.loc[lipid_df['LAR']=='true'].reset_index(drop=True)
        lipid_df = lipid_df['index'].astype(int).unique()
        dataframeset['LAR'] = dataframeset['resi'].apply(lambda x: span_multi(x, lipid_df))

    if deepTMHMM_file!='':
        deepTMHMM_df = pd.read_csv(deepTMHMM_file, header=0)
        dataframeset['deepTM'] = dataframeset['resi']
        for idx, row in deepTMHMM_df.iterrows():
            deepTM_start = row['start']
            deepTM_end = row['end']
            info_name = row['location']
            dataframeset['deepTM'] = dataframeset['deepTM'].apply(lambda x: deepTMHMM_multi(x, deepTM_start, deepTM_end, info_name))

    sha = sha_tag.split('tag')[0]
    tag = sha_tag.split('tag')[1]

    if rosetta_info == None:
        rosetta_info = {
            "version": f"{tag} ({sha})",
        }

    if MP:
        if scale==1.0:
            units = f') REU'
        else:
            units = f'/{scale}) REU'
    else:
        units = f'/{scale}) kcal/mol'

    metadata = {
        # new version always when the data is edited or major changes are made
        # to metadata
        "version": version,
        # extracted from uniprot
        "protein": {
            "name": sys_name,
            "organism": "",
            "uniprot": "",
            "sequence": sequence,
        },
        "rosetta": rosetta_info,
        # columns is dependent on the data. different conditions go into
        # different PRISM_files
        "columns": {
            "mean_ddG": f"mean Rosetta ddG values (mean((MUT-mean(WT)){units})",
            "std_ddG": f"std Rosetta ddG values (std((MUT-mean(WT)){units})",
        },
    }
    if span_file!='':
        metadata['columns']['TMspan'] = 'Residue within the TM region defined by the Rosetta span file'
    if lipid_file!='':
        metadata['columns']['LAR'] = 'Lipid accessible residue defined by Rosetta'
    if deepTMHMM_file!='':
        metadata['columns']['deepTM'] = 'Protein regions defined by deepTMHMM'

    if first_residue_number != 1:
        metadata['protein']['first_residue_number'] = first_residue_number

    comment = [
        f"version {version} - {datetime.date(datetime.now())} - ddG pipeline",
    ]

    variant_dataset = VariantData(metadata, dataframeset)
    parser = PrismParser()
    parser.write(prism_file, variant_dataset, comment_lines=comment)


def write_prism(metadata, dataframe, prism_file, comment=''):
    variant_dataset = VariantData(metadata, dataframe)
    parser = PrismParser()
    parser.write(prism_file, variant_dataset, comment_lines=comment)


def read_prism(prism_file):
    parser = PrismParser()
    data = parser.read(prism_file)
    return data
