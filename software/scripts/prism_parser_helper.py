#!/groups/sbinlab/software/miniconda3/bin/python3

# Copyright (C) 2021 Johanna K. S. Tiemann <johanna.tiemann@gmail.com>

"""Script to create a prism file using the prism parser from a pdb file 

"""

# Standard library imports
import logging as log
import os
import shutil
import sys


# Third party imports


# Local application imports
log_message="verbose"

if log_message=="verbose":
    log.basicConfig(
        format='%(levelname)s:%(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=log.INFO
    )
elif log_message=="debug":
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
    
try:
    if os.uname()[1].startswith('fend'):
        base_script_path = '/groups/sbinlab/tiemann/repos/PRISM/'
    else:
        base_script_path = '/storage1/tiemann/dev/repos/'

    sys.path.insert(1, os.path.join(base_script_path, 'prism/scripts/'))
    from PrismData import PrismParser, VariantData

except (ModuleNotFoundError, ImportError) as e:
    logger.error("{} fileure".format(type(e)))

else:
    logger.info("Import succeeded")



def write_prism(metadata, dataframe, prism_file, comment=''):
    variant_dataset = VariantData(metadata, dataframe)
    parser = PrismParser()
    parser.write(prism_file, variant_dataset, comment_lines=comment)


def read_from_prism(primsfile):
    parser = PrismParser()
    dataframe = parser.read(primsfile).dataframe
    meta_data = parser.read_header(primsfile)
    return meta_data, dataframe


def read_prism(primsfile):
    parser = PrismParser()
    data = parser.read(primsfile)
    return data


def merge_prism(filenames, output_dir=None, identity=0.9, merge='outer', verbose=False):
    parser = PrismParser()
    data_list = []
    first_resn = 1
    for file in filenames:
        data = parser.read(file)
        data_list.append(data)
    target_seq = (data_list[0]).metadata['protein']['sequence']
    if 'first_residue_number' in (data_list[0]).metadata['protein']:
        first_resn = int((data_list[0]).metadata['protein']['first_residue_number'])
        
    merged_data = data_list[0].merge(data_list[1:], target_seq, first_resn, merge=merge, 
                                     min_identity=identity, min_coverage=.01, mismatch="remove", 
                                     allow_inserts=True, allow_deletions=True, verbose=verbose)

    prism_file = None
    if output_dir:
        prism_file = os.path.join(output_dir, 'prism_merged_files.txt')
        write_prism(merged_data.metadata, merged_data.dataframe, prism_file, comment='')

    return merged_data, prism_file


def merge_prism_right(filenames, output_dir=None, identity=0.9, verbose=False):
    parser = PrismParser()
    data_list = []
    first_resn = 1
    for file in filenames:
        data = parser.read(file)
        data_list.append(data)
    data_list[0].dataframe = data_list[0].dataframe[['variant', 'aa_ref', 'resi', 'aa_var', 'n_mut']]#[0:0]
    data_list[0].metadata['columns'] = {}
    target_seq = (data_list[0]).metadata['protein']['sequence']
    if 'first_residue_number' in (data_list[0]).metadata['protein']:
        first_resn = int((data_list[0]).metadata['protein']['first_residue_number'])

    merged_data = data_list[0].merge(data_list[1:], target_seq, first_resn, merge='outer', 
                                     min_identity=identity, min_coverage=.01, mismatch="remove", 
                                     allow_inserts=True, allow_deletions=True, verbose=verbose)
    merged_data.dataframe = merged_data.dataframe.dropna(subset=merged_data.metadata['columns'].keys(), how='all').reset_index(drop=True)    
    prism_file = None
    if output_dir:
        prism_file = os.path.join(output_dir, 'prism_merged_files.txt')
        write_prism(merged_data.metadata, merged_data.dataframe, prism_file, comment='')
        
    return merged_data, prism_file
        

def merge_prism_df(data_list, identity=0.9):
    target_seq = (data_list[0]).metadata['protein']['sequence']
    if 'first_residue_number' in (data_list[0]).metadata['protein']:
        first_resn = int((data_list[0]).metadata['protein']['first_residue_number'])
    else:
        first_resn = 1
        
    merged_data = data_list[0].merge(data_list[1:], target_seq, first_resn, merge='outer', 
                                     min_identity=identity, min_coverage=.01, mismatch="remove", 
                                     allow_inserts=True, allow_deletions=True)
    return merged_data


def stats_prism(out_csv):
    shell_command = ['python', os.path.join(path_to_prism_parser, 'scripts/PrismData.py'), 
                     os.path.join(path_to_prism_parser, 'prism_ddG/prism_ddG_*'), 
                     '--dump_header_csv', out_csv]
    
    subprocess.run(shell_command,
                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,
                   check=True, text=True)
    pandas_df = pd.read_csv(out_csv, header=1, sep=';')
    return pandas_df
