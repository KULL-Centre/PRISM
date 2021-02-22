#!/groups/sbinlab/software/miniconda3/bin/python3

# Copyright (C) 2021 Johanna K. S. Tiemann <johanna.tiemann@gmail.com>

"""Script to create a prism file using the prism parser from a pdb file 

"""

# Standard library imports
from argparse import ArgumentParser
from datetime import datetime
from functools import reduce
import logging as log
import os
import shutil
import sys


# Third party imports
from Bio.PDB import PDBParser
from Bio.PDB.PDBList import PDBList 
from Bio.PDB.DSSP import DSSP
import pandas as pd

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

    from prism_parser_helper import write_prism

    
except (ModuleNotFoundError, ImportError) as e:
    logger.error("{} fileure".format(type(e)))
else:
    logger.info("Import succeeded")


def parse_args():
    """
    Argument parser function
    """

    parser = ArgumentParser( description="" )

    # Initiating/setup command line arguments
    parser.add_argument( 'pdbID',
        type=str,
        help="Input pdb id"
        )
    parser.add_argument( '--pdb_file', '-i',
        type=str,
        default='None',
        help="Optional input pdb file. Otherwise file will be downloaded."
        )
    parser.add_argument('--output_dir', '-o',
        type=str,
        default='.',
        help="Directory where files will be written. Default is the execution directory."
        )
    parser.add_argument('--chain', '-c',
        type=str,
        default='all',
        help="Select specific chains. Default is 'all'. Multiple chains are separated by a comma: A,B"
        )
    
    args = parser.parse_args()

    if args.pdb_file == 'None':
        args.pdb_file = None

    try:
        os.makedirs(args.output_dir, exist_ok = True)
        logger.info(f"Directory {args.output_dir} created successfully")
    except OSError as error:
        logger.warning(f"Directory {args.output_dir} can not be created")
    return args


def download_pdb(pdb_id, output_dir='.'):
    pdbl = PDBList()
    pdbl.retrieve_pdb_file(pdb_id, file_format='pdb', pdir=output_dir)
    pdb_path = os.path.join(output_dir, f'{pdb_id}.pdb')
    shutil.move(os.path.join(output_dir, f'pdb{pdb_id.lower()}.ent'), pdb_path)
    return pdb_path


def make_dssp_df(pdb_file, pdbID=None, chain='all'):

    pdb_p = PDBParser()
    if not pdbID:
        pdbID = os.path.basename(pdb_file)
    structure = pdb_p.get_structure(pdbID, pdb_file)
    model = structure[0]
    dssp = DSSP(model, pdb_file)
    
    variant_list = [['variant', 'SS', 'ASA', 'chain']]
    if chain == 'all':
        chains = list(set([key[0] for key in dssp.keys() ]))
    else:
        chains = chain.split(',')

    for chain in chains:
        for key in dssp.keys():
            if key[0] == chain:
                arr = [None] * len(variant_list[0])
                arr[0] = f'{dssp[key][1]}{key[1][1]}='
                # secondary structure
                arr[1] = dssp[key][2] 
                # accessible surface area
                arr[2] = dssp[key][3]
                arr[3] = chain
                variant_list.append(arr)
    dssp_df = pd.DataFrame(data = variant_list[1:], columns=variant_list[0])

    variant_dic = { 'SS': 'Secondary structure with H=alpha helix (4-12), B=isolated beta-bridge residue, E=Strand, G=3-10 helix, i=Pi helix, T=Turn, S=Bend, -=None', 
                    'ASA':'Accessible surface area',}

    return dssp_df, variant_dic


def merge_dfs(df_array):
    return reduce(lambda  left,right: pd.merge(left,right,on=['variant', 'chain'],
                                            how='outer'), df_array).fillna('void')


def merge_dics(meta_dict_array):
    return dict(j for i in meta_dict_array for j in i.items())


def dfs_to_prism(df, meta, pdbID, output_dir='.', organism=None, uniprot_id=None, version=1, chain='all'):
    #create prism file with all available info
    os.makedirs(output_dir, exist_ok = True)
    prism_file_list = []
    if chain == 'all':
        chains = df['chain'].unique()
    else:
        chains = chain.split(',')
    for chain in chains:
        output_df = df.copy()
        output_df = output_df[output_df['chain']==chain]
        output_df = output_df.drop(columns=['chain'])
        output_df = output_df.dropna(axis=1, how='all')
        output_df = output_df.reset_index(drop=True)
        
        prism_file = os.path.join(output_dir, f'prism_pdb_XXX_{pdbID}_{chain}.txt')            
        
        logger.info('Get sequence')
        tmp_df = output_df.copy()
        tmp_df['resnum'] = tmp_df['variant'].apply(lambda x: int(x[1:-1]))
        tmp_df['var'] = tmp_df['variant'].apply(lambda x: x[0])
        tmp_df['resnum'].min()
        if tmp_df['resnum'].min() <= 0:
            first_residue_number = tmp_df['resnum'].min()
            add = tmp_df['resnum'].min()
        else:
            add = -1
            first_residue_number = 1

        seq_array = ['-']*(tmp_df['resnum'].max() + add + 1)
        for i, row in tmp_df.iterrows():
            seq_array[row['resnum']+add] = row['var']
        sequence = "".join(seq_array)
        
        if first_residue_number == 1:
            for ind, elem in enumerate(sequence):
                first_residue_number = ind+1
                if elem != '-':
                    sequence=sequence[ind:]
                    break

        logger.info('Generate metadata')
        metadata = {
            "version": str(version),
            "protein": {
                "name": f'{pdbID}_{chain}',
                "organism": organism,
                "uniprot": uniprot_id,
                "sequence": sequence,
            },
            "pdb": {
                'chain': chain,
                'pdbID': pdbID
            },
            "columns": meta,
        }
        
        if first_residue_number != 1:
            metadata['protein']['first_residue_number'] = first_residue_number
        comment = [ f"version {version} - {datetime.date(datetime.now())} - pdb info script",] 

        try:
            write_prism(metadata, output_df, prism_file, comment=comment)
            prism_file_list.append(prism_file)
            logger.info(f'chain {chain} pdb prism file written!')
        except Exception as e:
            logger.info(f'Problem while writing prism file for chain {chain}: {e}')

    return prism_file_list

def pdb_to_prism(pdbID, pdb_file=None, output_dir='.', chain='all'):

    if not pdb_file:
        logger.info(f'PDB {pdbID} will be downloaded')
        pdb_file = download_pdb(pdbID, output_dir=output_dir)
    else:
        pdb_file = os.path.abspath(pdb_file)

    logger.info('Extract information from DSSP')
    dssp_df, dssp_meta_dict = make_dssp_df(pdb_file, pdbID, chain=chain)

    # for future options, if other databases are consulted, you can use this merge option
    merged_dic = merge_dics([dssp_meta_dict])
    merged_df = merge_dfs([dssp_df])

    logger.info('Generate prism file')
    file_list = dfs_to_prism(merged_df, merged_dic, pdbID, output_dir=output_dir, chain=chain)

    return file_list



def main():
    """
    Main function called as default at the end.
    """
    # get user input arguments
    args = parse_args()
    
    # generate pdb prism files
    file_list = pdb_to_prism(args.pdbID, pdb_file=args.pdb_file, output_dir=args.output_dir, chain=args.chain)



if __name__ == '__main__':
    main()
