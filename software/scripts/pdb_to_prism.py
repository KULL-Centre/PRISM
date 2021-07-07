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
import subprocess
import sys
import time


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
    from FillVariants import copy_wt_variants

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
    parser.add_argument('--fill', '-f',
        type=str,
        default='False',
        help="Fill all variants according to WT. Default is 'False'."
        )
    
    
    args = parser.parse_args()

    if args.pdb_file == 'None':
        args.pdb_file = None
    args.fill = bool(args.fill)
    try:
        os.makedirs(args.output_dir, exist_ok = True)
        logger.info(f"Directory {args.output_dir} created successfully")
    except OSError as error:
        logger.warning(f"Directory {args.output_dir} can not be created")
    return args


def calc_contacts(infile, pdbID, chain, tmp_base_dir, uniprot_id='', exec_path='/groups/sbinlab/tiemann/repos/getcontacts/'):
    time_stamp = time.time()
    tmp_dir_sub = os.path.join(tmp_base_dir, f'tmp_deleted_{time_stamp}')
    os.makedirs(tmp_base_dir, exist_ok = True)
    os.makedirs(tmp_dir_sub, exist_ok = True)

    contact_tsv = os.path.join(tmp_dir_sub, f"contacts_{uniprot_id}_{pdbID}_{chain}.tsv")
    submitter = f'python3 {os.path.join(exec_path, "get_static_contacts.py")} --structure {infile} --output {contact_tsv} --itypes all'
    pipes = subprocess.Popen(submitter, shell=True, cwd=tmp_dir_sub,stdout=subprocess.PIPE,stderr=subprocess.PIPE,)
    std_out, std_err = pipes.communicate()

    resfreq_tsv = os.path.join(tmp_dir_sub, f"resfreq_{uniprot_id}_{pdbID}_{chain}.tsv")
    submitter = f'python3 {os.path.join(exec_path, "get_contact_frequencies.py")} --input_files {contact_tsv} --output {resfreq_tsv} --itypes all'
    pipes = subprocess.Popen(submitter, shell=True, cwd=tmp_dir_sub, stdout=subprocess.PIPE,stderr=subprocess.PIPE,)
    std_out, std_err = pipes.communicate()

    df = pd.read_csv(resfreq_tsv, sep='\t', skiprows=2, names=['resi_1', 'resi_2', 'frame_count','contact_freq'])
    def convert_resi(x):
        d3d1={'ALA':'A', 'ARG':'R', 'ASP':'D', 'ASN':'N', 'CYS':'C', 'GLY':'G', 'GLU':'E', 'GLN':'Q', 'HIS': 'H', 'ILE': 'I',
             'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}
        x = x.split(':')
        wt = x[1].split('_')[0]
        resi = x[2].split('_')[-1]
        return f"{d3d1[wt]}{resi}="
    df['variant'] = df['resi_1'].apply(lambda x: convert_resi(x))
    df['variant_2'] = df['resi_2'].apply(lambda x: convert_resi(x)[:-1])
    df = df.drop(columns=['resi_1', 'resi_2', 'frame_count', 'contact_freq']).sort_values(by=['variant']).reset_index(drop=True)
    df_pvt = df.pivot_table(index = 'variant', aggfunc='count')
    df_pvt = df_pvt.reset_index(drop=False)
    df_pvt = df_pvt.rename(columns={'variant_2': 'contact;count'})

    df = df.pivot_table(index = 'variant', aggfunc={'variant_2':lambda x: "|".join(x)})
    df = df.reset_index(drop=False)
    df = df.rename(columns={'variant_2': 'contact;res'})
    df = pd.merge(df_pvt, df, on=['variant'])

    df2 = pd.read_csv(contact_tsv, sep='\t', skiprows=2, names=['frame','interaction_type', 'resi_1' 'resi_2']).reset_index(drop=False)
    df2 = df2.rename(columns={'frame': 'interaction_type', 'interaction_type': 'resi_1', 'resi_1resi_2': 'resi_2'}).drop(columns='index')
    df2['variant'] = df2['resi_1'].apply(lambda x: convert_resi(x))
    df2['variant_2'] = df2['resi_2'].apply(lambda x: convert_resi(x))
    df2 = df2.drop(columns=['resi_1', 'resi_2']).sort_values(by=['variant']).reset_index(drop=True)
    df2_pvt = df2.pivot_table(index = 'variant', aggfunc='count', columns = 'interaction_type')
    df2_pvt.columns = df2_pvt.columns.droplevel()
    df2_pvt = df2_pvt.reset_index(drop=False)
    df2_pvt.columns.name = None
    columns_types_raw = list(df2_pvt.columns)
    columns_types = list(df2_pvt.columns)
    columns_types.remove('vdw')
    df2_pvt['total_no_vdw'] = df2_pvt[columns_types].sum(axis=1)
    df2_pvt['total'] = df2_pvt[columns_types_raw].sum(axis=1)
    df2 = df2.pivot_table(index = 'variant', aggfunc={'variant_2':lambda x: "|".join(x)})
    df2 = df2.reset_index(drop=False)
    df2 = df2.rename(columns={'variant_2': 'contact;res'})
    df2 = pd.merge(df2_pvt, df2, on=['variant'])

    shutil.rmtree(tmp_dir_sub)

    return df


def rosetta_energy_to_prism(infile, prism_file, pdbID, chain, tmp_base_dir, uniprot_id='', organism='', version=1, count=True):
    df_list = []
    with open(infile, 'r') as fp:
        started=False
        for line in fp:
            if not started:
                if line.startswith('#BEGIN_POSE_ENERGIES_TABLE'):
                    started=True
            else:
                if (line.startswith('#END_POSE_ENERGIES_TABLE') or line.startswith('MEM')):
                    break
                elif not (line.startswith('pose') or line.startswith('weights')):
                    df_list.append(line.split())
                else:
                    pass
    df = pd.DataFrame(data=df_list[1:], columns=df_list[0])
    
    def convert_resi(x):
        d3d1={'ALA':'A', 'ARG':'R', 'ASP':'D', 'ASN':'N', 'CYS':'C', 'GLY':'G', 'GLU':'E', 'GLN':'Q', 'HIS': 'H', 'ILE': 'I',
             'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}
        x = x.split(':')
        if len(x)>1:
            wt = x[0].split('_')[0]
            resi = x[-1].split('_')[-1]
        else:
            wt = x[0].split('_')[0]
            resi = x[0].split('_')[-1]
        return f"{d3d1[wt]}{resi}="
        
    
    df['variant'] = df['label'].apply(lambda x: convert_resi(x))
    df['resi'] = df['variant'].apply(lambda x: int(x[1:-1]))
    df['wt'] = df['variant'].apply(lambda x: x[0])
    max_val = df['resi'].max()
    len_seq = len(df['resi'])
    if max_val >=len_seq:
        sequence = ['X']*(max_val+1)
    else:
        sequence = ['X']*(len_seq+1)
    varis = df['variant'].unique()
    for var in varis:
        sequence[int(var[1:-1])-1] = var[0]
    sequence = "".join(sequence)
    df = df.drop(columns=['label', 'resi', 'wt'])
    df = df.dropna(axis=1, how='all')
    df = df.reset_index(drop=True)
    df = df.add_prefix('energy;')
    df = df.rename(columns={'energy;variant': 'variant'})
    
    if count:
        df_count = calc_contacts(infile, pdbID, chain, tmp_base_dir, uniprot_id=uniprot_id, exec_path='/groups/sbinlab/tiemann/repos/getcontacts/')
        df = pd.merge(df, df_count, on=['variant']).reset_index(drop=True)
    
    
    meta = {}
    for column in df.columns:
        if not column in ['variant', 'n_mut', 'aa_ref', 'resi']:
            meta[column] = column
    df = df[['variant']+list(meta.keys())]

    metadata = {
        "version": str(version),
        "protein": {
            "name": f'{pdbID}_{chain}',
            "organism": organism,
            "uniprot": uniprot_id,
            "sequence": sequence,
        },
        "rosettapdb": {
            'chain': chain,
            'pdbID': pdbID
        },
        "columns": meta,
    }
    comment = [ f"version {version} - {datetime.date(datetime.now())} - rosetta_pdb info script",] 
    write_prism(metadata, df, prism_file, comment=comment)


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

        seq_array = ['X']*(tmp_df['resnum'].max() + add + 1)
        for i, row in tmp_df.iterrows():
            seq_array[row['resnum']+add] = row['var']
        sequence = "".join(seq_array)
        
        if first_residue_number == 1:
            for ind, elem in enumerate(sequence):
                first_residue_number = ind+1
                if elem != 'X':
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

def pdb_to_prism(pdbID, pdb_file=None, output_dir='.', chain='all', fill=False):

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

    if fill:
        logger.info('Filling the dataframes')
        merged_df = copy_wt_variants(merged_df)
    
    logger.info('Generate prism file')
    file_list = dfs_to_prism(merged_df, merged_dic, pdbID, output_dir=output_dir, chain=chain)

    return file_list


def pdb_renumb(pdb_input, output_dir=None, keepchain='all', chainorder=None, keep_ligand=None):
    def run_through_lines(pdb_input, keepchain, fp2, resnum=1, lastresstring="", keep_ligand=None, atom_num=1):
        with open(pdb_input, 'r') as fp:
            if keep_ligand:
                for line in fp:
                    if (line[17:20] == keep_ligand) and (line[:6] in ['HETATM', 'ATOM  ']):
                        resstring = line[22:27]
                        if lastresstring == "" or resstring != lastresstring :
                            if lastresstring != "" : 
                                resnum += 1
                            resnum = resnum
                            lastresstring = resstring
                        newresstring = str(resnum) + " "
                        if len(newresstring) == 2: 
                            newresstring = "   " + newresstring
                        elif len(newresstring) == 3: 
                            newresstring = "  " + newresstring
                        elif len(newresstring) == 4: 
                            newresstring = " " + newresstring
                        new_line = 'HETATM' + ' '*(5-len(str(atom_num))) + str(atom_num) +line[11:22] + newresstring + line[27:]
                        atom_num += 1
                        fp2.write(new_line)
            else:
                for line in fp:
                    if line[0:4] == "ATOM":
                        chain = line[21]
                        if chain == keepchain or keepchain=='all':
                            resstring = line[22:27]
                            if lastresstring == "" or resstring != lastresstring :
                                if lastresstring != "" : 
                                    resnum += 1
                                resnum = resnum
                                lastresstring = resstring
                            newresstring = str(resnum) + " "
                            if len(newresstring) == 2: 
                                newresstring = "   " + newresstring
                            elif len(newresstring) == 3: 
                                newresstring = "  " + newresstring
                            elif len(newresstring) == 4: 
                                newresstring = " " + newresstring
                            new_line = 'ATOM  ' + ' '*(5-len(str(atom_num))) + str(atom_num) +line[11:22] + newresstring + line[27:]
                            atom_num += 1
                            fp2.write(new_line)
        return resnum, lastresstring, atom_num
    
    resnum = 1
    lastresstring = ""
    atom_num = 1
    if output_dir:
        out_pdb = os.path.join(output_dir, f'{os.path.basename(pdb_input)[:-4]}_renum.pdb')
    else:
        out_pdb = os.path.join(os.path.dirname(pdb_input), f'{os.path.basename(pdb_input)[:-4]}_renum.pdb')
    with open(out_pdb, 'w') as fp2:
        if chainorder:
            for ref_chain in chainorder.split(','):
                if (ref_chain in keepchain) or (keepchain == 'all'):
                    resnum, lastresstring, atom_num = run_through_lines(pdb_input, ref_chain, fp2, resnum=resnum, lastresstring=lastresstring)
        else:
            resnum, lastresstring, atom_num = run_through_lines(pdb_input, keepchain, fp2, resnum=resnum, lastresstring=lastresstring, atom_num=atom_num)
        if keep_ligand:
            resnum, lastresstring, atom_num = run_through_lines(pdb_input, keepchain, fp2, resnum=resnum, lastresstring=lastresstring, atom_num=atom_num, keep_ligand=keep_ligand)
    return out_pdb


def main():
    """
    Main function called as default at the end.
    """
    # get user input arguments
    args = parse_args()
    
    # generate pdb prism files
    file_list = pdb_to_prism(args.pdbID, pdb_file=args.pdb_file, output_dir=args.output_dir, chain=args.chain, fill=args.fill)



if __name__ == '__main__':
    main()
