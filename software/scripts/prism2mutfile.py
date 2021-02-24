#!/groups/sbinlab/software/miniconda3/bin/python3

# Copyright (C) 2021 Johanna K. S. Tiemann <johanna.tiemann@gmail.com>

"""Script to create rosetta mut-files from a prism file, optional via a pdb 

"""

# Standard library imports
import logging as log
import os
import pandas as pd
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
        base_script_path = '/groups/sbinlab/tiemann/repos/PRISM/'

    sys.path.insert(1, os.path.join(base_script_path, 'prism/scripts/'))
    from pdb_to_prism import pdb_to_prism, download_pdb, dfs_to_prism, pdb_renumb
    from prism_parser_helper import merge_prism, read_from_prism

except (ModuleNotFoundError, ImportError) as e:
    logger.error("{} fileure".format(type(e)))

else:
    logger.info("Import succeeded")



def prism2mut(input_file, output_dir='.', pipeline_combined_mut=True, combined_mut=True, residue_files=False, inclWT=True, drop_multi_mut=False):
    AA=['A','C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
       'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    
    os.makedirs(output_dir, exist_ok = True) 
    meta, df = read_from_prism(input_file)
    df = df[['variant', 'aa_ref', 'resi', 'aa_var', 'n_mut']]
    if drop_multi_mut:
        df = df[df['n_mut']==1]
    # drop WT variants and deletions
    df['drop'] = df['aa_var'].apply(lambda x: "".join(x) )
    drop = df[df['drop'].str.contains('~')].index
    df= df.drop(drop)

    df['resi_comb'] = df['resi'].apply(lambda x: ";".join([str(elem) for elem in x]))
    df['aa_wt_comb'] = df['aa_ref'].apply(lambda x: ";".join(x))
    df['count'] = 1

    
    if not inclWT:
        drop = df[df['drop'].str.contains('=')].index
        df= df.drop(drop)
        drop = df[df['aa_var'] == df['aa_ref']].index
        df= df.drop(drop)
        drops = df.copy()
    else:
        drops = df.copy()
        drop = drops[drops['drop'].str.contains('=')].index
        drops= drops.drop(drop)
        drop = drops[drops['aa_var'] == drops['aa_ref']].index
        drops= drops.drop(drop)
    if len(drops) !=0:
        drops = drops[['variant', 'resi_comb', 'aa_wt_comb', 'count']]
        drops = drops.groupby(['resi_comb', 'aa_wt_comb'],as_index=False).agg(lambda x : x.sum() if x.dtype in ['float64', 'int'] else ';'.join(x))
        drops['sum_muts'] = drops['variant'].apply(lambda x: len([item for subl in [i.split(':') for i in x.split(';')] for item in subl]))
        drops['sum_incl_wt'] = drops['sum_muts'] + drops['sum_muts']/drops['count']

    if len(df) == 0:
        logger.warning('No variants to process. Does it contain only WT and function switched on?')
        return
    df = df[['variant', 'resi_comb', 'aa_wt_comb', 'count']]
    df = df.groupby(['resi_comb', 'aa_wt_comb'],as_index=False).agg(lambda x : x.sum() if x.dtype in ['float64', 'int'] else ';'.join(x))
    df['sum_muts'] = df['variant'].apply(lambda x: len([item for subl in [i.split(':') for i in x.split(';')] for item in subl]))
    df['sum_incl_wt'] = df['sum_muts'] + df['sum_muts']/df['count']


    if residue_files:
        mut_dir = os.path.join(output_dir, 'mut_files')
        os.makedirs(mut_dir, exist_ok = True)
        for i, row in df.iterrows():
            if len(drops) !=0:
                row2 = drops[drops['variant']==row['variant']]
            if len(row['resi_comb']) > 5:
                naming = "_".join(row['resi_comb'].split(';'))
            else:
                naming = "0"*(5-len(row['resi_comb'])) + row['resi_comb']
            mutfile = os.path.join(mut_dir, f'mutfile{naming}')
            logger.info(mutfile)
            with open(mutfile, 'w') as fp:
                if inclWT:
                    if len(drops) != 0:
                        total_numb = int(row2['sum_incl_wt'])
                    else:
                        total_numb = int(row['sum_muts'])
                else:
                    total_numb = int(row['sum_muts'])
                fp.write(f"total {total_numb}\n")
                if inclWT:
                    fp.write(f"{len(row['aa_wt_comb'].split(';'))}\n")
                    for variant in row['variant'].split(';')[0].split(':'):
                        mutant_line = f"{variant[0]} {variant[1:-1]} {variant[0]}\n"
                        fp.write(mutant_line)
                if len(drops) !=0 and len(row2)!=0:
                    for variants in row['variant'].split(';'):
                        fp.write(f"{len(row['aa_wt_comb'].split(';'))}\n")
                        for variant in variants.split(':'):
                            if variant[-1] in AA:
                                mutant_line = f"{variant[0]} {variant[1:-1]} {variant[-1]}\n"
                                fp.write(mutant_line)
    if combined_mut:
        mutfile_all = os.path.join(output_dir, 'mutfile_all')
        with open(mutfile_all, 'w') as fp:
            if inclWT:
                total_numb = int(df['sum_incl_wt'].sum())
            else:
                total_numb = df['sum_muts'].sum()
            fp.write(f"total {total_numb}\n")
            for i, row in df.iterrows():
                if inclWT:
                    fp.write(f"{len(row['aa_wt_comb'].split(';'))}\n")
                    for variant in row['variant'].split(';')[0].split(':'):
                        mutant_line = f"{variant[0]} {variant[1:-1]} {variant[0]}\n"
                        fp.write(mutant_line)
                for variants in row['variant'].split(';'):
                    mutant_line = []
                    for variant in variants.split(':'):
                        if variant[-1] in AA:
                            mutant_line.append(f"{variant[0]} {variant[1:-1]} {variant[-1]}\n")
                    if len(mutant_line)>0:
                        fp.write(f"{len(mutant_line)}\n")
                        fp.write("".join(mutant_line))


    if pipeline_combined_mut:
        combined_mut_file = os.path.join(output_dir, 'pipeline_mutfile.txt')
        with open(combined_mut_file, 'w') as fp:
            for i, row in df.iterrows():
                outstring = row['variant']
                mutant_line = [variants.split(':') for variants in row['variant'].split(';')]
                resstring = []
                for resind, residue in enumerate(mutant_line[0]):
                    aa_ref = residue[0]
                    res = residue[1:-1]
                    muts = "".join([ mutant_line[l][resind][-1] for l in range(len(mutant_line)) if mutant_line[l][resind][-1] in AA])
                    if inclWT:
                        muts = aa_ref+muts
                    resstring.append(f'{aa_ref} {res} {muts}')
                fp.write(f"{' '.join(resstring)}\n")

                
                
