#!/groups/sbinlab/software/miniconda3/bin/python3

# Copyright (C) 2021 Johanna K. S. Tiemann <johanna.tiemann@gmail.com>

"""Script to create rosetta mut-files from a prism file, optional via a pdb 

"""

# Standard library imports
from argparse import ArgumentParser
import logging as log
import os
import sys


# Third party imports
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
        base_script_path = '/storage1/tiemann/dev/repos'

    sys.path.insert(1, os.path.join(base_script_path, 'prism/scripts/'))
    from pdb_to_prism import pdb_to_prism, pdb_renumb
    from prism_parser_helper import merge_prism, read_from_prism, merge_prism_right

except (ModuleNotFoundError, ImportError) as e:
    logger.error("{} fileure".format(type(e)))

else:
    logger.info("Import succeeded")


AA=['A','C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
       'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    
def parse_args():
    """
    Argument parser function
    """

    parser = ArgumentParser( description="" )

    # Initiating/setup command line arguments
    parser.add_argument( 'prism',
        type=str,
        help="Input prism file"
        )
    parser.add_argument( '--pdb_file', '-i',
        type=str,
        default='None',
        help="Optional input pdb file. If provided, the input prism file will be aligned and only overlapping AA are selected. Multi mutants are only kept if all residues are present in pdb. First chain will be used for alignment."
        )
    parser.add_argument('--pdb2rosetta', '-r',
        default=False,
        type=lambda s: s.lower() in ['true', 't', 'yes', '1'],
        help="Converts the provided PDB to Rosetta numbering before alignment"
        )
    parser.add_argument('--output_dir', '-o',
        type=str,
        default='.',
        help="Directory where files will be written. Default is the execution directory."
        )
    parser.add_argument('--drop_multi_mut', '-d',
        default=False,
        type=lambda s: s.lower() in ['true', 't', 'yes', '1'],
        help="Removes all multi mutants"
        )
    parser.add_argument('--inclWT', '-w',
        default=True,
        type=lambda s: s.lower() in ['false', 'f', 'no', '0'],
        help="Includes/adds WT variants"
        )
    parser.add_argument('--resisatu', '-s',
        default=False,
        type=lambda s: s.lower() in ['true', 't', 'yes', '1'],
        help="Residual saturation for all residues which have known variants"
        )
    parser.add_argument('--mut_output', '-m',
        choices=['all', 'pipeline_combined_mut', 'combined_mut',
                 'residue_files'],
        default='all',
        help=('Define what mutfile output should be created:\n'
              '\tall: all 3 options are created \n'
              '\tpipeline_combined_mut: combined mutation input file (no Rosetta input) for current stability pipeline \n'
              '\tcombined_mut: one single Rosetta mutfile including all mutations \n'
              '\tresidue_files: one Rosetta mutfile for each residue saved in a single file in a created directory \n'
              'Default value: all'
              )
        )
    
    args = parser.parse_args()

    if args.pdb_file == 'None':
        args.pdb_file = None
    args.inclWT = bool(args.inclWT)
    args.resisatu = bool(args.resisatu)
    args.drop_multi_mut = bool(args.drop_multi_mut)
    args.pdb2rosetta = bool(args.pdb2rosetta)
    if args.pdb2rosetta==True and args.pdb_file==None:
        logger.warning(f'pdb2rosetta set to True but no PDB provided. This step will be skiped')
        args.pdb2rosetta = False
    try:
        os.makedirs(args.output_dir, exist_ok = True)
        logger.info(f"Directory {args.output_dir} created successfully")
    except OSError as error:
        logger.warning(f"Directory {args.output_dir} can not be created")
    return args


def make_satu(aa_wt_comb, resi_comb, inclWT=False):
    print(aa_wt_comb, resi_comb)
    AA_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    return_list = []
    for AA in AA_list:
        if not ((inclWT==False) and (AA==aa_wt_comb)):
            return_list.append(f'{aa_wt_comb}{resi_comb}{AA}')
    return ";".join(return_list)


def prepare_df(in_df, inclWT=True, drop_multi_mut=False, resisatu=False):
    df = in_df.copy()
    if drop_multi_mut:
        df = df[df['n_mut']==1]
    # drop WT variants and deletions
    df['drop'] = df['aa_var'].apply(lambda x: "".join(x) )
    drop = df[df['drop'].str.contains('~')].index
    df= df.drop(drop)

    df['resi_comb'] = df['resi'].apply(lambda x: ";".join([str(elem) for elem in x]))
    df['aa_wt_comb'] = df['aa_ref'].apply(lambda x: ";".join(x))
    df['aa_var_comb'] = df['aa_var'].apply(lambda x: ";".join(x))
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
        df = drops.copy()
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

    if resisatu:
        print(df)
        df['variant'] = df.apply(lambda x: make_satu(x['aa_wt_comb'], x['resi_comb'], inclWT=inclWT), axis=1)
        if inclWT:
            df['count'] = 20
            df['sum_muts'] = 20
            df['sum_incl_wt'] = 20.0
        else:
            df['count'] = 19
            df['sum_muts'] = 19
            df['sum_incl_wt'] = 20.0

    return df, drops


def residue_files_func(output_dir, in_df, inclWT=True, drop_multi_mut=False, resisatu=False):

    df, drops = prepare_df(in_df, inclWT=inclWT, drop_multi_mut=drop_multi_mut, resisatu=resisatu)

    mut_dir = os.path.join(output_dir, 'mut_files')
    os.makedirs(mut_dir, exist_ok = True)
    for i, row in df.iterrows():
        if len(drops) !=0:
            row2 = drops[drops['variant']==row['variant']]
        if len(row['resi_comb'].split(';')) > 1:
            naming = "_".join(row['resi_comb'].split(';'))
        else:
            naming = "0"*(5-len(row['resi_comb'])) + row['resi_comb']
        mutfile = os.path.join(mut_dir, f'mutfile{naming}')
        logger.info(mutfile)
        with open(mutfile, 'w') as fp:
            if inclWT:
                if len(drops) != 0:
                    total_numb = row2['sum_incl_wt'].astype(int)
                else:
                    total_numb = row['sum_muts'].astype(int)
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


def combined_mut_func(output_dir, in_df, inclWT=True, drop_multi_mut=False, resisatu=False):

    df, drops = prepare_df(in_df, inclWT=inclWT, drop_multi_mut=drop_multi_mut, resisatu=resisatu)

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


def pipeline_combined_mut_func(output_dir, in_df, inclWT=True, drop_multi_mut=False, resisatu=False):

    df, drops = prepare_df(in_df, inclWT=inclWT, drop_multi_mut=drop_multi_mut, resisatu=resisatu)
    
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


def prism2mut(input_file, output_dir='.', pipeline_combined_mut=True, combined_mut=True, residue_files=False, inclWT=True, resisatu=False, drop_multi_mut=False):

    os.makedirs(output_dir, exist_ok = True) 
    meta, df = read_from_prism(input_file)
    df = df[['variant', 'aa_ref', 'resi', 'aa_var', 'n_mut']]

    if residue_files:
        residue_files_func(output_dir, df, inclWT=inclWT, drop_multi_mut=drop_multi_mut, resisatu=resisatu)
    if combined_mut:
        combined_mut_func(output_dir, df, inclWT=inclWT, drop_multi_mut=drop_multi_mut, resisatu=resisatu)
    if pipeline_combined_mut:
        pipeline_combined_mut_func(output_dir, df, inclWT=inclWT, drop_multi_mut=drop_multi_mut, resisatu=resisatu)
        
                
def main():
    """
    Main function called as default at the end.
    """
    # get user input arguments
    args = parse_args()
    
    
    if args.pdb_file:
        if args.pdb2rosetta:
            args.pdb_file = pdb_renumb(args.pdb_file, output_dir=args.output_dir)
        # generate pdb prism files
        pdb_prism = pdb_to_prism('infile', pdb_file=args.pdb_file, output_dir=args.output_dir, fill=True)[0]
        
        meta_data, dataframe = read_from_prism(args.prism)
        if meta_data['variants']['width'] == 'single mutants':
            variantdata, args.prism = merge_prism([pdb_prism, args.prism], output_dir=args.output_dir, 
                                   identity=0.5, merge='inner', verbose=True)
        else:
            variantdata, args.prism = merge_prism_right([pdb_prism, args.prism], output_dir=args.output_dir, 
                                   identity=0.5, verbose=True)
    
    residue_files = False
    pipeline_combined_mut = False
    combined_mut = False
    if args.mut_output == 'pipeline_combined_mut':
        pipeline_combined_mut = True
    elif args.mut_output == 'combined_mut':
        combined_mut = True
    elif args.mut_output == 'residue_files':
        residue_files = True
    elif args.mut_output == 'all':
        residue_files = True
        pipeline_combined_mut = True
        combined_mut = True

    prism2mut(args.prism, output_dir=args.output_dir, pipeline_combined_mut=pipeline_combined_mut, 
              combined_mut=combined_mut, residue_files=residue_files, inclWT=args.inclWT, resisatu=args.resisatu, drop_multi_mut=args.drop_multi_mut)
    logger.info('Conversion sucessful!')

if __name__ == '__main__':
    main()
