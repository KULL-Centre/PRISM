#!/groups/sbinlab/software/miniconda3/bin/python3

# Copyright (C) 2022 Johanna K. S. Tiemann <johanna.tiemann@gmail.com>

"""Script to create input for homo-oligomers 

"""

# Standard library imports
from argparse import ArgumentParser
import logging as log
import os
import sys


# Third party imports
from Bio.PDB import PDBParser
from Bio import pairwise2
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
    from prism_parser_helper import read_from_prism
    from pdb_to_prism import pdb_renumb
    from prism2mutfile import residue_files_func, combined_mut_func, pipeline_combined_mut_func

except (ModuleNotFoundError, ImportError) as e:
    logger.error("{} fileure".format(type(e)))

else:
    logger.info("Import succeeded")


d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'BNG': 'X'}


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
    parser.add_argument( 'pdb_file',
        type=str,
        help="Dimer pdb files with each monomer having a unique chain ID"
        )
    parser.add_argument('--output_dir', '-o',
        type=str,
        default='.',
        help="Directory where files will be written. Default is the execution directory."
        )
    parser.add_argument('--additional', '-d',
        default=False,
        type=lambda s: s.lower() in ['true', 't', 'yes', '1'],
        help=('Will create additional mutfiles:\n'
              '\tpipeline_combined_mut: combined mutation input file (no Rosetta input) for current stability pipeline \n'
              '\tresidue_files: one Rosetta mutfile for each residue saved in a single file in a created directory \n'
              'Default value: False'
              )
        )
    parser.add_argument('--keep_ligand', '-l',
        default=True,
        type=lambda s: s.lower() in ['false', 'f', 'no', '0'],
        help="keep ligands/HETATM"
        )
    parser.add_argument('--inclWT', '-w',
        default=True,
        type=lambda s: s.lower() in ['false', 'f', 'no', '0'],
        help="Includes/adds WT variants"
        )

    
    args = parser.parse_args()

    args.inclWT = bool(args.inclWT)
    args.keep_ligand = bool(args.keep_ligand)
    try:
        os.makedirs(args.output_dir, exist_ok = True)
        logger.info(f"Directory {args.output_dir} created successfully")
    except OSError as error:
        logger.warning(f"Directory {args.output_dir} can not be created")
    return args


def rechain(pdb_file, output_dir):
    rechain_pdb = os.path.join(output_dir, f'{os.path.basename(pdb_file)[:-4]}_uniquechain.pdb')
    with open(pdb_file, 'r') as fp, open(rechain_pdb, 'w') as fp2:
        for line in fp:
            if line[0:4] == "ATOM":
                chain = line[21]
                new_line = f"{line[:21]}A{line[22:]}"
                fp2.write(new_line)
            else:
                fp2.write(line)
    return rechain_pdb


def create_homo_oligomer_input(pdb_file, prism_file, output_dir, additional=False, inclWT=True, keep_ligand=None):
    #PDB to rosetta numbering
    if keep_ligand == True:
        ligs = []
        with open(pdb_file, 'r') as fp:
            for line in fp:
                if line.startswith('HETATM'):
                    ligs.append(line[17:20].strip())
        ligs = list(set(ligs))
        renumb_pdb = pdb_renumb(pdb_file, output_dir=output_dir, keepchain='all', keep_ligand=ligs)
    else:
        renumb_pdb = pdb_renumb(pdb_file, output_dir=output_dir, keepchain='all')

    #get sequence
    pdb_p = PDBParser()
    structure = pdb_p.get_structure('test', renumb_pdb)
    model = structure[0]
    all_data = {}
    nchains = []
    for chain in model:
        nchain = str(chain).split('=')[1][0]
        nchains.append(nchain)
        sequence = []
        residues = []
        i = 0
        for res in chain:
            residue = str(res).split()
            if residue[1] in d3to1.keys():
                sequence.append(d3to1[residue[1]])
                residues.append(int(residue[3].split('=')[1]))
                i += 1
        all_data[nchain] = [sequence, residues]    
    
    meta_data, dataframe = read_from_prism(prism_file)
    prism_seq = meta_data['protein']['sequence']
    prism_var = {}
    for var in dataframe['variant'].tolist():
        if not var in ['WT']:
            resi = int(var[1:-1])
            if resi in prism_var.keys():
                prism_var[resi].append(var)
            else:
                prism_var[resi] = [var]
    if 'first_res_num' in meta_data['protein']:
        sys.exit('We have a first res num case! Needs to be handled - let Johanna know ;-)')

    for key in prism_var.keys():
        native = f"{prism_var[key][0][0]}{key}{prism_var[key][0][0]}"
        if not native in prism_var[key]:
            prism_var[key].append(native)

    alignments_all = []
    alignments_all_nresi = []
    for chain_n in nchains:
        seq_X = "".join(all_data[chain_n][0])

        alignmentsX = pairwise2.align.globalms(prism_seq, seq_X, 2, -1, -0.5, -0.1)[0]
        alignments_all.append(alignmentsX)

        alignmentsX_nresi = [[-1, '-']]*len(alignmentsX[1])
        i = 0
        for index, elem in enumerate(alignmentsX[1]):
            if elem !='-':
                alignmentsX_nresi[index] = [all_data[chain_n][1][i], all_data[chain_n][0][i]]
                i+=1
        alignments_all_nresi.append(alignmentsX_nresi)
    
    def check_align(alignments_all, index):
        for aligns in alignments_all[1:]:
            if (len(alignments_all[0][1]) >= index) and (len(aligns[1]) >= index):
                if (alignments_all[0][1][index-1] == aligns[1][index-1]):
                    pass
                else:
                    return False
            else:
                return False
        return True

    def get_muts(alignments_all_nresi, index, variant):
        muts = []
        for alignmentsX_nresi in alignments_all_nresi:
            mutX = f"{alignmentsX_nresi[index-1][1]}{alignmentsX_nresi[index-1][0]}{variant[-1]}"
            muts.append(mutX)
        return muts

    vari_list = []
    for index, res in enumerate(alignments_all[0][1]):
        if (alignments_all[0][1][index-1] != '-') & check_align(alignments_all, index):
            if index in prism_var.keys():
                for variant in prism_var[index]:
                    muts = get_muts(alignments_all_nresi, index, variant)
                    merged_variant = ":".join(muts)
                    vari_list.append(merged_variant)


    df = pd.DataFrame(data=vari_list, columns=['variant'])
    df['n_mut'] = 2
    df['aa_ref'] = df['variant'].apply(lambda x: [resi[0] for resi in x.split(':')])
    df['resi'] = df['variant'].apply(lambda x: [int(resi[1:-1]) for resi in x.split(':')])
    df['aa_var'] = df['variant'].apply(lambda x: [resi[-1] for resi in x.split(':')])
    df = df[['variant', 'aa_ref', 'resi', 'aa_var', 'n_mut']]


    if additional:
        residue_files_func(output_dir, df, inclWT=inclWT, drop_multi_mut=False)
        combined_mut_func(output_dir, df, inclWT=inclWT, drop_multi_mut=False)
    pipeline_combined_mut_func(output_dir, df, inclWT=inclWT, drop_multi_mut=False)
    
    renamed_chain_pdb = rechain(renumb_pdb, output_dir)

    return renamed_chain_pdb, os.path.join(output_dir, 'mutfile_all')


def main():
    """
    Main function called as default at the end.
    """
    # get user input arguments
    args = parse_args()
    
    create_homo_oligomer_input(args.pdb_file, args.prism, args.output_dir, additional=args.additional, inclWT=args.inclWT, keep_ligand=args.keep_ligand)
    
    logger.info('Conversion sucessful!')

if __name__ == '__main__':
    main()