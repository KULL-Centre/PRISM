"""relax_parse_results.py 
This function will parse the results of a rosetta pre_relax run. it will be callable from an 
sbatch script, that can wait for the relaxation to finish, and then select the best one. The 
point is to implement the 20x pre relaxation for the stability pipeline, that Amelie requested.

Author: Anders Frederiksen

Date of last major changes: 2020-04

"""

# Standard library imports
import logging as logger
import glob
import os
import io
import sys
from os.path import join
import subprocess

# 3rd party library imports
from Bio import pairwise2
from Bio import PDB
from Bio.Seq import Seq 
import numpy as np
import pandas as pd

# Local application imports
from helper import AttrDict, create_symlinks, create_copy, find_copy
from get_memory_stats import check_memory
from mp_helper import extract_from_opm, get_res_in_all, getnums, get_seq, get_res_in_mem


d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def parse_relax_results(folder, pdb_id='', sc_name='score_bn15_calibrated', logger_mode='info', keep=1, is_MP=False, do_checking=True):
    '''This function parses the scorefile from a rosetta
    pre-relaxation, and selects the lowest scoring one'''
    print(sc_name, keep, logger_mode, folder)
    try:
        job_id_pos= join(folder.relax_run, 'job_id_relax.txt')
        with open(job_id_pos, 'r') as job_id_file:
            relax_process_id=str(job_id_file.readlines()[-1])
            print(relax_process_id)
         
        #Get stats 
        shell_command = f'sacct -j {relax_process_id} --format="JobID,Start,End,CPUTime,NCPUS,Elapsed,ReqMem,MaxRSS,MaxVMSize,AveVMSize,JobName" > memory_usage_{relax_process_id}.log'
        subprocess.call(shell_command, cwd=folder.relax_run, shell=True)
                                                
        check_memory(relax_process_id,folder.relax_run) 
        
    except: 
        print('no memory file')
    
    
    relax_scores = {}
    relax_models = {}
    for n in range(20):
        if sc_name=='score_bn15_calibrated':
            path_to_scorefile = os.path.join(
                folder.relax_run, f'{n}-{sc_name}.sc')
        else:
            path_to_scorefile = os.path.join(
                folder.relax_run, f'{sc_name}.sc')
        with open(path_to_scorefile) as scorefile:
            scorelines = scorefile.readlines()

        # Put the relaxation scores in this dict
        for line in scorelines:
            score_fields = line.split()
            name = score_fields[-1].strip()
            if score_fields[0].strip() == 'SCORE:' and name != 'description':
                score = float(score_fields[1])
                relax_scores[score] = name
                relax_models[name] = score
    print(relax_scores, keep)
    most_relaxed_array = sorted(list(relax_scores.keys()), reverse=False)[:int(keep)]
    most_relaxed = relax_scores[most_relaxed_array[0]]
    print(most_relaxed_array)
    most_relaxed_models = [relax_scores[score] for score in most_relaxed_array]
    logger.info(f'most relaxed structure is {most_relaxed}.')

    # checking for consitency
    if do_checking:
        input_pdb = glob.glob(os.path.join(folder.relax_input, '*.pdb'))[0]
        output_pdb = os.path.join(folder.relax_run, f'{most_relaxed}.pdb')
        chain = 'A'
        with open(output_pdb, 'r') as fp:
            for line in fp:
                if line.startswith('ATOM'):
                    chain = line[21]
                    break
        consistency = check_consitency(input_pdb, output_pdb, chain=chain, is_MP=is_MP, pdb_id=pdb_id)
        if consistency == False:
            print('input and relaxed structures are not consitent or too different')
            sys.stderr()
            sys.exit()
            return

    for key in relax_models:
        if not key in most_relaxed_models:
            path_to_tense = os.path.join(folder.relax_run, f'{key}.pdb')
            print('deleting', path_to_tense)
            os.remove(path_to_tense)
    if int(keep) > 1:
        for ind, ddg_subfolder in enumerate(folder.ddG_input):
            relax_output_strucfile = find_copy(
                folder.relax_run, f'{relax_scores[most_relaxed_array[ind]]}.pdb', folder.relax_output, f'output_{ind}.pdb')
            create_copy(
                os.path.join(folder.relax_output, f'output_{ind}.pdb'), ddg_subfolder, name=f'input.pdb')
    else:
        relax_output_strucfile = find_copy(
            folder.relax_run, f'{most_relaxed}.pdb', folder.relax_output, 'output.pdb')
        create_copy(
            os.path.join(folder.relax_output, 'output.pdb'), folder.ddG_input, name='input.pdb')

    return os.path.join(folder.relax_output, f'{most_relaxed}.pdb')


def check_seq_completeness(input_pdb, output_pdb):
    #get sequence info about TM region and where best to align

    ref_align_atoms = get_res_in_all(input_pdb, isfile=True)

    seq1, seq1num, maxi1 = get_seq(output_pdb, isfile=True)
    seq2, seq2num, maxi2 = get_seq(input_pdb, isfile=True)
    seq1 = Seq(seq1)
    seq2 = Seq(seq2)
    alignments = pairwise2.align.globalxx(seq1, seq2)
    maxx = maxi2 if maxi1 < maxi2 else maxi1
    for align in alignments:
        if maxx == align[-1]:
            break
    align
    align[0]
    seq1 = [align[0][i:i+1] for i in range(0, len(align[0]), 1)]
    seqnum1 = getnums(seq1, seq1num)
    seq2 = [align[1][i:i+1] for i in range(0, len(align[1]), 1)]
    seqnum2 = getnums(seq2, seq2num)

    df = pd.DataFrame(np.array([seq1, seqnum1, seq2, seqnum2]).T, columns=['final', 'final_num', 'original', 'original_num'])
    df2 = df.dropna(axis=0, how='all', subset=['final_num', 'original_num'])
    if '-' in df2['final'].tolist()+df2['original'].tolist():
        print('Sequences are different')
        return False
    else:
        print('Sequences are the same')
        return True

def check_struc_alignment(ref_pdb, target, target_chain='A', ref_model_id=0, target_model_id=0):
    # Adapted from https://gist.github.com/andersx/6354971
    # Copyright (c) 2010-2016 Anders S. Christensen

    # Parse reference (from file) and target structure (from file)
    parser = PDB.PDBParser(QUIET=True)
    bio_ref_struc_raw = parser.get_structure("reference", ref_pdb)
    bio_target_struc_raw = parser.get_structure("target", target)

    ref_align_atoms = get_res_in_all(ref_pdb, isfile=True)

    seq1, seq1num, maxi1 = get_seq(target, isfile=True)
    seq2, seq2num, maxi2 = get_seq(ref_pdb, isfile=True)
    seq1 = Seq(seq1)
    seq2 = Seq(seq2)
    alignments = pairwise2.align.globalxx(seq1, seq2)
    maxx = maxi2 if maxi1 < maxi2 else maxi1
    for align in alignments:
        if maxx == align[-1]:
            break

    seq1 = [align[0][i:i+1] for i in range(0, len(align[0]), 1)]
    seqnum1 = getnums(seq1, seq1num)
    seq2 = [align[1][i:i+1] for i in range(0, len(align[1]), 1)]
    seqnum2 = getnums(seq2, seq2num)

    df = pd.DataFrame(np.array([seq1, seqnum1, seq2, seqnum2]).T, columns=['infile', 'infile_num', 'opm', 'opm_num'])
    df = df.replace(to_replace='None', value=np.nan).dropna(subset=['infile_num', 'opm_num'], how='all').reset_index(drop=True)
    df = df.replace(to_replace='None', value=np.nan).dropna(subset=['infile_num', 'opm_num'], how='any').reset_index(drop=True)

    ref_align_atoms = df['opm_num'].astype(int).unique().tolist()
    target_align_atoms = df['infile_num'].astype(int).unique().tolist()

    # Select the model number - normally always the first
    bio_ref_struc = bio_ref_struc_raw[ref_model_id]
    bio_target_struc = bio_target_struc_raw[target_model_id]

    # List of residues to align
    align_ref_atoms = []
    for ind, chain in enumerate(bio_ref_struc_raw.get_chains()):
        if chain.id == target_chain:
            for res in chain.get_residues():  # bio_ref_struc.get_residues():
                if ref_align_atoms == [] or res.get_id()[1] in ref_align_atoms:
                    for atom in res:
                        if atom.get_name() == 'CA':
                            align_ref_atoms.append(atom)
    align_target_atoms = []
    for ind, chain in enumerate(bio_target_struc.get_chains()):
        if chain.id == target_chain:
            for res in chain.get_residues():  # bio_target_struc.get_residues():
                if target_align_atoms == [] or res.get_id()[1] in target_align_atoms:
                    for atom in res:
                        if atom.get_name() == 'CA':
                            align_target_atoms.append(atom)

    # Superposer
    super_imposer = PDB.Superimposer()
    super_imposer.set_atoms(align_ref_atoms, align_target_atoms)
    print(f"RMSD of not yet superimposed structures of all residues: {super_imposer.rms}")
    #super_imposer.apply(bio_target_struc)
    #print(f"RMSD of superimposed structures of all residues: {super_imposer.rms}")

    if super_imposer.rms > 7.5:
        print(f'Structural alignment too bad ({super_imposer.rms} > 75.) - align structure manually')
        return False
    else:
        return True

def check_consitency(input_pdb, output_pdb, chain='A', is_MP=True, pdb_id='-'):
    complete = check_seq_completeness(input_pdb, output_pdb)
    if is_MP:
        aligned = check_struc_alignment(input_pdb, output_pdb, target_chain=chain)
    else:
        aligned = True
    if (complete == True) and (aligned == True):
        return True
    else:
        return False

if __name__ == '__main__':
    folder = AttrDict()
    print(sys.argv)
    if sys.argv[1]=='True':
        do_checking = True
    else:
        do_checking = False
    if sys.argv[2]=='True':
        is_MP = True
    else:
        is_MP = False
    if sys.argv[3]=='-':
        pdb_id = ''
    else:
        pdb_id = sys.argv[3]

    if len(sys.argv) <= 9:
        folder.update({'relax_input': sys.argv[4], 'relax_run': sys.argv[5], 'relax_output': sys.argv[6], 
            'ddG_input': sys.argv[7]})
    else:
        folder.update({'relax_input': sys.argv[4], 'relax_run': sys.argv[5], 'relax_output': sys.argv[6], 
            'ddG_input': [sys.argv[x] for x in range(7, len(sys.argv)-2)]})
    print(sys.argv)
    print(folder)
    if len(sys.argv) == 9:
        parse_relax_results(folder, sc_name=sys.argv[8], is_MP=is_MP, pdb_id=pdb_id, do_checking=do_checking)
    elif len(sys.argv) > 9:
        parse_relax_results(folder, sc_name=sys.argv[-2], keep=int(sys.argv[-1]), is_MP=is_MP, pdb_id=pdb_id, do_checking=do_checking)
    else:
        parse_relax_results(folder)
