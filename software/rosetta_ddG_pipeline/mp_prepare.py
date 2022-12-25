"""mp_prepare.py includes several preprocessing and script generating functions to prepare membrane proteins.

Author: Johanna K.S. Tiemann

Date of last major changes: 2022-07-26

"""

# Standard library imports
import io
import logging as logger
import os
import random
import re
import shutil
import subprocess
import sys
import time
import glob
import json

# Third party imports
from Bio import PDB
from Bio.Seq import Seq 
from Bio import pairwise2
import biolib
import numpy as np
import pandas as pd


# Local application imports
import pdb_to_fasta_seq
import rosetta_paths
from mp_helper import extract_from_opm, get_res_in_all, get_res_in_mem, get_seq, getnums 

d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def mp_superpose_opm(reference_chain, target, filename, target_chain='A',
                     ref_model_id=0, target_model_id=0, write_opm=True, inTM=True):
    # Adapted from https://gist.github.com/andersx/6354971
    # Copyright (c) 2010-2016 Anders S. Christensen

    # Get reference structure e.g. from OPM
    def get_ref_struc(keyword):
        try:
            ref_struc = extract_from_opm(keyword)
            if write_opm:
                ref_file_name = os.path.join(os.path.dirname(filename), 'ref_opm.pdb')
                with open(ref_file_name, 'w') as fp:
                    fp.write(ref_struc)
                print(f"write_opm set to true - OPM structure saved")
        except:
            print("no OPM - id found - add a workaround, e.g. PDBTM")
            return
        else:
            print("Obtain structure from OPM: successful")
        return ref_struc
    reference = reference_chain.split('_')[0]
    ref_chain = reference_chain.split('_')[1]
    ref_struc = get_ref_struc(reference)
    ref_file_name = os.path.join(os.path.dirname(filename), 'ref_opm.pdb')

    sucess = superpose_MMLigner(target, target_chain[0], 
        ref_file_name, ref_chain, filename, 
        exec_dir = os.path.dirname(filename))

    if sucess == False:
        # Parse reference (from string) and target structure (from file)
        parser = PDB.PDBParser(QUIET=True)
        bio_ref_struc_raw = parser.get_structure(
            "reference", io.StringIO(ref_struc))
        bio_target_struc_raw = parser.get_structure("target", target)
        #get sequence info about TM region and where best to align

        if inTM:
            ref_align_atoms = get_res_in_mem(ref_struc, isfile=False)
        else:
            ref_align_atoms = get_res_in_all(ref_struc, isfile=False)

        seq1, seq1num, maxi1 = get_seq(target, isfile=True)
        seq2, seq2num, maxi2 = get_seq(ref_struc, isfile=False)

        superpose_struc(bio_ref_struc_raw, bio_target_struc_raw, ref_align_atoms, 
                         seq1, seq1num, maxi1, seq2, seq2num, maxi2, filename, target_chain=target_chain[0], ref_chain=ref_chain,
                         ref_model_id=ref_model_id, target_model_id=target_model_id)

def superpose_struc(bio_ref_struc_raw, bio_target_struc_raw, ref_align_atoms, 
                     seq1, seq1num, maxi1, seq2, seq2num, maxi2, filename, target_chain='A',ref_chain='A',
                     ref_model_id=0, target_model_id=0):
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
    #df = df.dropna(subset=["infile_num","opm_num"], how='any').reset_index(drop=True)
    df['infile_num'] = df['infile_num'].astype('int')
    df['opm_num'] = df['opm_num'].astype('int')

    target_align_atoms = []
    for res in ref_align_atoms:
        target_align_atoms.append(df.loc[(df['opm_num']==int(res)), 'infile_num'].tolist()[0])


    # Select the model number - normally always the first
    bio_ref_struc = bio_ref_struc_raw[ref_model_id]
    bio_target_struc = bio_target_struc_raw[target_model_id]

    # List of residues to align
    align_ref_atoms = []
    for ind, chain in enumerate(bio_ref_struc_raw.get_chains()):
        if chain.id == ref_chain[0]:
            for res in chain.get_residues():  # bio_ref_struc.get_residues():
                if ref_align_atoms == [] or res.get_id()[1] in ref_align_atoms:
                    for atom in res:
                        if atom.get_name() == 'CA':
                            align_ref_atoms.append(atom)
    align_target_atoms = []
    for ind, chain in enumerate(bio_target_struc.get_chains()):
        if chain.id == target_chain[0]:
            for res in chain.get_residues():  # bio_target_struc.get_residues():
                if target_align_atoms == [] or res.get_id()[1] in target_align_atoms:
                    for atom in res:
                        if atom.get_name() == 'CA':
                            align_target_atoms.append(atom)
    # Superposer
    super_imposer = PDB.Superimposer()
    super_imposer.set_atoms(align_ref_atoms, align_target_atoms)
    super_imposer.apply(bio_target_struc)

    print(f"RMSD of superimposed structures: {super_imposer.rms}")
    if super_imposer.rms > 10:
        print(f'Structural alignment too bad ({super_imposer.rms} > 10) - align structure manually')
        sys.exit()

    bioio = PDB.PDBIO()
    bioio.set_structure(bio_target_struc)
    bioio.save(filename)
    print(f"Aligned structure saved")


def superpose_MMLigner(target_pdb, target_chain, reference_pdb, reference_chain, output_pdb, exec_dir = ''):
    
    # make clean version of reference file"
    reference_pdb_cleaned = os.path.join(os.path.dirname(reference_pdb), 'ref_cleaned.pdb')
    with open(reference_pdb, 'r') as fp, open(reference_pdb_cleaned, 'w') as fp2:
        prev = ''
        for line in fp:
            if line.startswith('# All scores below are weighted scores, not raw scores.'):
                break
            elif line.startswith('CONECT'):
                pass
            elif line[17:20].strip() in ['MEM', 'EMB']:
                pass
            elif (prev == 'TER') and (line[:3] == prev):
                pass
            else:
                prev = line[:3]
                fp2.write(line)

    exect_mmlinger = (f"{rosetta_paths.MMLigner_exec} "
        f"{reference_pdb_cleaned}:{reference_chain[0]} "
        f"{target_pdb}:{target_chain[0]} "
        "--superpose"
        "")
    print(f"Superpose pdbs using MMLigner: {exect_mmlinger} - here: {exec_dir}")
    mmlinger_call = subprocess.Popen(
                exect_mmlinger, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, cwd=exec_dir)
    mmlinger_process_id_info = mmlinger_call.communicate()

    if len(glob.glob(os.path.join(exec_dir, '*superposed*.pdb')))==0:
        return False
    else:
        shutil.copy(glob.glob(os.path.join(exec_dir, '*superposed*.pdb'))[0], output_pdb)
        print(glob.glob(os.path.join(exec_dir, '*superposed*.pdb'))[0], output_pdb)
        print(f"Aligned structure saved")
        return True


def mp_TMalign_opm(reference_chain, target, filename, target_chain='A',
                   ref_model_id=0, target_model_id=0,
                   ref_align_atoms=[], target_align_atoms=[], write_opm=False):
    # Get reference structure e.g. from OPM
    def get_ref_struc(keyword):
        try:
            ref_struc = extract_from_opm(keyword)
        except:
            print("no OPM - id found - add a workaround, e.g. PDBTM")
            return
        else:
            print("Obtain structure from OPM: successful")
        return ref_struc
    reference = reference_chain.split('_')[0]
    ref_chain = reference_chain.split('_')[1]
    ref_struc = get_ref_struc(reference)

    ref_struc_dest = os.path.join(os.path.dirname(filename), 'ref_opm.pdb')
    with open(ref_struc_dest, 'w') as fp:
        fp.write(ref_struc)
    print(f"OPM structure saved")


    path_to_TMalign = rosetta_paths.path_to_TMalign
    shell_call = f'{path_to_TMalign} {target} {ref_struc_dest} -o {filename.split(".pdb")[0]}'
    subprocess.call(shell_call, shell=True)

    print(f"Aligned structure saved")


def mp_superpose_span(pdbinput, outdir_path, span_file_list, filename, SLURM=False, chain='A'):
    """
    Calculates the residues which are accessible by lipids. 
    Rosetta scripts used: mp_lipid_acc
    Cite:   Koehler Leman, Lyskov, & Bonneau (2017) BMC Bioinformatics (https://doi.org/10.1186/s12859-017-1541-z)
    Documentation: https://www.rosettacommons.org/docs/latest/application_documentation/membrane_proteins/RosettaMP-App-MPSpanFromPDB
    """
    spanfile = span_file_list[0]
    Rosetta_lipacc_exec = os.path.join(
        rosetta_paths.path_to_rosetta, f'bin/mp_transform.{rosetta_paths.Rosetta_extension}')
    lipacc_command = (f'{Rosetta_lipacc_exec} '
                    f'-in:file:s {pdbinput} '
                    f'-database {rosetta_paths.Rosetta_database_path} '
                    f'-mp:setup:spanfiles {spanfile} '
                    '-restore_talaris_behavior -mp:scoring:hbond true -mp:restore_lazaridis_imm_behavior 1'
                    '')
    print(f"MP_transform call function: {lipacc_command}")

    if SLURM:
        logger.warn("need to write the slurm script!")
        sys.exit()
        lipacc_call = subprocess.Popen('sbatch rosetta_lipacc.sbatch', stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT, shell=True, cwd=outdir_path)
        lipacc_process_id_info = lipacc_call.communicate()
        lipacc_process_id = str(lipacc_process_id_info[0]).split()[3][0:-3]
        print(lipacc_process_id_info)
    else:
        lipacc_call = subprocess.Popen(
            lipacc_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, cwd=outdir_path)
        lipacc_process_id_info = lipacc_call.communicate()

        if len(glob.glob(os.path.join(outdir_path, '*.pdb')))==0:
            print('Missing heavy atom. run fixbb')
            fixbb_output_dir = os.path.join(outdir_path, 'fixbb')
            os.makedirs(fixbb_output_dir, exist_ok=True)

            Rosetta_fixbb_exec = os.path.join(
                rosetta_paths.path_to_rosetta, f'bin/fixbb.{rosetta_paths.Rosetta_extension}')
            pre_align_cmd = (f'{Rosetta_fixbb_exec} '
                f'-in:file:s {pdbinput} '
                f'-in:file:fullatom -nstruct 1 '
                '-minimize_sidechains false -min_pack true '
                '-off_rotamer_pack true -ex1 -ex2 ' 
                '-packing:ndruns 1 -packing:repack_only True '
                '')
            print(pre_align_cmd)
            fixbb_call = subprocess.Popen(
                pre_align_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, cwd=fixbb_output_dir)
            fixbb_process_id_info = fixbb_call.communicate()

            newpdbinput = glob.glob(os.path.join(fixbb_output_dir, '*.pdb'))[0]
            lipacc_command = (f'{Rosetta_lipacc_exec} '
                    f'-in:file:s {newpdbinput} '
                    f'-database {rosetta_paths.Rosetta_database_path} '
                    f'-mp:setup:spanfiles {spanfile} '
                    '-restore_talaris_behavior -mp:scoring:hbond true -mp:restore_lazaridis_imm_behavior 1'
                    '')

            print(lipacc_command)
            lipacc_call = subprocess.Popen(
                lipacc_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, cwd=outdir_path)
            lipacc_process_id_info = lipacc_call.communicate()
    
    ref_struc = glob.glob(os.path.join(outdir_path, '*.pdb'))[0]

    chains = []
    with open(ref_struc, 'r') as fp:
        for line in fp:
            if line.startswith('ATOM'):
                chains.append(line[21])
    chains = list(set(chains))
    for chain_re in chains:
        sucess = superpose_MMLigner(pdbinput, chain[0], 
            ref_struc, chain_re, filename, 
            exec_dir = os.path.dirname(filename))
        if sucess == True:
            break

    if sucess == False:
        for chain_re in chains:
            try:
                # Parse reference (from string) and target structure (from file)
                parser = PDB.PDBParser(QUIET=True)
                bio_ref_struc_raw = parser.get_structure("reference", ref_struc)
                bio_target_struc_raw = parser.get_structure("target", pdbinput)
                #get sequence info about TM region and where best to align

                ref_align_atoms = get_res_in_all(ref_struc, isfile=True)

                seq1, seq1num, maxi1 = get_seq(pdbinput, isfile=True)
                seq2, seq2num, maxi2 = get_seq(ref_struc, isfile=True)

                superpose_struc(bio_ref_struc_raw, bio_target_struc_raw, ref_align_atoms, 
                                seq1, seq1num, maxi1, seq2, seq2num, maxi2, filename, target_chain=chain[0], ref_chain=chain_re,
                                ref_model_id=0, target_model_id=0)
            except:
                pass


def check_for_repeats(seq):
    REPEATER = re.compile(r"(.+?)\1+$")

    def repeated(s):
        match = REPEATER.match(s)
        return match.group(1) if match else None

    sub = repeated(seq)
    if sub:
        ari = seq.split(sub)
        for ar in ari:
            if ar != '':
                return False, 1, len(seq)
        return True, len(ari), len(sub)
    else:
        return False, 1, len(seq)


def mp_span_from_deepTMHMM(pdbinput, outdir_path, signal_TM=False):

    spanfiles = []

    # get input sequence
    input_secs = {}
    with open(pdbinput, 'r') as fp:
        for line in fp:
            if line.startswith('ATOM'):
                resseq = int(line[22:26])
                resname = line[17:20].strip()
                chain = line[21]
                if line[12:16].strip() == 'CA':
                    if line[16] in [' ', 'A']:
                        resname = d3to1[resname]
                        if chain in input_secs.keys():
                            inputseq = input_secs[chain]
                        else:
                            inputseq = []
                        inputseq.append(resname)
                        input_secs[chain] = inputseq
    for chain in input_secs.keys():
        tmp_output_path = os.path.join(outdir_path, chain)
        os.makedirs(tmp_output_path, exist_ok=True)

        inputseq = input_secs[chain]
        inputseq = "".join(inputseq)

        # check for seqs:
        does_repeat, num_repeats, len_repeats = check_for_repeats(inputseq)
        if does_repeat == True:
            inputseq_back = inputseq
            inputseq = inputseq_back[0:len_repeats]

        # make fasta_file
        fasta_file = os.path.join(tmp_output_path, 'query.fasta')
        with open(fasta_file, 'w') as fp:
            fp.write('> input\n')
            fp.write(inputseq)

        def calc_deepTMHMM(fasta_file, tmp_output_path):
            # calculate TM regions
            deeptmhmm = biolib.load('DTU/DeepTMHMM')
            deeptmhmm_res = deeptmhmm.cli(args=f'--fasta {fasta_file}')
            deeptmhmm_res.save_files(tmp_output_path)

            result_file = os.path.join(tmp_output_path, 'TMRs.gff3')
            return result_file

        run_count = 0
        while run_count < 5 :
            try:
                result_file = calc_deepTMHMM(fasta_file, tmp_output_path)
            except Exception as e:
                print("The error raised is: ", e)
                print("Maybe consider updating the library: pip3 install -U pybiolib")
                result_file = ''
            if os.path.isfile(result_file):
                run_count = 9
            else:
                run_count += 1
                sleeptime = random.randint(1*60, 5*60)
                time.sleep(sleeptime) # wait between 1 to 5 minutes before trying to reach the server again

        # get total length
        total_length = len(inputseq)
        with open(result_file, 'r') as fp:
            for line in fp:
                if line.startswith('#'):
                    if 'Length' in line:
                        line = line.split()
                        total_length = line[3]
                        break

        # get TM regions
        df = pd.read_csv(result_file, skiprows=3, delimiter='\t', header=None)
        df.rename(columns={0: 'ID', 1: 'location', 2:'start', 3: 'end'}, inplace=True)
        df.to_csv(os.path.join(outdir_path, f"{os.path.splitext(os.path.basename(pdbinput))[0]}_{chain}_deepTMHMM.csv"))
        print(df['location'].unique())
        if signal_TM:
            non_tm = ['outside', 'inside', 'periplasm']
            tm = ['Beta sheet', 'TMhelix', 'signal']
        else:
            non_tm = ['outside', 'inside', 'periplasm', 'signal']
            tm = ['Beta sheet', 'TMhelix']
        TM_df = df.loc[~df['location'].isin(non_tm)]
        num_span = len(TM_df['start'])

        # get info about structure
        order = 'antiparallel' # parallel possible - not sure how to calculate - only possible for n1c (horizontal)

        # get pdb resnumber info
        res_dic = []
        with open(pdbinput, 'r') as fp:
            for line in fp:
                if line.startswith('ATOM'):
                    if line[21] == chain:
                        if line[12:16].strip() == 'CA':
                            if line[16] in [' ', 'A']:
                                res_dic.append(int(line[22:26]))
        res_dic.sort()

        # write span file
        if len(TM_df)>0:
            span_file = os.path.join(outdir_path, f"{os.path.splitext(os.path.basename(pdbinput))[0]}_{chain}.span")
            with open(span_file, 'w') as fp:
                fp.write('manual-generated spanfile from DeepTMHMM\n')
                if does_repeat == True:
                    fp.write(f'{num_span*num_repeats} {int(total_length)*(num_repeats-1)}\n')
                else:
                    fp.write(f'{num_span} {int(total_length)}\n')
                fp.write(f'{order}\n')
                fp.write('n2c\n')
                for index, row in TM_df.iterrows():
                    fp.write(f"\t\t{row['start']+res_dic[0]-1}\t{row['end']+res_dic[0]-1}\n")
                if does_repeat == True:
                    for reps in range(1, num_repeats-1):
                        for index, row in TM_df.iterrows():
                            fp.write(f"\t\t{row['start']+res_dic[0]-1+(len_repeats*reps)}\t{row['end']+res_dic[0]-1+(len_repeats*reps)}\n")

            spanfiles.append(span_file)
    print("Span process done")
    return spanfiles


def mp_span_from_pdb_octopus(pdbinput, outdir_path, SLURM=False):

    Rosetta_span_exec = os.path.join(
        rosetta_paths.path_to_rosetta, f'bin/spanfile_from_pdb.{rosetta_paths.Rosetta_extension}')
    span_command = f'{Rosetta_span_exec} -in:file:s {pdbinput}'
    print(f"Span call function: {span_command}")

    if SLURM:
        logger.warn("need to write the slurm script!")
        span_call = subprocess.Popen('sbatch rosetta_span.sbatch', stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT, shell=True, cwd=outdir_path)
        span_process_id_info = span_call.communicate()
        span_process_id = str(span_process_id_info[0]).split()[3][0:-3]
        print(span_process_id_info)
    else:
        span_call = subprocess.Popen(
            span_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, cwd=outdir_path)
        span_process_id_info = span_call.communicate()

    with open(os.path.join(outdir_path, 'span.log'), 'w') as fp:
        fp.writelines(span_process_id_info[0].decode())

    if not any(fname.endswith('.span') for fname in os.listdir(outdir_path)):
        print("Span file not created. Check log file.")
    else:
        spanfiles = []
        for fname in os.listdir(outdir_path):
            if (fname.endswith('.span')):
                new_filename = f"{os.path.splitext(os.path.basename(pdbinput))[0]}.span"
                os.rename(os.path.join(outdir_path, fname),
                          os.path.join(outdir_path, new_filename))
                spanfiles.append(os.path.join(outdir_path, new_filename))
        print("Span process done")
        return spanfiles


def mp_lipid_acc_resi(pdbinput, outdir_path, folder_spanfile, thickness=15, SLURM=False):
    """
    Calculates the residues which are accessible by lipids. 
    Rosetta scripts used: mp_lipid_acc
    Cite:   Koehler Leman, Lyskov, & Bonneau (2017) BMC Bioinformatics (https://doi.org/10.1186/s12859-017-1541-z)
    Documentation: https://www.rosettacommons.org/docs/latest/application_documentation/membrane_proteins/RosettaMP-App-MPSpanFromPDB
    """
    for root, dirs, files in os.walk(folder_spanfile):
        for file in files:
            if file.endswith('.span'):
                spanfile = os.path.join(root, file)
    Rosetta_lipacc_exec = os.path.join(
        rosetta_paths.path_to_rosetta, f'bin/mp_lipid_acc.{rosetta_paths.Rosetta_extension}')
    lipacc_command = (f'{Rosetta_lipacc_exec} '
                    f'-in:file:s {pdbinput} '
                    f'-database {rosetta_paths.Rosetta_database_path} '
                    f'-mp:setup:spanfiles {spanfile} '
                    '')
    print(f"Lipid_acc_resi call function: {lipacc_command}")

    if SLURM:
        logger.warn("need to write the slurm script!")
        sys.exit()
        lipacc_call = subprocess.Popen('sbatch rosetta_lipacc.sbatch', stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT, shell=True, cwd=outdir_path)
        lipacc_process_id_info = lipacc_call.communicate()
        lipacc_process_id = str(lipacc_process_id_info[0]).split()[3][0:-3]
        print(lipacc_process_id_info)
    else:
        lipacc_call = subprocess.Popen(
            lipacc_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, cwd=outdir_path)
        lipacc_process_id_info = lipacc_call.communicate()
    
    lipacc_dic = {}
    out_pdb_file = glob.glob(os.path.join(outdir_path, '*.pdb'))[0]
    output_json_file = os.path.join(outdir_path, "mp_lipid_acc_dic.json")
    with open (out_pdb_file, 'r') as fp, open(output_json_file, 'w') as fp2:
        for line in fp:
            if line.startswith('ATOM'):
                resid = int(line[22:26])
                if line[60:66]==' 50.00':
                    lipacc_dic[resid] = 'true'
                else:
                    lipacc_dic[resid] = 'false'
        json.dump(lipacc_dic, fp2, indent=4) 
    return lipacc_dic, output_json_file


def mp_span_from_pdb_dssp(pdbinput, outdir_path, thickness=15, SLURM=False):
    """
    Calculates the membrane spanning Rosetta input file for membrane proteins using the pdb orientation and dssp. 
    Rosetta scripts used: mp_span_from_pdb
    Cite:   Alford et al. Gray (2015) PLoS Comput. Biol (https://doi.org/10.1371/journal.pcbi.1004398)
            Koehler Leman et al. Gray (2017) Bioinformatics (https://doi.org/10.1093/bioinformatics/btw716)
    Documentation: https://www.rosettacommons.org/docs/latest/application_documentation/membrane_proteins/RosettaMP-App-MPSpanFromPDB
    """

    Rosetta_span_exec = os.path.join(
        rosetta_paths.path_to_rosetta, f'bin/mp_span_from_pdb.{rosetta_paths.Rosetta_extension}')
    span_command = (f'{Rosetta_span_exec} '
                    f'-in:file:s {pdbinput} '
                    f'-mp::thickness {thickness} '
                    f'-database {rosetta_paths.Rosetta_database_path} '
                    '-ignore_unrecognized_res true '
                    '-ignore_zero_occupancy false '
                    '')
    print(f"Span call function: {span_command}")

    if SLURM:
        logger.warn("need to write the slurm script!")
        sys.exit()
        span_call = subprocess.Popen('sbatch rosetta_span.sbatch', stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT, shell=True, cwd=outdir_path)
        span_process_id_info = span_call.communicate()
        span_process_id = str(span_process_id_info[0]).split()[3][0:-3]
        print(span_process_id_info)
    else:
        span_call = subprocess.Popen(
            span_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, cwd=outdir_path)
        span_process_id_info = span_call.communicate()

    with open(os.path.join(outdir_path, 'span.log'), 'w') as fp:
        fp.writelines(span_process_id_info[0].decode())
    pdbinput_dir = os.path.dirname(pdbinput)
    if not any(fname.endswith('.span') for fname in os.listdir(pdbinput_dir)):
        print("Span file not created. Check log file.")
    else:
        spanfiles = []
        for fname in os.listdir(pdbinput_dir):
            if (fname.endswith('.span')):
                new_filename = f"{os.path.splitext(os.path.basename(pdbinput))[0]}.span"
                os.rename(os.path.join(pdbinput_dir, fname),
                          os.path.join(outdir_path, new_filename))
                spanfiles.append(os.path.join(outdir_path, new_filename))
        print("Span process done")
        return spanfiles


def rosetta_relax_mp(folder, SLURM=False, num_struc=20, sys_name='mp', partition='sbinlab', repeats=2, lipid_type='DLPC', 
  mp_thickness=15, mp_switch_off=False, score_function='franklin2019'):
    Rosetta_script_exec = os.path.join(
        rosetta_paths.path_to_rosetta, f'bin/rosetta_scripts.{rosetta_paths.Rosetta_extension}')
    for root, dirs, files in os.walk(folder.relax_input):
        for file in files:
            if file.endswith('.span'):
                spanfile = os.path.join(root, file)
    if score_function=='franklin2019':
      energy_fawtb=1.5
    else:
      energy_fawtb=0
    if mp_switch_off:
        relax_command = (f'{Rosetta_script_exec} '
                      # Use the membrane relax protocol Rosetta script
                      f'-parser:protocol {os.path.join(folder.relax_input, "relax.xml")} '
                      # Repeatition of FastRelax
                      f'-parser:script_vars repeats={repeats} energy_func={score_function} '
                      # Input PDB Structure: PDB file for protein structure
                      f' -in:file:s {os.path.join(folder.relax_input, "input.pdb")} '
                      # Spanfile describing trans-membrane spans of the
                      # starting structure
#                      f'-mp:setup:spanfiles {spanfile} '
#                      '-mp:scoring:hbond '  # Turn on membrane depth-dependent hydrogen bonding weight
#                      f'-mp:lipids:composition {lipid_type} '
#                      f'-mp::thickness {mp_thickness} '
                      # Use the FastRelax mode of Rosetta Relax (uses 5-8
                      # repeat cycles)
      #                '-relax:fast '
#                      '-relax:jump_move true '  # Allow the MEM and other jumps to move during refinement
      #                '-relax:bb_move false ' # Set all backbone torsion angles to unmovable during minimization. --> not working
                      # Number of structures to generate
                      f'-nstruct 1 '
                      # Wait to pack until the membrane mode is turned on
                      '-packing:pack_missing_sidechains 0 '
                      '-out:pdb '  # Output all PDB structures of refined models
                      # Specify destination for score file
                      f'-out:file:scorefile {os.path.join(folder.relax_run, "relax_scores.sc")} '
                      '-out:prefix $SLURM_ARRAY_TASK_ID- '
                      f'-database {rosetta_paths.Rosetta_database_path} '
                      '-ignore_unrecognized_res true '
                      f'-score:weights {score_function} '
                      #                        '-ignore_zero_occupancy false '
                      '')
    else:
        relax_command = (f'{Rosetta_script_exec} '
                      # Use the membrane relax protocol Rosetta script
                      f'-parser:protocol {os.path.join(folder.relax_input, "relax.xml")} '
                      # Repeatition of FastRelax
                      f'-parser:script_vars repeats={repeats} energy_func={score_function} energy_fawtb={energy_fawtb} '
                      # Input PDB Structure: PDB file for protein structure
                      f'-in:file:s {os.path.join(folder.relax_input, "input.pdb")} '
                      # Spanfile describing trans-membrane spans of the
                      # starting structure
                      f'-mp:setup:spanfiles {spanfile} '
                      '-mp:scoring:hbond '  # Turn on membrane depth-dependent hydrogen bonding weight
                      f'-mp:lipids:composition {lipid_type} '
                      f'-mp::thickness {mp_thickness} '
                      # Use the FastRelax mode of Rosetta Relax (uses 5-8
                      # repeat cycles)
      #                '-relax:fast '
                      '-relax:jump_move true '  # Allow the MEM and other jumps to move during refinement
      #                '-relax:bb_move false ' # Set all backbone torsion angles to unmovable during minimization. --> not working
                      # Number of structures to generate
                      f'-nstruct 1 '
                      # Wait to pack until the membrane mode is turned on
                      '-packing:pack_missing_sidechains 0 '
                      '-out:pdb '  # Output all PDB structures of refined models
                      # Specify destination for score file
                      f'-out:file:scorefile {os.path.join(folder.relax_run, "relax_scores.sc")} '
                      '-out:prefix $SLURM_ARRAY_TASK_ID- '
                      f'-database {rosetta_paths.Rosetta_database_path} '
                      '-ignore_unrecognized_res true '
                      f'-score:weights {score_function} '
                      '-fa_max_dis 9 '
                      '-ex1 '
                      '-ex2 '
                      '-flip_HNQ '
                      '-missing_density_to_jump '
                      '-relax:coord_constrain_sidechains '
                      '-relax:constrain_relax_to_start_coords '
                      #                        '-ignore_zero_occupancy false '
                      '')
    print(f'MP relax call function: {relax_command}')

    if SLURM:
        path_to_sbatch = os.path.join(
            folder.relax_input, 'rosetta_relax.sbatch')
        with open(path_to_sbatch, 'w') as fp:
            fp.write(f'''#!/bin/sh
#SBATCH --job-name=relax_{sys_name}
#SBATCH --time=24:00:00
#SBATCH --mem 5000
#SBATCH --array=0-{num_struc-1}
#SBATCH --partition={partition}

# launching rosetta relax 
''')
            fp.write(f'{relax_command} ')
        print(f'Location of relax sbatch file: {path_to_sbatch}')
        return(path_to_sbatch)
    else:
        relax_call = subprocess.Popen(relax_command,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.STDOUT,
                                      shell=True,
                                      cwd=folder.relax_run)
        relax_process_id_info = relax_call.communicate()
        with open(os.path.join(folder.relax_run, 'relax.log'), 'w') as fp:
            fp.writelines(relax_process_id_info[0].decode())
