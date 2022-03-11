"""mp_prepare.py includes several preprocessing and script generating functions to prepare membrane proteins.

Author: Johanna K.S. Tiemann

Date of last major changes: 2020-04-16

"""

# Standard library imports
import io
import logging as logger
import os
import subprocess
import sys
import glob
import json

# Third party imports
from Bio import PDB
from Bio.Seq import Seq 
from Bio import pairwise2
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
        except:
            print("no OPM - id found - add a workaround, e.g. PDBTM")
            return
        else:
            print("Obtain structure from OPM: successful")
        return ref_struc
    reference = reference_chain.split('_')[0]
    ref_chain = reference_chain.split('_')[1]
    ref_struc = get_ref_struc(reference)
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

    target_align_atoms = []
    for res in ref_align_atoms:
        target_align_atoms.append(df.loc[(df['opm_num']==res), 'infile_num'].tolist()[0])




    # Select the model number - normally always the first
    bio_ref_struc = bio_ref_struc_raw[ref_model_id]
    bio_target_struc = bio_target_struc_raw[target_model_id]

    # List of residues to align
    align_ref_atoms = []
    for ind, chain in enumerate(bio_ref_struc_raw.get_chains()):
        if chain.id == ref_chain:
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
    super_imposer.apply(bio_target_struc)

    print(f"RMSD of superimposed structures: {super_imposer.rms}")
    if super_imposer.rms > 10:
        print(f'Structural alignment too bad ({super_imposer.rms} > 10) - align structure manually')
        sys.exit()

    bioio = PDB.PDBIO()
    bioio.set_structure(bio_target_struc)
    bioio.save(filename)
    print(f"Aligned structure saved")

    if write_opm:
        with open(os.path.join(os.path.dirname(filename), 'ref_opm.pdb'), 'w') as fp:
            fp.write(ref_struc)
        print(f"write_opm set to true - OPM structure saved")


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
