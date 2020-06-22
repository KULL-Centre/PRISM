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

# Third party imports
from Bio import PDB
import numpy as np


# Local application imports
import pdb_to_fasta_seq
import rosetta_paths
from mp_helper import extract_from_opm


def mp_superpose_opm(reference_chain, target, filename, target_chain='A',
                     ref_model_id=0, target_model_id=0,
                     ref_align_atoms=[], target_align_atoms=[], write_opm=False):
    # Adapted from https://gist.github.com/andersx/6354971
    # Copyright (c) 2010-2016 Anders S. Christensen

    # Get reference structure e.g. from OPM
    def get_ref_struc(keyword):
        try:
            ref_struc = extract_from_opm(keyword)
        except:
            logger.error("no OPM - id found - add a workaround, e.g. PDBTM")
            return
        else:
            logger.info("Obtain structure from OPM: successful")
        return ref_struc
    reference = reference_chain.split('_')[0]
    ref_chain = reference_chain.split('_')[1]
    ref_struc = get_ref_struc(reference)
    # Parse reference (from string) and target structure (from file)
    parser = PDB.PDBParser(QUIET=True)
    bio_ref_struc_raw = parser.get_structure(
        "reference", io.StringIO(ref_struc))
    bio_target_struc_raw = parser.get_structure("target", target)
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

    logger.info(f"RMSD of superimposed structures: {super_imposer.rms}")

    bioio = PDB.PDBIO()
    bioio.set_structure(bio_target_struc)
    bioio.save(filename)
    logger.info(f"Aligned structure saved")

    if write_opm:
        with open(os.path.join(os.path.dirname(filename), 'ref_opm.pdb'), 'w') as fp:
            fp.write(ref_struc)
        logger.info(f"write_opm set to true - OPM structure saved")


def mp_TMalign_opm(reference_chain, target, filename, target_chain='A',
                   ref_model_id=0, target_model_id=0,
                   ref_align_atoms=[], target_align_atoms=[], write_opm=False):
    # Get reference structure e.g. from OPM
    def get_ref_struc(keyword):
        try:
            ref_struc = extract_from_opm(keyword)
        except:
            logger.error("no OPM - id found - add a workaround, e.g. PDBTM")
            return
        else:
            logger.info("Obtain structure from OPM: successful")
        return ref_struc
    reference = reference_chain.split('_')[0]
    ref_chain = reference_chain.split('_')[1]
    ref_struc = get_ref_struc(reference)

    ref_struc_dest = os.path.join(os.path.dirname(filename), 'ref_opm.pdb')
    with open(ref_struc_dest, 'w') as fp:
        fp.write(ref_struc)
    logger.info(f"OPM structure saved")


    path_to_TMalign = rosetta_paths.path_to_TMalign
    shell_call = f'{path_to_TMalign} {target} {ref_struc_dest} -o {filename.split(".pdb")[0]}'
    subprocess.call(shell_call, shell=True)

    logger.info(f"Aligned structure saved")



def mp_span_from_pdb_octopus(pdbinput, outdir_path, SLURM=False):

    Rosetta_span_exec = os.path.join(
        rosetta_paths.path_to_rosetta, 'bin/spanfile_from_pdb.{rosetta_paths.Rosetta_extension}')
    span_command = f'{Rosetta_span_exec} -in:file:s {pdbinput}'
    logger.info(f"Span call function: {span_command}")

    if SLURM:
        logger.warn("need to write the slurm script!")
        span_call = subprocess.Popen('sbatch rosetta_span.sbatch', stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT, shell=True, cwd=outdir_path)
        span_process_id_info = span_call.communicate()
        span_process_id = str(span_process_id_info[0]).split()[3][0:-3]
        logger.info(span_process_id_info)
    else:
        span_call = subprocess.Popen(
            span_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, cwd=outdir_path)
        span_process_id_info = span_call.communicate()

    with open(os.path.join(outdir_path, 'span.log'), 'w') as fp:
        fp.writelines(span_process_id_info[0].decode())

    if not any(fname.endswith('.span') for fname in os.listdir(outdir_path)):
        logger.error("Span file not created. Check log file.")
    else:
        spanfiles = []
        for fname in os.listdir(outdir_path):
            if (fname.endswith('.span')):
                new_filename = f"{os.path.splitext(os.path.basename(pdbinput))[0]}.span"
                os.rename(os.path.join(outdir_path, fname),
                          os.path.join(outdir_path, new_filename))
                spanfiles.append(os.path.join(outdir_path, new_filename))
        logger.info("Span process done")
        return spanfiles


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
    logger.info(f"Span call function: {span_command}")

    if SLURM:
        logger.warn("need to write the slurm script!")
        sys.exit()
        span_call = subprocess.Popen('sbatch rosetta_span.sbatch', stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT, shell=True, cwd=outdir_path)
        span_process_id_info = span_call.communicate()
        span_process_id = str(span_process_id_info[0]).split()[3][0:-3]
        logger.info(span_process_id_info)
    else:
        span_call = subprocess.Popen(
            span_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, cwd=outdir_path)
        span_process_id_info = span_call.communicate()

    with open(os.path.join(outdir_path, 'span.log'), 'w') as fp:
        fp.writelines(span_process_id_info[0].decode())
    pdbinput_dir = os.path.dirname(pdbinput)
    if not any(fname.endswith('.span') for fname in os.listdir(pdbinput_dir)):
        logger.error("Span file not created. Check log file.")
    else:
        spanfiles = []
        for fname in os.listdir(pdbinput_dir):
            if (fname.endswith('.span')):
                new_filename = f"{os.path.splitext(os.path.basename(pdbinput))[0]}.span"
                os.rename(os.path.join(pdbinput_dir, fname),
                          os.path.join(outdir_path, new_filename))
                spanfiles.append(os.path.join(outdir_path, new_filename))
        logger.info("Span process done")
        return spanfiles


def rosetta_relax_mp(folder, SLURM=False, num_struc=3, sys_name='mp', partition='sbinlab', repeats=2, lipid_type='DLPC'):
    Rosetta_script_exec = os.path.join(
        rosetta_paths.path_to_rosetta, f'bin/rosetta_scripts.{rosetta_paths.Rosetta_extension}')
    for root, dirs, files in os.walk(folder.relax_input):
        for file in files:
            if file.endswith('.span'):
                spanfile = os.path.join(root, file)
    relax_command = (f'{Rosetta_script_exec} '
                      # Use the membrane relax protocol Rosetta script
                      f'-parser:protocol {os.path.join(folder.relax_input, "relax.xml")} '
                      # Repeatition of FastRelax
                      f'-parser:script_vars repeats={repeats} '
                      # Input PDB Structure: PDB file for protein structure
                      f'-in:file:s {os.path.join(folder.relax_input, "input.pdb")} '
                      # Spanfile describing trans-membrane spans of the
                      # starting structure
                      f'-mp:setup:spanfiles {spanfile} '
                      '-mp:scoring:hbond '  # Turn on membrane depth-dependent hydrogen bonding weight
                      f'-mp:lipids:composition {lipid_type} '
                      # Use the FastRelax mode of Rosetta Relax (uses 5-8
                      # repeat cycles)
                      '-relax:fast '
                      '-relax:jump_move true '  # Allow the MEM and other jumps to move during refinement
                      # Number of structures to generate
                      f'-nstruct {num_struc} '
                      # Wait to pack until the membrane mode is turned on
                      '-packing:pack_missing_sidechains 0 '
                      '-out:pdb '  # Output all PDB structures of refined models
                      # Specify destination for score file
                      f'-out:file:scorefile {os.path.join(folder.relax_run, "relax_scores.sc")} '
                      # f'-mp::thickness {thickness} '
                     f'-database {rosetta_paths.Rosetta_database_path} '
                      '-ignore_unrecognized_res true '
                      '-score:weights franklin2019 '
                      #                        '-ignore_zero_occupancy false '
                      '')
    logger.info(f'MP relax call function: {relax_command}')

    if SLURM:
        path_to_sbatch = os.path.join(
            folder.relax_input, 'rosetta_relax.sbatch')
        with open(path_to_sbatch, 'w') as fp:
            fp.write(f'''#!/bin/sh
#SBATCH --job-name=relax_{sys_name}
#SBATCH --time=10:00:00
#SBATCH --mem 5000
#SBATCH --partition={partition}

# launching rosetta relax 
''')
            fp.write(f'{relax_command} ')
        logger.info(f'Location of relax sbatch file: {path_to_sbatch}')
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
