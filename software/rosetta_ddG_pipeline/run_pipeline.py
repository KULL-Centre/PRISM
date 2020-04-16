"""run_pipeline.py main executable to run the Rosetta ddG pipeline.

Author: Anders Frederiksen
Contribution: Johanna K.S. Tiemann

Date of last major changes: 2020-04-15

"""

# Standard library imports
import logging as logger
import os
import re
import shutil
import sys

# Third party imports
# import getopt

# Local application imports
from AnalyseStruc import get_structure_parameters
from args_pipeline import parse_args2
from checks import compare_mutfile, pdbxmut
from folders import folder2
from helper import create_symlinks, create_copy, find_copy
from pdb_to_fasta_seq import pdb_to_fasta_seq
import rosetta_paths
import run_modes
import storeinputs
from structure_input import structure


def predict_stability(argv):
    logger.info("Starting pipeline")
    os.chdir(os.getcwd())
    logger.info(f'Current working directory: {os.getcwd()}')

    # Obtain, redirect and adapt user arguments
    args = parse_args2()

    chain_id = args.CHAIN
    ddgfile = args.DDG_FLAG_FILE
    mode = args.MODE
    mutation_input = args.MUTATION_INPUT
    outpath = args.OUTPUT_FILE
    overwrite_path = args.OVERWRITE_PATH
    relaxfile = args.RELAX_FLAG_FILE
    structure_list = args.STRUC_FILE
    uniprot_accesion = args.UNIPROT_FILE
    # System name
    name = os.path.splitext(os.path.basename(structure_list))[0]

    # Initiate folder structure
    folder = folder2(outpath, overwrite_path)

    # Store input files
    input_dict = storeinputs.storeinputfuc(name, args, folder)

    if mode == "proceed" or mode == "relax" or mode == "ddg_calculation":
        mutation_input == "proceed"
        logger.info(f'No preparation, proceeding to execution')

##########################################################################
#                                         HANDLING THE STRUCTURES
##########################################################################

    # Preprocessing
    if mode == 'create' or mode == 'fullrun':
        logger.info(f'Preparation started')
        # Get input files
        prep_struc = create_copy(
            input_dict['STRUC_FILE'], folder.prepare_input, name='input.pdb')

        # Defining structure parameters

        # Create structure instance
        logger.info(f'Creating structure instance')
        structure_instance = structure()
        structure_instance.sys_name = name
        structure_instance.chain_id = chain_id
        structure_instance.path = prep_struc
        structure_dic = get_structure_parameters(
            folder.prepare_checking, prep_struc,chain_id)

        # Cleaning pdb and making fasta based on pdb or uniprot-id if provided
        logger.info(f'Prepare the pdb and extract fasta file')
        structure_instance.path_to_cleaned_pdb = structure_instance.clean_up_and_isolate(
            folder.prepare_cleaning, prep_struc, chain_id, name='input')
        structure_instance.fasta_seq = pdb_to_fasta_seq(
            structure_instance.path_to_cleaned_pdb)
        if uniprot_accesion != "":
            structure_instance.uniprot_seq = structure_instance.read_fasta(
                uniprot_accesion)
            structure_instance.muscle_align_to_uniprot(
                folder.prepare_checking, structure_instance.uniprot_seq)
        else:
            structure_instance.muscle_align_to_uniprot(
                folder.prepare_checking, structure_instance.fasta_seq)

        # Making mutfiles and checks
        logger.info(f'Generate mutfiles.')
        check2 = structure_instance.make_mutfiles(
            input_dict['MUTATION_INPUT'], folder.prepare_mutfiles, structure_dic)
        # check1 = compare_mutfile(structure_instance.fasta_seq,structure_instance.path_to_run_folder,mutation_input)
        #check3, errors = pdbxmut(folder.prepare_mutfiles, resdata)
        check1 = False
        check2 = False
        check3 = False
        if check1 == True or check2 == True or check3 == True:
            logger.error(
                "ERROR: STOPPING SCRIPT DUE TO RESIDUE MISMATCH BETWEEN MUTFILE AND PDB SEQUENCE")
            sys.exit()

        # Create hard link to mutfile directory and to output structure
        prepare_output_ddg_mutfile_dir = create_copy(
            folder.prepare_mutfiles, folder.prepare_output, name='mutfiles', directory=True)
        prepare_output_struc = create_copy(
            structure_instance.path_to_cleaned_pdb, folder.prepare_output, name='output.pdb')

        # Generate sbatch files
        logger.info(f'Generate sbatch files')
        # Parse sbatch relax file
        relax_input_struc = create_copy(
            prepare_output_struc, folder.relax_input, name='input.pdb')
        relax_input_relaxfile = create_copy(
            input_dict['RELAX_FLAG_FILE'], folder.relax_input, name='relax_flagfile')
        path_to_relax_sbatch = structure_instance.rosetta_sbatch_relax(
            folder, relaxfile=relax_input_relaxfile, sys_name=name)
        # Parse sbatch relax parser
        path_to_parse_relax_results_sbatch = structure_instance.parse_relax_sbatch(
            folder, sys_name=name)

        # Parse sbatch ddg file
        ddg_input_mutfile_dir = create_copy(
            prepare_output_ddg_mutfile_dir, folder.ddG_input, name='mutfiles', directory=True)
        ddg_input_ddgfile = create_copy(
            input_dict['DDG_FLAG_FILE'], folder.ddG_input, name='ddg_flagfile')
        path_to_ddg_calc_sbatch = structure_instance.write_rosetta_cartesian_ddg_sbatch(
            folder, ddg_input_mutfile_dir, ddgfile=ddg_input_ddgfile, sys_name=name)
        # Parse sbatch ddg parser
        path_to_parse_ddg_sbatch = structure_instance.write_parse_cartesian_ddg_sbatch(
            folder)

    # Execution
    # Single SLURM execution
    if mode == 'relax':
        parse_relax_process_id = run_modes.relaxation(folder)
        relax_output_strucfile = find_copy(
            folder.relax_run, '.pdb', folder.relax_output, 'output.pdb')

    if mode == 'ddg_calculation':
        run_modes.ddg_calculation(folder)

    # Full SLURM execution
    if mode == 'proceed' or mode == 'fullrun':
        # Start relax calculation
        parse_relax_process_id = run_modes.relaxation(folder)
        relax_output_strucfile = find_copy(
            folder.relax_run, '.pdb', folder.relax_output, 'output.pdb')
        # Start ddG calculation
        ddg_input_struc = create_copy(
            os.path.join(folder.relax_output, 'output.pdb'), folder.ddG_input, name='input.pdb')
        run_modes.ddg_calculation(folder, parse_relax_process_id)
        ddg_output_score = find_copy(
            folder.ddG_run, '.sc', folder.ddG_output, 'output.sc')

   # if mode == "print":
   #     open("relax_flag_file_copy", "w").writelines(open(relax_flag_file).readlines())
   #     open("ddg_flag_file_copy", "w").writelines(open(ddg_flag_file).readlines())

##########################################################################
#                                          Initiate scripts
##########################################################################
if __name__ == '__main__':
    predict_stability(sys.argv[1:])
