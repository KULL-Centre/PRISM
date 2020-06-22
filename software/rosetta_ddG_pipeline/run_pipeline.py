"""run_pipeline.py main executable to run the Rosetta ddG pipeline.

Author: Anders Frederiksen
Contribution: Johanna K.S. Tiemann

Date of last major changes: 2020-04-29

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
from analysis import calc_all
from args_pipeline import parse_args2
from checks import compare_mutfile, pdbxmut
from folders import folder2
from helper import create_symlinks, create_copy, find_copy, get_mut_dict
import mp_prepare
import mp_ddG
from pdb_to_fasta_seq import pdb_to_fasta_seq
from plotting import plot_all
from prism_rosetta_parser import prism_to_mut, read_from_prism
import rosetta_paths
import run_modes
import storeinputs
from structure_input import structure


def predict_stability(args):
    logger.info("Starting pipeline")
    os.chdir(os.getcwd())
    logger.info(f'Current working directory: {os.getcwd()}')

    # Obtain, redirect and adapt user arguments
    chain_id = args.CHAIN
    ddgfile = args.DDG_FLAG_FILE
    mode = args.MODE
    mutation_input = args.MUTATION_INPUT
    outpath = args.OUTPUT_FILE
    overwrite_path = args.OVERWRITE_PATH
    relaxfile = args.RELAX_FLAG_FILE
    structure_list = args.STRUC_FILE
    uniprot_accesion = ''  # args.UNIPROT_ID
    run_struc = args.RUN_STRUC
    ligand = args.LIGAND
    mp_span = args.MP_SPAN_INPUT

    if run_struc == None:
        run_struc = chain_id
    # System name
    name = os.path.splitext(os.path.basename(structure_list))[0]

    # Initiate folder structure
    folder = folder2(outpath, overwrite_path, is_mp=args.IS_MP)

    # Store input files
    input_dict = storeinputs.storeinputfuc(name, args, folder)

    if mode == "proceed" or mode == "relax" or mode == "ddg_calculation":
        mutation_input == "proceed"
        logger.info(f'No preparation, proceeding to execution')

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
        run_name = 'input'

        # adjust mp structure if MP_ALIGN_MODE is selected
        if args.IS_MP == True and args.MP_ALIGN_MODE != 'False':
            logger.info(f'Align the structure along the membrane using {args.MP_CALC_SPAN_MODE}')
            if args.MP_ALIGN_MODE == 'OPM':
                if args.MP_ALIGN_REF != '':
                    run_name = 'input_mp_aligned'
                    structure_instance.path = os.path.join(
                        folder.prepare_mp_superpose, f'{run_name}.pdb')
                    try:
                        mp_prepare.mp_superpose_opm(
                            args.MP_ALIGN_REF, prep_struc, structure_instance.path, target_chain=structure_instance.chain_id, write_opm=True)
                    except:
                        mp_prepare.mp_TMalign_opm(
                            args.MP_ALIGN_REF, prep_struc, structure_instance.path, target_chain=structure_instance.chain_id, write_opm=True)                        
                elif args.UNIPROT_ID != '':
                    logger.error('Uniprot-ID to ref pdb not implemented yet')
                    sys.exit()
                else:
                    logger.error(
                        'No reference or Uniprot-ID provided. Automatic extraction via sequence not yet implemented.')
                    sys.exit()
            else:
                logger.error(
                    'Other modes (PDBTM, TMDET, MemProtMD) not yet implemented.')
                sys.exit()

        structure_dic = get_structure_parameters(
            folder.prepare_checking, structure_instance.path, structure_instance.chain_id)

        # Cleaning pdb and making fasta based on pdb or uniprot-id if provided
        logger.info(f'Prepare the pdb and extract fasta file')
        structure_instance.path_to_cleaned_pdb = structure_instance.clean_up_and_isolate(
            folder.prepare_cleaning, structure_instance.path, structure_instance.chain_id, run_struc, name=run_name, ligand=ligand)
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

        # Get span file for mp from cleaned file if not provided
        if args.IS_MP == True:
            if input_dict['MP_SPAN_INPUT'] == None:
                logger.info(f'Calculate span file with option {args.MP_CALC_SPAN_MODE}')
                if args.MP_CALC_SPAN_MODE == 'DSSP':
                    structure_instance.span = mp_prepare.mp_span_from_pdb_dssp(
                        structure_instance.path_to_cleaned_pdb, folder.prepare_mp_span, thickness=args.MP_THICKNESS, SLURM=False)
                elif args.MP_CALC_SPAN_MODE == 'octopus':
                    structure_instance.span = mp_prepare.mp_span_from_pdb_octopus(
                        structure_instance.path_to_cleaned_pdb, folder.prepare_mp_span, thickness=args.MP_THICKNESS, SLURM=False)
                elif args.MP_CALC_SPAN_MODE == 'False':
                    logger.warn(
                        'No span file provided and no calculation method selected.')
                else:
                    logger.error(
                        'Other modes (struc, bcl, Boctopus) not yet implemented.')
                    sys.exit()
            elif input_dict['MP_SPAN_INPUT'] != None:
                structure_instance.span = create_copy(
                    input_dict['MP_SPAN_INPUT'], folder.prepare_mp_span, name='input.span')

        # Making mutfiles and checks
        print(f'Convert prism file if present: {input_dict["PRISM_INPUT"]}')
        if input_dict['PRISM_INPUT'] == None:
            new_mut_input = input_dict['MUTATION_INPUT']
            mut_dic = get_mut_dict(input_dict['MUTATION_INPUT'])
        else:
            new_mut_input = os.path.join(folder.prepare_input, 'input_mutfile')
            mut_dic = prism_to_mut(input_dict['PRISM_INPUT'], new_mut_input)

        logger.info(f'Generate mutfiles.')
        print(input_dict['MUTATION_INPUT'])
        check2 = structure_instance.make_mutfiles(
            new_mut_input, folder.prepare_mutfiles, structure_dic, structure_instance.chain_id)
        check1 = compare_mutfile(structure_instance.fasta_seq,
                                 folder.prepare_mutfiles, folder.prepare_checking, new_mut_input)
        check3, errors = pdbxmut(folder.prepare_mutfiles, structure_dic)
        check2 = False

        if check1 == True or check2 == True or check3 == True:
            print("check1:", check1, "check2:", check2, "check3:", check3)
            logger.error(
                "ERROR: STOPPING SCRIPT DUE TO RESIDUE MISMATCH BETWEEN MUTFILE AND PDB SEQUENCE")
            sys.exit()

        # Create hard link to mutfile directory and to output structure
        prepare_output_struc = create_copy(
            structure_instance.path_to_cleaned_pdb, folder.prepare_output, name='output.pdb')
        if args.IS_MP == True:
            prepare_output_span_dir = create_copy(folder.prepare_mp_span, f'{folder.prepare_output}', name='spanfiles', directory=True)
        else:
            prepare_output_ddg_mutfile_dir = create_copy(
                folder.prepare_mutfiles, folder.prepare_output, name='mutfiles', directory=True)

        # Copy files for relax & run
        relax_input_struc = create_copy(
            prepare_output_struc, folder.relax_input, name='input.pdb')

        # Generate sbatch files
        logger.info(f'Generate sbatch files')
        if args.IS_MP == True:

            # copy MP relax input files
            logger.info('Copy MP relax input files')
            relax_input_xml = create_copy(
                input_dict['RELAX_XML_INPUT'], folder.relax_input, name='relax.xml')
            relax_input_span_dir = create_copy(
                prepare_output_span_dir, folder.relax_input, name='spanfiles', directory=True)

            # Parse sbatch relax file
            logger.info('Create MP relax sbatch files.')
            path_to_relax_sbatch = mp_prepare.rosetta_relax_mp(
                folder, SLURM=True, num_struc=3, sys_name=name, partition=args.SLURM_PARTITION)

            # Parse sbatch relax parser
            path_to_parse_relax_results_sbatch = structure_instance.parse_relax_sbatch(
                folder, sys_name=f'{name}_relax', sc_name='relax_scores', partition=args.SLURM_PARTITION)

            # Parse sbatch ddg file
            ddg_input_ddgfile = create_copy(
                input_dict['DDG_FLAG_FILE'], folder.ddG_input, name='ddg_flagfile')
            ddg_input_span_dir = create_copy(
                prepare_output_span_dir, folder.ddG_input, name='spanfiles', directory=True)

            if args.MP_PH == -1:
                is_pH = 0
                pH_value = 7
            else:
                is_pH = 1
                pH_value = args.MP_PH
            path_to_ddg_calc_sbatch = mp_ddG.rosetta_ddg_mp_pyrosetta(
                folder, mut_dic, SLURM=True, sys_name=name, partition=args.SLURM_PARTITION,
                repack_radius=args.BENCH_MP_REPACK, lipids=args.MP_LIPIDS,
                temperature=args.MP_TEMPERATURE, repeats=args.BENCH_MP_REPEAT,
                is_pH=is_pH, pH_value=pH_value)
            # Parse sbatch ddg parser
            path_to_parse_ddg_sbatch = mp_ddG.write_parse_rosetta_ddg_mp_pyrosetta_sbatch(
                folder, uniprot=args.UNIPROT_ID, sys_name=name, output_name='ddG.out', partition=args.SLURM_PARTITION)
        else:
            # Parse sbatch relax file
            relax_input_relaxfile = create_copy(
                input_dict['RELAX_FLAG_FILE'], folder.relax_input, name='relax_flagfile')
            path_to_relax_sbatch = structure_instance.rosetta_sbatch_relax(
                folder, relaxfile=relax_input_relaxfile, sys_name=name, partition=args.SLURM_PARTITION)
            # Parse sbatch relax parser
            path_to_parse_relax_results_sbatch = structure_instance.parse_relax_sbatch(
                folder, partition=args.SLURM_PARTITION)

            # Parse sbatch ddg file
            ddg_input_ddgfile = create_copy(
                input_dict['DDG_FLAG_FILE'], folder.ddG_input, name='ddg_flagfile')
            ddg_input_mutfile_dir = create_copy(
                prepare_output_ddg_mutfile_dir, folder.ddG_input, name='mutfiles', directory=True)
            path_to_ddg_calc_sbatch = structure_instance.write_rosetta_cartesian_ddg_sbatch(
                folder, ddg_input_mutfile_dir, ddgfile=ddg_input_ddgfile, sys_name=name, partition=args.SLURM_PARTITION)
            # Parse sbatch ddg parser
            path_to_parse_ddg_sbatch = structure_instance.write_parse_cartesian_ddg_sbatch(
                folder, structure_instance.fasta_seq, structure_instance.chain_id, sys_name=name, partition=args.SLURM_PARTITION)

    # Execution
    # Single SLURM execution
    if mode == 'relax':
        parse_relax_process_id = run_modes.relaxation(folder)
        relax_output_strucfile = find_copy(
            folder.relax_run, '.pdb', folder.relax_output, 'output.pdb')


# if SLURM == False:
#    path_to_scorefile = os.path.join(structure_instance.path_to_run_folder + '/relax_scores.sc')
#    relax_pdb_out = relax_parse_results.parse_relax_results(path_to_scorefile, path_to_run_folder)
# else:
#    path_to_parse_relax_results_sbatch = structure_instance.parse_relax_sbatch(os.path.join(structure_instance.path_to_run_folder + '/relax_scores.sc'), structure_instance.path_to_run_folder)
#    relax_pdb_out = parse_relax_process_id = run_modes.relaxation(structure_instance.path_to_run_folder)
# logger.info(f"Relaxed structure for ddG calculations: {relax_pdb_out}")

    if mode == 'ddg_calculation':
        run_modes.ddg_calculation(folder)
#        ddg_output_score = find_copy(
#            folder.ddG_run, '.sc', folder.ddG_output, 'output.sc')

    if mode == 'analysis':
        calc_all(folder, sys_name=name)
        plot_all(folder, sys_name=name)

    # Full SLURM execution
    if mode == 'proceed' or mode == 'fullrun':
        # Start relax calculation
        parse_relax_process_id = run_modes.relaxation(folder)
        # relax_output_strucfile = find_copy(
        # folder.relax_run, '.pdb', folder.relax_output, 'output.pdb')
        # Start ddG calculation
        # ddg_input_struc = create_copy(
        # os.path.join(folder.relax_output, 'output.pdb'), folder.ddG_input,
        # name='input.pdb')
        run_modes.ddg_calculation(folder, parse_relax_process_id)
#        ddg_output_score = find_copy(
#            folder.ddG_run, '.sc', folder.ddG_output, 'output.sc')


##########################################################################
#                                          Initiate scripts
##########################################################################
if __name__ == '__main__':

    args = parse_args2()
    predict_stability(args)
