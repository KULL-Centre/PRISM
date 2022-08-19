"""run_pipeline.py main executable to run the Rosetta ddG pipeline.

Author: Anders Frederiksen
Contribution: Johanna K.S. Tiemann

Date of last major changes: 2020-04-29

"""

# Standard library imports
import logging 
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
from helper import create_symlinks, create_copy, find_copy, get_mut_dict, read_fasta, check_path, make_no_hetatm_file
import mp_prepare
import mp_ddG
from pdb_to_fasta_seq import pdb_to_fasta_seq
import rosetta_paths
import run_modes
import storeinputs
from structure_input import structure
from make_logs import make_log


def predict_stability(args):
    logging.info("Starting pipeline")
    os.chdir(os.getcwd())
    logging.info(f'Current working directory: {os.getcwd()}')

    # Obtain, redirect and adapt user arguments
    chain_id = args.CHAIN
    ddgfile = check_path(args.DDG_FLAG_FILE)
    mode = args.MODE
    mutation_input = check_path(args.MUTATION_INPUT)
    outpath = check_path(args.OUTPUT_FILE)
    overwrite_path = args.OVERWRITE_PATH
    relaxfile = check_path(args.RELAX_FLAG_FILE)
    structure_list = check_path(args.STRUC_FILE)
    uniprot_accesion = check_path(args.UNIPROT_ID)
    run_struc = args.RUN_STRUC
    ligand = args.LIGAND
    mp_span = args.MP_SPAN_INPUT
    verbose = args.VERBOSE
    partition=args.SLURM_PARTITION

    SHA_TAG = f'{rosetta_paths.sha}tag{rosetta_paths.tag}'

    if run_struc == None:
        run_struc = chain_id
        
    # System name
    name = os.path.splitext(os.path.basename(structure_list))[0]

    # Initiate folder structure
    folder = folder2(outpath, overwrite_path, is_mp=args.IS_MP, mp_multistruc=args.MP_MULTISTRUC_PROTOCOL)
    logger = make_log(folder,verbose)

    # Store input files
    input_dict = storeinputs.storeinputfuc(name, args, folder)

    if mode == "proceed" or mode == "relax" or mode == "ddg_calculation":
        mutation_input == "proceed"
        logger.info(f'No preparation, proceeding to execution')

    # Preprocessing
    if mode == 'create' or mode == 'fullrun':
        logger.info(f'Preparation started')
        # Get input files
        prep_struc = check_path(create_copy(
            input_dict['STRUC_FILE'], folder.prepare_input, name='input.pdb'))

        # Defining structure parameters

        # Create structure instance
        logger.info(f'Creating structure instance')
        structure_instance = structure(chain_id,name,folder,prep_struc,run_struc,logger,input_dict,uniprot_accesion=uniprot_accesion)
        run_name = 'input'

        # adjust mp structure if MP_ALIGN_MODE is selected
        if args.IS_MP == True and (not args.MP_ALIGN_MODE in ['False', 'span']):
            logger.info(f'Align the structure along the membrane using {args.MP_CALC_SPAN_MODE}')
            if args.MP_ALIGN_MODE == 'OPM':
                if args.MP_ALIGN_REF != '-':
                    run_name = 'input_mp_aligned'
                    structure_instance.path = os.path.join(
                        folder.prepare_mp_superpose, f'{run_name}.pdb')
                    try:
                        logger.info('Run OPM alignment with superpose')
                        mp_prepare.mp_superpose_opm(
                            args.MP_ALIGN_REF, prep_struc, structure_instance.path, target_chain=structure_instance.chain_id, write_opm=True, inTM=args.SUPERPOSE_ONTM)
                        logger.info('Superpose successful')
                    except:
                        logger.info('Superpose failed - Run OPM alignment with TM align')
                        mp_prepare.mp_TMalign_opm(
                            args.MP_ALIGN_REF, prep_struc, structure_instance.path, target_chain=structure_instance.chain_id, write_opm=True)                        
                    prep_struc = create_copy(
                        structure_instance.path, folder.prepare_input, name='input.pdb')
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
            folder.prepare_checking, prep_struc)

        # Cleaning pdb and making fasta based on pdb or uniprot-id if provided
        logger.info(f'Prepare the pdb and extract fasta file')
        structure_instance.path_to_cleaned_pdb, struc_dic_cleaned = structure_instance.clean_up_and_isolate(ligand=args.LIGAND)
        structure_instance.fasta_seq_full,structure_instance.fasta_seq = pdb_to_fasta_seq(
            structure_instance.path_to_cleaned_pdb,chain_id)
        if os.path.isfile(uniprot_accesion):
            structure_instance.muscle_align_to_uniprot(structure_instance.uniprot_seq)
        else:
            structure_instance.muscle_align_to_uniprot(structure_instance.fasta_seq)

        # Get span file for mp from cleaned file if not provided
        if args.IS_MP == True:
            if input_dict['MP_SPAN_INPUT'] == None:
                logger.info(f'Calculate span file with option {args.MP_CALC_SPAN_MODE}')
                if args.MP_CALC_SPAN_MODE == 'deepTMHMM':
                    structure_instance.span = mp_prepare.mp_span_from_deepTMHMM(
                        structure_instance.path_to_cleaned_pdb, folder.prepare_mp_span)
                elif args.MP_CALC_SPAN_MODE == 'DSSP':
                    structure_instance.span = mp_prepare.mp_span_from_pdb_dssp(
                        structure_instance.path_to_cleaned_pdb, folder.prepare_mp_span, thickness=args.MP_THICKNESS, SLURM=False)
                elif args.MP_CALC_SPAN_MODE == 'octopus':
                    structure_instance.span = mp_prepare.mp_span_from_pdb_octopus(
                        structure_instance.path_to_cleaned_pdb, folder.prepare_mp_span, thickness=args.MP_THICKNESS, SLURM=False)
                elif args.MP_CALC_SPAN_MODE == 'False':
                    logger.warning(
                        'No span file provided and no calculation method selected.')
                else:
                    logger.error(
                        'Other modes (struc, bcl, Boctopus) not yet implemented.')
                    sys.exit()
            elif input_dict['MP_SPAN_INPUT']:
                structure_instance.span = create_copy(
                    input_dict['MP_SPAN_INPUT'], folder.prepare_mp_span, name='input.span')


            # superpose if span mode on
            if args.MP_ALIGN_MODE == 'span':
                logger.info(f'Superpositioning based on span file')
                run_name = 'input_mp_aligned'
                prep_struc = create_copy(
                        structure_instance.path_to_cleaned_pdb, folder.prepare_cleaning, name='backup.pdb')
                structure_instance.path_to_cleaned_pdb = os.path.join(folder.prepare_mp_superpose, 'superposed_struc.pdb')
                mp_prepare.mp_superpose_span(prep_struc, folder.prepare_mp_superpose, structure_instance.span, 
                    structure_instance.path_to_cleaned_pdb, chain=args.CHAIN)


            logger.info(f'Calculate lipid accessible residues')
            if args.LIGAND:
                no_lig_file = make_no_hetatm_file(structure_instance.path_to_cleaned_pdb)
                lipacc_dic, lipacc_file = mp_prepare.mp_lipid_acc_resi(no_lig_file, folder.prepare_mp_lipacc, folder.prepare_mp_span, thickness=args.MP_THICKNESS, SLURM=False)
            else:
                lipacc_dic, lipacc_file = mp_prepare.mp_lipid_acc_resi(structure_instance.path_to_cleaned_pdb, folder.prepare_mp_lipacc, folder.prepare_mp_span, thickness=args.MP_THICKNESS, SLURM=False)

        # Making mutfiles and checks
        if args.MUT_MODE == 'mut_file':
            new_mut_input = input_dict['MUTATION_INPUT']
        else:
            new_mut_input = None
        #    mut_dic = get_mut_dict(input_dict['MUTATION_INPUT'])

        logger.info(f'Generate mutfiles.')
        logger.info(input_dict['MUTATION_INPUT'])
        
        check2, mut_dic = structure_instance.make_mutfiles(
            new_mut_input, args.MUT_MODE)

        new_mut_input = os.path.join(folder.prepare_cleaning, 'mutation_clean.txt')

        check1 = compare_mutfile(structure_instance.fasta_seq,
                                 folder.prepare_mutfiles, folder.prepare_checking, 
                                 structure_instance.struc_dic_cleaned["resdata"], new_mut_input, chainid=structure_instance.chain_id)
        check3, errors = pdbxmut(folder.prepare_mutfiles, struc_dic_cleaned)
        #check3= False

        if check1 == True or check2 == True or check3 == True:
            logger.info(f"check1: {check1}, check2: {check2}, check3: {check3}")
            logger.error(
                "ERROR: STOPPING SCRIPT DUE TO RESIDUE MISMATCH BETWEEN MUTFILE AND PDB SEQUENCE")
            sys.exit()

        # Create hard link to mutfile directory and to output structure
        prepare_output_struc = check_path(create_copy(
            structure_instance.path_to_cleaned_pdb, folder.prepare_output, name='output.pdb'))
        if args.IS_MP == True:
            prepare_output_span_dir = create_copy(folder.prepare_mp_span, f'{folder.prepare_output}', name='spanfiles', directory=True)
        prepare_output_ddg_mutfile_dir = create_copy(
            folder.prepare_mutfiles, folder.prepare_output, name='mutfiles', directory=True)

        # Copy files for relax & run
        relax_input_struc = check_path(create_copy(
            prepare_output_struc, folder.relax_input, name='input.pdb'))

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
            if (args.MP_CART_DDG != 0) and (args.MP_ENERGY_FUNC == 'franklin2019'):
                relax_weights = args.MP_ENERGY_FUNC_WEIGHTS
            else:
                relax_weights = args.MP_ENERGY_FUNC
                
            path_to_relax_sbatch = mp_prepare.rosetta_relax_mp(
                folder, SLURM=True, repeats=args.BENCH_MP_RELAX_REPEAT, num_struc=args.BENCH_MP_RELAX_STRUCS, 
                lipid_type=args.MP_LIPIDS, sys_name=name, partition=partition, mp_thickness=args.MP_THICKNESS, 
                mp_switch_off=args.MP_IGNORE_RELAX_MP_FLAGS, score_function=relax_weights)

            # Parse sbatch relax parser
            path_to_parse_relax_results_sbatch = structure_instance.parse_relax_sbatch(
                folder, sys_name=f'{name}_relax', partition=args.SLURM_PARTITION, sc_name='relax_scores', 
                mp_multistruc=args.MP_MULTISTRUC_PROTOCOL, is_MP=args.IS_MP, ref_pdb=args.MP_ALIGN_REF, 
                do_checking=args.DO_CHECKING)

            if args.MP_PH == -1:
                is_pH = 0
                pH_value = 7
            else:
                is_pH = 1
                pH_value = args.MP_PH

            # Parse sbatch ddg file
            if  args.MP_MULTISTRUC_PROTOCOL == 0:
                ddg_input_ddgfile = create_copy(
                    input_dict['DDG_FLAG_FILE'], folder.ddG_input, name='ddg_flagfile')
                ddg_input_span_dir = create_copy(
                    prepare_output_span_dir, folder.ddG_input, name='spanfiles', directory=True)
                ddg_input_mutfile_dir = create_copy(
                    prepare_output_ddg_mutfile_dir, folder.ddG_input, name='mutfiles', directory=True)
                path_to_ddg_calc_sbatch = mp_ddG.rosetta_ddg_mp_pyrosetta(
                    folder.ddG_input, folder.ddG_run, mut_dic, SLURM=True, sys_name=name, partition=args.SLURM_PARTITION,
                    repack_radius=args.BENCH_MP_REPACK, lipids=args.MP_LIPIDS,
                    temperature=args.MP_TEMPERATURE, repeats=args.BENCH_MP_REPEAT, dump_pdb = args.DUMP_PDB,
                    is_pH=is_pH, pH_value=pH_value, lipacc_dic=lipacc_file, score_function=args.MP_ENERGY_FUNC,
                    repack_protocol=args.MP_REPACK_PROTOCOL, mutfiles=ddg_input_mutfile_dir, 
                    cartesian=args.MP_CART_DDG, ddgfile=ddg_input_ddgfile, score_function_file=args.MP_ENERGY_FUNC_WEIGHTS)
                if args.MP_CART_DDG:
                    path_to_parse_ddg_sbatch = structure_instance.write_parse_cartesian_ddg_sbatch(
                    folder,  partition=partition, output_gaps=args.GAPS_OUTPUT, zip_files=args.ZIP_FILES, sha_tag=SHA_TAG, is_MP=args.IS_MP, scale_factor=args.SCALE_FACTOR)
                else:
                    path_to_parse_ddg_sbatch = mp_ddG.write_parse_rosetta_ddg_mp_pyrosetta_sbatch(
                        folder, chain_id=args.CHAIN, sys_name=name, output_name='ddG.out', add_output_name='ddG_additional.out', partition=partition, 
                        output_gaps=args.GAPS_OUTPUT, zip_files=args.ZIP_FILES, sha_tag=SHA_TAG)
            else:
                for indi, sub_ddg_folder in enumerate(folder.ddG_input):
                    ddg_input_ddgfile = create_copy(
                        input_dict['DDG_FLAG_FILE'], sub_ddg_folder, name='ddg_flagfile')
                    ddg_input_span_dir = create_copy(
                        prepare_output_span_dir, sub_ddg_folder, name='spanfiles', directory=True)
                    ddg_input_mutfile_dir = create_copy(
                        prepare_output_ddg_mutfile_dir, sub_ddg_folder, name='mutfiles', directory=True)


                    path_to_ddg_calc_sbatch = mp_ddG.rosetta_ddg_mp_pyrosetta(
                        folder.ddG_input[indi], folder.ddG_run[indi], mut_dic, SLURM=True, sys_name=name, partition=args.SLURM_PARTITION,
                        repack_radius=args.BENCH_MP_REPACK, lipids=args.MP_LIPIDS, lowest=0,
                        temperature=args.MP_TEMPERATURE, repeats=1, dump_pdb = args.DUMP_PDB,
                        is_pH=is_pH, pH_value=pH_value, lipacc_dic=lipacc_file, score_function=args.MP_ENERGY_FUNC,
                        repack_protocol=args.MP_REPACK_PROTOCOL, mutfiles=ddg_input_mutfile_dir, 
                        cartesian=args.MP_CART_DDG, ddgfile=ddg_input_ddgfile, score_function_file=args.MP_ENERGY_FUNC_WEIGHTS)

                # Parse sbatch ddg parser
                #folds = [folder.prepare_checking, folder.ddG_run, folder.ddG_output, folder.ddG_input, folder.output]
                if args.MP_CART_DDG:
                    path_to_parse_ddg_sbatch = structure_instance.write_parse_cartesian_ddg_sbatch(
                    folder,  partition=partition, output_gaps=args.GAPS_OUTPUT, zip_files=args.ZIP_FILES, sha_tag=SHA_TAG, is_MP=args.IS_MP, scale_factor=args.SCALE_FACTOR)
                else:
                    path_to_parse_ddg_sbatch = mp_ddG.write_parse_rosetta_ddg_mp_pyrosetta_sbatch(
                        folder, chain_id=args.CHAIN, sys_name=name, output_name='ddG.out', add_output_name='ddG_additional.out', partition=partition, 
                        output_gaps=args.GAPS_OUTPUT, mp_multistruc=args.MP_MULTISTRUC_PROTOCOL, zip_files=args.ZIP_FILES, sha_tag=SHA_TAG)
        else:
            # Parse sbatch relax file
            relax_input_relaxfile = check_path(create_copy(
                input_dict['RELAX_FLAG_FILE'], folder.relax_input, name='relax_flagfile'))
            path_to_relax_sbatch = structure_instance.rosetta_sbatch_relax(
                folder, relaxfile=relax_input_relaxfile, sys_name=name,  partition=partition)
            # Parse sbatch relax parser
            path_to_parse_relax_results_sbatch = structure_instance.parse_relax_sbatch(
                folder, partition=args.SLURM_PARTITION, is_MP=args.IS_MP, ref_pdb='-', 
                do_checking=args.DO_CHECKING)

            # Parse sbatch ddg file
            ddg_input_ddgfile = check_path(create_copy(
                input_dict['DDG_FLAG_FILE'], folder.ddG_input, name='ddg_flagfile'))
            ddg_input_mutfile_dir = create_copy(
                prepare_output_ddg_mutfile_dir, folder.ddG_input, name='mutfiles', directory=True)
            path_to_ddg_calc_sbatch = structure_instance.write_rosetta_cartesian_ddg_sbatch(
                folder, ddg_input_mutfile_dir, ddgfile=ddg_input_ddgfile, sys_name=name,  
                partition=partition, dump_pdb=args.DUMP_PDB)
            # Parse sbatch ddg parser
            path_to_parse_ddg_sbatch = structure_instance.write_parse_cartesian_ddg_sbatch(
                folder,  partition=partition, output_gaps=args.GAPS_OUTPUT, zip_files=args.ZIP_FILES, sha_tag=SHA_TAG, is_MP=args.IS_MP, scale_factor=args.SCALE_FACTOR)


    # Execution
    # Single SLURM execution
    if mode == 'relax':
        parse_relax_process_id = run_modes.relaxation(folder)
        relax_output_strucfile = find_copy(
            folder.relax_run, '.pdb', folder.relax_output, 'output.pdb')

    if mode == 'ddg_calculation':
        run_modes.ddg_calculation(folder,parse_relax_process_id=None, mp_multistruc=args.MP_MULTISTRUC_PROTOCOL)
#        ddg_output_score = find_copy(
#            folder.ddG_run, '.sc', folder.ddG_output, 'output.sc')

    if mode == 'proceed':
        # Check if relax calculation is finished
        parse_relax_process_id = run_modes.check_relax_done_launch_rest(folder)
        # Start ddG calculations (by checking if it did run already)
        run_modes.ddg_calculation(folder, parse_relax_process_id=parse_relax_process_id, mp_multistruc=args.MP_MULTISTRUC_PROTOCOL)

    # Full SLURM execution
    if mode == 'fullrun':
        # Start relax calculation
        parse_relax_process_id = run_modes.relaxation(folder)
        # Start ddG calculations (by checking if it did run already)
        run_modes.ddg_calculation(folder, parse_relax_process_id=parse_relax_process_id, mp_multistruc=args.MP_MULTISTRUC_PROTOCOL)



##########################################################################
#                                          Initiate scripts
##########################################################################
if __name__ == '__main__':

    args = parse_args2()
    print(args)
    predict_stability(args)
