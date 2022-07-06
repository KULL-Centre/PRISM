
"""folder.py creates and stores all relevant folders.

Author: Anders Frederiksen
Contribution: Johanna K.S. Tiemann

Date of last major changes: 2020-04-15

"""

# Standard library imports
from argparse import ArgumentParser, RawTextHelpFormatter, SUPPRESS
import os
import re

# Local application imports
import rosetta_paths

def parse_args2():
    """
    Argument parser function
    """

    parser = ArgumentParser(description=(
        'Consult the README for further information\n'), formatter_class=RawTextHelpFormatter
    )

    parser.add_argument('--structure', '-s',
                        # type=lambda s: s.lower() in ['true', 't', 'yes',
                        # '1'],
                        #default=False,
                        dest='STRUC_FILE',
                        help='Structure file'
                        )
    parser.add_argument('--outputpath', '-o',
                        default=os.path.join(os.getcwd(), 'Run'),
                        dest='OUTPUT_FILE',
                        help='Output path. default is the current working directory + Run '
                        )
    parser.add_argument('--mode', '-i',
                        choices=['print', 'create', 'proceed',
                                 'fullrun', 'relax', 'ddg_calculation'],
                        default='create',
                        dest='MODE',
                        help=('Mode to run:\n'
                              '\tprint: prints default flag files \n'
                              '\tcreate: Creates all run files \n'
                              '\tproceed: Starts calculations with created run files (incl. relax and ddG calculation) \n'
                              '\trelax: Starts relax calculations with created run files\n'
                              '\tddg_calculation: Starts ddg_calculation calculations with created run files\n'
                              '\tfullrun: runs full pipeline\n'
                              'Default value: create'
                              )
                        )
    parser.add_argument('--mutate_mode', '-mm',
                        choices=['all', 'mut_file'],
                        default='all',
                        dest='MUT_MODE',
                        help=('Mutation modes:\n'
                              '\tall: mutate residues in pdb \n'
                              '\tmut_file: mutate variants present in pipeline mutation file, rosetta mut-file or directory with rosetta mut-files \n'
                              'Default value: all'
                              )
                        )
    parser.add_argument('--mutations', '-m',
                        default=None,
                        dest='MUTATION_INPUT',
                        help='mutation input file'
                        )
    parser.add_argument('--chainid', '-c',
                        default='A',
                        dest='CHAIN',
                        help='chain ID for ddG mutagenesis'
                        )
    parser.add_argument('--run_struc',
                        default=None,
                        dest='RUN_STRUC',
                        help=('Insert what chains you want to be part of the full structure format: ABC \n'
                              'ignorechain for full structure'
                              )
                        )
    parser.add_argument('--ligand',
                        default=None,
                        dest='LIGAND',
                        help='Set to true if you want to keep ligand'
                        )

    parser.add_argument('--overwrite_path',
                        default=False,
                        type=lambda s: s.lower() in ['true', 't', 'yes', '1'],
                        dest='OVERWRITE_PATH',
                        help='Overwrites paths when creating folders'
                        )
    parser.add_argument('--slurm_partition',
                        default='sbinlab',
                        dest='SLURM_PARTITION',
                        help='Partition to run the SLURM jobs'
                        )
    parser.add_argument('--gapped_output',
                        default=False,
                        type=lambda s: s.lower() in ['true', 't', 'yes', '1'],
                        dest='GAPS_OUTPUT',
                        help='Generates prism and pdb files which include gaps and starts with the residue-numbering from the original pdb'
                        )
    parser.add_argument('--dump_pdb', '-dp',
                        default=0,
                        type=lambda s: s.lower() in ['true', 't', 'yes', '1'],
                        dest='DUMP_PDB',
                        help='Dumps all mutant pdbs, default False.'
                        )
    parser.add_argument('--skip_zip', '-nzip',
                        default=False,
                        type=lambda s: s.lower() in ['true', 't', 'yes', '1'],
                        dest='NO_ZIP',
                        help='Skip folder zipping. default off'
                        )
    parser.add_argument('--do_checks',
                        default=True,
                        type=lambda s: s.lower() in ['false', 'f', 'no', '0'],
                        dest='DO_CHECKING',
                        help='Does consistency checkings after relax and before (might be good to set of for HM structures). default on'
                        )
    parser.add_argument('--verbose',
                        default=False,
                        dest='VERBOSE',
                        help='Make pipeline more verbose'
                        )


    parser.add_argument('--is_membrane', '-mp',
                        default=False,
                        type=lambda s: s.lower() in ['true', 't', 'yes', '1'],
                        dest='IS_MP',
                        help='Set to true if you want to run the membrane protein ddG pipeline'
                        )
    parser.add_argument('--mp_span',
                        default=None,
                        dest='MP_SPAN_INPUT',
                        help=('If span input file (defining the membrane spanning region) is prodived \n'
                              'coordinates will be used from there. Otherwise (default) it will be calculated \n'
                              '(see --mp_calc_span_mode).')
                        )
    parser.add_argument('--mp_calc_span_mode',
                        choices=['False', 'deepTMHMM', 'struc', 'DSSP', 'octopus',
                                 'bcl', 'Boctopus'],
                        default='DSSP',
                        dest='MP_CALC_SPAN_MODE',
                        help=('Function/mode to calculate the membrane spanning region/file:\n'
                              '\tFalse: file will not be calculated \n'
                              '\tstruc: uses the information provided the structure \n'
                              '\tDSSP: uses DSSP & pdb orientation to calculate the span region (mp_span_from_pdb) \n'
                              '\toctopus: uses octopus \n'
                              '\tbcl: should be used for helix & beta sheets \n'
                              '\tBoctopus: should be used for beta sheets \n'
                              'Default value: DSSP'
                              )
                        )
    parser.add_argument('--mp_align_ref',
                        default='-',
                        dest='MP_ALIGN_REF',
                        help=('Reference PDB-id to membrane protein alignment.'
                              'Required for --mp_prep_align_mode options [OPM]'
                              'Format: {PDBid}_{chain}'
                              )
                        )
    parser.add_argument('--mp_prep_align_mode',
                        choices=['False', 'OPM', 'PDBTM',
                                 'TMDET', 'MemProtMD'],
                        default='OPM',
                        dest='MP_ALIGN_MODE',
                        help=('Function/mode to align the membrane protein structure:\n'
                              '\tFalse: structure will not be rearranged \n'
                              '\tOPM: uses the information provided in OPM \n'
                              '\tPDBTM: uses the information provided in PDBTM (better than OPM but not tested) \n'
                              '\tTMDET: uses the information provided in TMDET \n'
                              '\tMemProtMD: uses the information provided in MemProtMD \n'
                              'Default value: OPM'
                              )
                        )
    parser.add_argument('--superpose_onTM', '-spTM',
                        default=True,
                        type=lambda s: s.lower() in ['false', 'f', 'no', '0'],
                        dest='SUPERPOSE_ONTM',
                        help='Set to false if you want to align on the complete structure'
                        )
    #additional flags for benchmark or changes to protocol
    parser.add_argument('--ddgflags', '-d',
                        default=os.path.join(
                            rosetta_paths.path_to_data, 'sp', 'cartesian_ddg_flagfile'),
                        dest='DDG_FLAG_FILE',
                        help='ddG flag file'
                        )
    parser.add_argument('--relaxflags', '-r',
                        default=os.path.join(
                            rosetta_paths.path_to_data, 'sp', 'relax_flagfile'),
                        dest='RELAX_FLAG_FILE',
                        help='Relaxation flag file'
                        )
    parser.add_argument('--mp_relax_xml',
                        default=os.path.join(
                            rosetta_paths.path_to_data, 'mp', 'mp_relax.xml'),
                        dest='RELAX_XML_INPUT',
                        help='Relaxation xml file for membrane pipeline'
                        )
    parser.add_argument('--uniprot', '-u',
                        default='',
                        dest='UNIPROT_ID',
                        help=SUPPRESS,
                        #help='Uniprot accession ID'
                        )
    parser.add_argument('--scale', '-sc',
                        default=-999,
                        type=float,
                        dest='SCALE_FACTOR',
                        help=SUPPRESS,
                        #help='scaling of ∆∆G to kcal/mol (default for soluble pipeline: 2.9, for membrane pipeline 1.0)'
                        )


    #influence of features not tested for ddG calculation
    parser.add_argument('--mp_thickness',
                        default=15,
                        type=int,
                        dest='MP_THICKNESS',
                        help=SUPPRESS,
                        #help='Half thickness of membrane.'
                        )
    parser.add_argument('--mp_lipids',
                        choices=['DLPC', 'DMPC', 'DOPC', 'DPPC', 'POPC', 'DLPE', 
                        'DMPE', 'DOPE', 'DPPE', 'POPE', 'DLPG', 'DMPG', 'DOPG', 'DPPG', 'POPG'],
                        default='DLPC',
                        dest='MP_LIPIDS',
                        help=SUPPRESS,
                        #help=('Lipid composition choices:\n'
                        #      '\tDLPC: 1,2-dilauroyl-sn-glycero-3-phosphocholine \n'
                        #      '\tDMPC: 1,2-dimyristoyl-sn-glycero-3-phosphocholine \n'
                        #      '\tDOPC: 1,2-dioleoyl-sn-glycero-3-phosphocholine \n'
                        #      '\tDPPC: 1,2-dipalmitoyl-sn-glycero-3-phosphocholine \n'
                        #      '\tPOPC: 1-palmitoyl-2-oleoyl-glycero-3-phosphocholine \n'
                        #      '\tDLPE: 1,2-dilauroyl-sn-glycero-3-phosphoethanolamine \n'
                        #      '\tDMPE: 1,2-dimyristoyl-sn-glycero-3-phosphoethanolamine \n'
                        #      '\tDOPE: 1,2-dioleoyl-sn-glycero-3-phosphoethanolamine \n'
                        #      '\tDPPE: 1,2-dipalmitoyl-sn-glycero-3-phosphoethanolamine \n'
                        #      '\tPOPE: 1-palmitoyl-2-oleoyl-sn-glycero-3-phosphoethanolamine \n'
                        #      '\tDLPG: 1,2-dilauroyl-sn-glycero-3-phospho-(1’-rac-glycerol) \n'
                        #      '\tDMPG: 1,2-dimyristoyl-sn-glycero-3-phospho-(1’-rac-glycerol) \n'
                        #      '\tDOPG: 1,2-dioleoyl-sn-glycero-3-phospho-(1’-rac-glycerol) \n'
                        #      '\tDPPG: 1,2-dipalmitoyl-sn-glycero-3-phospho-(1’-rac-glycerol) \n'
                        #      '\tPOPG: 1-palmitoyl-2-oleoyl-sn-glycero-3-phospho-(1’-rac-glycerol) \n'
                        #      'Default value: DLPC'
                        #      )
                        )
    parser.add_argument('--mp_temperature',
                        default=20.0,
                        type=float,
                        dest='MP_TEMPERATURE',
                        help=SUPPRESS,
                        #help='Experimental temperature. Default value: 37.0'
                        )
    parser.add_argument('--mp_pH',
                        default=-1.0,
                        type=float,
                        dest='MP_PH',
                        help=SUPPRESS,
                        #help=('Experimental pH value between 0-14. -1=off.\n'
                        #      'Default value: -1')
                        )

    # for MP benchmarking only - otherwise defaults set
    parser.add_argument('--benchmark_mp_repack',
                        default=8.0,
                        type=float,
                        dest='BENCH_MP_REPACK',
                        help=SUPPRESS,
                        #help='For benchmark purpose: repack value'
                        )
    parser.add_argument('--benchmark_mp_repeat',
                        default=5,
                        type=int,
                        dest='BENCH_MP_REPEAT',
                        help=SUPPRESS,
                        #help='For benchmark purpose: repeat value'
                        )
    parser.add_argument('--benchmark_mp_relax_repeat',
                        default=5,
                        type=int,
                        dest='BENCH_MP_RELAX_REPEAT',
                        help=SUPPRESS,
                        #help='For benchmark purpose: relax repeat value'
                        )
    parser.add_argument('--benchmark_mp_relax_strucs',
                        default=20,
                        type=int,
                        dest='BENCH_MP_RELAX_STRUCS',
                        help=SUPPRESS,
                        #help='For benchmark purpose: relax structure value output'
                        )
    parser.add_argument('--mp_ignore_relax_mp_flags',
                        default=False,
                        type=lambda s: s.lower() in ['true', 't', 'yes', '1'],
                        dest='MP_IGNORE_RELAX_MP_FLAGS',
                        help=SUPPRESS,
                        #help='For relax checking'
                        )
    parser.add_argument('--mp_energy_func',
                        default='franklin2019',
                        dest='MP_ENERGY_FUNC',
                        help=SUPPRESS,
                        #help='MP Energy function (mainly for benchmarking). Examples: franklin2019, mpframework_smooth_fa_2012, ref2015_memb.'
                        )
    parser.add_argument('--mp_energy_func_weights',
                        default=os.path.join(
                            rosetta_paths.path_to_data, 'mp', 'f19_cart_1.5.wts'),
                        dest='MP_ENERGY_FUNC_WEIGHTS',
                        help=SUPPRESS,
                        #help='MP Energy function (mainly for benchmarking). Examples: franklin2019, mpframework_smooth_fa_2012, ref2015_memb.'
                        )
    parser.add_argument('--mp_repack_protocol',
                        default='MP_flex_relax_ddG',
                        dest='MP_REPACK_PROTOCOL',
                        choices=['MP_repack', 'MP_flex_relax_ddG'],
                        help=SUPPRESS,
                        #help="MP repacking algorithm (mainly for benchmarking). Default=MP_repack, other options are 'MP_flex_relax_ddG' "
                        )
    parser.add_argument('--mp_multistruc_protocol',
                        default=0,
                        type=int,
                        dest='MP_MULTISTRUC_PROTOCOL', 
                        help=SUPPRESS,
                        #help="MP generates x relaxed structures and calculates exactly 1 ddG from each structure  (mainly for benchmarking). Default=False "
                        )
    parser.add_argument('--mp_cart_ddg',
                        default=1,
                        type=int,
                        dest='MP_CART_DDG', 
                        help=SUPPRESS,
                        #help="MP cartesian space. Default=1 "
                        )
    

    args = parser.parse_args()

    # Handle user input errors
    if args.MUT_MODE == 'mut_file' and args.MUTATION_INPUT == None:
        parser.error("Please specify a mutation input file or change the mutation mode.")
    if args.MUTATION_INPUT != None:
        print('Mutation mode changed to mut_file')
        args.MUT_MODE = 'mut_file'
    if args.DUMP_PDB != 0:
        args.DUMP_PDB = 1
    if args.NO_ZIP != False:
        args.ZIP_FILES = False
    else:
        args.ZIP_FILES = True
    if args.LIGAND != None:
        args.LIGAND = True
    if args.IS_MP == False:
        args.MP_MULTISTRUC_PROTOCOL == 0
    if args.IS_MP:
        if args.SUPERPOSE_ONTM != True:
            args.SUPERPOSE_ONTM = False
        if (args.MP_ALIGN_MODE=='OPM') and (args.MP_ALIGN_REF == '-'):
            parser.error('Please specify a reference PDBid and chain for alginment into membrane plane or switch "MP_ALIGN_MODE" to false')
            sys.exit()
        if args.MP_ALIGN_REF != '-':
            args.MP_ALIGN_MODE='OPM'
        if args.MP_CART_DDG == 1:
            args.RELAX_XML_INPUT = os.path.join(rosetta_paths.path_to_data, 'mp', 'mp_cart_relax.xml')
        if args.MP_SPAN_INPUT != None:
            args.MP_ALIGN_REF = '-'
    if args.SCALE_FACTOR == -999:
        if args.IS_MP:
            args.SCALE_FACTOR = 1
        else:
            args.SCALE_FACTOR = 2.9
    return args
