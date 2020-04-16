
"""folder.py creates and stores all relevant folders.

Author: Anders Frederiksen
Contribution: Johanna K.S. Tiemann

Date of last major changes: 2020-04-15

"""

# Standard library imports
from argparse import ArgumentParser, RawTextHelpFormatter
import os
import re

# Local application imports
import rosetta_paths


def parse_args2():
    """
    Argument parser function
    """

    parser = ArgumentParser(description=(
        '***Available options***\n'
        '-o output_path; \n'
        '-s Structure.pdb; \n'
        '-m Mutation_Input (optional); \n'
        '-r Relax_flags (optional); \n'
        '-c cartesian_ddg_flags (optional); \n'
        '-i modes [print,create,proceed,fullrun]; \n'
        '\tprint: prints default flag files \n'
        '\tcreate: Creates all run files \n'
        '\tproceed: Starts calculations with created run files \n'
        '\tfullrun: runs full pipeline'), formatter_class=RawTextHelpFormatter
    )

    parser.add_argument('--structure', '-s',
                        # type=lambda s: s.lower() in ['true', 't', 'yes',
                        # '1'],
                        default=False,
                        dest='STRUC_FILE',
                        help='Structure file'
                        )
    parser.add_argument('--uniprot', '-u',
                        default='',
                        dest='UNIPROT_FILE',
                        help='Uniprot sequence'
                        )
    parser.add_argument('--mutations', '-m',
                        default=None,
                        dest='MUTATION_INPUT',
                        help='mutation input file'
                        )
    parser.add_argument('--outputpath', '-o',
                        default=os.getcwd() + '/Run',
                        dest='OUTPUT_FILE',
                        help='Output path'
                        )
    parser.add_argument('--ddgflags', '-d',
                        default=rosetta_paths.path_to_parameters + '/cartesian_ddg_flagfile',
                        dest='DDG_FLAG_FILE',
                        help='ddG flag file'
                        )
    parser.add_argument('--relaxflags', '-r',
                        default=rosetta_paths.path_to_parameters + '/relax_flagfile',
                        dest='RELAX_FLAG_FILE',
                        help='Relaxation flag file'
                        )
    parser.add_argument('--mode', '-i',
                        choices=['print', 'create', 'proceed',
                                 'fullrun', 'relax', 'ddg_calculation'],
                        default='create',
                        dest='MODE',
                        help=(
                            'Mode to run:\n'
                            '\tprint: prints default flag files \n'
                            '\tcreate: Creates all run files \n'
                            '\tproceed: Starts calculations with created run files (incl. relax and ddG calculation) \n'
                            '\trelax: Starts relax calculations with created run files\n'
                            '\tddg_calculation: Starts ddg_calculation calculations with created run files\n'
                            '\tfullrun: runs full pipeline'
                            'Default value: create'
                        )

                        )
    parser.add_argument('--chainid',
                        default='ignorechain',
                        dest='CHAIN',
                        help='chain ID'
                        )
    parser.add_argument('--overwrite_path',
                        default=False,
                        dest='OVERWRITE_PATH',
                        help='Overwrites paths when creating folders'
                        )
    args = parser.parse_args()

    return args
