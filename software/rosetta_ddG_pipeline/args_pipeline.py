from argparse import ArgumentParser
import os
import re
import rosetta_paths

def parse_args2():
    """
    Argument parser function
    """
    
    parser = ArgumentParser(description=' -o output_path \n -s Structure.pdb \n -m Mutation_Input (optional) \n -r Relax_flags (optional) \n -c cartesian_ddg_flags (optional) \n -i modes [print,create,proceed,fullrun] \n        print: prints default flag files \n        create: Creates all run files \n        proceed: Starts calculations with created run files \n        fullrun: runs full pipeline')
    
    parser.add_argument( '--structure', '-s',
#    type=lambda s: s.lower() in ['true', 't', 'yes', '1'],
                    default=False,
                    dest="STRUC_FILE",
                    help="Structure file"
                    )
    parser.add_argument( '--uniprot', '-u',
                    default="",
                    dest="UNIPROT_FILE",
                    help="Uniprot sequence"
                    )
    parser.add_argument( '--mutations', '-m',
                    default=None,
                    dest="MUTATION_INPUT",
                    help="mutation input file"
                    )    
    parser.add_argument( '--outputpath', '-o',
                    default=os.getcwd()+"/Run",
                    dest="OUTPUT_FILE",
                    help="Output path"
                    )   
    parser.add_argument( '--cartesianflags', '-c',
                    default=rosetta_paths.path_to_parameters + '/cartesian_ddg_flagfile',
                    dest="CART_FLAG_FILE",
                    help="Cartesian flag file"
                    )
    parser.add_argument( '--relaxflags', '-r',
                    default=rosetta_paths.path_to_parameters + '/relax_flagfile',
                    dest="RELAX_FLAG_FILE",
                    help="Relaxation flag file"
                    )
    parser.add_argument( '--mode', '-i',
                    default="create",
                    dest="MODE",
                    help="Mode to run"
                    )    
    parser.add_argument( '--chainid',
                    default="ignorechain",
                    dest="CHAIN",
                    help="chain ID"
                    )
    parser.add_argument( '--overwrite_path',
                    default=False,
                    dest="OVERWRITE_PATH",
                    help="Overwrites paths when creating folders"
                    ) 
    args = parser.parse_args()
    
    return args