"""storeinputs.py stores all input files.

Author: Anders Frederiksen
Contribution: Johanna K.S. Tiemann

Date of last major changes: 2020-04-15

"""

# Standard library imports
from json import dump
import logging as logger
from os.path import join, basename
from shutil import copyfile


def storeinputfuc(name, args, folder):

    input_dict = {}
    logger.info(f'Input files can be found at {folder.input}')
    # Create customized info file
    with open(join(folder.input, 'inputs.info'), 'w') as file:
        file.write(
            (f'NAME = {name}\n'
             f'CHAIN = {args.CHAIN}\n'
             f'STRUCTURE = {args.STRUC_FILE}\n'
             f'OUTPATH = {folder.output_path}\n'
             f'{args.UNIPROT_FILE}\n')
        )

    # Create json file from input args
    with open(join(folder.input, 'args.info'), 'w') as fp:
        dump(vars(args), fp, indent=4)

    # Copy input files, incl rosetta flag files (relax, ddg), structure and
    # mutation file (if specified)
    copyfile(args.RELAX_FLAG_FILE, join(
        folder.input, basename(args.RELAX_FLAG_FILE)))
    input_dict['RELAX_FLAG_FILE'] = join(
        folder.input, basename(args.RELAX_FLAG_FILE))

    copyfile(args.DDG_FLAG_FILE, join(
        folder.input, basename(args.DDG_FLAG_FILE)))
    input_dict['DDG_FLAG_FILE'] = join(
        folder.input, basename(args.DDG_FLAG_FILE))

    copyfile(args.STRUC_FILE, join(
        folder.input, basename(args.STRUC_FILE)))
    input_dict['STRUC_FILE'] = join(
        folder.input, basename(args.STRUC_FILE))

    if args.MUTATION_INPUT:
        copyfile(args.MUTATION_INPUT, join(folder.input, 'mutations'))
        input_dict['MUTATION_INPUT'] = join(folder.input, 'mutations')
    return input_dict
