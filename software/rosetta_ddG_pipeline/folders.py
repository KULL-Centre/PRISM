"""folder.py creates and stores all relevant folders.

Author: Anders Frederiksen
Contribution: Johanna K.S. Tiemann

Date of last major changes: 2020-04-15

"""

# Standard library imports
import logging as logger
from os import makedirs
from os.path import join
from sys import exit

# Local application imports
from args_pipeline import parse_args2
from helper import check_path


def check_paths(dir_path, overwrite_path=True):
    try:
        makedirs(dir_path)
        logger.info(f'Global working directory ({dir_path}) created.')
    except FileExistsError:
        if overwrite_path:
            logger.warn(f'Directory {dir_path} already exists. overwrite_path set to {overwrite_path}, so we will use this path.')
        else:
            logger.error(f'Directory {dir_path} already exists. overwrite_path set to {overwrite_path}, so we stop the execution. Please provide a different output directory.')
            exit()  # will terminate the complete script
    dir_path = check_path(dir_path)
    return dir_path


class folder2:

    def __init__(self, output_path, overwrite_path, is_mp=False):
        # Create
        # Global folder
        self.output_path = check_paths(
            output_path, overwrite_path=overwrite_path)

        # Main folders
        self.logs = check_paths(join(self.output_path, 'logs'))
        self.input = check_paths(join(self.output_path, 'input'))
        self.prepare = check_paths(join(self.output_path, 'prepare'))
        self.relax = check_paths(join(self.output_path, 'relax'))
        self.ddG = check_paths(join(self.output_path, 'ddG'))
        self.output = check_paths(join(self.output_path, 'output'))

        # Subfolders
        self.prepare_input = check_paths(join(self.prepare, 'input'))
        self.prepare_mutfiles = check_paths(join(self.prepare, 'mutfiles'))
        self.prepare_cleaning = check_paths(join(self.prepare, 'cleaning'))
        self.prepare_checking = check_paths(join(self.prepare, 'checking'))
        self.prepare_output = check_paths(join(self.prepare, 'output'))

        if is_mp:
            # Membrane Protein sub-subfolders
            self.prepare_mp_files = check_paths(join(self.prepare, 'mp_files'))
            self.prepare_mp_superpose = check_paths(
                join(self.prepare_mp_files, 'superpose'))
            self.prepare_mp_span = check_paths(
                join(self.prepare_mp_files, 'membrane_span'))
            self.prepare_mp_lipacc = check_paths(
                join(self.prepare_mp_files, 'mp_lipid_acc'))

        self.relax_input = check_paths(join(self.relax, 'input'))
        self.relax_run = check_paths(join(self.relax, 'run'))
        self.relax_output = check_paths(join(self.relax, 'output'))

        self.ddG_input = check_paths(join(self.ddG, 'input'))
        self.ddG_run = check_paths(join(self.ddG, 'run'))
        self.ddG_output = check_paths(join(self.ddG, 'output'))

        return
