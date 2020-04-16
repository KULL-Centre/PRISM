"""run_modes.py contains helpful functions.

Author: Anders Frederiksen
Contribution: Johanna K.S. Tiemann

Date of last major changes: 2020-04-15

"""

# Standard library imports
import logging as logger
from os.path import join
import subprocess


def relaxation(folder):

    relax_call = subprocess.Popen(f'sbatch {join(folder.relax_input, "rosetta_relax.sbatch")}', stdout=subprocess.PIPE, shell=True, cwd=folder.relax_run)

    relax_process_id_info = relax_call.communicate()
    relax_process_id = str(relax_process_id_info[0]).split()[3][0:-3]

    parse_relaxation_call = subprocess.Popen(f'sbatch --dependency=afterany:{relax_process_id} {join(folder.relax_input, "parse_relax.sbatch")}', stdout=subprocess.PIPE, shell=True, cwd=folder.relax_run)

    parse_relax_process_id_info = parse_relaxation_call.communicate()
    parse_relax_process_id = str(
        parse_relax_process_id_info[0]).split()[3][0:-3]

    return parse_relax_process_id


def ddg_calculation(folder, parse_relax_process_id=None):
    if parse_relax_process_id == None:
        dependency = ''
        parse_relax_process_id = ''
    else:
        dependency = '--dependency=afterany:'

    ddg_call = subprocess.Popen(f'sbatch {dependency}{parse_relax_process_id} {join(folder.ddG_input, "rosetta_ddg.sbatch")}', stdout=subprocess.PIPE, shell=True, cwd=folder.ddG_run)

    ddg_process_id_info = ddg_call.communicate()

    logger.info(f'ddG process ID info: {ddg_process_id_info}')
    ddg_process_id = str(ddg_process_id_info[0]).split()[3][0:-3]

    parse_results_call = subprocess.Popen(f'sbatch --dependency=afterany:{ddg_process_id} {join(folder.ddG_input, "parse_ddgs.sbatch")}', stdout=subprocess.PIPE, shell=True, cwd=folder.ddG_run)

    return
