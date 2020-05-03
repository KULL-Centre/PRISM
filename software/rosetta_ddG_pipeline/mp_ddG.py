"""mp_ddG.py includes several scriptsrelated to the ddG calculation of membrane proteins.

Author: Johanna K.S. Tiemann

Date of last major changes: 2020-05-01

"""

# Standard library imports
import logging as logger
import os
import sys


# Local application imports
import rosetta_paths


def rosetta_ddg_mp_pyrosetta(folder, mut_dict, SLURM=True, sys_name='',
                             partition='sbinlab', output_name='ddG.out', repack_radius=0,
                             score_file_name='scores', is_pH=0, pH_value=7):
    ddg_script_exec = os.path.join(
        rosetta_paths.path_to_stability_pipeline, 'rosetta_mp_ddG_adapted.py')
    input_struc = os.path.join(folder.ddG_input, 'input.pdb')
    input_span = os.path.join(folder.ddG_input, 'spanfiles', 'input.span')
    output_file = os.path.join(folder.ddG_run, output_name)
    score_file = os.path.join(folder.ddG_run, f'{score_file_name}.sc')

    ddG_command = (f'python3 {ddg_script_exec} '
                   f'--in_pdb {input_struc} '
                   f'--input_span {input_span} '
                   f'--repack_radius {repack_radius}'
                   f'--out {output_file}'
                   f'--output_breakdown {score_file}'
                   f'--include_pH {is_pH} '
                   f'--pH_value {pH_value} '
                   '')

    if SLURM:
        path_to_sbatch = os.path.join(folder.ddG_input, 'rosetta_ddg.sbatch')
        with open(path_to_sbatch, 'w') as fp:
            fp.write(f'''#!/bin/sh
# SBATCH --job-name=mp_ddG_{sys_name}
# SBATCH --array=0-{len(mut_dict.keys())}
# SBATCH --time=32:00:00
# SBATCH --mem 5000
# SBATCH --partition={partition}
# SBATCH --nice
RESIS={tuple(mut_dict.keys())}
MUTS={tuple(mut_dict.values())}
OFFSET=0
INDEX=$((OFFSET+SLURM_ARRAY_TASK_ID))
echo $INDEX

# launching rosetta
''')
            new_ddG_command = ddG_command + \
                '--res ${{RESIS[$INDEX]}} --mut ${{MUTS[$INDEX]}} '
            fp.write(new_ddG_command)
        logger.info(path_to_sbatch)

    else:
        logger.warn("need to write the script!")
        for resid in mut_dict.keys():
            new_ddG_command = ddG_command + (f'--res {resid} '
                                             f'--mut {mut_dict[resid]} ')
            logger.info(new_ddG_command)
        sys.exit()
