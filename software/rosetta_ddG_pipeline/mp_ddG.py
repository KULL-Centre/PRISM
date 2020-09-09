"""mp_ddG.py includes several scriptsrelated to the ddG calculation of membrane proteins.

Author: Johanna K.S. Tiemann

Date of last major changes: 2020-05-01

"""

# Standard library imports
import logging as logger
import json
import os
import sys


# Local application imports
from helper import AttrDict, create_copy, generate_output, runtime_memory_stats
import rosetta_paths


def rosetta_ddg_mp_pyrosetta(folder, mut_dict, SLURM=True, sys_name='',
                             partition='sbinlab', output_name='ddG.out', 
                             add_output_name='ddG_additional.out', repack_radius=0,
                             lipids='DLPC', temperature=37.0, repeats=3,
                             score_file_name='scores', is_pH=0, pH_value=7):
    ddg_script_exec = os.path.join(
        rosetta_paths.path_to_stability_pipeline, 'rosetta_mp_ddG_adapted.py')
    input_struc = os.path.join(folder.ddG_input, 'input.pdb')
    output_file = os.path.join(folder.ddG_run, output_name)
    score_file = os.path.join(folder.ddG_run, f'{score_file_name}.sc')
    for root, dirs, files in os.walk(folder.ddG_input):
        for file in files:
            if file.endswith('.span'):
                input_span = os.path.join(root, file)

    ddG_command = (f'python3 {ddg_script_exec}'
                   f' --in_pdb {input_struc}'
                   f' --in_span {input_span}'
                   f' --out {output_file}'
                   f' --out_add {add_output_name}'
                   f' --repack_radius {repack_radius}'
                   f' --output_breakdown {score_file}'
                   f' --include_pH {is_pH}'
                   f' --pH_value {pH_value}'
                   f' --repeats {repeats}'
                   f' --lipids {lipids}'
                   f' --temperature {temperature}'
                   '')

    if SLURM:
        path_to_sbatch = os.path.join(folder.ddG_input, 'rosetta_ddg.sbatch')
        with open(path_to_sbatch, 'w') as fp:
            fp.write(f'''#!/bin/sh
#SBATCH --job-name=mp_ddG_{sys_name}
#SBATCH --array=0-{len(mut_dict.keys())}
#SBATCH --time=32:00:00
#SBATCH --mem 5000
#SBATCH --partition={partition}
#SBATCH --nice
RESIS=({' '.join(mut_dict.keys())})
MUTS=({' '.join(mut_dict.values())})
OFFSET=0
INDEX=$((OFFSET+SLURM_ARRAY_TASK_ID))
echo $INDEX

# launching rosetta
''')
            new_ddG_command = ddG_command + \
                ' --res ${RESIS[$INDEX]} --mut ${MUTS[$INDEX]} '
            fp.write(new_ddG_command)
        logger.info(path_to_sbatch)

    else:
        logger.warn("need to write the script!")
        for resid in mut_dict.keys():
            new_ddG_command = ddG_command + (f'--res {resid} '
                                             f'--mut {mut_dict[resid]} ')
            logger.info(new_ddG_command)
        sys.exit()


def postprocess_rosetta_ddg_mp_pyrosetta(folder, output_name='ddG.out', sys_name='', version=1, prims_nr='XXX', chain_id='A', output_gaps=False):
  runtime_memory_stats(folder.ddG_run)
  generate_output(folder, output_name=output_name, sys_name=sys_name, version=version, prims_nr=prims_nr, chain_id=chain_id, output_gaps=output_gaps)
    # The ddg_file should only contain the data, looking like this:
    # M1T,-0.52452 # first value=variant, second=mean([var1-WT1,var2-WT2, ...]) - comma separated
    # M1Y,0.2352,0.2342,.... # it may contain more values like the var-mut of
    # each run

    #create_copy(os.path.join(folder.prepare_checking, 'fasta_file.fasta'), folder.output, name=f'{sys_name}_seq.fasta')
    #create_copy(prims_file_pdb_nmbr, folder.output)
    #create_copy(prims_file, folder.output)



def write_parse_rosetta_ddg_mp_pyrosetta_sbatch(folder, chain_id='A', sys_name='input', output_name='ddG.out', partition='sbinlab', output_gaps=False):
    score_sbatch_path = os.path.join(folder.ddG_input, 'parse_ddgs.sbatch')
    with open(score_sbatch_path, 'w') as fp:
        fp.write(f'''#!/bin/sh 
#SBATCH --job-name=collect_rosetta_ddgs_{sys_name} 
#SBATCH --array=1 
#SBATCH --nodes=1 
#SBATCH --time=0:10:00 
#SBATCH --partition={partition} 

#This sbatch script launches the parse postprocess_rosetta_ddg_mp_pyrosetta function, from the mp_ddG 
''')
        fp.write((f'python3 {os.path.join(rosetta_paths.path_to_stability_pipeline, "mp_ddG.py")} '
                  f'{folder.prepare_checking} {folder.ddG_run} {folder.ddG_output} '
                  f'{folder.ddG_input} {folder.output} '
                  f'{chain_id} {sys_name} {output_name} {output_gaps}'
                  ))
    return score_sbatch_path


if __name__ == '__main__':
    folder = AttrDict()
    folder.update({'prepare_checking': sys.argv[1], 'ddG_run': sys.argv[2],
                   'ddG_output': sys.argv[3], 'ddG_input': sys.argv[4], 'output': sys.argv[5]})

    if len(sys.argv) > 9:
        postprocess_rosetta_ddg_mp_pyrosetta(
            folder, chain_id=sys.argv[6], sys_name=sys.argv[7], output_name=sys.argv[8], output_gaps=sys.argv[9])
    if len(sys.argv) > 8:
        postprocess_rosetta_ddg_mp_pyrosetta(
            folder, chain_id=sys.argv[6], sys_name=sys.argv[7], output_name=sys.argv[8])
    elif len(sys.argv) > 7:
        postprocess_rosetta_ddg_mp_pyrosetta(
            folder, chain_id=sys.argv[6], sys_name=sys.argv[7])
    elif len(sys.argv) > 6:
        postprocess_rosetta_ddg_mp_pyrosetta(folder, chain_id=sys.argv[6])
    else:
        postprocess_rosetta_ddg_mp_pyrosetta(folder)
