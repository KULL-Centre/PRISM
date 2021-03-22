"""mp_ddG.py includes several scriptsrelated to the ddG calculation of membrane proteins.

Author: Johanna K.S. Tiemann

Date of last major changes: 2020-05-01

"""

# Standard library imports
import logging as logger
import json
import os
import sys

import pandas as pd


# Local application imports
from helper import AttrDict, create_copy, generate_output, runtime_memory_stats
import rosetta_paths


def rosetta_ddg_mp_pyrosetta(folder_ddG_input, folder_ddG_run, mut_dict, SLURM=True, sys_name='',
                             partition='sbinlab', output_name='ddG.out', 
                             add_output_name='ddG_additional.out', repack_radius=0,
                             lipids='DLPC', temperature=37.0, repeats=5, lowest=1,
                             score_file_name='scores', is_pH=0, pH_value=7, dump_pdb=0,
                             score_function='franklin2019', repack_protocol='MP_repack', 
                             lipacc_dic='mp_lipid_acc_dic.json', mutfiles='mutfiles'):
    ddg_script_exec = os.path.join(
        rosetta_paths.path_to_stability_pipeline, 'pyrosetta_ddG.py')
    input_struc = os.path.join(folder_ddG_input, 'input.pdb')
    for root, dirs, files in os.walk(folder_ddG_input):
        for file in files:
            if file.endswith('.span'):
                input_span = os.path.join(root, file)

    ddG_command = (f'python3 {ddg_script_exec}'
                   f' --in_pdb {input_struc}'
                   f' --in_span {input_span}'
                   f' --outdir {folder_ddG_run}'
                   f' --outfile {output_name}'
                   f' --out_add {add_output_name}'
                   f' --repack_radius {repack_radius}'
                   f' --include_pH {is_pH}'
                   f' --pH_value {pH_value}'
                   f' --repeats {repeats}'
                   f' --lowest {lowest}'
                   f' --lipids {lipids}'
                   f' --temperature {temperature}'
                   f' --score_function {score_function}'
                   f' --repack_protocol {repack_protocol}'
                   f' --lip_has_pore {lipacc_dic}'
                   f' --dump_pdb {dump_pdb}'
                   '')


    if SLURM:
        path_to_sbatch = os.path.join(folder_ddG_input, 'rosetta_ddg.sbatch')
        if mutfiles == '':
            mutfiles = os.path.join(folder_ddG_input, 'mutfiles')
        muts = os.listdir(mutfiles)

        with open(path_to_sbatch, 'w') as fp:
            fp.write(f'''#!/bin/sh 
#SBATCH --job-name={sys_name}_MPddg
#SBATCH --array=0-{len(muts)-1}
#SBATCH --time=48:00:00
#SBATCH --mem 5000
#SBATCH --partition={partition}
#SBATCH --nice 
LST=(`ls {mutfiles}/mutfile*`)
OFFSET=0 
INDEX=$((OFFSET+SLURM_ARRAY_TASK_ID))
echo $INDEX

# launching rosetta 
''')
            new_ddG_command = ddG_command + \
                ' --mutfile ${LST[$INDEX]} '
            fp.write(new_ddG_command)
        logger.info(path_to_sbatch)
    else:
        logger.warn("need to write the script!")
        sys.exit()


def postprocess_rosetta_ddg_mp_pyrosetta(folder, output_name='ddG.out', sys_name='', version=1, prism_nr='XXX', chain_id='A', output_gaps=False, mp_multistruc=0):
    print(mp_multistruc, folder)
    if mp_multistruc == 0:
        runtime_memory_stats(folder.ddG_run)
        generate_output(folder, output_name=output_name, sys_name=sys_name, version=version, prism_nr=prism_nr, chain_id=chain_id, output_gaps=output_gaps)
    else:
        dfs = [['variant', 'ddG', 'std']]
        for sub_ddG_folder in folder.ddG_run.split(":"):
            with open(os.path.join(sub_ddG_folder, output_name), 'r') as fp2:
                for line in fp2:
                    dfs.append(line.split(','))
        dfs = pd.DataFrame(data=dfs[1:], columns=dfs[0])
        dfs = dfs[['variant', 'ddG']]
        dfs['ddG'] = dfs['ddG'].astype(float)
        dfs = dfs.groupby('variant').agg({'ddG': ['mean', 'std']}) 
        dfs.to_csv(os.path.join(folder.ddG_postparse_run, output_name), sep=',', na_rep='', header=False, index=True)
        folder.ddG_run = folder.ddG_postparse_run
        folder.ddG_output = folder.ddG_postparse_output


        generate_output(folder, output_name=output_name, sys_name=sys_name, version=version, prism_nr=prism_nr, chain_id=chain_id, output_gaps=output_gaps)


    # The ddg_file should only contain the data, looking like this:
    # M1T,-0.52452 # first value=variant, second=mean([var1-WT1,var2-WT2, ...]) - comma separated
    # M1Y,0.2352,0.2342,.... # it may contain more values like the var-mut of
    # each run

    #create_copy(os.path.join(folder.prepare_checking, 'fasta_file.fasta'), folder.output, name=f'{sys_name}_seq.fasta')
    #create_copy(prism_file_pdb_nmbr, folder.output)
    #create_copy(prism_file, folder.output)



def write_parse_rosetta_ddg_mp_pyrosetta_sbatch(folder, chain_id='A', sys_name='input', output_name='ddG.out', 
        partition='sbinlab', output_gaps=False, mp_multistruc=0):
    if mp_multistruc != 0:
        ddG_run = ":".join(folder.ddG_run)
        ddG_output = ":".join(folder.ddG_output )
        ddG_input = folder.ddG_input[0]
        ddG_input_script = folder.ddG_postparse_input
    else:
        ddG_run = folder.ddG_run
        ddG_output = folder.ddG_output
        ddG_input = folder.ddG_input
        ddG_input_script = folder.ddG_input
    score_sbatch_path = os.path.join(ddG_input_script, 'parse_ddgs.sbatch')
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
                  f'{folder.prepare_checking} {ddG_run} {ddG_output} '
                  f'{ddG_input} {folder.output} '
                  f'{chain_id} {sys_name} {output_name} {output_gaps}'
                  ))
        if mp_multistruc != 0:
            fp.write(f' {mp_multistruc} {folder.ddG_postparse_run} {folder.ddG_postparse_output}')
    return score_sbatch_path


if __name__ == '__main__':
    print(sys.argv)
    folder = AttrDict()
    folder.update({'prepare_checking': sys.argv[1], 'ddG_run': sys.argv[2],
                   'ddG_output': sys.argv[3], 'ddG_input': sys.argv[4], 'output': sys.argv[5]})

    if len(sys.argv) > 10:
        folder.update({'ddG_postparse_run': sys.argv[11], 'ddG_postparse_output': sys.argv[12]})
        print(folder)
        postprocess_rosetta_ddg_mp_pyrosetta(
            folder, chain_id=sys.argv[6], sys_name=sys.argv[7], output_name=sys.argv[8], output_gaps=sys.argv[9], mp_multistruc=sys.argv[10] )
    elif len(sys.argv) > 8:
        postprocess_rosetta_ddg_mp_pyrosetta(
            folder, chain_id=sys.argv[6], sys_name=sys.argv[7], output_name=sys.argv[8], output_gaps=sys.argv[9])
    elif len(sys.argv) > 7:
        postprocess_rosetta_ddg_mp_pyrosetta(
            folder, chain_id=sys.argv[6], sys_name=sys.argv[7], output_name=sys.argv[8])
    elif len(sys.argv) > 6:
        postprocess_rosetta_ddg_mp_pyrosetta(
            folder, chain_id=sys.argv[6], sys_name=sys.argv[7])
    elif len(sys.argv) > 5:
        postprocess_rosetta_ddg_mp_pyrosetta(folder, chain_id=sys.argv[6])
    else:
        postprocess_rosetta_ddg_mp_pyrosetta(folder)
