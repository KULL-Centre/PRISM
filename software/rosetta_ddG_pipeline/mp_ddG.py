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
import run_modes


def rosetta_ddg_mp_pyrosetta(folder_ddG_input, folder_ddG_run, mut_dict, SLURM=True, sys_name='',
                             partition='sbinlab', output_name='ddG.out', 
                             add_output_name='ddG_additional.out', repack_radius=0,
                             lipids='DLPC', temperature=37.0, repeats=5, lowest=1,
                             score_file_name='scores', is_pH=0, pH_value=7, dump_pdb=0,
                             score_function='franklin2019', repack_protocol='MP_repack', 
                             lipacc_dic='mp_lipid_acc_dic.json', mutfiles='mutfiles', 
                             cartesian=0, ddgfile='', score_function_file='f19_cart_1.5.wts'):

    input_struc = os.path.join(folder_ddG_input, 'input.pdb')
    for root, dirs, files in os.walk(folder_ddG_input):
        for file in files:
            if file.endswith('.span'):
                input_span = os.path.join(root, file)
    if cartesian==0:
        ddg_script_exec = os.path.join(
            rosetta_paths.path_to_stability_pipeline, 'pyrosetta_ddG.py')

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
    
# cartesian mp ddg calculation
    else:
        path_to_sbatch = os.path.join(folder_ddG_input, 'rosetta_ddg.sbatch')
        if ddgfile == '':
            # path_to_ddgflags = os.path.join(
            # rosetta_paths.path_to_data, 'sp', 'cartesian_ddg_flagfile')
            path_to_ddgflags = os.path.join(folder_ddG_input, 'ddg_flagfile')
        else:
            path_to_ddgflags = ddgfile

        if mutfiles == '':
            mutfiles = os.path.join(folder_ddG_input, 'mutfiles')
        muts = os.listdir(mutfiles)

        with open(path_to_sbatch, 'w') as fp:
            fp.write(f'''#!/bin/sh 
#SBATCH --job-name={sys_name}_MPcartddg
#SBATCH --array=0-{len(muts)-1}
#SBATCH --time=48:00:00
#SBATCH --mem 2000
#SBATCH --partition={partition}
#SBATCH --nice 
LST=(`ls {mutfiles}/mutfile*`)
OFFSET=0 
INDEX=$((OFFSET+SLURM_ARRAY_TASK_ID))
echo $INDEX

# launching rosetta 
if test -f "{input_struc}"; then
''')
            fp.write((f'\t{os.path.join(rosetta_paths.path_to_rosetta, f"bin/cartesian_ddg.{rosetta_paths.Rosetta_extension}")} '
                      f'-database {rosetta_paths.Rosetta_database_path} '
                      f'-s {input_struc} '
                      f'-ddg:mut_file ${{LST[$INDEX]}} '
                      f'-score:weights {score_function_file} '
                      f'-mp:lipids:composition {lipids} '
                      f'-ddg:iterations {repeats} '
                      f'-ddg::dump_pdbs {dump_pdb} '
                      f'-mp:setup:spanfiles {input_span} '
                      f'-in:membrane '
                      f'-fa_max_dis 9 '
                      f'-missing_density_to_jump '
                      f'-has_pore 0 '
                      f'-ddg:legacy true '
                      f'-ddg:optimize_proline 1 '
                      f'-ddg:frag_nbrs 4 '
                      f'-ddg:bbnbrs 1 '
                      f'-ddg::cartesian '
                      f'-out:prefix ddg-$SLURM_ARRAY_JOB_ID-$SLURM_ARRAY_TASK_ID '))#@{path_to_ddgflags} '))
            fp.write(f'''\nelse
\techo "{input_struc} does not exist - exiting the call"
fi
''')
        logger.info(path_to_sbatch)
        return path_to_sbatch

def postprocess_rosetta_ddg_mp_pyrosetta(folder, output_name='ddG.out', sys_name='', version=1, prism_nr='XXX', chain_id='A', zip_files=True, output_gaps=False, mp_multistruc=0, sha_tag='', scale=2.9):
    print(mp_multistruc, folder)
    if mp_multistruc == 0:
        runtime_memory_stats(folder.ddG_run)
        generate_output(folder, output_name=output_name, sys_name=sys_name, version=version, prism_nr=prism_nr, chain_id=chain_id, output_gaps=output_gaps, zip_files=zip_files, sha_tag=sha_tag, MP=True, scale=scale)
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


        generate_output(folder, output_name=output_name, sys_name=sys_name, version=version, prism_nr=prism_nr, chain_id=chain_id, output_gaps=output_gaps, zip_files=zip_files, sha_tag=sha_tag, MP=True)


    # The ddg_file should only contain the data, looking like this:
    # M1T,-0.52452 # first value=variant, second=mean([var1-WT1,var2-WT2, ...]) - comma separated
    # M1Y,0.2352,0.2342,.... # it may contain more values like the var-mut of
    # each run

    #create_copy(os.path.join(folder.prepare_checking, 'fasta_file.fasta'), folder.output, name=f'{sys_name}_seq.fasta')
    #create_copy(prism_file_pdb_nmbr, folder.output)
    #create_copy(prism_file, folder.output)



def write_parse_rosetta_ddg_mp_pyrosetta_sbatch(folder, chain_id='A', sys_name='input', output_name='ddG.out', add_output_name='ddG_additional.out',
        partition='sbinlab', output_gaps=False, mp_multistruc=0, zip_files=True, sha_tag=''):
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
                  f'{chain_id} {zip_files} {sys_name} {output_name} {output_gaps} '
                  f'{sha_tag} '
                  f'{add_output_name} {folder.prepare_cleaning}'
                  ))
        if mp_multistruc != 0:
            fp.write(f' {mp_multistruc} {folder.ddG_postparse_run} {folder.ddG_postparse_output}')
    return score_sbatch_path


def quickcheck( out, out_add, base_mut ):
    # the following variants should be calculated
    arr = []
    with open(base_mut, 'r') as fp:
        for line in fp:
            line = line.split()
            wt = line[0]
            resi = line[1]
            varrs = line[2].strip()
            for var in varrs:
                if var != wt:
                    arr.append(f"{wt}{resi}{var}")
    df_goal = pd.DataFrame(data={'variant':arr})
    df_goal['target'] = True
    df_goal['resi'] = df_goal['variant'].apply(lambda x: int(x[1:-1]))

    # the following variants are already calculated 5/max-calc times
    if os.path.getsize(out) > 0:
        df = pd.read_csv(out, header=None)
        df = df.rename(columns={0:'variant', 1:'mean_ddG', 2:'std_ddG'})
    else:
        df = pd.DataFrame(columns=['variant', 'mean_ddG', 'std_ddG'])

    # the following variants are already calculated but not written out in the final file (only upon competen of all runs)
    if os.path.getsize(out_add) > 0:
        df_add = pd.read_csv(out_add, header=None)
        df_add = df_add.rename(columns={0:'variant', 1:'dG', 2:'rmsd', 3: 'resi'})
    else:
        df_add = pd.DataFrame(columns=['variant', 'dG', 'rmsd', 'resi'])
    df_add['base'] = df_add['variant'].apply(lambda x: x.split('_')[0])
    df_add['variant'] = df_add['variant'].apply(lambda x: x.split('_')[1])
    df_add['resi'] = df_add['variant'].apply(lambda x: int(x[1:-1]))
    # calculating the counts
    df_add_count = df_add.copy()
    df_add_count = df_add_count.groupby(['variant']).count()['dG']#.unique()
    df_add_count = pd.DataFrame(df_add_count).reset_index(drop=False)
    df_add_count = df_add_count.rename(columns={'dG': 'count'})
    df_add = pd.merge(df_add, df_add_count, on=['variant'], how='outer')
    df_add

    # merging
    df_all = pd.merge(df_goal, df_add, on=['variant', 'resi'], how='outer')
    df_all['target'] = df_all['target'].fillna(False)
    df_all = pd.merge(df_all, df, on=['variant'], how='outer')
    df_all

    # do some stats
    max_written_out = len(df_all.loc[(~df_all.duplicated(subset=['variant'], keep='first')) & (df_all['base']=='MUT')]['variant'])
    current_written_out = len(df_all.loc[(~df_all.duplicated(subset=['variant'], keep='first')) & (df_all['base']=='MUT') & (~df_all['mean_ddG'].isna())]['variant'])
    missing_written_out = len(df_all.loc[(~df_all.duplicated(subset=['variant'], keep='first')) & (df_all['base']=='MUT') & (df_all['mean_ddG'].isna())]['variant'])
    logger.info(f"written out in final file: {current_written_out} / {max_written_out} = {(current_written_out/float(max_written_out)):.1%}")
    logger.info(f"not yet written out in final file: {missing_written_out} / {max_written_out} = {(missing_written_out/float(max_written_out)):.1%}")

    max_calculated = len(df_all['variant'])
    current_calculated = len(df_all.loc[(~df_all['dG'].isna())]['variant'])
    missing_calculated = len(df_all.loc[(df_all['dG'].isna())]['variant'])
    logger.info(f"currently calculated: {current_calculated} / {max_calculated} = {(current_calculated/float(max_calculated)):.1%}")
    logger.info(f"not yet calculated: {missing_calculated} / {max_calculated} = {(missing_calculated/float(max_calculated)):.1%}")
    logger.info(f"not yet calculated variants: {df_all.loc[(df_all['dG'].isna())]['variant'].unique()}")

    if (current_written_out == max_written_out) & (missing_written_out == 0):
        all_written_out = True
    else:
        all_written_out = False
    if (current_calculated == max_calculated) & (missing_calculated == 0):
        all_calculated = True
    else:
        all_calculated = False
    
    return all_written_out, all_calculated, df_all






if __name__ == '__main__':

    out = os.path.join(sys.argv[2], sys.argv[9])#ddG.out
    out_add = os.path.join(sys.argv[2], sys.argv[12])#ddG_additional.out
    base_mut = os.path.join(sys.argv[13], 'mutation_clean.txt')
    all_written_out, all_calculated, df_all = quickcheck( out, out_add, base_mut )

    if all_written_out:
        print(sys.argv)
        folder = AttrDict()
        folder.update({'prepare_checking': sys.argv[1], 'ddG_run': sys.argv[2],
                       'ddG_output': sys.argv[3], 'ddG_input': sys.argv[4], 'output': sys.argv[5]})

        if len(sys.argv) > 14:
            folder.update({'ddG_postparse_run': sys.argv[15], 'ddG_postparse_output': sys.argv[16]})
            print(folder)
            postprocess_rosetta_ddg_mp_pyrosetta(
                folder, chain_id=sys.argv[6], zip_files=sys.argv[7], sys_name=sys.argv[8], output_name=sys.argv[9], output_gaps=sys.argv[10], mp_multistruc=sys.argv[14], sha_tag=sys.argv[11] )
        elif len(sys.argv) > 9:
            postprocess_rosetta_ddg_mp_pyrosetta(
                folder, chain_id=sys.argv[6], zip_files=sys.argv[7], sys_name=sys.argv[8], output_name=sys.argv[9], output_gaps=sys.argv[10], sha_tag=sys.argv[11])
        elif len(sys.argv) > 8:
            postprocess_rosetta_ddg_mp_pyrosetta(
                folder, chain_id=sys.argv[6], zip_files=sys.argv[7], sys_name=sys.argv[8], output_name=sys.argv[9], sha_tag=sys.argv[10])
        elif len(sys.argv) > 7:
            postprocess_rosetta_ddg_mp_pyrosetta(
                folder, chain_id=sys.argv[6], zip_files=sys.argv[7], sys_name=sys.argv[8], sha_tag=sys.argv[9])
        elif len(sys.argv) > 6:
            postprocess_rosetta_ddg_mp_pyrosetta(folder, chain_id=sys.argv[6], zip_files=sys.argv[7], sha_tag=sys.argv[8])
        else:
            postprocess_rosetta_ddg_mp_pyrosetta(folder)
    else:
        print(sys.argv)
        folder = AttrDict()
        folder.update({'ddG_input': sys.argv[4], 'ddG_run': sys.argv[2]})
        if len(sys.argv) > 14:
            folder.update({'ddG_postparse_run': sys.argv[15], 'ddG_postparse_output': sys.argv[16]})
            print(folder)
            run_modes.ddg_calculation(folder, parse_relax_process_id=None, mp_multistruc=sys.argv[14])
        else:
            run_modes.ddg_calculation(folder, parse_relax_process_id=None)
