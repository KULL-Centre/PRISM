import json
import sys
from os.path import join
import subprocess
from parse_cartesian_functions import rosetta_cartesian_read, ddgs_from_dg
import numpy as np
import re

from argparse import ArgumentParser
from helper import read_slurms, runtime_memory_stats
import os

import pandas as pd

from helper import AttrDict, generate_output
import run_modes


def csv_to_prism(data,structure_input,chain_id):
    """This script convert the regular data file, to a prism-like file. This file contains more information than regular files and are ideal for sharing data"""
    
    with open(structure_input) as json_file:
        strucdata = json.load(json_file)
    path_to_data=data
    data=pd.read_csv(path_to_data)
    path_to_output= str(path_to_data).split('/')
    path_to_output='/'.join(path_to_output[:-1])
    p=100
    
    relative_seq = strucdata["strucdata"][str(chain_id)][2]
    relative_seq
    new_variants=[]
    resdata=strucdata['resdata']
    for n in data['variant']:
        k = n[1:-1]
        pos=(resdata[str(k)][1])
        new_variant=n[0]+str(pos)+n[-1]
        new_variants=np.hstack((new_variants,new_variant))
        
    numpy_array = np.zeros(shape = (len(data), 2))    
    output_df = pd.DataFrame(numpy_array, columns = ["variant", "Rosetta_ddg_score"])
    output_df["variant"] = new_variants
    output_df["Rosetta_ddg_score"] = data['Rosetta_ddg_score']
    
    for m in range(len(data['variant'])):
        
        wt_old=data['variant'][m][0]
        var_old= data['variant'][m][-1]
        pos_new=output_df['variant'][m][1:-1]
        wt_new=output_df['variant'][m][0]
        var_new= output_df['variant'][m][-1]
        if (data['Rosetta_ddg_score'][m] == output_df['Rosetta_ddg_score'][m]) and (wt_old == wt_new and var_old == var_new) and (str(relative_seq[int(pos_new)-1])==str(wt_new)):
            continue
        else:
            break
            print('error')
            print('Check-1   ',str(relative_seq[int(pos_new)-1])==str(wt_new))
            print('Check-2   ',data['Rosetta_ddg_score'][m]==output_df['Rosetta_ddg_score'][m])
            print('Check-3   ',wt_old == wt_new and var_old == var_new)
        
    with open(join(path_to_output, 'data.prism'), 'w') as prism_data:       
        for n in strucdata["DBREF"]:
            if strucdata["DBREF"][str(n)][2] == str(chain_id):
                p=n
        if p != 100:        
            ref_seq=strucdata["strucdata"][str(chain_id)][1]
            coord_seq = strucdata["strucdata"][str(chain_id)][0]
            relative_seq = strucdata["strucdata"][str(chain_id)][2]
            prism_data.write('# ------------------------------------\n')
            prism_data.write('# protein:\n')
            prism_data.write('#\t reference sequence: '+ ref_seq +"\n")
            prism_data.write('#\t coordinate sequence: '+ coord_seq+"\n")
            prism_data.write('#\t Relative sequence: '+ relative_seq+"\n")
            prism_data.write('#\t pdb: '+ strucdata["DBREF"][p][1]+"\n")
            prism_data.write('#\t chain_id: '+ chain_id+"\n")
            output_df.to_csv(prism_data, index = False,sep=' ')

def parse_rosetta_ddgs(sys_name, chain_id, fasta_seq, ddG_run, ddG_output, structure_input, ddG_input, output, prepare_checking, output_gaps=False, zip_files=True, sha_tag='', is_MP=False, scale_factor=2.9):
    """This script parses the results from the ddG calculations into two files. A regular data file containing only the data and the variants and a prism-like file with data variants and additional information"""
    
    runtime_memory_stats(ddG_run)
    
    path_to_run_folder = ddG_run
    
    # Collects results from result files and run the parse cartesian functions script
    rosetta_summary_file = f'{sys_name}_{chain_id}.rosetta_cartesian.dgs'
    shell_command = f'cat *.ddg | grep -v WT > {rosetta_summary_file}'
    subprocess.call(shell_command, cwd=path_to_run_folder, shell=True)

    rosetta_cartesian_ddgs_dict, rosetta_cartesian_ddgs_array = ddgs_from_dg(rosetta_cartesian_read(
        join(path_to_run_folder, rosetta_summary_file), fasta_seq), scale_factor=scale_factor)
    
    
    protein_sequence=fasta_seq 
    aa_order_alphabetical = pd.Series(["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M",
           "N", "P", "Q", "R", "S", "T", "V", "W", "Y"])
    aa_order_group = pd.Series(["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P",
            "A", "V", "I", "L", "M", "F", "Y", "W"])[::-1]

    convert_aa_order_alph_to_group = {}
    for residue, pos_alph in zip(aa_order_alphabetical, range(len(aa_order_alphabetical))):
        pos_group = np.where(aa_order_group == residue)[0]
        convert_aa_order_alph_to_group[pos_alph] = pos_group
    
    #Getting information from structure json file
    with open(structure_input) as json_file:
        strucdata = json.load(json_file)

    accession=''
    protein_name=''
    for n in strucdata['DBREF']:
        print(strucdata['DBREF'][str(n)])
        if str(strucdata['DBREF'][str(n)][2]) == str(chain_id):
            accession=strucdata['DBREF'][str(n)][6]
            protein_name=strucdata['DBREF'][str(n)][7]
    uniprot_accession_name=accession
    pdb_name=sys_name
    PROTEIN_NAME=protein_name
    ddg_sequence=fasta_seq
    
    #Saving files
    ddgs=rosetta_cartesian_ddgs_dict
    numpy_array = np.zeros(shape = (len(ddgs), 2))
    output_df = pd.DataFrame(numpy_array, columns = ["variant", "Rosetta_ddg_score"])
    output_df["variant"] = ddgs.keys()
    output_df["Rosetta_ddg_score"] = ddgs.values()
    output_filename = PROTEIN_NAME+"_"+uniprot_accession_name+"_" + pdb_name +"_rosetta"
    stats = "Positions: " + str(len(fasta_seq)) + ", variants: " + str(len(output_df))
    print(uniprot_accession_name, pdb_name, stats)
    output_file_path = join(ddG_run,output_filename)
    output_df.to_csv(output_file_path + ".csv", index = False)
    #Trying to make .prism file
    with open (join(ddG_run,'ddg.out'), 'w') as fp:
        for elem in rosetta_cartesian_ddgs_array:
            fp.write(f'{elem[0]},{elem[1]},{elem[2]}\n')
    folder = AttrDict()
    folder.update({'prepare_checking': prepare_checking, 'ddG_run': ddG_run,
                   'ddG_output': ddG_output, 'ddG_input': ddG_input, 'output': output})
    generate_output(folder, output_name='ddg.out', sys_name=sys_name, 
        chain_id=chain_id, output_gaps=output_gaps, zip_files=zip_files, sha_tag=sha_tag, MP=is_MP, scale=scale_factor)

#    try:
#        data=output_file_path + ".csv"
#        csv_to_prism(data,structure_input,chain_id)
#    except:
#        print('Cannot create .prism file due to error')
        
    #Checking for warnings and cancels
    try: 
        read_slurms(ddG_run,printing=True)
    except:
        print('cant print warnings due to error. Please check slurm files manually')


def quickcheck( run_dir, base_mut ):
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

    df_add_all = []
    for files in os.listdir (run_dir):
        if files.startswith("mutfile"):
            out_add = os.path.join(run_dir, files)
            if os.path.getsize(out_add) > 0:
                df_add = pd.read_csv(out_add, header=None)
                df_add[4] = df_add[0].apply(lambda x: x.split()[2])
                df_add[5] = df_add[0].apply(lambda x: x.split()[3])
                df_add = df_add.loc[df_add[4].str.startswith('MUT')]
                d1t3 = {'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L', 
                       'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R', 'SER':'S', 'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR': 'Y'}
                def get_vars(x):
                    res = x.split('_')
                    if len(res)>2:
                        res = res[1]
                    else:
                        res = res[1][:-1]
                    mut_var = res[-3:]
                    mut_var = d1t3[mut_var]
                    resi = int(res[:-3])
                    wt_var = df_goal.loc[df_goal['resi']==int(resi), 'variant'].to_list()[0][0]
                    return f"{wt_var}{resi}{mut_var}"
                df_add['variant'] = df_add[4].apply(lambda x: get_vars(x))
                df_add['dG'] = df_add[5]
                df_add = df_add[['variant', 'dG']]
                df_add['rmsd'] = None
            else:
                df_add = pd.DataFrame(columns=['variant', 'dG', 'rmsd', 'resi'])
            df_add['resi'] = df_add['variant'].apply(lambda x: int(x[1:-1]))
            df_add_all.append(df_add)
        else:
            df_add = pd.DataFrame(columns=['variant', 'dG', 'rmsd', 'resi'])
            df_add['resi'] = df_add['variant'].apply(lambda x: int(x[1:-1]))
            df_add_all.append(df_add)
    df_add = pd.concat(df_add_all)
    df_add['base'] = 'MUT'
    df_add['var'] = df_add['variant'].apply(lambda x: x[-1])
    df_add['wt'] = df_add['variant'].apply(lambda x: x[0])
    df_add.loc[df_add['var']==df_add['wt'], 'base'] = 'WT'
    df_add = df_add.drop(columns=['var', 'wt'])
    
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
    df_all

    # do some stats
    max_calculated = len(df_all['variant'])
    current_calculated = len(df_all.loc[(~df_all['dG'].isna())]['variant'])
    missing_calculated = len(df_all.loc[(df_all['dG'].isna())]['variant'])
    print(f"currently calculated: {current_calculated} / {max_calculated} = {(current_calculated/float(max_calculated)):.1%}")
    print(f"not yet calculated: {missing_calculated} / {max_calculated} = {(missing_calculated/float(max_calculated)):.1%}")
    print(f"not yet calculated variants: {df_all.loc[(df_all['dG'].isna())]['variant'].unique()}")

    if (current_calculated == max_calculated) & (missing_calculated == 0):
        all_calculated = True
    else:
        all_calculated = False
    
    return all_calculated, df_all


if __name__ == '__main__':
    print(sys.argv)
    base_mut = os.path.join(sys.argv[13], 'mutation_clean.txt')
    all_calculated, df_all = quickcheck( sys.argv[4], base_mut )
    print(f'all calculated: {all_calculated}')

    if all_calculated:
        parse_rosetta_ddgs(sys_name=sys.argv[1], chain_id=sys.argv[2], fasta_seq=sys.argv[3], 
            ddG_run=sys.argv[4], ddG_output=sys.argv[5], structure_input=sys.argv[6], 
            ddG_input=sys.argv[7], output=sys.argv[8], prepare_checking=sys.argv[9], output_gaps=eval(sys.argv[10]), 
            zip_files=eval(sys.argv[11]), sha_tag=sys.argv[12], is_MP=eval(sys.argv[14]), scale_factor=float(sys.argv[15]))
    else:
        print(f"sys args before rerun: {sys.argv}")
        folder = AttrDict()
        folder.update({'ddG_input': sys.argv[7], 'ddG_run': sys.argv[4]})
        run_modes.ddg_calculation(folder, parse_relax_process_id=None)
