import json
import sys
from os.path import join
import subprocess
from parse_cartesian_functions import rosetta_cartesian_read, ddgs_from_dg
import numpy as np
import re

from argparse import ArgumentParser
from helper import read_slurms
import os

import pandas as pd
from get_memory_stats import check_memory

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

def parse_rosetta_ddgs(sys_name, chain_id, fasta_seq, ddG_run, ddG_output, structure_input):
    """This script parses the results from the ddG calculations into two files. A regular data file containing only the data and the variants and a prism-like file with data variants and additional information"""
    
    #This part collects informations about runtime and memory use
    try:
        job_id_pos= join(ddG_run, 'job_id_ddg.txt')
        with open(job_id_pos, 'r') as job_id_file:
            ddg_process_id=str(job_id_file.readlines()[-1])
            print(ddg_process_id)
         
        #Get stats 
        shell_command = f'sacct --format="JobID,Start,End,CPUTime,ReqMem,MaxRSS,MaxVMSize,AveVMSize,JobName" > memory_usage_{ddg_process_id}.log'
        subprocess.call(shell_command, cwd=ddG_run, shell=True)
                                                
        check_memory(ddg_process_id,ddG_run) 
        
    except: 
        job_id_pos= join(ddG_run, 'job_id_ddg.txt')
        with open(job_id_pos, 'r') as job_id_file:
            ddg_process_id=str(job_id_file.readlines()[-1])
            print(ddg_process_id)
         
        #Get stats 
        shell_command = f'sacct --format="JobID,Start,End,CPUTime,ReqMem,MaxRSS,MaxVMSize,AveVMSize,JobName" > memory_usage_{ddg_process_id}.log'
        subprocess.call(shell_command, cwd=ddG_run, shell=True)
                                                
        check_memory(ddg_process_id,ddG_run)         
        print('no memory file')
    
    
    path_to_run_folder = ddG_run
    
    # Collects results from result files and run the parse cartesian functions script
    rosetta_summary_file = f'{sys_name}_{chain_id}.rosetta_cartesian.dgs'
    shell_command = f'cat *.ddg | grep -v WT > {rosetta_summary_file}'
    subprocess.call(shell_command, cwd=path_to_run_folder, shell=True)

    rosetta_cartesian_ddgs_dict = ddgs_from_dg(rosetta_cartesian_read(
        join(path_to_run_folder, rosetta_summary_file), fasta_seq))
    
    
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
    output_file_path = join(ddG_output,output_filename)
    output_df.to_csv(output_file_path + ".csv", index = False)
    #Trying to make .prism file
    try:
        data=output_file_path + ".csv"
        csv_to_prism(data,structure_input,chain_id)
    except:
        print('Cannot create .prism file due to error')
        
    #Checking for warnings and cancels
    try: 
        read_slurms(ddG_run,printing=True)
    except:
        print('cant print warnings due to error. Please check slurm files manually')
    
if __name__ == '__main__':
    parse_rosetta_ddgs(sys_name=sys.argv[1], chain_id=sys.argv[2], fasta_seq=sys.argv[3], ddG_run=sys.argv[4], ddG_output=sys.argv[5], structure_input=sys.argv[6])
