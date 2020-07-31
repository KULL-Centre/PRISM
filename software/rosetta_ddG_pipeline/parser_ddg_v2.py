import json
import sys
from os.path import join
import subprocess
from parse_cartesian_functions import rosetta_cartesian_read, ddgs_from_dg
import numpy as np
import re
import seaborn as sns
from argparse import ArgumentParser
from helper import read_slurms
import os
import matplotlib
matplotlib.use('Agg')
import pandas as pd


def csv_to_prism(data,structure_input,chain_id):
    with open(structure_input) as json_file:
        strucdata = json.load(json_file)
    resdata=strucdata['resdata']
    path_to_hm3=data
    data=pd.read_csv(path_to_hm3)
    path_to_output= str(path_to_hm3).split('/')
    path_to_output='/'.join(path_to_output[:-1])
    p=100
    with open(join(path_to_output, 'data.prism'), 'w') as prism_data:       
        for n in strucdata["DBREF"]:
            if strucdata["DBREF"][str(n)][2] == str(chain_id):
                p=n
        if p != 100:        
            ref_seq=strucdata["strucdata"][str(chain_id)][1]
            coord_seq = strucdata["strucdata"][str(chain_id)][0]
            prism_data.write('# ------------------------------------\n')
            prism_data.write('# protein:\n')
            prism_data.write('#\t reference sequence: '+ ref_seq +"\n")
            prism_data.write('#\t coordinate sequence: '+ coord_seq+"\n")
            prism_data.write('#\t pdb: '+ strucdata["DBREF"][p][1]+"\n")
            prism_data.write('#\t chain_id: '+ chain_id+"\n")
            data.to_csv(prism_data, index = False,sep=' ')

def parse_rosetta_ddgs(sys_name, chain_id, fasta_seq, ddG_run, ddG_output, structure_input):

    path_to_run_folder = ddG_run
    print('the path to run folder is')
    print(path_to_run_folder)

    rosetta_summary_file = f'{sys_name}_{chain_id}.rosetta_cartesian.dgs'
    shell_command = f'cat *.ddg | grep -v WT > {rosetta_summary_file}'
    print('calling to the shell:')
    print(shell_command)
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

    ddgs=rosetta_cartesian_ddgs_dict
    numpy_array = np.zeros(shape = (len(ddgs), 2))
    output_df = pd.DataFrame(numpy_array, columns = ["variant", "Rosetta_ddg_score"])
    output_df["variant"] = ddgs.keys()
    output_df["Rosetta_ddg_score"] = ddgs.values()
    output_filename = PROTEIN_NAME+"_"+uniprot_accession_name+"_" + pdb_name +"_rosetta"
    print("Saving Rosetta variants in .prism format:\n"+ output_filename+".csv\n")
    stats = "Positions: " + str(len(fasta_seq)) + ", variants: " + str(len(output_df))
    print(uniprot_accession_name, pdb_name, stats)
    output_file_path = join(ddG_output,output_filename)
    output_df.to_csv(output_file_path + ".csv", index = False)
    try:
        data=output_file_path + ".csv"
        csv_to_prism(data,structure_input,chain_id)
    except:
        print('exception')
    print(ddG_run)   
    try: 
        read_slurms(ddG_run,printing=True)
    except:
        print('cant print warnings')
    
if __name__ == '__main__':
    parse_rosetta_ddgs(sys_name=sys.argv[1], chain_id=sys.argv[2], fasta_seq=sys.argv[3], ddG_run=sys.argv[4], ddG_output=sys.argv[5], structure_input=sys.argv[6])
