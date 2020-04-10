#! /usr/bin/env python3
import sys
import os
import numpy as np

def compare_mutfile(fasta_seq, path_to_run_folder,mutation_input=None):
    mutfiles_folder = path_to_run_folder + '/mutfiles/'
    error=False
    path_to_alignment = '{}/uniprot_index_list.txt'.format(output_path)
    alignment = np.loadtxt(path_to_alignment)   
    alignment_dic={}
    for n in enumerate(alignment):
        alignment_dic[int(n[1])] = int(n[0]+1)
         
    
    if mutation_input == None:
        for residue_number in range(1, len(fasta_seq)+1):
            mutfile = open(mutfiles_folder+'mutfile{:0>5}'.format(str(residue_number)), 'r')
            f= mutfile.readlines()
            fasta_seq_list=list(fasta_seq)
            print(residue_number," ",f[2][0]," ",fasta_seq_list[residue_number-1][0])
            if f[2][0] != fasta_seq_list[residue_number-1][0]:
                error=True
                print("ERROR: RESIDUE MISMATCH")
                break
        return(error)
    
    else:
        mut_input = open(mutation_input,'r')
        s = mut_input.readlines()
        listToStr = ''.join([str(elem) for elem in s]) 
        a= listToStr.split()
        print(a)
        ran = list(range(0,len(a),3))
        print("ran",ran)
        for residue_number in ran:
            res = a[residue_number+1]

            residue_number_ros = alignment_dic[res] 
            mutfile = open(mutfiles_folder+'mutfile{:0>5}'.format(str(residue_number_ros)), 'r')
            f= mutfile.readlines()
            fasta_seq_list=list(fasta_seq)
            print(res," ",f[2][0]," ",fasta_seq_list[int(residue_number_ros)-1][0])
            if f[2][0] != fasta_seq_list[int(residue_number_ros)-1][0]:
                error=True
                print("ERROR: RESIDUE MISMATCH")
                break
        return(error)

def pdbxmut(input_mutfiles,resdata):
    
    p=os.listdir(input_mutfiles);errors =0;check3 = False;
    for files in p:
        decs=6
        with open(input_mutfiles+files,'r') as mutfiles:
            mutfile_lines = mutfiles.readlines()
            mutations = (mutfile_lines[2:40:2])
        
        mutfiles.close
        while decs > 2:
            try:
                mut_res=(mutations[1]); m_id=int((mut_res[2:decs])); x,y,z=resdata[m_id]
                decs=0
            except:
                decs+=-1
        
        if x != str(mut_res)[0]:
            print('ERROR',(x,str(mut_res)[0]))
            errors +=1
            check3 = True


    print("Amount of errors = ",errors)
    return(check3,errors)
    
    
    
    
    
    
    
    
    
    
    
