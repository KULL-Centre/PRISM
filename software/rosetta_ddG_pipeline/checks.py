#! /usr/bin/env python3
import sys
import os
import numpy as np

def compare_mutfile(fasta_seq, path_to_run_folder,prepare_checking,mutation_input=None):
    """This function checks the created mutfiles and compares them to the fasta_sequence"""
    
    mutfiles_folder = path_to_run_folder +'/'
    error=False
    path_to_alignment = os.path.join(prepare_checking, 'uniprot_index_list.txt')
    alignment = np.loadtxt(path_to_alignment)   
    alignment_dic={}
    for n in enumerate(alignment):
        alignment_dic[int(n[1])] = int(n[0]+1)
         
    
    if mutation_input == None:
        for residue_number in range(1, len(fasta_seq)+1):
            mutfile = open(os.path.join(mutfiles_folder,f'mutfile{str(residue_number):0>5}'), 'r')
            f= mutfile.readlines()
            fasta_seq_list=list(fasta_seq)
            if f[2][0] != fasta_seq_list[residue_number-1][0]:
                error=True
                print("ERROR: RESIDUE MISMATCH")
                break
        return(error)
    
    else:
        a = []
        with open(mutation_input,'r') as mut_input:
            for line in mut_input:
                lines = line.split()
                if len(lines) ==3:
                    a.append(lines)
                else:
                    key=[]
                    wt=[]
                    val=[]
                    for multi_muts in range(int(len(lines)/3)):
                        wt_single = lines[0+(3*multi_muts)]
                        var_single = lines[2+(3*multi_muts)]
                        if not wt_single in var_single:
                            var_single += wt_single
                        var_single = ''.join(set(var_single))
                        wt.append(wt_single)
                        key.append(lines[1+(3*multi_muts)])
                        val.append(var_single)
                    a.append(["_".join(wt), "_".join(key), "_".join(val)])
        a = [elem for ind in a for elem in ind]
        ran = list(range(0,len(a),3))
        for residue_number in ran:
            res = a[residue_number+1]
            if len(res.split("_")) == 1:
                residue_number_ros = alignment_dic[int(res)]
                mutfile = open(os.path.join(mutfiles_folder,f'mutfile{str(residue_number_ros):0>5}'), 'r')
            else:

                mutfile = open(os.path.join(mutfiles_folder,f'mutfile{res}'), 'r')
            
            f= mutfile.readlines()
            fasta_seq_list=list(fasta_seq)
            for indi, resi in enumerate(res.split("_")):
                if f[2+indi].split()[0] != fasta_seq_list[int(resi)-1][0]:
                    print(f[2+indi].split()[0],fasta_seq_list[int(resi)-1][0])
                    error=True
                    print("ERROR: RESIDUE MISMATCH")
                    break
        return(error)

def pdbxmut(input_mutfiles,structure_dic):
    """This function checks the created mutfiles and compares them to the structure directory"""
    
    resdata = structure_dic["resdata"]
    p=os.listdir(input_mutfiles);errors =0;check3 = False;
    for files in p:
        decs=6
        with open(os.path.join(input_mutfiles,files),'r') as mutfiles:
            mutfile_lines = mutfiles.readlines()
            mutations = (mutfile_lines[2:40:2])
        while decs > 2:
            try:
                mut_res=(mutations[1])
                m_id=int((mut_res[2:decs]))
                x,y,z=resdata[m_id]
                decs=0
            except:
                decs+=-1
        
        if x != str(mut_res)[0]:
            print('ERROR',(x,str(mut_res)[0]))
            errors +=1
            check3 = True


    print("Amount of errors = ",errors)
    return(check3,errors)
    
    
    
    
    
    
    
    
    
    
    
