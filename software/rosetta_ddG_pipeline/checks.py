#! /usr/bin/env python3
import sys


def compare_mutfile(fasta_seq, path_to_run_folder,mutation_input=None):
    mutfiles_folder = path_to_run_folder + '/mutfiles/'
    error=False
    
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
        
        for residue_number in ran:
            res = a[residue_number+1]
    
            mutfile = open(mutfiles_folder+'mutfile{:0>5}'.format(str(res)), 'r')
            f= mutfile.readlines()
            fasta_seq_list=list(fasta_seq)
            print(res," ",f[2][0]," ",fasta_seq_list[int(res)-1][0])
            if f[2][0] != fasta_seq_list[int(res)-1][0]:
                error=True
                print("ERROR: RESIDUE MISMATCH")
                break
        return(error)

