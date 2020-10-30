"""msa includes functions for multiple sequence alignment and parsing of the alignments.

Author: Anders Frederiksen

Date of last major changes: 2020-10-09

"""

# Standard library imports
import io
from os.path import join

# Third party imports
import numpy as np
import pandas as pd
import Bio
from Bio.PDB import *
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline

# Local settings
clustal='/groups/sbinlab/haagenb/software/clustal/bin/clustalo'



def run_alignment(uniprot_seq, df, output_location):
    """This function runs the alignment based on the input of the dataframe"""
    
    #Writing sequence file for Clustal
    df_coords = df
    seqs = [uniprot_seq]
    seqs_names = ['Uniprot']
    coord_data = {}
    for n in range(0, len(df_coords)):        
        try:
            if len(df_coords['coord_sequence'][n]) > 10:
                seqs = np.hstack((seqs, df_coords['coord_sequence'][n]))
                seqs_names = np.hstack((seqs_names, df_coords['pdb'][n]+'_'+df_coords['chainid'][n]))
                coord_data[df_coords['pdb'][n]+'-'+df_coords['chainid'][n]] = [
                    df_coords['coord_sequence'][n]]
        except:
            if df_coords['coord_sequence'][n] > 10:
                seqs = np.hstack((seqs, df_coords['coord_sequence'][n]))
                seqs_names = np.hstack((seqs_names, df_coords['pdb'][n]+'_'+df_coords['chainid'][n]))
                coord_data[df_coords['pdb'][n]+'-'+df_coords['chainid'][n]] = [
                    df_coords['coord_sequence'][n]]
    coord_data['Uniprot'] = [uniprot_seq]

    # Write to file
    seqss = io.StringIO()
    num = 0
    for n, m in zip(seqs, seqs_names):
        if num == 0:
            string = f'>{m}\n'+n+'\n'
            seqss.write(f'{string}')
            num += 1
        else:
            string = f'>{m}\n'+n+'\n'
            seqss.write(f'{string}')
            num += 1           
    x = seqss.getvalue()
          
    infile = join(output_location,'temp1')
    outfile= join(output_location,'clustal.aln')
    
    with open(infile, 'w') as f:
        f.write(x)
    
    # Run ClustalO
    cline = ClustalwCommandline(clustal, infile=infile, outfile=outfile)
    cline()

    return(outfile)


def parse_alignment(outfile):
    """This function parse the output of the alignment"""
    
    ## Load alignment file
    file_output=outfile   
    alignment = AlignIO.read(open(file_output), "clustal")
    alignment_output_dic={}
    new_alignment_dic={}
    
    for record in alignment:
        alignment_output_dic[str(record.id)]=[str(record.seq)]
    
    for m in alignment_output_dic:
        sequence=''
        for n in range(len(alignment_output_dic['Uniprot'][0])):
            if alignment_output_dic[m][0][n] != alignment_output_dic['Uniprot'][0][n] and alignment_output_dic['Uniprot'][0][n] != '-' :
                sequence=sequence+'-'
            if alignment_output_dic[m][0][n] == alignment_output_dic['Uniprot'][0][n]:
                sequence=sequence+(alignment_output_dic[m][0][n])
        new_alignment_dic[m]=sequence
        
    num = 0
    seq_numbers = []
    for n in new_alignment_dic['Uniprot']:
        if n != '-':
            seq_numbers = np.hstack((seq_numbers, int(1)))
            num += 1
        if n == '-':
            seq_numbers = np.hstack((seq_numbers, int(0)))
            num += 1
            
    new_seq = ''
    nums = 0
    aligned_data={}
    arr=[]
    for m in new_alignment_dic.keys():   
        if m != 'Uniprot':
            arr=[]
            for p,q in zip(new_alignment_dic[m],seq_numbers):
                x=p*int(q)
                arr=np.hstack((arr,x))
        aligned_data[m]=arr

    binary_data = {}
    dat = []
    for m in aligned_data.keys():
        dat = []
        seq_temp = aligned_data[m]
        for n in seq_temp:
            if n == '-':
                dat = np.hstack((dat, 0))
            else:
                dat = np.hstack((dat, 1))
        binary_data[m] = dat

    binary_data
    sum_list = np.zeros(len(binary_data.values()))
    num = 1
    for n in binary_data.keys():
        if num == 1:
            temp = np.zeros(len(binary_data[n]))
        if n != 'Uniprot':
            x=binary_data[n]
            if x == 'NaN':
                x=0
            num += 1
            temp = x+temp
    return(binary_data)

def bin_to_df(df_coords,binary_data):
    """This function incorporates the result in the DataFrame"""
    
    ## Load binary data into DataFrame
    for m, p in zip(df_coords['pdb'], range(len(df_coords['pdb']))):
        #if 1==1:
        try:
            if m == df_coords['pdb'][p]:
                #print(m,'-',df_coords['chainid'][p])
                #print(binary_data[str(m+'_'+df_coords['chainid'][p])])
                df_coords['binary_sequence'][p] = str(binary_data[str(m+'_'+df_coords['chainid'][p])])
               # df_coords['length'][p] = sum(binary_data[str(m+'_'+df_coords['chainid'][p])])
        except:
            print('error adding binary data to df = ',m)
            continue
    return(df_coords)