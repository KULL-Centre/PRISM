import sys
import os
import Bio
from Bio.PDB.PDBParser import PDBParser
parser = PDBParser(PERMISSIVE=1)
from Bio.PDB import *
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import json

def get_structure_parameters(outpath,structure_id,chain_id,printing=True):
    """This script creates the .json file containing all the structure and protein information. This is all extracted from the pdb file, so if the information is not in the pdb, it will not be in file.
    
    RESDATA: Information about pdb numbering, rosetta numbering and chain
    STRUCDATA: Information about sequence for the chains. Both coordinate seq and reference seq
    DBREF: Information about protein name, starting and end pos, uniprot ID and more
    ALIGNMENT: Contains alignments between coordinate sequence and reference sequence 
    """
    name = structure_id.split("/")
    name = name[-1].split(".")[-2] 
    structure = parser.get_structure(name, structure_id)
    
    aminocodes = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y"
    }
    resdata={};strucdata={};count=0;exceptions = 0
    resdata_reverse={}
    resdata_reverse2={}
    model = structure[0]
    
    #Making resdata
    true_count = 0
    for chain in structure[0]:
        print(chain)
        for residue in chain:
            if residue.get_id()[0] != " " and residue.resname in modres:
                print(residue.get_id())
            #if residue.get_id()[0] == " " or (residue.get_id()[0] != " " and residue.resname in modres):
            if residue.get_id()[0] == " ":
                print(residue.get_id())
                count += 1
                restore_res_id = f"{residue.get_id()[1]}{residue.get_id()[2]}".strip()
                try:
                    residue_letter = aminocodes[residue.get_resname()]
                    resdata_reverse[restore_res_id] = count
                    if chain.get_id() in [x for x in chain_id]:
                        true_count += 1
                        resdata[true_count] = residue_letter,restore_res_id,chain.get_id()
                        resdata_reverse2[f'{restore_res_id};{chain.get_id()}'] = true_count
                except: 
                    residue_letter = str(residue.get_resname())
                    exceptions += 1
                    resdata_reverse[restore_res_id] = count
                    if chain.get_id() in [x for x in chain_id]:
                        true_count += 1
                        resdata[true_count] = residue_letter,str(restore_res_id),chain.get_id()
                        resdata_reverse2[f'{restore_res_id};{chain.get_id()}'] = true_count
    #Counts special residues such as DNA 
    if printing == True:
        print("Special residues in structure = ",exceptions)
    
    #Making strucdata
    for chain in model:
        sequence = []
        for res in range(1,len(resdata)+1):
    
            if str(resdata[res][2]) == str(chain)[-2:-1]:
                sequence.append(resdata[res][0])
        sequence_chain = ''.join([str(elem) for elem in sequence]) 
        strucdata[str(chain)[-2:-1]]= [sequence_chain, '', '']
  
    for record in SeqIO.parse(structure_id, "pdb-seqres"):
        seq=str(record.seq)
        if str(record.annotations["chain"]) in strucdata.keys():
            strucdata[str(record.annotations["chain"])] = [strucdata[str(record.annotations["chain"])][0],seq,''] 
        
    with open(structure_id, 'r') as pdblines:
        fasta_seq_full = ''
        fasta_seq_chain = ''
        chainspec = 'NULL'
        residue_number = 0
        #Relative sequence
        for line in pdblines:
            if len(line) > 1:
                line_fields = line.split()
                if line_fields[0] == 'ATOM' and line[17:20] in aminocodes :
                    entry_type = line_fields[0]
                    atom_number = line_fields[1]
                    atom_type = line_fields[2]
                    residue_type = line[17:20]
                    residue_letter = aminocodes[residue_type]
                    previous_chainspec = chainspec
                    chainspec = line[21]              
                    previous_residue_number = residue_number
                    residue_number = line[22:26]# maybe keep at 22:26
                    residue_number = residue_number.lstrip()
                    
                    if (previous_chainspec == 'NULL' or (residue_number != previous_residue_number and chainspec == previous_chainspec)):
                        if ((int(residue_number) - int(previous_residue_number)) == 1 )or residue_number ==0:
                            fasta_seq_chain += residue_letter     
                        if ((int(residue_number) - int(previous_residue_number)) != 1) and residue_number != 0:
                            fasta_seq_chain += "-" * (int(residue_number) - int(previous_residue_number) -1 )         
                            fasta_seq_chain += residue_letter
                    elif chainspec != previous_chainspec:
                        strucdata[previous_chainspec] = [strucdata[previous_chainspec][0], strucdata[previous_chainspec][1],fasta_seq_chain]
                        fasta_seq_chain = residue_letter
        strucdata[previous_chainspec] = [strucdata[previous_chainspec][0], strucdata[previous_chainspec][1],fasta_seq_chain]

    #Making DBREF
    with open(structure_id, 'r') as pdb_file:
        pdblines = pdb_file.readlines()
    ref_count=0
    dbref={}
    for line in pdblines:
        if len(line) > 1:
            line_fields = line.split()
            if line_fields[0] == 'DBREF' :
                if str(line_fields[2]) in strucdata.keys():
                    refs=line_fields
                    ref_count+=1
                    dbref[ref_count] = refs
    
    #Making alignments using pairwise2
    align ={}
    for chain in strucdata:
        
        seq1=strucdata[str(chain)][0]
        seq2=strucdata[str(chain)][1]
        alignments = pairwise2.align.globalxx(seq1, seq2)
        align[str(chain)] = alignments

    #Saving in a dictionary
    structure_dic = {"resdata": resdata, "resdata_reverse":resdata_reverse, "resdata_reverse2":resdata_reverse2, "strucdata": strucdata, "DBREF":dbref, "alignment":align}

    ##Creating overview txt file
    #with open(outpath +"/structure_{}.txt".format(name),'w') as strucfile:
#
    #    strucfile.write('#Structure features \n')
    #    strucfile_line = '{}' + '\t {}'+ '\t {}' +  '\t {}' + '\n'
    #    strucfile.write(strucfile_line.format("Ros-nr","Pdb-nr","AA","Chain"))
    #    for rosnumber in resdata:    
    #        AA, pdbnumber, chainid = resdata[rosnumber]
    #        strucfile.write(strucfile_line.format(str(rosnumber),str(pdbnumber),str(AA),str(chainid)))
            
    #Making json file        
    with open(outpath +f"/structure_{name}.json", 'w') as outfile:
        json.dump(structure_dic, outfile,indent=0)        
    return(structure_dic)

if __name__ == '__main__':
    structure_id = sys.argv[2]
    outpath = sys.argv[1]
    get_structure_parameters(outpath,structure_id)
