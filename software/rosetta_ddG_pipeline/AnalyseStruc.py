import sys
import os
import Bio
from Bio.PDB.PDBParser import PDBParser
parser = PDBParser(PERMISSIVE=1)
from Bio.PDB import *
import json

def get_structure_parameters(outpath,structure_id,chain_id):
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
    model = structure[0]
   
    for chain in structure[0]:
        for residue in chain:
            if residue.get_id()[0] == " ":
                count += 1
                try:
                    residue_letter = aminocodes[residue.get_resname()]
                    resdata[count] = residue_letter,residue.get_id()[1],chain.get_id()
                except: 
                    residue_letter = str(residue.get_resname())
                    exceptions += 1
                    resdata[count] = residue_letter,str(residue.get_id()[1]),chain.get_id()
    print("Special residues in structure = ",exceptions)
    
    
    for chain in model:
        sequence = []
        for res in range(1,len(resdata)+1):
    
            if str(resdata[res][2]) == str(chain)[-2:-1]:
                sequence.append(resdata[res][0])
        sequence_chain = ''.join([str(elem) for elem in sequence]) 
        strucdata[str(chain)[-2:-1]]= sequence_chain
        
    structure_dic = {"resdata": resdata, "strucdata": strucdata}
    with open(outpath +"/structure_{}.txt".format(name),'w') as strucfile:

        strucfile.write('#Structure features \n')
        strucfile_line = '{}' + '\t {}'+ '\t {}' +  '\t {}' + '\n'
        strucfile.write(strucfile_line.format("Ros-nr","Pdb-nr","AA","Chain"))
        for rosnumber in resdata:    
            AA, pdbnumber, chainid = resdata[rosnumber]
            strucfile.write(strucfile_line.format(str(rosnumber),str(pdbnumber),str(AA),str(chainid)))
    with open(outpath +f"/structure_{name}.json", 'w') as outfile:
        json.dump(structure_dic, outfile)        
    return(structure_dic)

if __name__ == '__main__':
    structure_id = sys.argv[2]
    chain_id = sys.argv[3]
    outpath = sys.argv[1]
    get_structure_parameters(outpath,structure_id,chain_id)