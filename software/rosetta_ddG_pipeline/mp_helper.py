"""mp_prepare.py includes several preprocessing and script generating functions to prepare membrane proteins.

Author: Johanna K.S. Tiemann

Date of last major changes: 2020-04-16

"""

# Standard library imports
import json
import urllib.request
from ptm_dict import modres


d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def extract_from_opm(pdb_id):
    '''
    Extracts from OPM the pdb-structure as a string

    Examples
    --------
    >>> pdb_id = "3gp6"
    >>> extract_from_opm( pdb_id )
    '''

    url = f"https://opm-assets.storage.googleapis.com/pdb/{pdb_id.lower()}.pdb"
    r = urllib.request.urlopen(url)
    data = r.read()
    data = data.decode('utf-8')
    return data

def getnums(seq, seqnum):
    i = 0
    numseq = []
    for elem in seq:
        if elem != '-':
            numseq.append(seqnum[i])
            i = i + 1
        else:
            numseq.append(None)
    return numseq

def get_seq(file_pointer, isfile=True):
    global modres
    i = 0
    if isfile:
        fp = open(file_pointer)
    else:
        fp = file_pointer.split('\n')
    for line in fp:
        if line.startswith('ATOM'):
            resseq = int(line[22:26])
            if line[12:16].strip() == 'CA':
                if resseq > i:
                    i = resseq
    seq = ["-"]*i
    seq2 = []
    if isfile:
        fp.close()
        fp = open(file_pointer)
    else:
        fp = file_pointer.split('\n')
    for line in fp:
        if line.startswith('ATOM'):
            resseq = int(line[22:26])
            resname = line[17:20].strip()
            if line[12:16].strip() == 'CA':
                if line[16] in [' ', 'A']:
                    # seq[resseq-1] = d3to1[modres[resname]]
                    seq[resseq-1] = d3to1[resname]
                    seq2.append(resseq)
    return "".join(seq), seq2, i

def get_res_in_mem(file_pointer, isfile=False):
    all_z=[]
    if isfile:
        fp = open(file_pointer)
    else:
        fp = file_pointer.split('\n')
    for line in fp:
        if line.startswith('HETATM') and line[17:20]=='DUM':
            z = float(line[46:54])
            all_z.append(z)
    all_z = sorted(list(set(all_z)))

    residues_in_membrane = []
    ca_atoms_in_membrane = []
    if isfile:
        fp.close()
        fp = open(file_pointer)
    else:
        fp = file_pointer.split('\n')
    for line in fp:
        if line.startswith('ATOM'):
            z = float(line[46:54])
            resseq = int(line[22:26])
            atomnum = int(line[6:11])
            if (z >= all_z[0]) and (z <= all_z[1]):
                residues_in_membrane.append(resseq)
                if line[12:16].strip() == 'CA':
                    ca_atoms_in_membrane.append(atomnum)
    residues_in_membrane = sorted(list(set(residues_in_membrane)))

    return residues_in_membrane

def get_res_in_all(file_pointer, isfile=False):
    all_z=[]
    if isfile:
        fp = open(file_pointer)
    else:
        fp = file_pointer.split('\n')

    residues_in_all = []
    ca_atoms_in_all = []
    if isfile:
        fp.close()
        fp = open(file_pointer)
    else:
        fp = file_pointer.split('\n')
    for line in fp:
        if line.startswith('ATOM'):
            resseq = int(line[22:26])
            atomnum = int(line[6:11])
            residues_in_all.append(resseq)
            if line[12:16].strip() == 'CA':
                ca_atoms_in_all.append(atomnum)
    residues_in_all = sorted(list(set(residues_in_all)))

    return residues_in_all
