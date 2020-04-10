#! /usr/bin/env python3
import sys


def pdb_to_fasta_seq(path_to_pdb):
    
    with open(path_to_pdb, 'r') as pdb_file:
        pdblines = pdb_file.readlines()

    fasta_seq = ''

    chainspec = 'NULL'
    residue_number = 'NULL'

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
                residue_number = line[22:26]

                residue_number = residue_number.lstrip()
                if residue_number != previous_residue_number or chainspec != previous_chainspec:
                        fasta_seq += residue_letter

    return(fasta_seq)

if __name__ == '__main__':
    fasta_sequence = pdb_to_fasta_seq(sys.argv[1])
    print(fasta_sequence)
