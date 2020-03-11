#! /usr/bin/env python3
import sys


def pdb_to_fasta_seq(path_to_pdb):
    '''This function reads a pdb file and returns the fasta sequence.
    This function ignores chain specifications and simply returns a string of residues'''

    with open(path_to_pdb, 'r') as pdb_file:
        pdblines = pdb_file.readlines()

    # this string will hold the fasta sequence
    fasta_seq = ''
    # these variables hold some information about the atom
    # we need to keep track of the chain specification, so we are not
    # tricked in case a second chain has a residue that is identical in type and number
    # to what the previous chained ended with
    chainspec = 'NULL'
    residue_number = 'NULL'

    # this dictionary is just for translating from 3-letter codes
    # to 1-letter code
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
        # this just removes empty lines
        if len(line) > 1:
            line_fields = line.split()
            if line_fields[0] == 'ATOM':
                entry_type = line_fields[0]
                atom_number = line_fields[1]
                atom_type = line_fields[2]
                # pdbs have AA three letter codes, let us
                # translate that to single letter
                residue_type = line[17:20]
                residue_letter = aminocodes[residue_type]

                # the chainspecification always happens at the
                # 22nd character of the line. It is not always
                # seperated by spaces. And some pdbs use the
                # legacy <space> instead of a letter, if there
                # is only one chain
                previous_chainspec = chainspec
                chainspec = line[21]
                # since there is not always room for a space
                # in the pdb, we have to get a little creative
                # with how we select the residue index
                previous_residue_number = residue_number
                residue_number = line[22:26]
                # and remove the whitespace
                residue_number = residue_number.lstrip()
                # check to see if the residue is a new one
                if residue_number != previous_residue_number or chainspec != previous_chainspec:
                        fasta_seq += residue_letter

    return(fasta_seq)


# this last part is just so you can call it from the shell
if __name__ == '__main__':
    fasta_sequence = pdb_to_fasta_seq(sys.argv[1])
    print(fasta_sequence)
