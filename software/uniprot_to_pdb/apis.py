"""apis have functions related to using apis, such as SIFTS and retrieving pdb information for rcsb.

Author: Anders Frederiksen

Date of last major changes: 2020-10-09

"""

# Standard library imports
import io
import requests
import json
import sys



def req_sifts(uniprot_accession):
    """This function retrieves the data froms SIFTS for uniprotID"""
    
    ## Request SIFTS
    requestURL = f'https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{uniprot_accession}'
    r = requests.get(requestURL)
    if not r.ok:
        r.raise_for_status()
        sys.exit()
        
    ## Make response dictionary
    response_dict = json.loads(r.text)
    
    return(response_dict)
    
    
def get_pdb_coord_seq(structure_id, chain_id='A'):
    """Get the specific coordinate sequence from a pdb file"""
    
    ## Get pdb from rcsb
    requestURL = f'https://files.rcsb.org/view/{structure_id}.pdb'
    r = requests.get(requestURL)
    chain_id = chain_id
    
    ## Get the coordinate sequence from the pdb
    sequence = pdb_to_fasta_seq(r, chain_id=chain_id)
    sequence = sequence[1]
    return(sequence)


def pdb_to_fasta_seq(r, chain_id='NULL'):
    """Converts pdb file residue coordinate sequence to a fasta sequence"""
    
    ## Create variable and dictionaries
    seq = io.StringIO(r.text)
    pdblines = seq.readlines()
    fasta_seq_full = ''
    fasta_seq_chain = ''
    chainspec = 'NULL'
    residue_number = 0
    chain_id = chain_id
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
    
    ## Readlines in pdb file
    for line in pdblines:
        if len(line) > 1:
            line_fields = line.split()
            if line_fields[0] == 'ENDMDL' and len(fasta_seq_chain) > 10:
                break
            if line_fields[0] == 'ATOM' and line[17:20] in aminocodes and str(line[21]) == str(chain_id):
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
                if (residue_number != previous_residue_number or chainspec != previous_chainspec):
                    if ((int(residue_number) - int(previous_residue_number)) == 1) and int(previous_residue_number) != 0:
                        fasta_seq_chain += residue_letter
                    if ((int(residue_number) - int(previous_residue_number)) != 1) or int(previous_residue_number) == 0:
                        fasta_seq_chain += "-" * \
                            ((int(residue_number) - int(previous_residue_number)) - 1)
                        fasta_seq_chain += residue_letter

    return(fasta_seq_full, fasta_seq_chain)


def read_fasta(uniprot_accession):
    """"Read a fasta file. Input uniprot accession. Outputs sequence"""
    
    ## Request fasta from Uniprot
    requestURL = f'https://www.uniprot.org/uniprot/{uniprot_accession}.fasta'
    r = requests.get(requestURL)
    fastau_file = io.StringIO(r.text)
    fastau_lines = fastau_file.readlines()
    fastau_file.close()
    uniprot_seq = ''
    ## readlines
    for line in fastau_lines[1:]:
        uniprot_seq = uniprot_seq + line.strip()
    return(uniprot_seq)    