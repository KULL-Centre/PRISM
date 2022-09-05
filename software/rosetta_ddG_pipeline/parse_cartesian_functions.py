import numpy as np
import scipy
import os
import json

def rosetta_cartesian_read(pathtofile, protein_seq='abcd', struc_dat=''):
    """This script takes the individual score files in the run folder and outputs a dictionary of dGs"""
    
    score_file = open(pathtofile, "r")
    score_data = score_file.readlines()
    score_file.close()

    if struc_dat!='':
        with open(struc_dat) as json_file:
            strucdata = json.load(json_file)

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

    cartesian_scores = {}
    for line in score_data:
        score_fields = line.split()
        descriptions = score_fields[2][4:-1]
        key = []
        for description in descriptions.split('_'):
            three_letter_code = description[-3:]
            one_letter = aminocodes[three_letter_code]
            res_number = description[:-3]
            dg = float(score_fields[3])
            if struc_dat=='':
                key.append(protein_seq[int(res_number) - 1] + res_number + one_letter)
            else:
                key.append(strucdata['resdata'][str(res_number)][0] + res_number + one_letter)
        key = ":".join(key)
        if key in cartesian_scores:
            cartesian_scores[key].append(dg)

        else:
            cartesian_scores[key] = [dg]
    return cartesian_scores


def ddgs_from_dg(dictionary_of_dGs, scale_factor=2.9):
    """This scripts take a dictionaru of dGs and first creates a dictionary of WT dGS and then substract the variants dGs. After all substractions the ddG score is divided by 2.9 to convert to kcal/mol """
    
    #Creating dictionary of WT dGs
    wt_dGs = {}
    for entrys in dictionary_of_dGs:
        residue_numbers = []
        for entry in entrys.split(':'):
            if entry[0] == entry[-1]:
                residue_numbers.append(entry[1:-1])
            else:
                residue_numbers = []
        if len(residue_numbers) > 0:
            residue_number = ":".join(residue_numbers)
            for item in dictionary_of_dGs[entrys]:
                if residue_number in wt_dGs:
                    wt_dGs[residue_number].append(float(item))
                else:
                    wt_dGs[residue_number] = [float(item)]
    
    #Creating dictionary of variant dGs
    dgs_as_floats = {}
    for mutation in dictionary_of_dGs:
        dgs_as_floats[mutation] = []
        for value in dictionary_of_dGs[mutation]:
            dgs_as_floats[mutation].append(float(value))

    ddgs = {}
    ddgs_array = []
    # (variant - WT) / 2.9
    for mutation in dictionary_of_dGs:
        residue_numbers = []
        for mutations in mutation.split(':'):
            residue_numbers.append(mutations[1:-1])
        residue_number = ":".join(residue_numbers)
        ddg = []
        for indi in range(len(dgs_as_floats[mutation])):
            ddg.append((dgs_as_floats[mutation][indi]-np.mean(wt_dGs[residue_number]))/scale_factor)
        ddgs_array.append([mutation, np.mean(ddg), np.std(ddg)])
        ddgs[mutation] = np.divide((np.mean(dgs_as_floats[mutation]) - np.mean(wt_dGs[residue_number])), scale_factor)

    return ddgs, ddgs_array


def postprocess_rosetta_ddg_prism_copy(folder, output_name='ddG.out', sys_name='', uniprot='', version=1, prism_nr='XXX', scale=2.9):
    # The ddg_file should only contain the data, looking like this:
    # M1T,-0.52452 # first value=variant, second=mean([var1-WT1,var2-WT2, ...]) - comma separated
    # M1Y,0.2352,0.2342,.... # it may contain more values like the var-mut of
    # each run
    ddg_file = os.path.join(folder.ddG_run, output_name)
    prism_file = os.path.join(folder.ddG_output, f'prism_rosetta_{prism_nr}_{sys_name}.txt')
    with open(os.path.join(folder.prepare_checking, 'fasta_file.fasta'), 'r') as fp:
        fp.readline()
        sequence = fp.readline().strip()
    rosetta_to_prism(ddg_file, prism_file, sequence, rosetta_info=None,
                     version=version, uniprot=uniprot, sys_name=sys_name, scale=scale)
    create_copy(os.path.join(folder.ddG_input, 'input.pdb'), folder.output, name=f'{sys_name}_final.pdb')
    create_copy(os.path.join(folder.prepare_checking, 'fasta_file.fasta'), folder.output, name=f'{sys_name}_seq.fasta')
    create_copy(prism_file, folder.output)
