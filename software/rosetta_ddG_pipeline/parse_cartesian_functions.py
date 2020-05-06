import numpy as np
import scipy
import os


def rosetta_cartesian_read(pathtofile, protein_seq='abcd'):
    score_file = open(pathtofile, "r")
    score_data = score_file.readlines()
    score_file.close()

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
        description = score_fields[2]
        three_letter_code = description[-4:-1]
        one_letter = aminocodes[three_letter_code]
        res_number = description[4:-4]
        dg = float(score_fields[3])

        key = protein_seq[int(res_number) - 1] + res_number + one_letter

        if key in cartesian_scores:
            cartesian_scores[key].append(dg)

        else:
            cartesian_scores[protein_seq[int(res_number) - 1] + res_number
                             + one_letter] = [dg]

    return cartesian_scores


def ddgs_from_dg(dictionary_of_dGs):
    wt_dGs = {}
    for entry in dictionary_of_dGs:

        if entry[0] == entry[-1]:

            residue_number = entry[1:-1]
            for item in dictionary_of_dGs[entry]:
                if residue_number in wt_dGs:
                    wt_dGs[residue_number].append(float(item))

                else:
                    wt_dGs[residue_number] = [float(item)]

    ddgs = {}
    dgs_as_floats = {}
    for mutation in dictionary_of_dGs:
        dgs_as_floats[mutation] = []
        for value in dictionary_of_dGs[mutation]:
            dgs_as_floats[mutation].append(float(value))

    for mutation in dictionary_of_dGs:

        residue_number = mutation[1:-1]
        ddgs[mutation] = np.divide(
            (np.mean(dgs_as_floats[mutation]) - np.mean(wt_dGs[residue_number])), 2.9)

    return ddgs


def postprocess_rosetta_ddg_prism_copy(folder, output_name='ddG.out', sys_name='', uniprot='', version=1, prims_nr='XXX'):
    # The ddg_file should only contain the data, looking like this:
    # M1T,-0.52452 # first value=variant, second=mean([var1-WT1,var2-WT2, ...]) - comma separated
    # M1Y,0.2352,0.2342,.... # it may contain more values like the var-mut of
    # each run
    ddg_file = os.path.join(folder.ddG_run, output_name)
    prims_file = os.path.join(folder.ddG_output, f'prims_rosetta_{prims_nr}_{sys_name}.txt')
    with open(os.path.join(folder.prepare_checking, 'fasta_file.fasta'), 'r') as fp:
        fp.readline()
        sequence = fp.readline().strip()
    rosetta_to_prism(ddg_file, prims_file, sequence, rosetta_info=None,
                     version=version, uniprot=uniprot, sys_name=sys_name)
    create_copy(os.path.join(folder.ddG_input, 'input.pdb'), folder.output, name=f'{sys_name}_final.pdb')
    create_copy(os.path.join(folder.prepare_checking, 'fasta_file.fasta'), folder.output, name=f'{sys_name}_seq.fasta')
    create_copy(prims_file, folder.output)
