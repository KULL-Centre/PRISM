import numpy as np
import scipy

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

        key = protein_seq[int(res_number)-1]+res_number+one_letter

        if key in cartesian_scores:
            cartesian_scores[key].append(dg)

        else:
            cartesian_scores[protein_seq[int(res_number)-1]+res_number
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
        ddgs[mutation] = np.divide((np.mean(dgs_as_floats[mutation])-np.mean(wt_dGs[residue_number])),2.9)
        

    return ddgs


