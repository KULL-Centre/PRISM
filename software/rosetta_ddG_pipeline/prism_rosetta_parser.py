"""prims_rosetta_parser.py contains functions to convert prism to mut & ddG to prims

Author: Johanna K.S. Tiemann
Date of last major changes: 2020-05-01

"""

# Standard library imports
from datetime import datetime
import logging as logger
import re
import sys


# Third party imports
import pandas as pd


# Local application imports
import rosetta_paths
sys.path.insert(1, rosetta_paths.prims_parser)
from PrismData import PrismParser, VariantData


def read_from_prism(primsfile):
    logger.info('Reads the prism file')
    parser = PrismParser()
    dataframe = parser.read(primsfile).dataframe
    meta_data = parser.read_header(primsfile)
    return meta_data, dataframe


def prism_to_mut(primsfile, mutfile):
    # extracts the mutations from prims
    logger.info(
        'Extract information from prismfile and converti it into dic & mutfile')
    parser = PrismParser()
    data = parser.read(primsfile)
    data_frame1 = data.dataframe
    with open(mutfile, 'w') as fp:
        searchlist = data_frame1["resi"].explode().unique()
        cleanedList = [x for x in searchlist if str(x) != 'nan']
        for resid in cleanedList:
            data_frame2 = data.get_var_at_pos(resid)
            native = data_frame2['aa_ref'].explode().unique()[0]
            variants = ''.join([''.join(map(str, l))
                                for l in data_frame2['aa_var']]) + native
            regex = re.compile('[^a-zA-Z]')
            final_variants = ''.join(set(regex.sub('', variants)))
            fp.write(f'{native} {resid} {final_variants} \n')
    return 


def rosetta_to_prism(ddg_file, prism_file, sequence, rosetta_info=None, version=1, sys_name='', first_residue_number=1):
    # create prism file with rosetta values
    logger.info('Create prism file with rosetta ddG values')
    variant = []
    norm_ddG_value = []
    std_ddG_value = []
    ddG_value = []
    with open(ddg_file, 'r') as fp:
        for line in fp:
            split_line = line.split(',')
            variant.append(split_line[0])
            norm_ddG_value.append(split_line[1])
            std_ddG_value.append(split_line[2].strip())
     #       ddG_value.append(split_line[1])

    data = {
        'variant': pd.Series(variant),
        'norm_ddG': pd.Series(norm_ddG_value),
        'std_ddG': pd.Series(std_ddG_value),
        #      'ddG': pd.Series(ddG_value),
        'n_mut': 1,  # pd.Series([1 for x in range(len(variant))]),
        'aa_ref': [[seg[0]] for seg in variant],
        "resi": [[seg[1:-1]] for seg in variant],
        "aa_var": [[seg[-1]] for seg in variant],
    }
    dataframeset = pd.DataFrame(data)

    if rosetta_info == None:
        rosetta_info = {
            "version": "XXX",
        }


    metadata = {
        # new version always when the data is edited or major changes are made
        # to metadata
        "version": version,
        # extracted from uniprot
        "protein": {
            "name": sys_name,
            "organism": "",
            "uniprot": "",
            "sequence": sequence,
        },
        "rosetta": rosetta_info,
        "variants": {
            #      "number": df["pos"].count(),
            #      "coverage": df["pos"].nunique() / len(seq[1][1]),
            "number": data["variant"].count(),
            "coverage": len(set(sum(data["resi"], []))) / len(sequence),
            "width": "single",
            #        "depth":,
        },
        # columns is dependent on the data. different conditions go into
        # different PRISM_files
        "columns": {
            "norm_ddG": "mean Rosetta ddG values normalized to WT",
            "std_ddG": "std Rosetta ddG values normalized to WT",
            #            "ddG": "Rosetta ddG values",
        },
    }
    if first_residue_number != 1:
        metadata['protein']['first_residue_number'] = first_residue_number

    comment = [
        f"version {version} - {datetime.date(datetime.now())} - ddG pipeline",
    ]

    variant_dataset = VariantData(metadata, dataframeset)
    parser = PrismParser()
    parser.write(prism_file, variant_dataset, comment_lines=comment)


def write_prism(metadata, dataframe, prism_file, comment=''):
    variant_dataset = VariantData(metadata, dataframe)
    parser = PrismParser()
    parser.write(prism_file, variant_dataset, comment_lines=comment)
