"""prims_rosetta_parser.py contains functions to convert prism to mut & ddG to prims

Author: Johanna K.S. Tiemann
Date of last major changes: 2020-05-01

"""

# Standard library imports
import logging as logger
import re
import sys


# Local application imports
import rosetta_paths
sys.path.insert(1, rosetta_paths.prims_parser)
from PrismData import PrismParser


def prism_to_mut(primsfile, mutfile):
    # extracts the mutations from prims
    logger.info(
        'Extract information from prismfile and converti it into dic & mutfile')
    parser = PrismParser()
    data = parser.read(primsfile)
    data_frame1 = data.dataframe
    mut_dic = {}
    with open(mutfile, 'w') as fp:
        for resid in data_frame1["resi"].explode().unique():
            data_frame2 = data.get_var_at_pos(resid)
            native = data_frame2['aa_ref'].explode().unique()[0]
            variants = ''.join([''.join(map(str, l))
                                for l in data_frame2['aa_var']]) + native
            regex = re.compile('[^a-zA-Z]')
            final_variants = ''.join(set(regex.sub('', variants)))
            mut_dic[resid] = final_variants
            fp.write(f'{native} {resid} {final_variants} \n')
    return mut_dic
