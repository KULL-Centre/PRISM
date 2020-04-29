"""mp_prepare.py includes several preprocessing and script generating functions to prepare membrane proteins.

Author: Johanna K.S. Tiemann

Date of last major changes: 2020-04-16

"""

# Standard library imports
import json
import urllib.request


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
