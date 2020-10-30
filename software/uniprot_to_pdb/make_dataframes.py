"""make_dataframes includes functions to generate dataframes.

Author: Anders Frederiksen

Date of last major changes: 2020-10-09

"""

# Standard library imports 
import warnings

# Third party imports
import pandas as pd
import numpy as np

# Local application imports
from apis import get_pdb_coord_seq



def make_df(data, uniprot_seq):
    """This function makes the dataframe from sifts"""
    
    struc_info = ''
    num = 0
    struc = ''
    key=list(data.keys())[0]
    for m in range(len(data[key])):
        pdb = data[key][m]
        struc = [pdb['pdb_id'], pdb['chain_id'], pdb['start'], pdb['end'],
                 pdb['unp_start'], pdb['unp_end'], pdb['experimental_method'], pdb['coverage']]
        try:
            struc = np.append(
                struc, pdb['resolution'])

        except:
            struc = np.append(struc, 'NAN')
        if num == 0:
            struc_info = struc
        else:
            struc_info = np.vstack((struc_info, struc))
        num += 1

    df_columns = ['pdb', 'chainid', 'pdb_start', 'pdb_end',
          'uniprot_start', 'uniprot_end', 'Method', 'coverage', 'resolution']
    df = pd.DataFrame(struc_info, columns=df_columns)

    return(df)


def get_coord_df_fix(df, uniprot_sequence):
    """This function makes the dataframe and includes coordinate information"""
    
    ## Define DataFrame
    full_protein_coverage = np.zeros(len(uniprot_sequence))
    full_protein_resolution = np.zeros(len(uniprot_sequence))
    df_coords_columns = ['pdb', 'chainid', 'pdb_start', 'pdb_end',
                         'uniprot_start', 'uniprot_end', 'coverage','length', 'Method', 'resolution','coord_sequence','binary_sequence']
    numpy_array = np.zeros(shape=(len(df['pdb']), len(df_coords_columns)))
    df_coords = pd.DataFrame(numpy_array, columns=df_coords_columns)

    ## Start reading items in DataFrame
    for p in range(len(df['pdb'])):     
        pdb = df['pdb'][p]
        method = df['Method'][p]
        chainid = df['chainid'][p]
        pdb_start = int(df['pdb_start'][p])
        pdb_end = int(df['pdb_end'][p])
        uniprot_start = int(df['uniprot_start'][p])
        uniprot_end = int(df['uniprot_end'][p])
        
        
        all_methods=True
        if all_methods == True:
            try:
            #if 1==1:
                individual_protein_coverage = np.zeros(len(uniprot_sequence))
                individual_protein_resolution = np.zeros(len(uniprot_sequence))
                
                # Get coordinate sequence from the pdb file
                coord_sequence = get_pdb_coord_seq(pdb, chainid)        
                coord_start = pdb_start
                coord_end = pdb_end 
                cover_seq = np.zeros(len(uniprot_sequence))
                
                # Crude cut of the sequence
                if coord_start < uniprot_start:
                    seq = uniprot_sequence[uniprot_start-1:uniprot_end]
                    for n in enumerate(coord_sequence[(uniprot_start-1):coord_end]):
                        if str(n[1]) == seq[n[0]]:
                            cover_seq[n[0]+uniprot_start] = 1      
                else:
                    seq = uniprot_sequence[coord_start-uniprot_start+1:uniprot_end]
                    for n in enumerate(coord_sequence[(coord_start-uniprot_start+1):coord_end]):
                        if str(n[1]) == seq[n[0]]:
                            cover_seq[n[0]+uniprot_start] = 1
                                        
                
                # Add to DataFrame
                df_coords['pdb'][p] = df['pdb'][p]
                df_coords['Method'][p] = df['Method'][p]
                df_coords['coverage'][p] = df['coverage'][p]
                df_coords['length'][p] = int(df['uniprot_end'][p])-int(df['uniprot_start'][p])
                df_coords['chainid'][p] = df['chainid'][p]
                df_coords['pdb_start'][p] = df['pdb_start'][p]
                df_coords['pdb_end'][p] = df['pdb_end'][p]
                df_coords['uniprot_start'][p] = df['uniprot_start'][p]
                df_coords['uniprot_end'][p] = df['uniprot_end'][p]
                df_coords['resolution'][p] = df['resolution'][p]
                df_coords['coord_sequence'][p] = str('NaN')
                df_coords['coord_sequence'][p] = str(coord_sequence)
                df_coords['binary_sequence'][p] = str('NaN')

            except:
                print('error for:  ', pdb, method)
                continue
            
    return(df_coords)

