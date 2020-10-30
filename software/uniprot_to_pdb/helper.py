"""helper includes functions to save, load and check for dataframes.

Author: Anders Frederiksen

Date of last major changes: 2020-10-09

"""

# Standard library imports
from os.path import isfile,join
from os import listdir

# Third party imports
import numpy as np
import pandas as pd




def check_for_file(uniprot_accession_list,output='/groups/sbinlab/haagenb/project/uniprot_to_pdb/protein_dictionaries/'):
    """Checks if a DataFrame for this protein is already in directory"""
    
    ## Check folder for files starting with uniprotID
    folder=output
    files = [f for f in listdir(folder) if isfile(join(folder, f))]
    
    
    file_list=[]
    for n in files:
        m=n.split('_')[0]
        file_list.append(m)
    
    
    if isinstance(uniprot_accession_list, str) == True:
        uniprot_accession_list = [uniprot_accession_list]
    
    ## Check what uniprotID that are in the Uniprot_accession list
    for uniprot_accession in uniprot_accession_list:
        if uniprot_accession in file_list:
            status=True
        else:
            status=False
    return(status)

def save_df(df, uniprot_accession, output='/groups/sbinlab/haagenb/project/uniprot_to_pdb/protein_dictionaries/'):
    """Saving the dataframe to the output folder"""
    
    ## Saving the dataframe as csv
    key_name = f'{uniprot_accession}_df.csv'
    output_key = join(output, key_name)
    df.to_csv(output_key)


def load_df(uniprot_accession, output):
    """Loading a dataframe"""
    
    ## Reading the .csv
    key_name = f'{uniprot_accession}_df.csv'
    output_key = join(output, key_name)
    df = pd.read_csv(output_key, index_col=['Unnamed: 0'])
    
    ## Fixing the binary sequence
    num=1
    temp=[]
    try:
        for p in range(len(df['binary_sequence'])):
            x = df['binary_sequence'][p]
            
            new_x=[]
    
            x=x.split(' ')
    
            for y in x[1:-1]:
                new_x=np.hstack((new_x,int(y[0])))
    
            if num == 1:
                temp = np.zeros(len(new_x))
            if p != 'Uniprot':
                if new_x == []:
                    new_x = np.zeros(len(temp))
                temp = new_x+temp
    
                num += 1
            df['binary_sequence'][p]=np.zeros(len(new_x))
            df['binary_sequence'][p]=new_x
    except:
        print('could not load binary')
    
    return(df)