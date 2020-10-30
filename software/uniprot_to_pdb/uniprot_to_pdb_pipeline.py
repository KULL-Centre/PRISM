"""uniprot_to_pdb_pipeline includes the framework that run all the other scripts related to UniprotToPDB.

Author: Anders Frederiksen

Date of last major changes: 2020-10-09

"""


# Standard library imports
import time
from os.path import isfile,join

# Third party imports
import pandas as pd
import numpy as np

# Local application imports
from apis import req_sifts, get_pdb_coord_seq, pdb_to_fasta_seq, read_fasta
from helper import check_for_file, save_df, load_df
from make_dataframes import make_df, get_coord_df_fix
from msa import run_alignment, parse_alignment, bin_to_df



def master_execute(uniprot_accession,output, overwrite=False):
    """This function executes all the functions below in the correct order"""

    print('Running protein:    ',uniprot_accession)
    #Check if we already have data for this protein
    status= check_for_file(uniprot_accession,output='/groups/sbinlab/haagenb/project/uniprot_to_pdb/protein_dictionaries/')       
    
    if status == False or overwrite == True:
        
        ## Prepare dataframe generation ##       
        #Get data from Sifts and uniprot sequence and make DataFrame
        data = req_sifts(uniprot_accession)
        uniprot_seq = read_fasta(uniprot_accession)
        df=make_df(data,uniprot_seq)  
        #Get coordinate sequence for all pdbs in DataFrame and Save DataFrame
        df_coords=get_coord_df_fix(df, uniprot_seq) 
        save_df(df_coords,uniprot_accession,output='/groups/sbinlab/haagenb/project/uniprot_to_pdb/protein_dictionaries_raw/')
        
        
        
        # Prepare MSA ##
        #Make sequence file and run alignment
        outfile=run_alignment(uniprot_seq,df_coords,output)
        #Parse Alignment and insert into DataFrame
        binary_data=parse_alignment(outfile)
        df_coords=bin_to_df(df_coords,binary_data)
        #Saving DataFrame to folder as .csv
        save_df(df_coords,uniprot_accession,output='/groups/sbinlab/haagenb/project/uniprot_to_pdb/protein_dictionaries/')        
        
    if status == True and overwrite == False:
        print('Found existing file: Loading')
        #Load DataFrame
        df_coords=load_df(uniprot_accession,output='/groups/sbinlab/haagenb/project/uniprot_to_pdb/protein_dictionaries/')
    
    time.sleep(10)
    return(df_coords)

#if __name__ == '__main__':
#        gene_list_folder='/groups/sbinlab/haagenb/project/uniprot_to_pdb/TorbensGenes.csv'
#    gene_list_df=pd.read_csv(gene_list_folder)
#    uniprot_id_list_df=gene_list_df['Entry']
#    uniprot_accession_list=[]
#    for n in uniprot_id_list_df:
#        uniprot_accession_list=np.hstack((uniprot_accession_list,n))
#        
#    
#    df_dic = {}
#    error_list=[]
#    for uniprot_accession in uniprot_accession_list:
#        try:
#            df_coords= master_execute(uniprot_accession,overwrite=True)
#            df_dic[uniprot_accession] = df_coords
#        except:
#            print('Error for protein:  ',uniprot_accession)
#            error_list.append(uniprot_accession)
#    
#    print('Completed') 
#    print(error_list)
#    run_all()



