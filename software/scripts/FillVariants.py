##############################################################
#
# Author: Johanna K. S. Tiemann
#
# Copyright: CC
#
# Last edited: 2020-11-12
##############################################################

"""
Function creating each possible variant with WT information

Commands:
=======


Example:
=======

python FillVariants.py prism_uniprot_XXX_uniprot_id.txt


"""

# Standard library imports
from argparse import ArgumentParser
from datetime import datetime
import os

# Third party imports
import pandas as pd 

# Local application imports
from PrismData import PrismParser, VariantData


def parse_args():
    """
    Argument parser function
    """

    parser = ArgumentParser( description="" )

    # Initiating/setup command line arguments
    parser.add_argument( 'file',
        type=str,
        help="Prism data file to process"
        )
    parser.add_argument('--output_dir', '-o',
        type=str,
        default='',
        help="Directory where files will be written. If not specified, the input directory will be used."
        )
    parser.add_argument('--incl_del', '-i',
        default=False,
        type=lambda s: s.lower() in ['true', 't', 'yes', '1'],
        help="generates values for ~ and * too. Default False"
        )
    
    args = parser.parse_args()

    # Generate and check saved outputs
    if args.output_dir!='':
        if not os.path.isdir(args.output_dir):
            os.makedirs(args.output_dir)
    else:
        args.output_dir = os.path.dirname(os.path.abspath(args.file))
    if args.incl_del != False:
        args.incl_del = True
    return args


def copy_wt_variants(dataframe, incl_del=False):
    #copy wt columns to all variants
    
    if 'aa_var' in dataframe.columns:
        dataframe = dataframe.drop(columns=['aa_var'])
    final_df = dataframe.copy()
    final_df['variant'] = final_df['variant'].str.replace(r'=', 'A')
    final_df['aa_var'] = ['A']*len(dataframe['variant'])
    if incl_del:
        resi_list = ['C','D','E', 'F','G','H', 'I','K','L','M','N','P','Q','R','S','T','V','W','Y', '~', '*']
    else:
        resi_list = ['C','D','E', 'F','G','H', 'I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    for variant in resi_list:
        df = dataframe.copy()
        df['variant'] = df['variant'].str.replace('=', variant)
        df['aa_var'] = [variant]*len(dataframe['variant'])
        final_df = pd.concat([final_df, df], axis=0, ignore_index=True)
    final_df['aa_var'] = final_df['variant'].apply(lambda x: x[-1])
    final_df['aa_ref'] = final_df['variant'].apply(lambda x: x[0])
    final_df['res'] = final_df['variant'].apply(lambda x: int(x[1:-1]))
    mask = final_df['aa_ref']==final_df['aa_var']
    final_df.loc[mask, 'aa_var'] = '='
    final_df['variant'] = final_df['aa_ref'] + final_df['res'].astype(str) + final_df['aa_var']
    final_df = final_df.sort_values(by=['res','aa_var'])
    for elem in ['n_mut', 'aa_ref', 'resi', 'aa_var', 'res']:
        if elem in final_df.columns:
            final_df = final_df.drop(columns=elem)
    final_df = final_df.drop_duplicates().reset_index(drop=True)
    return(final_df)


def write_prism(metadata, dataframe, prism_file, comment=''):
    variant_dataset = VariantData(metadata, dataframe)
    parser = PrismParser()
    parser.write(prism_file, variant_dataset, comment_lines=comment)


def read_from_prism(primsfile):
    parser = PrismParser()
    dataframe = parser.read(primsfile).dataframe
    meta_data = parser.read_header(primsfile)
    meta_data.pop('variants')
    return meta_data, dataframe


def main():
    """
    Main function called as default at the end.
    """
    # get user input arguments
    args = parse_args()
    
    # read metadata and dataframe
    meta_data, dataframe = read_from_prism(args.file)
    
    # create new dataframe
    new_df = copy_wt_variants(dataframe, incl_del=args.incl_del)
    
    # write dataframe
    file_name = os.path.basename(args.file)
    new_prism_file = os.path.join(args.output_dir, f'{file_name[:-4]}_filled{file_name[-4:]}')
    comment = [ f"{datetime.date(datetime.now())} - expanded variant data",] 
    write_prism(meta_data, new_df, new_prism_file, comment)


if __name__ == '__main__':
    main()
