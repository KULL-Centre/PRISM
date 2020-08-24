##############################################################
#
# Author: Johanna K. S. Tiemann
#
# Copyright: CC
#
# Last edited: 2020-08-24
##############################################################

"""
Extraction of domains, features, pdbs for a single uniprot ip or the whole human genome
information mainly extracted from uniprot and pfam

Commands:
=======

Example:
=======

python run.py -u P10912 -n 1


"""

# Standard library imports
from argparse import ArgumentParser
import logging as log
import os

# Local application imports
from get_domains import extract_single_protein_pfam, extract_pfam_nested, extract_pfam_pdb_mapping, extract_pfam_all_release
from get_human_proteome import extract_uniprot_human_proteome_ids
from get_uniprot_features import extract_uniprot_human_proteome_info, extract_uniprot_info
from map_domains_features import map_pfam_pdb, get_domain_features_human_proteome


log_message = "verbose"

if log_message == "verbose":
    log.basicConfig(
        format='%(levelname)s:%(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=log.INFO
    )
elif log_message == "debug":
    log.basicConfig(
        format='%(levelname)s:%(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=log.WARNING
    )
else:
    log.basicConfig(
        format='%(levelname)s:%(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=log.ERROR
    )

logger = log.getLogger(__name__)


def parse_args():
    """
    Argument parser function
    """

    parser = ArgumentParser( description="" )

    # Initiating/setup command line arguments
    parser.add_argument( '--human_proteome', '-hp',
        type=lambda s: s.lower() in ['true', 't', 'yes', '1'],
        default=False,
        help="Calculate domains and features for the human proteome"
        )
    parser.add_argument( '--uniprot_id', '-u',
        type=str,
        default='',
        help="Uniprot id for single runs (default: '')"
        )
    parser.add_argument( '--reviewed', '-r',
        type=str,
        default="yes",
        help="Obtain only reviewed entries from uniprot (options: *, yes [default]"
        )
    parser.add_argument( '--domainfamily_pdbs', '-dfp',
        type=lambda s: s.lower() in ['true', 't', 'yes', '1'],
        default=False,
        help="Get also all pdbs associated with each single domain family."
        )
    parser.add_argument( '--nested', '-n',
        type=lambda s: s.lower() in ['true', 't', 'yes', '1'],
        default=False,
        help="Get information about nested domains."
        )
    parser.add_argument('--output_dir', '-o',
        type=str,
        default='',
        help="Directory where files will be written. If not specified, no output files will be printed"
        )
    
    args = parser.parse_args()

    # Handle user input errors
    if args.human_proteome is False and args.uniprot_id=='':
        parser.error("Please specify a uniprot id or change to human proteome calulcation mode.")

    # Generate and check saved outputs
    if args.output_dir!='':
        if not os.path.isdir(args.output_dir):
            os.makedirs(args.output_dir)

    return args


def main():
    """
    Main function called as default at the end.
    """
    # get user input arguments
    args = parse_args()

    if args.human_proteome:
        human_proteome_uniprot_ids = extract_uniprot_human_proteome_ids(reviewed=args.reviewed)
        all_pfam_db_release_df = extract_pfam_all_release()
        if args.domainfamily_pdbs:
            pdb_pfam_df = extract_pfam_pdb_mapping()
        else:
            pdb_pfam_df = []
        uniprot_proteome_info_df = extract_uniprot_human_proteome_info(reviewed=args.reviewed)
        domain_feature_info_proteome = get_domain_features_human_proteome(
            human_proteome_uniprot_ids, all_pfam_db_release_df, pdb_pfam_df,
            uniprot_proteome_info_df, write_dir=args.output_dir)
        if args.output_dir=='':
            logger.info(domain_feature_info_proteome)

    if args.uniprot_id:
        uniprot_info_df = extract_uniprot_info(args.uniprot_id, reviewed=args.reviewed)
        all_pfam_db_release_df = extract_pfam_all_release()
        if args.domainfamily_pdbs:
            pdb_pfam_df = extract_pfam_pdb_mapping()
        else:
            pdb_pfam_df = []
        domain_feature_info = map_pfam_pdb(
            args.uniprot_id, all_pfam_db_release_df, pdb_pfam_df, uniprot_info_df, 
            write_dir=args.output_dir)
        if args.output_dir=='':
            logger.info(domain_feature_info)

    if args.nested:
        nested_pfam_df = extract_pfam_nested(write_dir=args.output_dir)
        if args.output_dir=='':
            logger.info(f'Nested domains in proteins: {nested_pfam_df}')


if __name__ == '__main__':
    main()
