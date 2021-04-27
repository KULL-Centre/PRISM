
#goes through results in /storage1/shared/data/netsurfp_preds/csv and turns them into prism_netsurfp files saved in /storage1/shared/data/prism_netsurfp/uniprotID[0:2]/uniprotID[2:4]/uniprotID[4:6]/ if those don't exist yet 

import argparse
import sys
import re
import os
from datetime import datetime
import numpy as np
import pandas as pd 

import_path_base = '/storage1/hezscha/src/'
sys.path.insert(1, import_path_base + 'PRISM/prism/scripts/')
from PrismData import PrismParser, VariantData#VariantParser, VariantData

version = '001'
base_prism_dir = '/storage1/shared/data/'

#parse arguments
################################################################################
parser = argparse.ArgumentParser()
parser.add_argument('-in_folder', dest='in_folder', help="Folder with the netsurfp csv output to go through. Default location is: "+base_prism_dir+"netsurfp_preds/csv")
parser.add_argument('-out_folder', dest='out_folder', help="Where the output files should be written. Default location is "+base_prism_dir+"prism_netsurfp/uniprot[0:2]/unipro[2:4]/uniprot[4:6]/")
parser.add_argument("-v","--verbose", action="count", default=0, help="Level of output, default zero is no output")
#parser.add_argument("--first_res_num", metavar="INT", help="In the output data, assign this number to the first given amino acid of the sequence")
args = parser.parse_args()
################################################################################



def write_prism(metadata, dataframe, prism_file, comment=''):
	variant_dataset = VariantData(metadata, dataframe)
	parser = PrismParser()
	parser.write(prism_file, variant_dataset, comment_lines=comment)

def read_from_prism(primsfile):
	parser = PrismParser()
	dataframe = parser.read(primsfile).dataframe
	meta_data = parser.read_header(primsfile)
	return meta_data, dataframe

def make_default_outfolder(uniprot_ID):
	#sometimes we cannot find the unirpot ID to a transcript. Put them in the none folder then
	if not uniprot_ID == 'None':
		if not os.path.exists(base_prism_dir+'prism_netsurfp/'+ uniprot_ID[0:2]):
			os.makedirs(base_prism_dir+'prism_netsurfp/'+ uniprot_ID[0:2])
		if not os.path.exists(base_prism_dir+'prism_netsurfp/'+ uniprot_ID[0:2]+ '/' + uniprot_ID[2:4]):
			os.makedirs(base_prism_dir+'prism_netsurfp/'+ uniprot_ID[0:2]+ '/' + uniprot_ID[2:4])
		if not os.path.exists(base_prism_dir+'prism_netsurfp/'+ uniprot_ID[0:2]+ '/' + uniprot_ID[2:4] + '/' + uniprot_ID[4:6]):
			os.makedirs(base_prism_dir+'prism_netsurfp/'+ uniprot_ID[0:2]+ '/' + uniprot_ID[2:4]+ '/' + uniprot_ID[4:6])
	
		return(base_prism_dir+'prism_netsurfp/'+ uniprot_ID[0:2]+ '/' + uniprot_ID[2:4]+ '/' + uniprot_ID[4:6])
	else:
		if not os.path.exists(base_prism_dir+'prism_netsurfp/None/'):
			os.makedirs(base_prism_dir+'prism_netsurfp/None/')
		return(base_prism_dir+'prism_netsurfp/None/')

def read_uniprot_datasets(db_path = '/storage1/shared/data/uniprot_datasets/'):
	'''
	Read in the datasets mentioned in the index of the db_path and return a dict that has one key for each short name and a python set of all uniprot IDs in that list as the value 
	'''
	dbs_d = {}
	with open(db_path+'index', 'r') as idx_IN:
		header = idx_IN.readline()
		for line in idx_IN:
			(file_name, short_name) = line.rstrip().split()
			dbs_d[short_name] = set()
			with open(os.path.join(db_path + file_name), 'r') as IN:
				for line in IN:
					dbs_d[short_name].add(line.rstrip())
	
	return dbs_d			
	

def check_membership(uniprot_id, dbs_d):
	return_d = {}
	for short_name in dbs_d:
		if uniprot_id in dbs_d[short_name]:
			return_d[short_name] = 'yes'
		else:
			return_d[short_name] = 'no'

	return return_d

#read in uniprot datasets
dbs_d = read_uniprot_datasets()

#go through all csv files an make prism files if those don't exit
if args.in_folder:
	csv_path = args.in_folder
else:	
	csv_path =  os.path.join(base_prism_dir, 'netsurfp_preds/csv')
#csv_files = os.listdir(csv_path)
#for csv in csv_files[:2]:
for csv in os.scandir(csv_path):
	uID = csv.name.split('/')[-1].split('.')[0]
	if args.out_folder:
		if args.out_folder == '.':
			args.out_folder = os.getcwd()
		prism_file = os.path.join(args.out_folder, 'prism_netsurfp_' + version + '_' + uID + '.txt')
	else:
		prism_path = make_default_outfolder(uID)
		prism_file = os.path.join(prism_path, 'prism_netsurfp_' + version + '_' + uID + '.txt')
	
	if os.path.exists(prism_file):
		continue
	else:
		print(uID)
		#read in as a df
		df = pd.read_csv(csv, header = 0)
		df.columns = ['name','aa','pos','rsa','asa','q3_ss_pred','p_q3_H','p_q3_E','p_q3_C',
		'q8_ss_pred','p_q8_G','p_q8_H','p_q8_I','p_q8_B','p_q8_E','p_q8_S','p_q8_T','p_q8_C',
		'phi','psi','disorder']
		#add burried/exposed and variant cols
		df["burried_exposed"] = np.where(df['rsa']>0.25, 'e', 'b')
		#https://stackoverflow.com/questions/11858472/string-concatenation-of-two-pandas-columns
		#df['baz'] = df.agg(lambda x: f"{x['bar']} is {x['foo']}", axis=1)
		df['variant'] = df.agg(lambda x: f"{x['aa']}{x['pos']}=", axis=1)
		
		#get the seq from the column		
		seq = ''.join(df['aa'])
		#print(seq)
		
		#remove some cols
		df.drop(columns=['name', 'aa', 'pos'])
		
		df=df.reindex(columns= ['variant', 'burried_exposed','rsa','asa','q3_ss_pred','p_q3_H','p_q3_E','p_q3_C',
		'q8_ss_pred','p_q8_G','p_q8_H','p_q8_I','p_q8_B','p_q8_E','p_q8_S','p_q8_T','p_q8_C',
		'phi','psi','disorder'])
		
		#check
		if args.verbose >= 1:
			print(df.columns)
			print(df.head())
		
		#print prism file
		columns_dic = {'burried_exposed': 'Classification if residue is burried (RSA <= 0.25) or exposed (RSA > 0.25)',
		'rsa': 'Relative Surface Accessibility', 
		'asa': 'Absolute Surface Accessibility', 
		'q3_ss_pred': 'Secondary structure prediction in 3 class scheme (q3)', 
		'p_q3_H': 'probability of class H (helix) in q3 assignment', 
		'p_q3_E': 'probability of class E (strand) in q3 assignment', 
		'p_q3_C':'probability of class C (coil) in q3 assignment',
		'q8_ss_pred': 'Secondary structure prediction in 8 class scheme (q8)', 
		'p_q8_G': 'probability of class G (3-turn helix) in q8 assignment', 
		'p_q8_H': 'probability of class H (4-turn helix) in q8 assignment', 
		'p_q8_I':'probability of class I (5-turn helix) in q8 assignment',
		'p_q8_B': 'probability of class B (isolated beta-bridge) in q8 assignment', 
		'p_q8_E': 'probability of class E (strand) in q8 assignment', 
		'p_q8_S':'probability of class S (bend) in q8 assignment',
		'p_q8_T': 'probability of class T (turn) in q8 assignment', 
		'p_q8_C': 'probability of class C (coil) in q8 assignment', 
		'phi': 'phi angle', 
		'psi': 'psi angle', 
		'disorder':'probability of residue being disordered',

		
		}
		
		#metadata dict
		metadata = {
			"version": version,
			"protein": {
				"name": uID,
				"organism": 'Homo sapiens',
				"sequence": seq,
				"uniprot": uID,
				#"Ensembl transcript ID": ','.join(neighbors[geneID][tID]),
				#"chromosome": args.chromosome,
				#"is canonical": meta[geneID]['is_canonical']
			},
			#"uniprot": {},
			"netsurfp": {
				"search method": 'hhblits'
			},
			"columns": columns_dic,
			}
			
		#evaluate membership in uniprot sets:
		membership = check_membership(uID, dbs_d)
		#if args.verbose >= 1:
		#	print(membership)
		for dataset in membership:
			metadata['protein'][dataset] = membership[dataset]	
		
		write_prism(metadata, df, prism_file)
		metadata, df_in = read_from_prism(prism_file)
	
	#break	#for debugging			

