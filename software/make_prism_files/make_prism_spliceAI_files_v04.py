
#v04: if fails to write or read in prism file, make a note and continue to next transcript in the gene. This is due to the error I'm getting for ENSG00000156239 N6AMT1 ENST00000303775 where mysteriously it seems spliceAI has a different ref nucleotide than ensembl causing an abbarent WT amino acid call

#v02: trying to get around several DNA SNVs causing the same protein var.
#the idea is to use a dict of dicts instead of a list of dicts during the data storing. the key will be the variant name and it will also be inside the sub dict. If we process a variant name we already have in the dod, we only keep the df line i.e. subdict with the more extreme delta scores.  
#Example: Q86YR6, POTED, ENSG00000166351, ENST00000299443
#51   E399D      0.93   0.00   0.00   0.03   0.93 13628452G>C
#52   E399D      0.93   0.00   0.00   0.02   0.93 13628452G>T
#since one of the delta scores is larger in G>C only keep that one. If all delta scores are equal just keep one


import argparse
#this is the data for relase 100. The current relase server is https://.rest.ensembl.org
#e100 doesn't work, need to try it in brower and see what it evaluates to
ensembl_server = "https://apr2020.rest.ensembl.org"
import requests #url requests
from requests import HTTPError
import sys
import re
import os
from datetime import datetime
import cyvcf2
import numpy as np
import pandas as pd 


import_path_base = '/storage1/hezscha/src/'
sys.path.insert(1, import_path_base + 'PRISM/prism/scripts/')
from PrismData import PrismParser, VariantData#VariantParser, VariantData
sys.path.insert(1, import_path_base + 'gnomad_to_prism/scripts/')
from vcf_parser_func_v09_revep import get_uniprot, get_prot_seq, parse_HGVSp
#from vcf_parser_func_v09_revep_dev import get_uniprot, get_prot_seq, parse_HGVSp

#parse arguments
################################################################################
parser = argparse.ArgumentParser()
parser.add_argument('-m', dest="mode", choices=['overwrite', 'leave'], default = 'leave', help="What do when the output file already exists. Leave (default) or overwrite")
parser.add_argument('-chrom', dest="chromosome", type = str, required=True, help="The chromosome for which to parse.")
parser.add_argument('-cutoff', dest="cutoff", choices=[0.2, 0.5], type = float, required=True, help="The delta cutoff for when spliceAI preds are considered splice altering.")
parser.add_argument('-v', dest="v", action='store_true', help="Verbose")
parser.add_argument('-out_folder', dest='out_folder', help="Use to change where the output files should be written. Default location is /storage1/shared/data/prism_spliceai/'+ uniprot_ID[0:2]+ '/' + uniprot_ID[2:4]+ '/' + uniprot_ID[4:6]")
args = parser.parse_args()
################################################################################

fail_out = open('/storage1/shared/data/prism_spliceai/failed_parsing.txt', 'a')

#functions
#PRISM parser functions
def write_prism(metadata, dataframe, prism_file, comment=''):
	variant_dataset = VariantData(metadata, dataframe)
	parser = PrismParser()
	parser.write(prism_file, variant_dataset, comment_lines=comment)

def read_from_prism(primsfile):
	parser = PrismParser()
	dataframe = parser.read(primsfile).dataframe
	meta_data = parser.read_header(primsfile)
	return meta_data, dataframe

def fill_t_matches(geneID):	
	#query ensembl with geneID, get all coding transcripts and compare their protein seqs to establish mapping of all t's with same protein product to one representative t_ID
	prot_seqs[geneID] = {}
	t_matches[geneID] = {}
	neighbors[geneID] = {}
	meta[geneID] = {}

	try:
		r = requests.get(ensembl_server+'/lookup/id/'+geneID+'?expand=1', headers={ "Content-Type" : "application/json"})
		#debug
		if args.v:
			print('Request: ' + ensembl_server+'/lookup/id/'+geneID+'?expand=1')
		#debug
		r.raise_for_status()
	except HTTPError as e:	
		status_code = e.r.status_code
		print('The following error occured', status_code)
	else:
		#if args.v:
			#print(r.text)
			
		dec = r.json()
		meta[geneID]['name'] = dec['display_name']
		meta[geneID]['organism'] = dec['species']
		
		#each transcript is a sublist
		for t in dec['Transcript']:
			#print(t['id'])
			#print(t['biotype'])
			#print(t['is_canonical'])
			if t['is_canonical']:
				meta[geneID]['is_canonical'] = t['id']
			
			if t['biotype'] == 'protein_coding':
				#get prot seq
				try:
					r_p = requests.get(ensembl_server+"/sequence/id/"+t['id']+"?type=protein", headers={ "Content-Type" : "text/plain"})
					r_p.raise_for_status()
				except HTTPError as e:	
					status_code = e.r.status_code
					print('The following error occured', status_code)
				else:
					curr_seq = r_p.text
					#does this match to any seqs already known for this gene? 
					for other_t in prot_seqs[geneID]:
						if prot_seqs[geneID][other_t] == curr_seq:
							#note that curr transcript matches this representative transcript
							t_matches[geneID][t['id']] = other_t
							neighbors[geneID][other_t].append(t['id'])
							break
					#Otherwise make new representative entry for curr transcript and match it to itself
					else:
						prot_seqs[geneID][t['id']] = curr_seq
						t_matches[geneID][t['id']] = t['id']
						neighbors[geneID][t['id']] = [t['id']]

def make_default_outfolder(uniprot_ID):
	#sometimes we cannot find the unirpot ID to a transcript. Put them in the none folder then
	base_dir = '/storage1/shared/data/prism_spliceai/'
	if not uniprot_ID == 'None':
		if not os.path.exists(base_dir+ uniprot_ID[0:2]):
			os.makedirs(base_dir+ uniprot_ID[0:2])
		if not os.path.exists(base_dir+ uniprot_ID[0:2]+ '/' + uniprot_ID[2:4]):
			os.makedirs(base_dir+ uniprot_ID[0:2]+ '/' + uniprot_ID[2:4])
		if not os.path.exists(base_dir+ uniprot_ID[0:2]+ '/' + uniprot_ID[2:4] + '/' + uniprot_ID[4:6]):
			os.makedirs(base_dir+ uniprot_ID[0:2]+ '/' + uniprot_ID[2:4]+ '/' + uniprot_ID[4:6])
	
		return(base_dir+ uniprot_ID[0:2]+ '/' + uniprot_ID[2:4]+ '/' + uniprot_ID[4:6] + '/')
	else:
		if not os.path.exists(base_dir+'None/'):
			os.makedirs(base_dir+'None/')
		return(base_dir+'None/')

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

#set up
version = '1'

#data structures:
dbs_d = read_uniprot_datasets()

#dict of transcript matches. Two transcripts match if they have the same protein product
#Format: t_matches[geneID] = {t1: t1 , t2: t1, t3: t1, t4: t4, t5: t4, ... } 
t_matches = {}
prot_seqs = {}
neighbors = {}
meta = {}
dfs = {}

#Step 1: open input vcf
data_dir = "/storage0/shared/data/spliceAI/vep/"
#output_dir = "/storage0/shared/data/prism_spliceAI/"
vcffile = data_dir+'spliceai_scores.masked.snv.hg38.chromosome'+args.chromosome+'.filter.spliceAI'+str(args.cutoff)+'.vcf.VEP.cache100.bgz'
#vcffile = data_dir + 'test.vcf'

print('Parsing vcf', vcffile, '...')
try:
	vcf_reader = cyvcf2.VCF(vcffile)
except:
	raise
	sys.exit()

#Step 2: find out positions of fields in vep anno
for line in vcf_reader.raw_header.split('\n'):
	if line.startswith('##INFO=<ID=CSQ'):
		#get this part:"... Format: Allele|..."
		fields = line.split('Format: ')[1].split('|')
		pos = {} #a dict storing the position of each field we search for, or False if the field is not available
		#find out where certain fields are or if they are not available
		#the fields I'm looking for:
		target_fields = ['Gene', 'SYMBOL', 'BIOTYPE', 'CANONICAL', 'Consequence', 'HGVSc','HGVSp', 'Codons', 'SWISSPROT', 'DOMAINS', 'Feature', 'Protein_position', 'Amino_acids']
		#the information I want printed per variant. F.x. if a transcript is canon does not change and neither does the uniprot ID so no reason to print them for every var 
		print_fields = ['Consequence', 'HGVSp', 'Codons', 'SWISSPROT', 'DOMAINS', 'Protein_position', 'Amino_acids']
		for target_field in target_fields:
			pos[target_field] = False
			for (position, field) in enumerate(fields):
				if field == target_field: 
					pos[target_field] = position
					break
		print('Found vep field information')
		#debug
		if args.v:
			print('Field positions:')
			for target_field in pos:
				print(target_field, pos[target_field]) 
		
		#debug
		break

	else:
		continue
else:
	#if we didn't break we didn't find the vep line
	print('Did not find annotation metadata line ##INFO=<ID=CSQ, please check vcf file')
	sys.exit()

#Step 3: Go through body
for record in vcf_reader:
	#first of, skip SNPs that did not PASS quality control. In vcf reader, PASS evaluates to nothing, so if there is something in the filter you want to skip the entry
	if record.FILTER:
		continue

	#remember the nucl change:
	#nuc_change = 

	#get spliceAI values
	sfields = record.INFO.get('SpliceAI').split('|')
	DS_AG = float(sfields[2])
	DS_AL = float(sfields[3])
	DS_DG = float(sfields[4])
	DS_DL = float(sfields[5])
	maxdelta = max(DS_AG, DS_AL, DS_DG, DS_DL)

	for transcript in record.INFO.get('CSQ').split(','): 
		fields = transcript.split('|')
		
		#we are only interested in coding transcripts and the variant types requested by argv (filter set is adjusted accordingly in step 1)
		if fields[pos['BIOTYPE']] == 'protein_coding':

			#get geneID
			geneID = fields[pos['Gene']]
			#check if geneID alread in t_matches
			#if not fill out t_matches for geneID
			if not geneID in t_matches:
				fill_t_matches(geneID)
				dfs[geneID] = {}
				
			#now, match the current transcript to the repr transcript
			#optimization: as discussed here it is more economical to create a list of dicts and then turn them into a df in the end than appending to a dataframe because then the data is copied every time (though dfs are supposedly mutable?) https://stackoverflow.com/questions/10715965/add-one-row-to-pandas-dataframe#24888331
			#the desired struct is something like: [{'A': 45, 'B': 57, 'C': 4, 'D': 68, 'E': 44}, {'A': 79, 'B': 38, 'C': 71, 'D': 55, 'E': 27}, {'A': 28, 'B': 7, 'C': 55, 'D': 21, 'E': 1}]. All dicts in the list must have the same keys obviously
			
			tID = t_matches[geneID][fields[pos['Feature']]]
			#is there already a dataframe for that transcript?
			if not tID in dfs[geneID]:
				#init df
				#for empty df, skip index and data arguments. 
				#add to empty dfs (no rows) the same way you set data in dfs already containing rows: df.loc[row,col] = 'foo'
				#dfs[geneID][tID] = pd.DataFrame(columns=['variant', 'maxdelta', 'DS_AG', 'DS_AL', 'DS_DG', 'DS_DL'])
				#lod
				#dfs[geneID][tID] = []
				#dod
				dfs[geneID][tID] = {}
				
			#save spliceAI info to that tID
			var_name = parse_HGVSp(transcript, HGSV_field_pos=pos['HGVSp'], prot_pos_field_pos=pos['Protein_position'], aa_field_pos=pos['Amino_acids'], v=False)
			#if there is no HGVSp, the parse_HGVSp will return None. This is how I control that I only get vars in that I can describe in the prism framework
			#add: skip stop loss mutations because they are not covered in parse_HGVSp and I'm not sure they would be useful in the prism framework because we can't model them in Rosetta and there will be no MAVE data either
			#debug
			#print(fields[pos['Consequence']])
			#debug
			if var_name is not None and not 'stop_lost' in fields[pos['Consequence']] and not 'stop_retained_variant' in fields[pos['Consequence']]:
				#debug
				#print(var_name, fields[pos['Consequence']])
				#print(record.POS, ', HGVSp:', fields[pos['HGVSp']], ', var_name:*', var_name, '*', sep = '')
				#debug
				nuc_change = fields[pos['HGVSc']].split(':')[-1][2:]
				temp_d = {'variant' : var_name, 'maxdelta': maxdelta, 'DS_AG': DS_AG, 'DS_AL':DS_AL, 'DS_DG':DS_DG, 'DS_DL':DS_DL, 'HGVSc': nuc_change}
				#dfs[geneID][tID].append(temp_d)
				#now check if aleady have a line for that protein var
				if var_name in dfs[geneID][tID]:
					#go through the delta values and compare. If any of the new ones are more extreme and all are at least as high the old ones, replace
					if ((temp_d['DS_AG'] > dfs[geneID][tID][var_name]['DS_AG']) or (temp_d['DS_AL'] > dfs[geneID][tID][var_name]['DS_AL']) or (temp_d['DS_DG'] > dfs[geneID][tID][var_name]['DS_DG']) or (temp_d['DS_DL'] > dfs[geneID][tID][var_name]['DS_DL'])) and ((temp_d['DS_AG'] >= dfs[geneID][tID][var_name]['DS_AG']) and (temp_d['DS_AL'] >= dfs[geneID][tID][var_name]['DS_AL']) and (temp_d['DS_DG'] >= dfs[geneID][tID][var_name]['DS_DG']) and (temp_d['DS_DL'] >= dfs[geneID][tID][var_name]['DS_DL'])):
						dfs[geneID][tID][var_name] = temp_d
					
					#else, if all deltas are the same, just add the DNA change
					elif (temp_d['DS_AG'] == dfs[geneID][tID][var_name]['DS_AG']) and (temp_d['DS_AL'] == dfs[geneID][tID][var_name]['DS_AL']) and  (temp_d['DS_DG'] == dfs[geneID][tID][var_name]['DS_DG']) and (temp_d['DS_DL'] == dfs[geneID][tID][var_name]['DS_DL']):
						dfs[geneID][tID][var_name]['HGVSc'] += ',' + nuc_change
					
				#else just create a subdict for that var_name	
				else:
					dfs[geneID][tID][var_name] = temp_d

#debug
#when a gene gets into meta it also gets into t_matches
if args.v:
	for geneID in meta:
		print(geneID, t_matches[geneID])


#Step 4: print prism files for each representative transcript in each gene
for geneID in meta:
	#the is one seq saved for each unique protein product and the keys are the representative transcripts
	for tID in prot_seqs[geneID]:
		
		#however, some transcript might not have any variants (prot seqs has all existing unique protein products)
		if not tID in dfs[geneID]:
			continue
		
		uniprot_ID = get_uniprot(tID, no_isoform = True)
		#check if output dir for the uniprot ID exists and otherwise make it
		if args.out_folder:
			if args.out_folder == '.':
				args.out_folder = os.getcwd()
			out_folder = args.out_folder
		else:
			out_folder = make_default_outfolder(uniprot_ID)
		
		prism_file = os.path.join(out_folder, 'prism_spliceai_XXX_' +uniprot_ID + '-' + tID + '-' + geneID +'.txt')
		
		if os.path.exists(prism_file):
			if args.mode == 'leave':
				print(prism_file, 'already exists. If you wish to overwrite it use -m overwrite.')
				continue #should go to next iteration of the closest still running for loop which should be the next file
		
		#make list of dicts from the dict of dicts
		temp_lod = []
		for var_name in dfs[geneID][tID]:
			temp_lod.append(dfs[geneID][tID][var_name])
		
		#turn list of dicts into df
		#output_df = pd.DataFrame(dfs[geneID][tID], columns = ['variant', 'maxdelta', 'DS_AG', 'DS_AL', 'DS_DG', 'DS_DL'])
		#lod from dod
		output_df = pd.DataFrame(temp_lod, columns = ['variant', 'maxdelta', 'DS_AG', 'DS_AL', 'DS_DG', 'DS_DL', 'HGVSc'])
		
		#if df is empty, skip
		if output_df.empty:
			continue
		
		#debug
		print(geneID, meta[geneID]['name'], tID)
		if args.v:
			print(prot_seqs[geneID][tID])
			print(output_df.head())
			output_df.to_csv(geneID+'_'+tID+'_df.csv', sep = '\t')
		#debug
		
		#columns dict for header
		columns_dic = {'maxdelta': 'naximum delta score', 'DS_AG': 'delta score for acceptor gain', 'DS_AL': 'delta score for acceptor loss', 'DS_DG': 'delta score for donor gain', 'DS_DL': 'delta score for donor loss', 'HGVSc':'change on nucleotide level'}
		#metadata dict
		metadata = {
			"version": version,
			"protein": {
				"name": meta[geneID]['name'],
				"organism": meta[geneID]['organism'],
				"sequence": prot_seqs[geneID][tID],
				"uniprot": uniprot_ID,
				"Ensembl transcript ID": ','.join(neighbors[geneID][tID]),
				"chromosome": args.chromosome,
				"is canonical": meta[geneID]['is_canonical']
			},
			#"uniprot": {},
			"spliceai": {
				"cutoff": str(args.cutoff)
			},
			"columns": columns_dic,
			}
		
		#evaluate membership in uniprot sets:
		membership = check_membership(uniprot_ID, dbs_d)
		#if args.verbose >= 1:
		#	print(membership)
		for dataset in membership:
			metadata['protein'][dataset] = membership[dataset]
		
		
		#make prism file
		comment = [ f"version {version} - {datetime.date(datetime.now())} - SpliceAI prediction",] 
		#logger.info('Writing prism file')
		print('Writing prism file')
		try:
			write_prism(metadata, output_df, prism_file, comment=comment)
		except Exception as e:
			print('Something went wrong parsing:', args.chromosome, geneID, meta[geneID]['name'], tID)
			print(e)
			print(args.chromosome, geneID, meta[geneID]['name'], tID, file = fail_out)
			continue
		#load in again to check
		try:
			meta_data, dataframe = read_from_prism(prism_file)
			print('Prism file passed check')
		except Exception as e:
			print('Something went wrong reading in prism file:', args.chromosome, geneID, meta[geneID]['name'], tID)
			print(e)
			print(args.chromosome, geneID, meta[geneID]['name'], tID, file = fail_out)
			continue


fail_out.close()
