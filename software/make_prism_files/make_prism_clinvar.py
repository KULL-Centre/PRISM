
#v06 another overhaul of how I save info until we reach the next symbol. This is assuming we completely process one symbol before going to the next one since the input file is sorted, and therefore only records info for 1 symbol at a time
#v05 I have now made a version of the input file that is sorted by HGNC (1 HGNC : 1 symbol), so when the symbol changes we can process the previous symbol because we've definitely see all the vars in that symbol and we can make a note that it's been processed
#v04 switch from going over the GI and elink to directly requesting the protein sequence from efetch with rettype=fasta_cds_aa. Because of this entries like NM_001708.2 that has been removed from refseq and entries that have been replaced with later versions can still get a protein seq even though elink doesn't work for those
#this also means I'm keeping the version number of the transcript for the protein seq lookup, but not fo uniprot lookup and the file name

'''
plan:

1. divide data into different HGNCs
2. within each HGNC, compare the protein sequences of all ncbi transcripts named in the clinvar entry names, Make a list of all transcripts that make the same protein product, these will share their variants since only look at vars on protein level
3. per unique protein product, make set of all uniprot IDs matching to all transcript IDs
4. make a file for each uniprotID in the above list with the uniprot ID in the file name and one representative of the transcript IDs. All assoc uniprot and transcr ID are named in the header


'''

import argparse
from datetime import datetime
import glob
import gzip
import json
import logging as log
import os
import requests
import shutil
import subprocess
import sys
import re
from requests import HTTPError

# Third party imports
#from Bio.PDB import PDBParser
#from Bio.PDB.PDBList import PDBList 
#from Bio.PDB.DSSP import DSSP
#from IPython.display import display
#import matplotlib.pyplot as plt
import numpy as np
import pandas as pd 
#import seaborn as sns
from scipy import stats
#import sklearn.model_selection
#import xmltodict
from time import sleep

log_message="verbose"

if log_message=="verbose":
	log.basicConfig(
		format='%(levelname)s:%(message)s',
		datefmt='%m/%d/%Y %I:%M:%S %p',
		level=log.INFO
	)
elif log_message=="debug":
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

#parse arguments
################################################################################
parser = argparse.ArgumentParser()
parser.add_argument('-outdir', dest="outdir", help="Output directory for prism_clinvar files. Will be /storage1/shared/data/prism_clinvar/[uniprotID[0:2]]/[uniprotID[2:4]]/[uniprotID[4:6]] if not given.")
args = parser.parse_args()
################################################################################

#server	
path_base = '/storage1/hezscha/'
import_path_base = path_base + 'src/'

#logOUT = open(output_dir+'clinvar_extraction.log', 'w')

#local
#path_base = '/home/henrike/Documents/PD_AS/'
#import_path_base = '/home/henrike/Documents/PD_AS/src/'

try:
	# The insertion index should be 1 because index 0 is this file
	#sys.path.insert(1, '/groups/sbinlab/tiemann/repos/PRISM/prism/scripts/')  # the type of path is string
	sys.path.insert(1, import_path_base + 'PRISM/prism/scripts/')
	# because the system path already have the absolute path to folder a
	# so it can recognize file_a.py while searching 
	from PrismData import PrismParser, VariantData#VariantParser, VariantData

except (ModuleNotFoundError, ImportError) as e:
	logger.error("{} fail".format(type(e)))
	print(e)
else:
	logger.info("Import succeeded")


version = '001'

three_to_one = {
	"Ala": "A",
	"Val": "V",
	"Ile": "I",
	"Leu": "L",
	"Met": "M",
	"Phe": "F",
	"Tyr": "Y",
	"Trp": "W",
	"Ser": "S",
	"Thr": "T",
	"Asn": "N",
	"Gln": "Q",
	"Cys": "C",
	"Gly": "G",
	"Pro": "P",
	"Arg": "R",
	"His": "H",
	"Lys": "K",
	"Asp": "D",
	"Glu": "E",
	"Asx": "B",
	"Sec": "U",
	"Xaa": "X",
	"Glx": "Z",
	"Ter": "*"
}

def make_default_outfolder(uniprot_ID, base_prism_dir = '/storage1/shared/data/'):
#def make_default_outfolder(uniprot_ID, base_prism_dir = '/storage1/hezscha/pos_spec_prism_files/results/'):
	#sometimes we cannot find the unirpot ID to a transcript. Put them in the none folder then
	if not uniprot_ID == 'None':
		if not os.path.exists(base_prism_dir+'prism_clinvar/'+ uniprot_ID[0:2]):
			os.makedirs(base_prism_dir+'prism_clinvar/'+ uniprot_ID[0:2])
		if not os.path.exists(base_prism_dir+'prism_clinvar/'+ uniprot_ID[0:2]+ '/' + uniprot_ID[2:4]):
			os.makedirs(base_prism_dir+'prism_clinvar/'+ uniprot_ID[0:2]+ '/' + uniprot_ID[2:4])
		if not os.path.exists(base_prism_dir+'prism_clinvar/'+ uniprot_ID[0:2]+ '/' + uniprot_ID[2:4] + '/' + uniprot_ID[4:6]):
			os.makedirs(base_prism_dir+'prism_clinvar/'+ uniprot_ID[0:2]+ '/' + uniprot_ID[2:4]+ '/' + uniprot_ID[4:6])
	
		return(base_prism_dir+'prism_clinvar/'+ uniprot_ID[0:2]+ '/' + uniprot_ID[2:4]+ '/' + uniprot_ID[4:6])
	else:
		if not os.path.exists(base_prism_dir+'prism_clinvar/NA/'):
			os.makedirs(base_prism_dir+'prism_clinvar/NA/')
		return(base_prism_dir+'prism_clinvar/NA/')

def get_var_name(input_string):
	res = re.search('[(]p[.]([A-Za-z]+)(\d+)([A-Za-z]+)[)]',input_string)
	if res:
		var_name = three_to_one[res.groups()[0]] + res.groups()[1] + three_to_one[res.groups()[2]]
		
		#there is a special case of a read-through (stop codon to coding) variant that I want to omit because they don't fit in the prism framework since they are on residue positions outside the sequence stated in the header. They look like this: Ter552Cys , i.e. *552C after parsing to 1 letter
		if var_name.startswith('*'):
			return('NA')
		else:
			return(var_name)
	#if the variant is not one we're interested in, i.e. not a substitution there will be no regex match and we will skip to the next line in the clinvar flat file in the code that evaluates whether the return from get_var_name == 'NA'
	else:
		#print('Unable to ')
		#sys.exit()
		return('NA')
	 
def get_transcript(input_string):
	res = re.search('([\w.]+)(?:[(][\w\-\_]+[)])*:',input_string)
	if res:
		#Are we interrested in early termination?
		
		# ~ return(res.groups()[0].split('.')[0])
		return(res.groups()[0])
	else:
		return('NA')

def get_ncbi_prot_seq(tID):
	#https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly
	#use efetch with rettype=fasta_cds_aa
	
	server = "https://eutils.ncbi.nlm.nih.gov"
	ext = "/entrez/eutils/efetch.fcgi?db=nuccore&id="+tID+"&rettype=fasta_cds_aa&retmode=text&api_key=fe42ab6b5e8d76a7de4d795a849a30c54f08"
	
	while True:
	
		try:
			r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})
			r.raise_for_status()
		except HTTPError as e:
			status_code = r.status_code
			print('While looking up protein sequence for', tID, 'the following error occured',status_code)
		
		else:
			prot_header = r.text.split('\n')[0]
			prot_seq = ''.join(r.text.split('\n')[1:])
			sleep(1)
			return(prot_seq)

def get_ncbi_organism(tID):
	server = "https://eutils.ncbi.nlm.nih.gov"
	ext = "/entrez/eutils/esummary.fcgi?db=nuccore&id="+tID+"&retmode=json&api_key=fe42ab6b5e8d76a7de4d795a849a30c54f08"
	
	while True:
	
		try:
			r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
			r.raise_for_status()
		except HTTPError as e:
			status_code = r.status_code
			print('While looking up', tID, 'on esummary the following error occured',status_code)
		
		else:
			if r.text == '':
				sleep(1)
				continue #this should return to the top of the while loop and therefore trigger a new request
			
			else:
				decoded = r.json()
			
				if 'error' in decoded.keys():
					print("While looking up transcript", tID, "the following error occurred: ", decoded['error'], '\nRe-trying after 5 seconds.')
					sleep(5)
					continue
				
				else:
					GI = decoded['result']['uids'][0]
					seq_name = decoded['result'][GI]['title']
					taxid = str(decoded['result'][GI]['taxid'])
					sleep(1)
					break
	
	#now lookup the species tothe taxid for the prism file header
	ext = "/entrez/eutils/esummary.fcgi?db=taxonomy&id="+taxid+"&retmode=json&api_key=fe42ab6b5e8d76a7de4d795a849a30c54f08"
	
	while True:
	
		try:
			r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
			r.raise_for_status()
		except HTTPError as e:
			status_code = r.status_code
			print('While looking up', taxid, 'on taxonomy the following error occured',status_code)
		
		else:
			if r.text == '':
				sleep(1)
				continue #this should return to the top of the while loop and therefore trigger a new request
			
			else:
				dec = r.json()
				
				if 'error' in dec.keys():
					print("While looking up taxID", taxid, "the following error occurred: ", dec['error'], '\nRe-trying after 5 seconds.')
					sleep(5)
					continue
					
				elif not taxid in dec['result'].keys():
					print("While looking up taxID", taxid, "something else went wrong.\nRe-trying after 5 seconds.")
					sleep(5)
					continue
				
				else:
					organism = dec['result'][taxid]['scientificname']
					return(organism)

def process_clinsig(input_string):
	#https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/
	
	clinsig_category = {
	'Benign': 'benign',
	'Likely benign': 'benign',
	'Benign/Likely benign': 'benign',
	'Pathogenic': 'pathogenic',
	'Likely pathogenic': 'pathogenic',
	'Pathogenic/Likely pathogenic': 'pathogenic',
	'Uncertain significance': 'VUS',
	'Affects':'other',
	'association' : 'VUS', #this is a hard one but the documentation says 'For variants identified in a GWAS study and further interpreted for their clinical significance.' and except for 6 entries there is no further interpretation 
	'association not found': 'NA',
	'no interpretation for the single variant': 'NA',
	'not provided': 'NA',
	'conflicting data from submitters': 'conflict',
	'Conflicting interpretations of pathogenicity': 'conflict',
	'drug response': 'other',
	'other': 'other',
	'protective': 'other',
	'risk factor': 'other',
	}
	
	#I see a lot of double assignments like this: Affects, association | Pathogenic/Likely pathogenic, drug response | Uncertain significance, association
	#but the first named seems to me the more relevant, and you never seem thse the other way around, i.e. drug response, pathogenic . So let's use the first part only 
	input_string = input_string.split(',')[0]
	
	return(clinsig_category[input_string])

def process_clinrev(input_string):
	#https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/
	clinrev_category = {
	'no assertion criteria provided': '0',
	'no assertion provided': '0',
	'no interpretation for the single variant': '0',
	'criteria provided, single submitter': '1',
	'criteria provided, conflicting interpretations': '1',
	'criteria provided, multiple submitters, no conflicts': '2',
	'reviewed by expert panel': '3',
	'practice guideline': '4'
	}
	return(clinrev_category[input_string])

def check_clinvar_varname(cID):
	server = "https://eutils.ncbi.nlm.nih.gov"
	ext = "/entrez/eutils/esummary.fcgi?db=clinvar&id="+cID+"&retmode=json&api_key=fe42ab6b5e8d76a7de4d795a849a30c54f08"
	
	#debug
	#print("query:", server+ext)
	#debug
	
	while True:
		try:
			r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
			r.raise_for_status()
		except HTTPError as e:
			status_code = r.status_code
			print('While looking up', cID, 'on clinvar the following error occured',status_code)
			if status_code == 500:
				print('Internal server error occured, trying again after 5 seconds.')
				sleep(5)
				continue
			
		else:
			#debug
			# ~ print(cID, '*', r.text, '*', file = sys.stderr)
			#print(cID, '*', repr(r), '*', file = sys.stderr)
			#debug
			
			#for whatever arcane reason the return of requests is sometimes randomly an empty string even though the raise_for_status returned code 200 = ok
			#but the decode fails since we're trying to get a json from an empty string
			
			if r.text == '':
				sleep(1)
				continue #this should return to the top of the while loop and therefore trigger a new request
			
			else:
				dec = r.json()
			
				if cID in dec['result']:
					if 'protein_change' in dec['result'][cID]:
						prot_change = dec['result'][cID]['protein_change']	
						sleep(1)
						return(prot_change)
					elif 'error' in dec['result'][cID]:
						print('the following error occured for clinvar var', cID, ':', dec['result'][cID]['error'], file = sys.stderr)
						return('')
					else:
						print('Something went wrong with clinvar var', cID, '. No protein change field is available.', file = sys.stderr)
						print(dec['result'][cID], file = sys.stderr)
						sys.exit()	
				else:
					#else we sleep 5 s and try again
					sleep(5)
					continue	
	

def write_prism(metadata, dataframe, prism_file, comment=''):
	variant_dataset = VariantData(metadata, dataframe)
	parser = PrismParser()
	parser.write(prism_file, variant_dataset, comment_lines=comment)

def read_from_prism(primsfile):
	parser = PrismParser()
	dataframe = parser.read(primsfile).dataframe
	meta_data = parser.read_header(primsfile)
	return meta_data, dataframe

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

def write_prism_wrapper(symbol, t_info, symbol_entry):
	print(symbol)
	#go through all main transcripts
	
	#some symbols are in the clinvar flat file but none of the entries is eligible so the symbol_entry dict should be empty and we don't need a prism_clinvar file for them.
	if not write_prism_wrapper:
		print(symbol, 'had no eligible clinvar variants so no prism file will be made')
		return()
	
	for transcript in symbol_entry:
		print(transcript)
	
		#get all uniprot IDs. Need to know all 'neighbors' to the main transcript, i.e all other transcripts with same protein product
		#some transcripts don't have a corresponding uniprot ID or at least it's not in the idmapping file
		uniprot_set = set()
		try:
			uniprot_set.update(transcript_to_uniprot[transcript.split('.')[0]])
		except KeyError:
			print('No uniprot ID to', transcript, 'of', symbol, file = sys.stdout)

		for tn in t_info[transcript]['neighbors']:
			try:
				uniprot_set.update(transcript_to_uniprot[tn.split('.')[0]])
			except KeyError:
				print('No uniprot ID to', tn, 'of', symbol, file = sys.stdout)
				continue

		#if there were no uniprot IDs to be found we can't do anything 
		if len(uniprot_set) == 0:
			print('No uniprot IDs found for transcript', transcript, 'of gene', symbol, 'so no prism file will be made', file = sys.stdout)
			continue
			#return()

		#make dataframe
		output_df = pd.DataFrame.from_dict(symbol_entry[transcript], columns=prism_cols, orient = 'index')
		#reset the index to numeric for compatibility with prismData.py 
		output_df.index = range(len(output_df.index))
		
		#get the other infos necessary for the prism file: uniprot ID and transcript protein sequence for the header
		organism = get_ncbi_organism(transcript)
		prot_seq = t_info[transcript]['prot_seq']
		

		#debug
		#print(output_df.head())
		#print(output_df.columns)
		#output_df.to_csv(output_dir+symbol+'_'+transcript+'_df.csv', sep = '\t')
		print(prot_seq)
		print(uniprot_set)
		#debug
		
		#now, make one file for each uniprot ID
		for uniprot in uniprot_set:
			
			#path to output file
			# ~ prism_file = os.path.join(output_dir, f'prism_clinvar_XXX_{uniprot}-{transcript.split('.')[0]}.txt')
			if not args.outdir:
				output_dir = make_default_outfolder(uniprot_ID=uniprot)
			prism_file = os.path.join(output_dir, 'prism_clinvar_XXX_'+uniprot+'-'+symbol+'-'+transcript.split('.')[0]+'.txt')
			#print(prism_file)
				
			metadata = {
			"version": version,
			"protein": {
				"name": symbol,
				"organism": organism,
				"uniprot": uniprot,
				"sequence": prot_seq,
			},
			#"uniprot": {},
			"clinvar": {
				"representative transcript" : transcript,
				"transcripts with same protein seq": ','.join(t_info[transcript]['neighbors']),
				"all associated uniprot IDs": ','.join(list(uniprot_set))
			},
			"columns": columns_dic,
			}

			#evaluate membership in uniprot sets:
	
			membership = check_membership(uniprot, dbs_d)
			
			for dataset in membership:
				metadata['protein'][dataset] = membership[dataset]

			write_prism(metadata, output_df, prism_file, comment=comment)
			
			#for checking compatibility
			meta_data, dataframe = read_from_prism(prism_file)
			#print('check: read in prism file again')
			#print(dataframe.head())
			#print(dataframe.columns)
			
	#the regular return
	return()

#in cases like this gene: https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000197381;r=21:45073853-45226560 
#different transcripts correspond to different uniprot entries, but they have the same HGNC. So even though we are using NCBI transcripts, we should base each file on one transcript (from the 'name' field in clinvar), not the HGNC since that can mix different things together 

#start doing things

#source file and output directory
base_prism_dir = '/storage1/shared/data/'
#base_prism_dir = '/storage1/hezscha/pos_spec_prism_files/results/'


#proper file
clinvar_file = '/storage1/hezscha/data/clinvar/variant_summary.sort.symbol'

#test files
#clinvar_file = '/storage1/hezscha/data/clinvar/variant_summary.NF1'
# ~ clinvar_file = '/storage1/hezscha/data/clinvar/variant_summary.test'
# ~ clinvar_file = '/storage1/hezscha/data/clinvar/PGM1_clinvar.HGNC.txt'
# ~ clinvar_file = '/storage1/hezscha/data/clinvar/variant_summary.CBS'
# ~ clinvar_file = '/storage1/hezscha/data/clinvar/var_sum.2'

#for real data
#output_dir = '/storage1/hezscha/src/PRISM/prism/prism_clinvar/'
#for testing
if args.outdir:
	if args.outdir == '.':
		output_base = os.getcwd()
	else:
		output_base = args.outdir
	#no substructure
	output_dir = output_base
		
else:
	output_base = os.path.join(base_prism_dir, 'prism_clinvar/')
	#output_dir = make_default_outfolder(uniprot_id)

#read in dbs
dbs_d = read_uniprot_datasets()

transcript_to_uniprot = {}
with open('/storage1/hezscha/data/uniprot_idmapping/HUMAN_9606_idmapping.dat','r') as IN:
	for line in IN:
		ls = line.rstrip().split('\t')
		if ls[1] == 'RefSeq_NT':
			#format: P78563-2	RefSeq_NT	NM_001112.3 . Omit isoforms for uniprot IDs
			transc = ls[2].split('.')[0]
			uniprot = ls[0].split('-')[0]
			if transc in transcript_to_uniprot:
				transcript_to_uniprot[transc].append(uniprot)
			else:
				transcript_to_uniprot[transc] = [uniprot]

#read in list of already processed symbols 
proc_set = set()
with open(os.path.join(output_base,'processed_symbols.txt'), 'r') as IN:
	for line in IN:
		proc_set.add(line.rstrip())

#debug
print(proc_set)
#debug

#dict of lists
#entries_per_gene = {}
#top level dict with one subdict per symbol. Structure: symbols[symbol][transcr][varname]
#subsequent, different transcript with the same prot product will be noted under the same entry as the first transc with that prot product
symbol_entry = {}
#t_info will have info assoc with the transcripts such as their protein seq and their 'neighbors', which are transcripts with the same prot seq
t_info = {}
#match will tell for each transcript we have seen what their 'representative' or matching transcript is, the one that represents all transc with the same prot seq
#the difference to t_info is that t_info only has one entry per unique protein seq and match has one entry per occuring transcript. transcripts can also match themselves, i.e. if t1 is representative for protseq1 then match[t1] = t1
match = {}
#transcripts without protein seqs
no_prot_seq_for = []

#for writing the prism files:
prism_cols = ['variant','Clinvar_ID', 'Clinvar_signifiance', 'Clinvar_review_status']
columns_dic = {}
for elem in prism_cols[1:]:
	columns_dic[elem]= elem
comment = [ f"version {version} - {datetime.date(datetime.now())} - clinvar extract script",]

#there may have been an issue that the file wasn't closed properly since the script crashed and therefore nothing was written to it in append mode. Changed to opening and closing the file every time I write to it
# ~ proc_OUT = open(output_dir+'processed_symbols.txt','a')
symbol = ''
with open(clinvar_file, 'r') as CV_in:
	#header = CV_in.readline()
	for line in CV_in:
		if line.startswith('#'):
			continue
		
		ls = line.rstrip().split('\t')
		curr_symbol = ls[4]
		
		#there are some entries that don't have a symbol or HGNC or GeneID like 626260. I would skip these?
		
		#first, check if this symbol has aleady been processes previously
		if curr_symbol in proc_set:
			print(curr_symbol, 'is already processed')
			continue
		
		#second, check if the symbol has changed, i.e. we are done with reading entire for the prev symbol and if yes, process the previous symbol into a prism_clinvar file
		#when we read the first line, the prev symbol=symbol is ''
		if (symbol) and (curr_symbol) and (curr_symbol != symbol):
			
			#write prism file
			write_prism_wrapper(symbol, t_info, symbol_entry)
			
			#add symbol to the processed file
			proc_OUT = open(os.path.join(output_base+'processed_symbols.txt'),'a')
			print(symbol, file = proc_OUT)
			proc_OUT.close()
			
			#re-init symbol_entry and t_info since we don't need the prev symbol anymore
			symbol_entry = {}
			t_info = {}
			match = {}
			
		#in the end, whether we processed a prev symbol or not, now switch symbol to be curr_symbol to process the current entry
		symbol = curr_symbol	
		
		#get the info necessary to accept or reject the current entry
		crev = process_clinrev(ls[24]) 
		assembly = ls[16]
		var_type = ls[1]
		cID = ls[30]
		curr_symbol = ls[4]
		
		if not assembly == 'GRCh38': #check assembly == GRCh38
			continue
		if not var_type == 'single nucleotide variant': #I think all substiutions are this type
			continue
		if crev == '0': #check review status is at least 1 star
			continue
		
		var_name = get_var_name(ls[2])
		#some names don't describe an amino acid substitution (on the canon transcript) such as: NM_001112.4(ADARB1):c.1397-354A>G , clinvar ID: 977164
		if var_name == 'NA': 
			continue
		
		#More fun: this entry https://www.ncbi.nlm.nih.gov/clinvar/variation/193343/ on CLN5 names a A26E protein change in its clinvar name but in actuality no change exists. I don't know if it's an old entry and the name was just not updated or a mistake. But the webentry clearly states that there is no protein change and def not in the transcript stated in the entry name. So I'm now implementing a check that the variant name extracted from the title is correct
		#some vars have several protein changes like CBS https://www.ncbi.nlm.nih.gov/clinvar/variation/117/ . Check that one of them matches the title
		#Idealy we'd want to only get the protein change that matches to the transcript named in the title
		check_varname = check_clinvar_varname(cID)
		for var in check_varname.split(', '):
			if var_name == var:
				#we found a match, go on with the extraction
				break
		#if we didn't break none of the protein changes named in the esummary matches the title. 
		else:
			print('Entry',ls[2], ', ID:', cID,'for symbol', symbol, 'states protein change', var_name, 'but this is incorrect according to esummary which returns', check_varname, '. The entry has therefore been skipped.', file = sys.stderr)
			continue
		
		#check done, get other infos and save the entry
		csig = process_clinsig(ls[6])
		#keeping version numbers because its necessary to identify the correct protein seq to the transcript
		transcript = get_transcript(ls[2])
		if transcript == 'NA':
			print('Could not extract transcript for', symbol, 'in', line)
			sys.exit()
		
		#start saving procedure
		#debug
		if not symbol_entry: 
			print('starting extraction for', symbol)
			#debug
		
		#sometimes we cannot get a protein sequence, i.e. when sequences have been removed. In that case, we can't process the entry, go to next line. A note should have been written in the log about this the first time we tried to get a prot_seq for this transcript
		if transcript in no_prot_seq_for:
			continue
		
		if not transcript in symbol_entry:
			
			#debug
			print(transcript, 'not in keys for this symbol')
			#debug
			
			#the first thing to do is find out if this transcript matches to another transcript that already has an entry in symbol_entry, and if so switch to that
			if transcript in match:
				#debug
				print('match for', transcript, 'is', match[transcript])
				#debug
				transcript = match[transcript]
			
			#else we have not seen this transcript before and we need to get the seq
			else:
				#debug
				print('getting protein seq for', transcript)
				#debug
				prot_seq = get_ncbi_prot_seq(transcript)
				
				#sometimes we cannot get a protein sequence, i.e. when sequences have been removed
				if prot_seq == 'NA':
					no_prot_seq_for.append(transcript)
					continue #this should go to the next line in the clinvar file because we're still inside the for line in cvIN
				
				#now we go through t_info which has one entry per unique prot seq and see if the seq of our current transcript matches to any of those
				for t in t_info:
					if prot_seq == t_info[t]['prot_seq']:
						#debug
						print('match between', transcript, 'and', t)
						#debug
						
						#register match of current transcript 'transcript' and known transcript 't'
						#also register transcript as a neighbor to 't'
						match[transcript] = t
						t_info[t]['neighbors'].append(transcript)
						#lastly, switch from current transcript to its match
						transcript = t
						break
				
				#else we have checked against all known unique protseq transcripts of this symbol and the sequence is different, so the current transcript has a different protein seq then and gets its own entry in t_info and matches itself
				else:
					#debug
					print('no match was found for transcript', transcript)
					#debug
					
					t_info[transcript] = {}
					t_info[transcript]['prot_seq'] = prot_seq
					t_info[transcript]['neighbors'] = []
					symbol_entry[transcript] = {}
					match[transcript] = transcript
		
		#debug
		print(transcript)
		#debug
		
		#do this in any case
		#check if this var_name already exists. several entries can refer to the same protein level change (since they are defined on DNA level changes). Accumulate these	
		if var_name in symbol_entry[transcript]:
			#check that the significance/interpreation is the same. If yes, just add the clinvar ID to the existing one. Otherwise add ID and change significance to 'conflict'
			#the entry must have at least one star, otherwise it wouldn't pass the check above
			if csig == symbol_entry[transcript][var_name][2]:
				symbol_entry[transcript][var_name][1] += ','+cID
				symbol_entry[transcript][var_name][3] += ','+crev
				
			else:
				symbol_entry[transcript][var_name][1] += ','+cID
				symbol_entry[transcript][var_name][2] = 'conflict'
				symbol_entry[transcript][var_name][3] += ','+crev
		else:	
			symbol_entry[transcript][var_name] = [var_name, cID, csig, crev]

#in the end, process the last symbol
if symbol:
	#write prism file
	write_prism_wrapper(symbol, t_info, symbol_entry)

	#add symbol to the processed file
	proc_OUT = open(os.path.join(output_base+'processed_symbols.txt'),'a')
	print(symbol, file = proc_OUT)
	proc_OUT.close()				

# ~ proc_OUT.close()

