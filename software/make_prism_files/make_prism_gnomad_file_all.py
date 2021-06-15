
#UNDER CONSTRUCTION!!!
#make prism_gnomad files by requesting a uniprot ID and make one file for each unique protein product of this uniprot ID by comparing the protein products of all Ensembl transcripts associated to the uniprot 

#v04: only use variant data with AN > 100 (both for exome and genome samples)
#v03: deal with server request errors from looking up Ensembl deprecated transcripts that are still listed on the uniprot page, and possibly finding the successor. Example: ENST00000620120, status: retired
#v0.2: new version including correct estimation of allel frequency

import argparse
#this is the data for relase 100. The current relase server is https://.rest.ensembl.org
#e100 doesn't work, need to try it in brower and see what it evaluates to
ensembl_server = "https://apr2020.rest.ensembl.org"
ensembl_release = '100'
import requests #url requests
from requests import HTTPError
import sys
import re
import os
from datetime import datetime
from numpy.random import default_rng #for sampling pseudocounts
from numpy import std, sqrt
import pandas as pd
# ~ from vcf_parser_func_v09_revep import clinvar_lookup, get_uniprot, get_prot_seq
from vcf_parser_func_v09_revep import get_uniprot, get_prot_seq
#from write_prism_header_full_chrom_v06 import write_prism_file

#local paths for HZ. Change this to your paths
path_base = '/storage1/hezscha/'
import_path_base = path_base + 'src/'
sys.path.insert(1, import_path_base + 'PRISM/prism/scripts/')
from PrismData import PrismParser, VariantData#VariantParser, VariantData

#parse arguments
################################################################################
parser = argparse.ArgumentParser()
parser.add_argument('-e', dest="extract", choices=['mis', 'syn', 'complex', 'all'], default = 'mis', help="Type of variants to write files for.")
parser.add_argument('-m', dest="mode", choices=['overwrite', 'leave'], default = 'leave', help="What do when the output file already exists. Leave (default) or overwrite")
parser.add_argument('-tmp_folder', dest='tmp_folder', help="Where the temporary step 1 files are stored. Default location is /storage1/hezscha/gnomad_to_prism_parser/step1_files/+args.chromosome/")
parser.add_argument('-out_folder', dest='out_folder', help="Where the output files should be written. Default location is /storage1/shared/data/prism_gnomad/uniprot[0:2]/unipro[2:4]/uniprot[4:6]/")
parser.add_argument('-d', dest="debug", action='store_true', help="Debug mode.")
parser.add_argument('-uniprot', dest='uniprot', required=True, help="Uniprot ID to get transcripts and make prism_gnomad files for.")
parser.add_argument("-v","--verbose", action="count", default=0, help="Level of output, default zero is no output")

args = parser.parse_args()
################################################################################

#0: setup
exome_version = '2.1.1 (exome data)'
genome_version = '3.0 (genome data)'
parser_version = '001' #increase this if you change the parser significantly

#ex_g_fields are specific to the exome and genome data. F.x. it matters if an allel count was in the exomes or genomes, but position of a variant in the protein and on the assembly, its consequence and clinvar entry are irrelevant of the data source
#common_fields are the same for the same var irrespective of whether data come from exomes or genomes.
#nucl_spec_fields are specific to a variant on DNa level. basepos because changes at a different baseposition (1 or 2 away) can still lead to the same prot variant
ex_g_fields = ['AC', 'AN', 'Homozygotes', 'Codons']
common_fields = []#['Protein_position'] #we don't need the variant's position in the protein anymore because the prism parser calculates depth ect. But we need a variant key for the df later
#the basepos is not included in the other fields because it is usually shared btw all instances of a protein var but might not be in specific cases like R coded by CGG and AGA, so one nucl change could be X->C(pos1) and one could be X->A(pos3) 
target_fields = ex_g_fields + common_fields #+ ['basepos']

#the fields/columns to print in the final prism file and their descriptions for the header
print_fields = ['AF_tot', 'std_AF_tot', 'AF_ex', 'std_AF_ex', 'AC_ex', 'AN_ex', 'AF_g', 'std_AF_g', 'AC_g', 'AN_g', 'Codons_ex', 'Codons_g', 'Homozygotes_ex', 'Homozygotes_g']
# ~ if args.extract == 'mis' or args.extract == 'syn':
	# ~ print_fields = ['AF_tot', 'std_AF_tot', 'AF_ex', 'std_AF_ex', 'AC_ex', 'AN_ex', 'AF_g', 'std_AF_g', 'AC_g', 'AN_g', 'Codons_ex', 'Codons_g', 'Homozygotes_ex', 'Homozygotes_g']
# ~ else:
	# ~ print('Sorry, right now it only works for missense and synonymous variants!')
	# ~ sys.exit()

descr = {
	'AF_tot' : 'Frequency estimate of the variant allele from exome and genome data combined (if present in both, otherwise just AF_ex or AF_g). Calculated as (AC*_ex + AC*_g) / (max(AN_ex) + max(AN_g)) where AC*_ex = sum(AF_ex_i) * max(AN_ex)', 
	'std_AF_tot' : 'Error estimate on AF_tot from error propagation',
	'AF_ex' : 'Variant allel frequency in exomes data',
	'std_AF_ex' : 'Error estimate on AF_ex derived from the standard deviation of AF when adding pseudocounts',
	'AC_ex' : 'Counts of variant allel in exomes data',
	'AN_ex' : 'Total allel counts at this site in exomes data',
	'AF_g' : 'Variant allel frequency in genomes data',
	'std_AF_g' : 'Error estimate on AF_g derived from the standard deviation of AF when adding pseudocounts',
	'AC_g' : 'Counts of variant allel in genomes data',
	'AN_g' : 'Total allel counts at this site in genomes data',
	'Codons_ex' : 'Change on codon level (VEP) in exome data',
	'Codons_g' : 'Change on codon level (VEP) in exome data',
	'Homozygotes_ex' : 'Number of individuals homozygous for the variant in exome data',
	'Homozygotes_g' : 'Number of individuals homozygous for the variant in genome data',
}

#write logfile

#logfile = '/storage1/hezscha/gnomad_to_prism_parser/log/'+ 'gnomad_parse_' + args.uniprot +  '_' + str(datetime.now()).split()[0] + '_' + re.sub(':|\.', '_',str(datetime.now()).split()[1]) + '.log'
#logOUT = open(logfile, 'w')
#print('Parsing files for Uniprot ID ', args.uniprot, '. A log will be written to ', logfile, ' that mentions issues like several uniprot IDs for a transcript ID.', sep = '', file = sys.stderr)

#functions
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


def make_default_outfolder(uniprot_ID):
	#sometimes we cannot find the unirpot ID to a transcript. Put them in the none folder then
	if not uniprot_ID == 'None':
		if not os.path.exists('/storage1/shared/data/prism_gnomad/'+ uniprot_ID[0:2]):
			os.makedirs('/storage1/shared/data/prism_gnomad/'+ uniprot_ID[0:2])
		if not os.path.exists('/storage1/shared/data/prism_gnomad/'+ uniprot_ID[0:2]+ '/' + uniprot_ID[2:4]):
			os.makedirs('/storage1/shared/data/prism_gnomad/'+ uniprot_ID[0:2]+ '/' + uniprot_ID[2:4])
		if not os.path.exists('/storage1/shared/data/prism_gnomad/'+ uniprot_ID[0:2]+ '/' + uniprot_ID[2:4] + '/' + uniprot_ID[4:6]):
			os.makedirs('/storage1/shared/data/prism_gnomad/'+ uniprot_ID[0:2]+ '/' + uniprot_ID[2:4]+ '/' + uniprot_ID[4:6])
	
		return('/storage1/shared/data/prism_gnomad/'+ uniprot_ID[0:2]+ '/' + uniprot_ID[2:4]+ '/' + uniprot_ID[4:6])
	else:
		if not os.path.exists('/storage1/shared/data/prism_gnomad/None/'):
			os.makedirs('/storage1/shared/data/prism_gnomad/None/')
		return('/storage1/shared/data/prism_gnomad/None/')

def est_AF_error(AC, AN, lam=5.0, n_iter=100):
	AC = int(AC)
	AN = int(AN)
	rng = default_rng()
	p_counts = rng.poisson(lam=lam,size=n_iter)
	
	AF_new = []
	for i in range(n_iter):
		AC_new = AC + p_counts[i]
		AN_new = AN + p_counts[i]
		AF_new.append(float(AC_new)/float(AN_new))
		
	return(std(AF_new))

#superseeded by check_transcript 
# ~ #v03:deal with server request errors from looking up Ensembl deprecated transcripts that are still listed on the uniprot page, and possibly finding the successor. Example: ENST00000620120, status: retired

def check_transcript(t_ID,ensembl_server=ensembl_server):
	ext = "/lookup/id/"+t_ID+"?" 
	
	#https://stackoverflow.com/questions/19342111/get-http-error-code-from-requests-exceptions-httperror
	try:
		response = requests.get(ensembl_server+ext, headers={ "Content-Type" : "application/json"})
		response.raise_for_status()
	except HTTPError as e:
		# Need to check its an 404, 503, 500, 403 etc.
		status_code = e.response.status_code
		#print(status_code)
		if status_code == 400:
			print('Bad Request for '+ensembl_server+'/lookup/id/'+t_ID+'. Possibly the transcript has been retired and is not in release' + ensembl_release +'for GRCh38.', file = sys.stderr)#, file = logOUT)
			return('NA')
	
	else:
	
		decoded = response.json()
		if decoded['seq_region_name'] == 'X' or decoded['seq_region_name'] == 'Y':
			decoded['chromosome'] = decoded['seq_region_name']
		else:
			try: 
				decoded['chromosome'] = str(int(decoded['seq_region_name'])) #we don't actually want an integer, just to check if it can be cast to one
			except:
				decoded['chromosome'] = 'not_main'
				
		return decoded		

#####################

#eval commandline args
if args.out_folder:
	if args.out_folder == '.':
		out_folder = os.getcwd()
	else:	
		out_folder = args.out_folder
else:
	out_folder = make_default_outfolder(args.uniprot)

var_type = args.extract
if args.verbose > 1:
	print('Extracting for:', var_type)

#translate between the var_type names used by args and reported in the step 1 files
short_name = {'missense_variant' : 'mis', 'synonymous_variant' : 'syn'}
if not (args.extract == 'mis' or args.extract == 'syn'):
	print('Sorry, right now it only works for missense and synonymous variants!')
	sys.exit()

args.uniprot = args.uniprot.split('-')[0]

#check if the folder ends in /, otherwise add one
if not out_folder.endswith('/'):
	out_folder += '/'
#also check if outfolder exists
if not os.path.exists(out_folder):
	print('folder', out_folder, "doesn't exist!", file = sys.stderr)
	sys.exit()

#tmp_folder is the folder that has the extracted data from the vcf files per transcript
tmp_folder_base = args.tmp_folder if args.tmp_folder else '/storage1/hezscha/gnomad_to_prism_parser/step1_files/'
dbs_d = read_uniprot_datasets()

#######################################################################
#1. make a list of coding transcripts with unique protein products. We will make one file for each

#get all transcripts to the uniprot
transcript_info = {} #dict of geneIDs and locations. The keys are the transcripts and serve as the transcript list to be iterated over
all_transcripts = []
r = requests.get('https://www.uniprot.org/uniprot/'+args.uniprot+'.txt')
if not r.ok:
	print('Lookup of uniprot ID has failed with:')
	r.raise_for_status()
	sys.exit()
else:
	r = r.text
	
for line_u in r.split('\n'):
	ls = line_u.rstrip().split()

	if ls and ls[0] == 'DR' and ls[1] == 'Ensembl;':
		t_ID = ls[2][:-1] #omit the ;
		ensg = ls[4][:-1]
		#look up the location and skip transcripts of genes not on the main assembly
		ens_ret = check_transcript(t_ID)
		if not ens_ret['chromosome'] == 'not_main' and ens_ret['biotype'] == 'protein_coding':
			if args.verbose >= 1:
				print('transcript:', t_ID, 'chromosome:', ens_ret['chromosome'])
			
			#record basic info about this transcript
			all_transcripts.append(t_ID)
			transcript_info[t_ID] = {}
			transcript_info[t_ID]['ensg'] = ensg
			transcript_info[t_ID]['chrom'] = ens_ret['chromosome']
			transcript_info[t_ID]['gene_name'] = ens_ret['display_name'].split('-')[0]
			transcript_info[t_ID]['neighbors'] = [t_ID] #transcripts with the same protein product are neighbors. By definition (made by me) each transcript is also it's own neighbor because this is the field I will use to write all associated transcript IDs into the prism file

if len(transcript_info) == 0:
	print('No eligible transcripts were found.')
	sys.exit()

#v03
#go through the protein sequences of all transcripts and only keep one transcript per unique seq
#all_transcripts = list(transcript_info.keys())
keep_list = {}
#the first transcript goes automatically into the dict of transcripts to keep. Every other transcript is compared to every transcript in the keep_list. 
#the keep_list has one key for each transc to keep and the value is all the transcripts with the same seq as the one named in the key. So minimally one, the key transcript itself
keep_list[all_transcripts[0]] = [all_transcripts[0]]
transcript_info[all_transcripts[0]]['prot_seq'] = get_prot_seq(all_transcripts[0], ensembl_server, prefix='Initial lookup')
#replace U (selenocysteine) with X in first transcript and make a note. Needs to be defined per transcript here, not globally!
transcript_info[all_transcripts[0]]['nsr'] = []
if 'U' in transcript_info[all_transcripts[0]]['prot_seq']:
	repl_seq = ''
	#need also the positions and all instances and I don't feel up for regex
	for pos,char in enumerate(transcript_info[all_transcripts[0]]['prot_seq']):
		if char == 'U':
			repl_seq += 'X'
			transcript_info[all_transcripts[0]]['nsr'].append('U'+str(pos+1))
		else:
			repl_seq += char
	
	transcript_info[all_transcripts[0]]['prot_seq'] = repl_seq

#then go through all other found transcript IDs and compare their sequences to the ones already in the keep_list. If we find a match, add the current tID to that entry in keep_list, otherwise if we go through the whole list and this sequence was not seen before make a new entry
for t_ID in all_transcripts[1:]:
	seq = get_prot_seq(t_ID, ensembl_server, prefix='Initial lookup')
	for key_transc in keep_list:
		if seq == transcript_info[key_transc]['prot_seq']:
			keep_list[key_transc].append(t_ID)
			transcript_info[key_transc]['neighbors'].append(t_ID)
			break
	#if we didn't break there is no transc in the keep_list with the same seq as the current transc, so add the current transc to the keep_list
	else:
		keep_list[t_ID] = [t_ID]
		transcript_info[t_ID]['prot_seq'] = seq
		#replace U (selenocysteine) with X and make a note. Needs to be defined per transcript here, not globally!
		transcript_info[t_ID]['nsr'] = []
		if 'U' in transcript_info[t_ID]['prot_seq']:
			repl_seq = ''
			#need also the positions and all instances and I don't feel up for regex
			for pos,char in enumerate(transcript_info[t_ID]['prot_seq']):
				if char == 'U':
					repl_seq += 'X'
					transcript_info[t_ID]['nsr'].append('U'+str(pos+1))
				else:
					repl_seq += char
			
			transcript_info[t_ID]['prot_seq'] = repl_seq

if args.verbose > 1:
	print('keep_list:', keep_list)

#now reduce the transcript_info dict to only the keys in the keep_list dict. I.e., we only need one entry per unique protein product
for item in set(all_transcripts) - set(keep_list.keys()):
	del transcript_info[item]
#print(transcript_info)

#print('Exiting for debug at line 224.')
#sys.exit()

#get file names of the step1 files. The thing is that the prefixes of the step1 files (raw gnomad data) look like this ENST00000315781_ENSG00000143740_Q5SQN1.mis but sometimes the uniprot ID is NA even though there is actually a uniprot ID since gnomad only contains uniprot info in the vep annotation for Swissprot. So we will identify the correct files in the common index from the ENST and store that prefix. Could also glob with the ENST but it might be more efficient to just read in the index since we have it.

#the way I understand biology all transcripts to a gene should be on the same chromosome, so I'm just getting the chrom of the first one in the keep list
idx_IN = open(tmp_folder_base+transcript_info[list(keep_list.keys())[0]]['chrom']+'/idx.txt', 'r')
for line_i in idx_IN:
	t_ID = line_i.split('_')[0]
	#this the current transcript one we are looking for?
	if t_ID in transcript_info:
		#transcript_info[t_ID]['tmp_file_prefix'] = line_i.rstrip()
		transcript_info[t_ID]['tmp_file_prefix'] = line_i.rstrip().split('.')[0]

idx_IN.close()

#sanity check: check that for each transcript if file was found in the index. some transcripts might have no variants, perhaps because they are strongly conserved and these will not be named in the index
rm_ls = []
for t_ID in transcript_info:
	if not 'tmp_file_prefix' in transcript_info[t_ID]:
		print(t_ID, 'has no gnomad variants.')
		rm_ls.append(t_ID)
		#print("Couldn't identify gnomad raw file for", t_ID)
		#sys.exit()

#remove transcripts that don't have vars from the dict
for item in rm_ls:
	del transcript_info[item]

if len(transcript_info) == 0:
	#print('has no coding transcripts or none on the main assembly or no step 1 file, i.e. no known variants.')
	print('No temp files were found for any of the transcript. Exiting.')
	sys.exit()

#last check before we start extracting:
if args.verbose > 1:
	print('transcript_info:')
	for t in transcript_info:
		print(t)
		print(transcript_info[t])

######################################################################
#2. open step1 files, two per transcript (exome and genome data, if exist), get data, write prism files

tmp_folder = tmp_folder_base + transcript_info[list(keep_list.keys())[0]]['chrom'] + '/'

#the below needs to be done once per transcript in transcript_info
for t_ID in transcript_info: 
	
	if args.verbose >= 1:
		print(t_ID)

	#process: since we don't know how many variants there are apriory, we can do as in make_prism_spliceai and create a dict of dicts while reading in, keyed to the var_name. Then conert to list of dict and cast that to df.
	#alternative, can try to just read in the step1 file as a df and assign var_names by parseHGVSp and then combine rows with same var names
	#"variants" is that dict of dicts, where each var has a subdict that has the var name as key: 
	#variants['M1A'] = {'AF_ex' = ???, 'AC_ex' = ???, ect}

	variants = {} #storing the properties of the variants like AC, AN, homozygotes ect
	pos = {} #a dict storing the position of each field we search for, or False if the field is not available
	pos['var_name'] = 0 #the variant's prism name is always the first field
	pos['AF'] = 1 #needs to be float and added instead of appended, so special treatment. 
	#mut_pos_d = {}
	exomes_exist = 0

	#go through all transcripts with the same protein product and find out if any of them is the canonical
	transcript_info[t_ID]['is_canonical'] = ''
	for neighbor_t in transcript_info[t_ID]['neighbors']:
		ext = "/lookup/id/"+neighbor_t+"?" 
		r_enst = requests.get(ensembl_server+ext, headers={ "Content-Type" : "application/json"})
		if not r_enst.ok:
			print(r_enst.raise_for_status())
			sys.exit()
		else:
			decoded = r_enst.json()
			if decoded['is_canonical'] == 1:
				transcript_info[t_ID]['is_canonical'] = neighbor_t
				break

	#check if output file exists and then act according to selected mode (leave, overwrite or ask)
	outpref = 'prism_gnomad_' + parser_version + '_' + args.uniprot + '-' + t_ID + '-'+ transcript_info[t_ID]['ensg']
	prism_file_name = outpref +'.'+var_type+'.txt'
	prism_file = os.path.join(out_folder, prism_file_name)
	#prism_file = os.path.join(out_folder, 'prism_gnomad_' + parser_version + '_' + args.uniprot + '-' + chosen_t['ID'] + '-'+ chosen_t['ensg'] + '.'+var_type+'.txt')

	if os.path.exists(prism_file):
		if args.mode == 'leave':
			print(prism_file, 'exists')
			continue
		elif args.mode == 'overwrite':
			pass	

	#exome file
	try:
		ex_IN = open(tmp_folder + transcript_info[t_ID]['tmp_file_prefix'] + '.' + var_type + '.exomes.tmp','r')
		exomes_exist = 1
		if args.verbose >= 1:
			print('Opened ', tmp_folder + transcript_info[t_ID]['tmp_file_prefix'] + '.' + var_type + '.exomes.tmp', sep = '')

		for line in ex_IN:
			#parse the header
			if line.startswith('#'):
				if line.startswith('#gnomad_version:'):
					continue #we assemble the gnomad version from exome and genome versions
				elif line.startswith('#transcript_ID:'):
					if not t_ID == line.rstrip().split()[-1]:
						print("Mismatch between transcript ID in step1 file:", line.rstrip().split()[-1], "and transcript ID currently being processed:", t_ID)
				elif line.startswith('#gene_ID:'):
					if not transcript_info[t_ID]['ensg'] == line.rstrip().split()[-1]:
						print("Mismatch between gene ID in step1 file:", line.rstrip().split()[-1], "and gene ID found of currently processed transcript:", transcript_info[t_ID]['ensg'])		
				elif line.startswith('#symbol:'):
					if not transcript_info[t_ID]['gene_name'] == line.rstrip().split()[-1]:
						print("Mismatch between gene name (symbol) in step1 file:", line.rstrip().split()[-1], "and gene name of currently processed transcript:", transcript_info[t_ID]['gene_name'])	
				elif line.startswith('#uniprot_ID:'):
					#vep only gives swissprot IDs so the id given her can be NA. Just continue, we use the one requested via args.uniprot
					#I.e. ENST00000432504_ENSG00000142188_NA.missense.txt actually matches C9K0I4 though it didn't have in the gnomad VEP annotation (which I think only queries SwissProt?): 
					#check if the stated uniprot ID matches the queried and if not, complain (unless stated is NA, then replace)
					continue
					
				elif line.startswith('#canonical:'):	
					#info_d['is_canonical'] = line.rstrip().split()[-1]
					continue
				elif line.startswith('#extracted_variants:'):	
					filter_set = line.rstrip().split()[-1] #set(line.rstrip().split()[-1].split(','))
					if not short_name[filter_set] == var_type:
						print("This file is for", filter_set, "but requested was", var_type, "so something is wrong. Exit for now.")
						sys.exit()
					
					
				else:
					#the only other # line is the column header. Now find the positions of fields that have been requested
					for target_field in target_fields:
						pos[target_field] = False
						for (position, field) in enumerate(line.rstrip().split()):
							if field == target_field: 
								pos[target_field] = position
								break
			
					if args.verbose > 1:
						print(pos)
			
			
			#else it's a variant line. Save in the dict if not prev. seen, update the info if prev. seen and if parsing missense vars add clinvar and domain lookup
			else:
				fields = line.rstrip().split()
				var_name = fields[0] #the variant's prism name is always the first field
				if var_name.startswith('U'):
					var_name = 'X'+ var_name[1:]
				
				#v04: don't consider data with AN < 100
				if int(fields[pos['AN']]) < 100:
					continue #go to the next variant line without processing of saving this one
				
				if var_name in variants:
					#if this var was seen before, we need to:
					#a) sum AF's
					#b) aggregate info from other columns
					#v08:
					#c) regarding the error estimate: I now know the only correct thing to do is error propagation as shown/discussed here: https://stats.stackexchange.com/questions/71419/average-over-two-variables-why-do-standard-error-of-mean-and-error-propagation
					#d) update max AN
					#debug
					#print(var_name, variants[var_name]['AF'], type(variants[var_name]['AF']))
					#debug
					
					variants[var_name]['AF_ex'] += float(fields[pos['AF']]) #a)
					variants[var_name]['max_AN_ex'] = int(fields[pos['AN']]) if int(fields[pos['AN']]) > variants[var_name]['max_AN_ex'] else variants[var_name]['max_AN_ex'] #d)
					
					#v08:
					#c: error propagation
					new_error = est_AF_error(AC=fields[pos['AC']], AN=fields[pos['AN']])
					variants[var_name]['std_AF_ex'] = sqrt((variants[var_name]['std_AF_ex'])**2 + (new_error)**2)

					for target_field in ex_g_fields:
						#skip unavailable fields. We don't need to add an extra NA
						if not pos[target_field]:
							continue
						else:
							variants[var_name][target_field+'_ex'] += ',' + fields[pos[target_field]]		
				
				#else this var was not seen before, just save
				else:
					
					variants[var_name] = {}
					variants[var_name]['variant'] = var_name #need a key for that for the df later
					variants[var_name]['AF_ex'] = float(fields[pos['AF']])
					variants[var_name]['max_AN_ex'] = int(fields[pos['AN']])
					
					#v04: add error estimate
					variants[var_name]['std_AF_ex'] = est_AF_error(AC=fields[pos['AC']], AN=fields[pos['AN']])
					#get all other requested target fields
					for target_field in ex_g_fields:
						variants[var_name][target_field+'_ex'] = fields[pos[target_field]] if (pos[target_field] and fields[pos[target_field]]) else 'NA'
					
					#common fields are shared by all instances of the same variant, so we only save the info the first time we see a variant
					for target_field in common_fields:
						variants[var_name][target_field] = fields[pos[target_field]] if (pos[target_field] and fields[pos[target_field]]) else 'NA'
					
					#the basepos is not included in the other fields because it is usually shared btw all instances of a protein var but might not be in specific cases like R coded by CGG and AGA, so one nucl change could be X->C(pos1) and one could be X->A(pos3) 
					#variants[var_name]['basepos'] = fields[pos['basepos']]			

	#if there is no exomes file for this transript + var type combo, try to open the genomes file instead		 	
	except FileNotFoundError:
		if args.verbose >= 1:
			print('There is no ',tmp_folder + transcript_info[t_ID]['tmp_file_prefix'] + '.' + var_type + '.exomes.tmp', sep = '')
			print('Proceeding to genomes tmp file')
		
		pass	

	#and now read the tmp for the same transcript for the genome data
	genomes_exist = 0

	try:
		gen_IN = open(tmp_folder + transcript_info[t_ID]['tmp_file_prefix'] + '.' + var_type + '.genomes.tmp','r')
		#ex_IN = open(tmp_folder + transcript_info[t_ID]['tmp_file_prefix'] + '.' + var_type + '.exomes.tmp','r')
		genomes_exist = 1
		if args.verbose >= 1:
			print('Opened ', tmp_folder + transcript_info[t_ID]['tmp_file_prefix'] + '.' + var_type + '.genomes.tmp', sep = '')
		
		for line in gen_IN:
			if line.startswith('#'):
				#if there was no exome file, parse the header
				if not exomes_exist:
					if line.startswith('#gnomad_version:'):
						continue #we assemble the gnomad version from exome and genome versions
					elif line.startswith('#transcript_ID:'):
						if not t_ID == line.rstrip().split()[-1]:
							print("Mismatch between transcript ID in step1 file:", line.rstrip().split()[-1], "and currently processed transcript ID:", t_ID)
					elif line.startswith('#gene_ID:'):
						if not transcript_info[t_ID]['ensg'] == line.rstrip().split()[-1]:
							print("Mismatch between gene ID in step1 file:", line.rstrip().split()[-1], "and gene ID of currently processed transcript:", transcript_info[t_ID]['ensg'])		
					elif line.startswith('#symbol:'):
						if not transcript_info[t_ID]['gene_name'] == line.rstrip().split()[-1]:
							print("Mismatch between gene name (symbol) in step1 file:", line.rstrip().split()[-1], "and gene name of currently processed transcript:", transcript_info[t_ID]['gene_name'])	
					elif line.startswith('#uniprot_ID:'):
						#vep only gives swissprot IDs so the id given her can be NA. Just continue, we use the one requested via args.uniprot
						#I.e. ENST00000432504_ENSG00000142188_NA.missense.txt actually matches C9K0I4 though it didn't have in the gnomad VEP annotation (which I think only queries SwissProt?): 
						#check if the stated uniprot ID matches the queried and if not, complain (unless stated is NA, then replace)
						continue
						
					elif line.startswith('#canonical:'):	
						#info_d['is_canonical'] = line.rstrip().split()[-1]
						continue
					elif line.startswith('#extracted_variants:'):	
						filter_set = line.rstrip().split()[-1] #set(line.rstrip().split()[-1].split(','))
						if not short_name[filter_set] == var_type:
							print("This file is for", filter_set, "but requested was", var_type, "so something is wrong. Exit for now.")
							sys.exit()
					
					else:
						#the only other # line is the column header. Now find the positions of fields that have been requested
						for target_field in target_fields:
							pos[target_field] = False
							for (position, field) in enumerate(line.rstrip().split()):
								if field == target_field: 
									pos[target_field] = position
									break
				
						if args.verbose > 1:
							print(pos)	
		
			#else it's a variant line. Save in the dict if not prev. seen, update the info if prev. seen and if parsing missense vars add clinvar and domain lookup
			else:
				fields = line.rstrip().split()
				var_name = fields[0] #the variant's prism name is always the first field
				if var_name.startswith('U'):
					var_name = 'X'+ var_name[1:]
				
				#!we need to differentiate whether this var was seen before in the genomes or the exomes when we calc the AF_g ect
				
				#v04: don't consider data with AN < 100
				if int(fields[pos['AN']]) < 100:
					continue #go to the next variant line without processing of saving this one
				
				if var_name in variants:
					
					#was this variant prev seen in genomes?
					if 'max_AN_g' in variants[var_name]:

						variants[var_name]['AF_g'] += float(fields[pos['AF']]) #a)
						variants[var_name]['max_AN_g'] = float(fields[pos['AN']]) if int(fields[pos['AN']]) > variants[var_name]['max_AN_g'] else variants[var_name]['max_AN_g'] #d)
						
						#error propagation
						new_error = est_AF_error(AC=fields[pos['AC']], AN=fields[pos['AN']])
						variants[var_name]['std_AF_g'] = sqrt((variants[var_name]['std_AF_g'])**2 + (new_error)**2)

						for target_field in ex_g_fields:
							#skip unavailable fields
							if not pos[target_field]:
								continue
							else:
								variants[var_name][target_field+'_g'] += ',' + fields[pos[target_field]]	
					
					#else this var was seen in exomes before but not in genomes, so create the genomes entries
					else:
						
						variants[var_name]['AF_g'] = float(fields[pos['AF']])
						variants[var_name]['max_AN_g'] = float(fields[pos['AN']])
						variants[var_name]['std_AF_g'] = est_AF_error(AC=fields[pos['AC']], AN=fields[pos['AN']])
						
						for target_field in ex_g_fields:
							variants[var_name][target_field+'_g'] = fields[pos[target_field]] if (pos[target_field] and fields[pos[target_field]]) else 'NA'
						
				#else this var was not seen before, just save
				else:
					
					variants[var_name] = {}
					variants[var_name]['variant'] = var_name #need a key for that for the df later
					variants[var_name]['AF_g'] = float(fields[pos['AF']])
					variants[var_name]['max_AN_g'] = float(fields[pos['AN']])
					variants[var_name]['std_AF_g'] = est_AF_error(AC=fields[pos['AC']], AN=fields[pos['AN']])

					#get all other requested target fields
					for target_field in ex_g_fields:
						variants[var_name][target_field+'_g'] = fields[pos[target_field]] if (pos[target_field] and fields[pos[target_field]]) else 'NA'
				
					#common fields are shared by all instances of the same variant, so we only save the info the first time we see a variant
					for target_field in common_fields:
						variants[var_name][target_field] = fields[pos[target_field]] if (pos[target_field] and fields[pos[target_field]]) else 'NA'

	except FileNotFoundError:
		if args.verbose >= 1:
			print('Could not open ', tmp_folder + transcript_info[t_ID]['tmp_file_prefix'] + '.' + var_type + '.genomes.tmp', sep = '')
		pass



	#now we go through all variants again and calc AF_tot
	for var_name in variants:
		#debug
		#print(var_name)
		#debug
		if 'max_AN_ex' in variants[var_name] and 'max_AN_g' in variants[var_name]:
			#debug
			#print('exome and genome data on this var')
			#debug
			
			#AF_ex is already the sum of all AF in the exome data for this var since if we have seen the var before AF = AF + new_AF
			AC_prime_ex = variants[var_name]['AF_ex'] * variants[var_name]['max_AN_ex']
			AC_prime_g = variants[var_name]['AF_g'] * variants[var_name]['max_AN_g']
			variants[var_name]['AF_tot'] = (AC_prime_ex + AC_prime_g) / (variants[var_name]['max_AN_ex'] + variants[var_name]['max_AN_g'])
			variants[var_name]['std_AF_tot'] = sqrt((variants[var_name]['std_AF_ex'])**2 + (variants[var_name]['std_AF_g'])**2)
		#if there is exome data on this var but no genome data
		elif 'max_AN_ex' in variants[var_name]:
			#debug
			#print('only exome data on this var')
			#debug
			
			#the calc simplifies to:
			variants[var_name]['AF_tot'] = variants[var_name]['AF_ex']
			variants[var_name]['std_AF_tot'] = variants[var_name]['std_AF_ex']
			#if there was genome data for this transcript (though not for this specific variant) the genome data fields will still be requested during printing, so fill them up with NAs
			for field in ['AF_g', 'std_AF_g', 'AC_g', 'AN_g', 'Codons_g', 'Homozygotes_g']:
				variants[var_name][field] = 'NA'
			
		#else there is only genome data
		elif 'max_AN_g' in variants[var_name]:
			#debug
			#print('only genome data on this var')
			#debug
			
			variants[var_name]['AF_tot'] = variants[var_name]['AF_g']
			variants[var_name]['std_AF_tot'] = variants[var_name]['std_AF_g']
			#fill up exome data fields for printing:
			for field in ['AF_ex', 'std_AF_ex', 'AC_ex', 'AN_ex', 'Codons_ex', 'Homozygotes_ex']:
				variants[var_name][field] = 'NA'
		
		#else something is weird
		else:
			print(variants[var_name], file = sys.stderr)
		
		#debug
		#print(variants[var_name].keys())
		#debug		

	#cast dict of dicts to list of dicts and then df
	temp_lod = []
	for var_name in variants:
		temp_lod.append(variants[var_name])

	#lod from dod
	df = pd.DataFrame(temp_lod, columns = ['variant', 'AF_tot', 'std_AF_tot', 'AF_ex', 'std_AF_ex', 'AC_ex', 'AN_ex', 'AF_g', 'std_AF_g', 'AC_g', 'AN_g', 'Codons_ex', 'Codons_g', 'Homozygotes_ex', 'Homozygotes_g'])
			

	if exomes_exist:
		if genomes_exist:
			transcript_info[t_ID]['data_source'] = 'exomes, genomes'
			
		else:
			transcript_info[t_ID]['data_source'] = 'exomes'
	else:
		transcript_info[t_ID]['data_source'] = 'genomes'

	metadata = {
			"version": parser_version,
			"protein": {
				"name": transcript_info[t_ID]['gene_name'],
				"organism": 'Human',
				"sequence": transcript_info[t_ID]['prot_seq'],
				"uniprot": args.uniprot,
				"Ensembl transcript ID": ','.join(transcript_info[t_ID]['neighbors']),
				"is canonical": transcript_info[t_ID]['is_canonical'],
				"chromosome": transcript_info[t_ID]['chrom'],
				"transcript has variants in": transcript_info[t_ID]['data_source'],
				"non-standard residues" : ','.join(transcript_info[t_ID]['nsr'])
				
			},
			"gnomad": {
			"version": exome_version+', '+genome_version,
			"build": 'GRCh38',
			"ensembl_gene_id": transcript_info[t_ID]['ensg'],
			"extracted variant type": filter_set,
			},
			"columns": descr,
		}

	comment = [ "containing only variants with AN > 100",]

	#add membership in uniprot sets:
	membership = check_membership(args.uniprot, dbs_d)
	for dataset in membership:
		metadata['protein'][dataset] = membership[dataset]	

	#check if the dataframe is empty. In the example of /storage1/shared/data/prism_gnomad/J3/KQ/82/prism_gnomad_001_J3KQ82-ENST00000454978-ENSG00000169800.mis.txt there is only one var in the exome tmp file but it has AN=7, i.e. super shallow seq depth so it is skipped
	if df.empty:
		if args.verbose >= 1:
			print(t_ID, 'has no variants.')
		continue

	print('Writing prism file')
	if args.verbose >= 1:
		print('Output file path:\n', prism_file, sep = '')
	write_prism(metadata, df, prism_file, comment=comment)
	metadata, df_in = read_from_prism(prism_file)
	print('Prism file passed check')
