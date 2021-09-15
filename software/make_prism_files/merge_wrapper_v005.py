
#v005: re-add checking alignment since I still think we can only do the mismacht handling (swapping to target_seq WT) is the alignment is ungapped or only has leading or trailing gaps, due to uncertainty of which residues are equivalent in target_seq and other_seq(gnomad, uniprot) around the gaps otherwise.
#v004: 
#example: prism_gnomad_001_O14786-ENST00000265371-ENSG00000099250.mis.txt V179A
#4.1. mismatches between gnomad and target seq are always resolved if possible. Resolution is possible if the WT in the target_seq at pos j is a var at that pos j in gnomad. In that case we swap the gnomad var from gnomad_WT+pos+var(==target_WT) to target_WT+pos+gnomad_WT. 
#the AF of the gnomad WT is estimated as 1 - sum(AF(all_other_gnomad_aa_at_pos_j))
#This does not imply that we belive in the target_WT more. It is merely done to deliver an estimate of a gnomad freq for the clinvar var
#4.2 other vars at pos j are also swapped from gnomad_WT+pos+var(!=target_WT) to target_WT+pos+var which keeping the same allel freq. See trello board for why
#4.3 for mismatches between position specific files such as prism_uniprot and prism_netsurfp we also swap from uniprot_WT+pos+= to target_WT+pos+= and add a column indicating that this swap was done. This is because in those pos spec files aa+pos_j actually implies that feature x was observed at pos j with aa as the WT and this is no longer true when we swap. The extra col allows to disregard any var lines that include that prism_uniprot pos from further analysis if we think we cannot trust that feature x still applies to target_WT != uniprot_WT
#see trello list merging, card SNP inversion https://trello.com/b/mhbusoya/proteome-wide-stuff


#v003: use same alignment cutoffs as in the parser (identity and coverage) instead of ngaps and nsubs. Allow SNP transfer to different isoforms so long as the cutoffs are not violated.
#v002: add SNP handling to merge wrapper instead of doing it inside the parser

import argparse
import requests #url requests
from requests import HTTPError
import sys
import re
import os
import glob
import subprocess
import pandas as pd
import numpy as np
import textwrap
from time import sleep
from datetime import datetime
from Bio import Seq,SeqRecord,SeqIO,pairwise2,SubsMat
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo
from collections import namedtuple
import math
import shutil
#for random letters
from string import ascii_letters, digits
from random import choice

import_path_base = '/storage1/hezscha/src/'
sys.path.insert(1, import_path_base + 'PRISM/prism/scripts/')
#from PrismData_HZ_allow_SNP_v4 import PrismParser, VariantData
from PrismData import PrismParser, VariantData
sys.path.insert(1, import_path_base + 'helper_python_functions/')
from helpers import do_align,get_uniprot_seq,fasta_print

#parse arguments
################################################################################
#https://stackoverflow.com/questions/50021282/python-argparse-how-can-i-add-text-to-the-default-help-message
parser = argparse.ArgumentParser(epilog=textwrap.dedent('''\
         additional information:
         The options -uniprot and -ff are mutually exlusive. Use only one of them.
         '''))
#parser.add_argument('-all', dest="all", action='store_true', help="Run on all gnomad files")
parser.add_argument('--typeslog', dest="typeslog", action='store_true', help="Write a log file listing the components of every merge file for stats later.")
#parser.add_argument("-list", dest='filelist', help="Run only on this list of files")
parser.add_argument('-uniprot', dest="uniprot", help="The uniprot ID for which to merge files. The authorative sequence for merging is selected with following hierarchy: 1. transcript with most clinvar vars 2. transcript with most gnomad vars 3. transcript with most spliceAI vars 4. uniprot isoform 1, unless otherwise specified with -seq")
#parser.add_argument('-ffs','--fromfileseq', dest="ffs", help="A file containing a header line specifying 'uniprot' and 'seq' columns. You can also run a bash loop over your list, calling this script once per line, and supply -uniprot and -seq as arguments.")
parser.add_argument('-ff','--fromfile', dest="ff", help="A file from which to read uniprot IDs (one per line).")
parser.add_argument('-out_folder', dest='out_folder', help="Where the output files should be written. Default location is /storage1/shared/data/prism_merge/uniprot[0:2]/uniprot[2:4]/uniprot[4:6]/")
parser.add_argument('-swiss', dest="swiss", action='store_true', help="Put a folder substructure into the designated output folder (only for -out_folder, it's done automatically for the default output folder).")
parser.add_argument('-m', dest="mode", choices=['overwrite', 'leave'], default = 'leave', help="What do when the output file already exists. Leave (default) or overwrite")
parser.add_argument("-v","--verbose", action="count", default=0, help="Level of output, default zero is no output")
parser.add_argument('-allow', dest="allow", action='store_true', help="Allow merging with files that do not have the exact same sequence (but still the same uniprot ID).")
parser.add_argument('-check_uniprot', dest="check_uniprot", action='store_true', help="If there is no prism_uniprot file, check if this ID is on the fail list (if it's not we can try to make a prism_uniprot file to this ID).")
parser.add_argument('-cleanup', dest="cleanup", action='store_true', help="Reduce the resulting prism file to only rows that have values in the gnomad/spliceAI/clinvar columns since those are actual variants we have data for. Remove rows, i.e. variants for which we only have uniprot feature or netsurfp data.")
parser.add_argument('-fail', dest="fail", action='store_true', help="Consult a file listing uniprot IDs for which we tried to make this type of merge file before and failed. Used to limit the amount of seq lookups we do to uniprot since this is failing a lot now")
parser.add_argument('-min_id', dest="min_id",default=0.8,help="Minimum identity (1.00 = 100 procent) to authoritative sequence needed to qualify file for the merge. Default 0.8.")
parser.add_argument('-min_cov', dest="min_cov",default=0.1,help="Minimum coverage (1.00 = 100 procent) of authoritative sequence needed to qualify file for the merge. Default 0.1.")
parser.add_argument('-snp_af', dest="snp_af",default=0.01,help="Variants with allel frequencies greater than or equal to this will be regarded as SNPs, i.e. common variants. Default 0.01.")
parser.add_argument("--fill_invert", action='store_true', default = False, help="If there are known SNP positions in the gnomad file this will attempt to resolve mismatches with the target seq by inverting those SNPs, also for uniprot and netsurfp files.")
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

def output_file_exists(uID, outsuffix='.all.txt'):
	#define output file name and folder and check if exists
	if args.out_folder:
		if args.out_folder == '.':
			args.out_folder = os.getcwd()
		
		#if no substruc req'd, put all files right in the output folder 
		if not args.swiss:	
			merged_prism_file = os.path.join(args.out_folder,'prism_merged_'+uID+outsuffix)
			#merged_prism_file = os.path.join(args.out_folder,file_name)
		#else, construct subfolder with the uniprot ID
		else:
			merged_prism_file = os.path.join(make_default_outfolder(uID, base=args.out_folder),'prism_merged_'+uID+outsuffix)
			
	else:	
		merged_prism_file = os.path.join(make_default_outfolder(uID),'prism_merged_'+uID+outsuffix)
		#merged_prism_file = os.path.join(make_default_outfolder(uID),file_name)
	
	if args.verbose >= 1:
		print('Output file:', merged_prism_file)
	
	if os.path.exists(merged_prism_file):
		if args.mode == 'leave':
			print(uID, ":" ,merged_prism_file, 'exists. Use -m overwrite to overwrite it.')
			return
		
		elif args.mode == 'overwrite':
		#else:
			return merged_prism_file
	else:
		return merged_prism_file		

def get_uniprot_seq(uID):
	while True:
		try:
			r = requests.get('https://www.uniprot.org/uniprot/'+uID+'.fasta', headers={ "Content-Type" : "text/plain"})
			seq = ''.join(r.text.split('\n')[1:])
			#in rare cases entries have become obsolete and those will have an empty string as seq. I guess we are not making merge files for those
			if not seq:
				print(uID, ': Returned sequence for isoform 1 of', uID, 'is empty. The entry may be obsolete. No merge file will be made since we have no authoritative sequence.')
				return
			else:
				return seq
		#except SSLError:
		except (requests.exceptions.SSLError, requests.exceptions.ConnectionError):
			sleep(10)

def fillVariants(prism_file, v=args.verbose):
	'''
	create a filled up prism_uniprot file by expanding each position to all 19 possible substitutions, i.e.
	S4A DHFR_1 NA NA NA NA NA NA
	S4C DHFR_1 NA NA NA NA NA NA
	S4D DHFR_1 NA NA NA NA NA NA
	...
	'''
	
	shell_command = ['python3', '/storage1/tiemann/dev/repos/prism/scripts/FillVariants.py', prism_file, '-i', '1']
	#shell_command = ['python3', '/storage1/hezscha/src/PRISM/software/domain_protein_features/scripts/FillVariants.py', prism_file]
	if v > 1:
		print('Running FillVariants.py')
		print(shell_command)
	subprocess.run(shell_command,stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True)
	#return the path to the filled up file
	return re.sub('.txt','_filled.txt', prism_file)		 


#new version: check seq identity since we align anyway and only qualify files that have at least x% identity with the target seq (right now 0.8) because at lower identities the parser will refuse. I can change that but it might not make sense to merge those. args.min_ID
def find_prism_file(uID, file_type, seq, v=args.verbose):
	#in contrast to the earlier version, find_prism_file_orig, this function only ever returns at max 1 file name
	
	if file_type not in prism_types:
		raise ValueError("Invalid prism file type. Chose one of: %s" % prism_types)
		
	path = os.path.join(base_prism_dir,'prism',uID[0:2], uID[2:4], uID[4:6])
	if v >= 1:
		print('checking ',path)

	if not os.path.exists(path):
		if v >= 1:
			print(path , "doesn't exist. If looking for clinvar or gnomad file perhaps there are just none.")
			return None

	#some proteins have longer than usual uniprot IDs, i.e. A0A087WTH5. Those can be together with other files in the same dir since their uniprot IDs have a common first 6 letters:
	#prism_uniprot_XXX_A0A087WW49.txt
	#prism_uniprot_XXX_A0A087WYE8.txt
	#prism_uniprot_XXX_A0A087X179.txt
	#therefore, check the uniprotID before adding the file to the list
	
	#divide into file_types that have only the uniprot ID and those that have others following to avoid that short uniprotIDs match in file names with longer uniprot IDs. I.e. A0A087 would match in prism_uniprot_XXX_A0A087WW49.txt if we didn't required the .txt after
	#switch to using v2 of the uniprot files
	if file_type == 'prism_uniprot':
		listing = glob.glob(os.path.join(path, file_type+'_002_'+uID+'.txt'))
		if v > 1:
			print('Looking for file', os.path.join(path, file_type+'_002_'+uID+'.txt'))
			print('Found:', listing)
		
	elif file_type == 'prism_netsurfp':
		listing = glob.glob(os.path.join(path, file_type+'_*_'+uID+'.txt'))
		if v > 1:
			print('Looking for file', os.path.join(path, file_type+'_*_'+uID+'.txt'))
			print('Found:', listing)
	else:
		#other file types have the transcript and geneIDs in their name too, separated by dash from the uniprot ID
		listing = glob.glob(os.path.join(path, file_type+'_*_'+uID+'-*'))
		if v > 1:
			print('Looking for file', os.path.join(path, file_type+'_*_'+uID+'-*'))
			print('Found:', listing)
	
	#if len(listing) > 1:
	if len(listing):
		if len(listing) > 1 and v >= 1:
			print('There is more than one' ,file_type, 'file for this uniprot ID:', uID)#, file = sys.stderr)
		if v > 1:
			print(listing)
		for prism_file in listing:
			meta, df = read_from_prism(prism_file)
			#debug
			#print(meta['protein']['sequence'])
			#debug
			if meta['protein']['sequence'] == seq:
				if v >= 1:
					print(prism_file, 'matches the authoritative sequence.')
				return prism_file
		
		else:
			if args.allow:
				max_score = float('-inf')
				best_match = None
				target_seq = seq
				#n_res_target = len(target_seq)
				#todo! If args.allow, go through all files again and find the closest seq to the requested
				for prism_file in listing:
					meta, df = read_from_prism(prism_file)
					probe_seq = meta['protein']['sequence']
					#n_res_data = len(probe_seq)
				
					#new version: switch to do_align
					align_obj = do_align(target_seq.upper(), probe_seq.upper(), verbose=v)
					if v > 1:
						print('Aligning', prism_file, 'to target seq.')
						print(align_obj.n_align,'alignments found.')
						print('Scores are:', align_obj.all_scores)

					if align_obj.best_score > max_score:
						max_score = align_obj.best_score
						best_match = prism_file
						n_res = align_obj.aln_len
						coverage = align_obj.cov
						identity = align_obj.id_a
						#same coverage and identity. calc from PrismData
						#n_res = max(n_res_target, n_res_data)
						#coverage = (n_res-align_obj.n_gaps)/n_res
						#identity = (n_res_target-len(align_obj.subs))/n_res_target
					
				#check cov and id are sufficient
				if coverage < args.min_cov or identity < args.min_id:
					if v >= 1:
						print("The best matching seq does not fullfill the requirements for minimum alignment coverage and identity. It has coverage", round(coverage,2), "and identity", round(identity,2), ". Required are at least coverage", args.min_cov, 'and identity', args.min_id)
					return None
				else:	
					#in the end, return the file that had the best alignment score, we'll try to merge that one
					if v >= 1:
						print("None of the files match the authoritative sequence but allow has been used so one of the files will be merged anyway.")
						print('Best matching sequence is', best_match)
					return best_match
			
			else:
				print("None of the files match the authoritative sequence. Use -allow to enable merging of prism files with not the same sequence (but the same uniprot ID). This feature doesn't work yet.")
				return None		
	else:
		if v >= 1:
			print('No', file_type,'file found for uniprot ID:', uID)#, file = sys.stderr)
		return None

def make_default_outfolder(uniprot_ID, base="prism_merge/"):
	#sometimes we cannot find the unirpot ID to a transcript. Put them in the none folder then
	if not uniprot_ID == 'None':
		if not os.path.exists(os.path.join(base_prism_dir, base, uniprot_ID[0:2])):
			os.makedirs(os.path.join(base_prism_dir, base, uniprot_ID[0:2]))
		if not os.path.exists(os.path.join(base_prism_dir, base, uniprot_ID[0:2], uniprot_ID[2:4])):
			os.makedirs(os.path.join(base_prism_dir, base, uniprot_ID[0:2], uniprot_ID[2:4]))
		if not os.path.exists(os.path.join(base_prism_dir, base, uniprot_ID[0:2], uniprot_ID[2:4], uniprot_ID[4:6])):
			os.makedirs(os.path.join(base_prism_dir, base, uniprot_ID[0:2], uniprot_ID[2:4], uniprot_ID[4:6]))
	
		return(os.path.join(base_prism_dir, base, uniprot_ID[0:2], uniprot_ID[2:4], uniprot_ID[4:6]))
	else:
		if not os.path.exists(os.path.join(base_prism_dir, base, 'None/')):
			os.makedirs(os.path.join(base_prism_dir, base, 'None/'))
		return(os.path.join(base_prism_dir, base, 'None/'))

def merge_prism(out, files, uID, target_seq, v=args.verbose):
	print('Merging files: ', ','.join(files))
	
	# ~ if args.fill_invert:
		# ~ more_args = ['--fill_invert']
	# ~ else:
		# ~ more_args =[]
		
	#https://stackoverflow.com/questions/4965159/how-to-redirect-output-with-subprocess-in-python (Ryan Thompson answer)
	if args.allow:
		with open('/'.join(out.split('/')[:-1])+'/'+uID+'.log', 'w') as outfile:
			#shell_command = ['python3', '/storage1/hezscha/src/PRISM/prism/scripts/PrismData_HZ_allow_SNP_v4.py',
			shell_command = ['python3', '/storage1/hezscha/src/PRISM/prism/scripts/PrismData.py',
					 '--merge', out] + files + ['-vv', '--target_seq', target_seq]
			if v > 1:
				print(shell_command)
			subprocess.run(shell_command,stdout=outfile, stderr=outfile, check=True, text=True)
		
	else:
		#shell_command = ['python3', '/storage1/hezscha/src/PRISM/prism/scripts/PrismData_HZ_allow_SNP_v4.py',
		shell_command = ['python3', '/storage1/hezscha/src/PRISM/prism/scripts/PrismData.py',
					 '--merge', out] + files + ['--target_seq', target_seq]
		if v > 1:
			print(shell_command)
		subprocess.run(shell_command,stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True)
				   
	meta_data, dataframe = read_from_prism(out)
	return meta_data, dataframe

def clean_filled_prism_file(prism_file):
	#filled_file = re.sub('.txt','_filled.txt', prism_file)
	if os.path.isfile(prism_file):
		os.remove(prism_file)

################################################################################
#global vars:
base_prism_dir = '/storage1/shared/data/'
prism_types = ['prism_uniprot', 'prism_clinvar', 'prism_gnomad', 'prism_netsurfp', 'prism_spliceai']
rand_str = ''

#parse arguments
#if (args.uniprot and args.ff) or (args.uniprot and args.ffs) or (args.ff and args.ffs):
if (args.uniprot and args.ff) :
	#print('The options -uniprot -ffs and -ff are mutually exlusive. Use only one of them.')
	print('The options -uniprot and -ff are mutually exlusive. Use only one of them.')
	sys.exit()

if args.verbose >= 1:
	print(str(datetime.now()))

if args.typeslog:
	if os.path.exists('prism_type_log'):
		rand_str = ''.join([choice(ascii_letters) for i in range(4)])
		shutil.move('prism_type_log','prism_type_log'+rand_str)
	
	typeOUT = open('prism_type_log', 'w')

if args.swiss and not args.out_folder:
	print('-swiss is meant for use with a user defined output folder, -out_folder, and ignored otherwise')

#summary file
if os.path.exists('merge_summary.txt'):
	if not rand_str:
		rand_str = ''.join([choice(ascii_letters) for i in range(4)])
	shutil.move('merge_summary.txt','merge_summary_'+rand_str+'.txt')
	
sumOUT = open('merge_summary.txt', 'w')

#fail list is a list of uniprot IDs that were tried before for this combination of file types and didn't yield results
fail_list = set()
if args.fail:
	fail_path = "/storage1/hezscha/genome_proteome_map/ect/merge_fail" + outsuffix
	#if file doesn't exist yet because it's the first time, do not attempt to open
	if os.path.exists(fail_path):	
		with open(fail_path, 'r') as IN:
			for line in IN:
				fail_list.add(line.rstrip())

if args.uniprot:
	uID = args.uniprot
	if args.fail:
		if uID in fail_list:
			print('Previously failed to make merge files for', uID, ' so no attempt will be made. If you think this is erronous try running without -fail.')
			if args.verbose >= 1:
				print(str(datetime.now()))
			sys.exit()	
	uniprot_ls = [uID]
	
elif args.ff:
	with open(args.ff, 'r') as IN:
		uniprot_ls = []
		for line in IN:
			if line.startswith('#'):
				continue
			uID = line.rstrip()
			#is this on the fail list for this cominbation of requested file types
			if args.fail:
				if uID in fail_list:
					print('Previously failed to make merge files for', uID, ' so no attempt will be made. If you think this is erronous try running without -fail.')
					if args.verbose >= 1:
						print(str(datetime.now()))
					continue
			
			uniprot_ls.append(uID)
	
for uID in uniprot_ls:	
	
	#step0: check if output file exists. The check for mode, overwrite or leave is done inside the funct, so if an output file name is returned we are supposed to make the file
	merged_prism_file = output_file_exists(uID)
	if merged_prism_file is None:
		continue
	
	#step1: find the clinvar file with the most non-VUS vars and set that sequence as the target_seq. 
	#If there are no clinvar files to this uniprot, I guess the gnomad file with the most vars is the next candidate for the target_seq since gnomad vars are also still data points and then the splice AI file if there is no gnomad file to this uniprot ID either. 
	#uniprot pos features are not actual vars. They apply to the WT protein at that position. They are mostly useful if we have additional info on actual vars found at that pos form either clinvar, gnomad or spliceAI. SpliceAI has the least weight of these three since those are predictions of variants that are potentially splice altering, not observed splice variants
	
	path = os.path.join(base_prism_dir,'prism',uID[0:2], uID[2:4], uID[4:6])
	file_type = 'prism_clinvar'
	listing = glob.glob(os.path.join(path, file_type+'_*_'+uID+'-*'))
	if args.verbose > 1:
		print('Looking for clinvar file at',  os.path.join(path, file_type+'_*_'+uID+'-*'))
		print('Found:', listing)
	
	if listing:
		#choose best clinvar file
		max_var = 0 #changed to -1 so that if there is 1 clinvar file and it has only VUS we still qualify it
		best_file = ''
		for prism_file in listing:
			metadata, df = read_from_prism(prism_file)
			n_var = sum(df['Clinvar_signifiance'] != 'VUS')
			if n_var > max_var:
				max_var = n_var
				best_file = prism_file
				target_seq = metadata['protein']['sequence']
				
		if max_var == 0:
			if args.verbose >= 1:
				print('Clinvar file only has VUS, looking for gnomad file.')
			#I don't see a good way to proceed here except for copying the entire code fragment from below :(
			#else choose best gnomad file, which is just the one with the most entries
			file_type = 'prism_gnomad'
			#path = os.path.join(base_prism_dir,file_type,uID[0:2], uID[2:4], uID[4:6])
			listing = glob.glob(os.path.join(path, file_type+'_*_'+uID+'-*'))	
			
			if listing:
				max_var = 0
				best_file = ''
				for prism_file in listing:
					metadata, df = read_from_prism(prism_file)
					n_var = len(df)
					if n_var > max_var:
						max_var = n_var
						best_file = prism_file
						target_seq = metadata['protein']['sequence']
						
			else:
				if args.verbose > 1:
					print('No gnomad files, looking for spliceai file.')
				#else, choose the best spliceAI file			
				file_type = 'prism_spliceai'
				#path = os.path.join(base_prism_dir,file_type,uID[0:2], uID[2:4], uID[4:6])
				listing = glob.glob(os.path.join(path, file_type+'_*_'+uID+'-*'))	
				
				if listing:
					max_var = 0
					best_file = ''
					for prism_file in listing:
						metadata, df = read_from_prism(prism_file)
						n_var = len(df)
						if n_var > max_var:
							max_var = n_var
							best_file = prism_file
							target_seq = metadata['protein']['sequence']
							
				else:
					if args.verbose > 1:
						print('No splice files, getting target seq from uniprot.')
					#else, get seq of uniprot isoform 1	
					best_file = 'uniprot'
					target_seq = get_uniprot_seq(uID)
					if target_seq is None:
						print(uID,": No clinvar, gnomad and spliceai files and no uniprot isoform 1 available, therefore no target sequence available.")
						continue
					if 'U' in target_seq:
						target_seq = re.sub('U', 'X', target_seq)		

	else:
		if args.verbose > 1:
			print('No clinvar files, looking for gnomad file.')
		#else choose best gnomad file, which is just the one with the most entries
		file_type = 'prism_gnomad'
		#path = os.path.join(base_prism_dir,file_type,uID[0:2], uID[2:4], uID[4:6])
		listing = glob.glob(os.path.join(path, file_type+'_*_'+uID+'-*'))	
		
		if listing:
			max_var = 0
			best_file = ''
			for prism_file in listing:
				metadata, df = read_from_prism(prism_file)
				n_var = len(df)
				if n_var > max_var:
					max_var = n_var
					best_file = prism_file
					target_seq = metadata['protein']['sequence']
					
		else:
			if args.verbose > 1:
				print('No gnomad files, looking for spliceai file.')
			#else, choose the best spliceAI file			
			file_type = 'prism_spliceai'
			#path = os.path.join(base_prism_dir,file_type,uID[0:2], uID[2:4], uID[4:6])
			listing = glob.glob(os.path.join(path, file_type+'_*_'+uID+'-*'))	
			
			if listing:
				max_var = 0
				best_file = ''
				for prism_file in listing:
					metadata, df = read_from_prism(prism_file)
					n_var = len(df)
					if n_var > max_var:
						max_var = n_var
						best_file = prism_file
						target_seq = metadata['protein']['sequence']
						
			else:
				if args.verbose > 1:
					print('No splice files, getting target seq from uniprot.')
				#else, get seq of uniprot isoform 1	
				best_file = 'uniprot'
				target_seq = get_uniprot_seq(uID)
				#in rare cases like A0A087X1A0 (which did end up having a spliceai file I didn't prev qualify as target seq due to another error that is now fixed) the uniprot entry is deprecated and get_target_seq returns None. Check for that. If so, we can't do anything since we need a target seq and looking for transcripts with clinvar, gnomad and spliceai data has failed
				if target_seq is None:
					print(uID,": No clinvar, gnomad and spliceai files and no uniprot isoform 1 available, therefore no target sequence available.")
					continue
				if 'U' in target_seq:
					target_seq = re.sub('U', 'X', target_seq)	
				

	if args.verbose >= 1:
		print('Getting authoritative transcript sequence from:', best_file)

	#now qualify at most 5 component files (one per type_file)
	
	uniprot_file = find_prism_file(uID=uID, seq=target_seq, file_type='prism_uniprot')
	if uniprot_file is None:
		if args.check_uniprot: 
			with open(os.path.join(base_prism_dir, 'prism_uniprot/fail_list.txt'), 'r') as IN:
				for line in IN:
					if uID == line.rstrip():
						print(uID, 'is either obsolete or had no features of interest and therefore no prism_uniprot file exists.')
						break
	#else:
	#	uniprot_file = fillVariants(uniprot_file)					

	gnomad_file = find_prism_file(uID=uID, seq=target_seq, file_type='prism_gnomad')
	clinvar_file = find_prism_file(uID=uID, seq=target_seq, file_type='prism_clinvar')
	netsurfp_file = find_prism_file(uID=uID, seq=target_seq, file_type='prism_netsurfp')
	spliceAI_file = find_prism_file(uID=uID, seq=target_seq, file_type='prism_spliceai')
	
	#if we want to try and match non-matching posititions to the target seq do this
	if args.fill_invert:
		if gnomad_file is not None:
			
			#load in gnomad file
			md, df = read_from_prism(gnomad_file)
			gnomad_seq = md['protein']['sequence']
			
			#find mismatches to the target_seq
			if gnomad_seq == target_seq:
				subs = []
				inner_gaps = 0
			else: 
				align_obj = do_align(gnomad_seq, target_seq, verbose = args.verbose)
				#there is no need to check the alignment quality since it was already checked when we picked the gnomad file
				subs = align_obj.subs
				inner_gaps = align_obj.inner_gaps
			
			#v005: re-add checking alignment since I still think we can only do the mismacht handling (swapping to target_seq WT) is the alignment is ungapped or only has leading or trailing gaps, due to uncertainty of which residues are equivalent in target_seq and other_seq(gnomad, uniprot) around the gaps otherwise.
			
			
			#unlist the residue number, ref and alt aa. We can do this since there are no multi mutants in gnomad, or rather we can't know if two mutations are in the same person, so by definition all variants are assumed to be single mutants
			df['resi'] = df['resi'].apply(lambda x: x[0])
			df['aa_ref'] = df['aa_ref'].apply(lambda x: x[0])
			df['aa_var'] = df['aa_var'].apply(lambda x: x[0])
			
			file_changes = 0
			#create column that tells whether we have matched the target WT (so the gnomad freq is an estimate, not an actual observation)
			if subs and inner_gaps == 0:
				df['match_target'] = False
				md['columns']['match_target'] = 'WT aa has been changed to match the aa in the target_seq at this position'
				
				#record the mismatches because I want to check them later
				mismatchOUT = open('merge_mismatches_'+uID+'.txt', 'w')
				#if args.verbose >= 1:
				print(uID, ': Mismatches between', gnomad_file, 'and target_seq:', subs, file = mismatchOUT)
			
				#these are positions we want to try and fix if the target_aa has also been observed at that pos in gnomad, as a variant aa
				for (pos,gnomad_aa,target_aa) in subs:
					#consider possible offset due to leading gap
					pos = pos - align_obj.a_offset
					this_pos = df.loc[df['resi'] == pos]
					#has the target_aa (aa in the target_seq) been observed in gnomad at this position?
					if any(this_pos.aa_var == target_aa):
						file_changes = 1
						#fix the seq in the metadata
						md['protein']['sequence'] = md['protein']['sequence'][:pos-1] + target_aa + md['protein']['sequence'][pos:]
						if args.verbose >= 1:
							print('Found', target_aa, 'in gnomad at', pos)
							print(uID, ': vars at mismatch pos', pos, ":", this_pos.variant.squeeze(), file = mismatchOUT)	
						#all variants in this_pos MUST start with gnomad_aa + pos since this is a prism file and only one WT is allowed. this means we want to convert all variants in this_pos to match the target_aa
						
						#In terms of AF, we need to recalc only for the one aa that is the gnomad_aa, since that one is the purpoted WT and does therefore not have an AF since we only have AFs for purpoted variants.
						#They are all gnomad_aa+pos as we established earlier. So better to say we need to recalc for the var gnomad_aa+pos+target_aa since we will turn that one into target_aa+pos+gnomad_aa
						#all other vars at the same pos are gnomad_aa+pos+other_aa(not target) and will be changed into target_aa+pos+other_aa
						#so we need to recalc for A179V (since we don't know the AF of V but we do know they should sum to 1) 
						#and NOT for any other A179*, including A179E
						sum_af = 0
						sum_af_g = 0
						sum_af_ex = 0
						for index, row in this_pos.iterrows():
							#print('looking at', row.variant)
							#variants refering to amino acid that are NOT the target_aa
							#while we're not sure what to do with these, swap the WT to target but put NA everywhere else
							if row.variant != gnomad_aa+str(pos)+target_aa:
								#print(row.variant)
								sum_af += row.AF_tot
								sum_af_g += row.AF_g if not math.isnan(row.AF_g) else 0
								sum_af_ex += row.AF_ex if not math.isnan(row.AF_ex) else 0
								
								new_var = target_aa+str(pos)+row.aa_var
								#safe the AFs
								#AF(target_aa+str(pos)+row.aa_var) ~ AF(gnomad_aa+pos+row.aa_var). For reasoning check Trello board https://trello.com/b/mhbusoya/proteome-wide-stuff, list merging, card SNP inversion
								AF_tot = df.loc[df[df['variant'] == row.variant].index[0],'AF_tot']
								AF_ex = df.loc[df[df['variant'] == row.variant].index[0],'AF_ex']
								AF_g = df.loc[df[df['variant'] == row.variant].index[0],'AF_g']
								if 'SNP' in df.columns:
									SNP = df.loc[df[df['variant'] == row.variant].index[0],'SNP']
								#print(AF_tot,AF_ex,AF_g)
								
								#replace variant gnomad_aa+pos+other_aa with target_aa+pos+other_aa
								#I also want to replace the other fields with NA to signify that this was not the orig observed substitution 
								rep_obj = (('NA',)*(df.shape[1]-5))
								df.loc[df[df['variant'] == row.variant].index[0]] = (new_var,) + rep_obj + ( 1, target_aa, pos, row.aa_var)
								#re-add the allel frequencies: AF(target_aa+str(pos)+row.aa_var) ~ AF(gnomad_aa+pos+row.aa_var . See above
								df.loc[df[df['variant'] == new_var].index[0], 'AF_tot'] = AF_tot
								df.loc[df[df['variant'] == new_var].index[0], 'AF_ex'] = AF_ex
								df.loc[df[df['variant'] == new_var].index[0], 'AF_g'] = AF_g
								df.loc[df[df['variant'] == new_var].index[0], 'match_target'] = True
								if 'SNP' in df.columns:
									df.loc[df[df['variant'] == new_var].index[0], 'SNP'] = SNP
								
								
						#now that we've seen all vars that are gnomad_aa+pos+* (* not target_aa), calc the freq for the gnomad_aa from 1 - sum(other_aas_at_that_pos)   
						#need to make the row into a series: https://datatofish.com/pandas-dataframe-to-series/
						row = this_pos.loc[this_pos['variant'] == gnomad_aa+str(pos)+target_aa].squeeze()
						#print('Now looking at', row.variant)
						sum_af += row.AF_tot
						sum_af_g += row.AF_g
						sum_af_ex += row.AF_ex
						
						inv_af = 1 - sum_af
						inv_af_g = 1 - sum_af_g
						inv_af_ex = 1 - sum_af_ex
						#print(sum_af, sum_af_g, sum_af_ex)
						#print(inv_af, inv_af_g, inv_af_ex)
						
						new_var = target_aa+str(pos)+gnomad_aa
						#print(row.variant, inv_af)
						#add new row with inverted AF and remove old row
						if 'SNP' in df.columns:
							SNP = df.loc[df[df['variant'] == row.variant].index[0], 'SNP']
					
						rep_obj = (('NA',)*(df.shape[1]-5))
						df.loc[df[df['variant'] == row.variant].index[0]] = (new_var,) + rep_obj + ( 1, target_aa, pos, row.aa_var)
						df.loc[df[df['variant'] == new_var].index[0], 'AF_tot'] = inv_af
						df.loc[df[df['variant'] == new_var].index[0], 'AF_ex'] = inv_af_ex
						df.loc[df[df['variant'] == new_var].index[0], 'AF_g'] = inv_af_g
						df.loc[df[df['variant'] == new_var].index[0], 'match_target'] = True
						if 'SNP' in df.columns:
							df.loc[df[df['variant'] == new_var].index[0], 'SNP'] = SNP
						
				#close here the log file
				mismatchOUT.close()
						
			#write out the changed gnomad file as SNPinv file
			if file_changes:
				new_gnomad = re.sub('.txt','_SNPinv.txt',gnomad_file)
				write_prism(md,df,new_gnomad)
				gnomad_file = new_gnomad
			
		#load in uniprot file if alignment quality was sufficient that it got selected
		#in case of mismatch, match the target_aa and put True in the match_target column
		for filename in [uniprot_file, netsurfp_file]:
			if filename is not None:
				md, df = read_from_prism(filename)
				uniprot_seq = md['protein']['sequence']
				
				#first, establish if SNP inv necessary, i.e. if there are mismatches to target seq
				if uniprot_seq == target_seq:
					continue #if they are identical we have nothing to do
				
				else:
					align_obj = do_align(uniprot_seq, target_seq, verbose = args.verbose)
					subs = align_obj.subs
					inner_gaps = align_obj.inner_gaps

				file_changes = 0
				if subs and inner_gaps == 0:
					df['match_target'] = False
					file_changes = 1
					md['columns']['match_target'] = 'WT aa has been changed to match the aa in the target_seq at this position'
					mismatchOUT = open('merge_mismatches_'+uID+'.txt', 'a')
					print(uID, ': Mismatches between', filename, 'and target_seq:', subs, file = mismatchOUT)
					
					for (pos,seq_aa,target_aa) in subs:		
						pos = pos - align_obj.a_offset			
						#swap seq_aa+pos+= with target_aa+pos+=
						new_var = target_aa+str(pos)+'='
						old_var = seq_aa+str(pos)+'='
						if args.verbose > 0:
							print('swap', old_var, 'with', new_var)
						df.loc[df['variant'] == old_var, 'variant'] = new_var
						df.loc[df['variant'] == new_var, 'aa_ref'] = [target_aa]
						df.loc[df['variant'] == new_var, 'match_target'] = True
						
						#fix the seq in the metadata
						md['protein']['sequence'] = md['protein']['sequence'][:pos-1] + target_aa + md['protein']['sequence'][pos:]
						#debug
						if args.verbose > 0:
							print('metadata:', md['protein']['sequence'][pos-1], file = mismatchOUT)
							print('new var line:', file = mismatchOUT)
							print(df.loc[df['variant'] == new_var], file = mismatchOUT)
						#debug
					mismatchOUT.close()
					
				if file_changes:
					new_filename = re.sub('.txt','_SNPinv.txt',filename)
					write_prism(md,df,new_filename)
					if 'uniprot' in filename:
						uniprot_file = new_filename
					else:
						netsurfp_file = new_filename
	
	#fill up uniprot and netsurfp files. We do this after the SNP inverting business
	if uniprot_file is not None:
		uniprot_file = fillVariants(uniprot_file)
	if netsurfp_file is not None:
		netsurfp_file = fillVariants(netsurfp_file)
	
	#merge the files that exist
	files = []
	for item in [uniprot_file, clinvar_file, gnomad_file, netsurfp_file, spliceAI_file]:
		if item is not None:
			files.append(item) 
	
	if len(files) < 2:
		print(uID, ': Not enough files to merge, need at least 2.', file = sumOUT)
		print(uID, ': Not enough files to merge, need at least 2.')
		continue
	#if only the netsurfp and uniprot files exist but we have asked for more types so they will be filled up and we have used -cleanup, the resulting df will always be empty since we remove all lines having only netsurfp and uniprot data. So just stop at this point is this is the case	
	elif len(files) == 2 and uniprot_file in files and netsurfp_file in files and uniprot_file.endswith('_filled.txt') and netsurfp_file.endswith('_filled.txt') and args.cleanup:
		print(uID ,': Only prism_uniprot and prism_netsurfp files are available and cleanup was requested which will result in an empty dataframe so no merging is done.', file = sumOUT)
		if args.fail:
			with open(fail_path, 'a') as fOUT:
				print(uID, file = fOUT)
		continue
		
	else:
		try:	
			(all_merged_metadata, all_merged_df) = merge_prism(out=merged_prism_file,files=files,uID=uID,target_seq=target_seq)
		#I've had one file fail to merge because of more than one var existing at a SNP site. 
		#Would like to see how many of those there are
		except:
			print(uID, ': Merge failed', file = sumOUT)
			continue	
		
	#added check that this is the filled up file since we now allow merging of unfilled files if only pos spec types are requested (uniprot and netsurfp for now) 
	if uniprot_file is not None and uniprot_file.endswith('_filled.txt'):
		clean_filled_prism_file(uniprot_file)
	if netsurfp_file is not None and netsurfp_file.endswith('_filled.txt'):
		clean_filled_prism_file(netsurfp_file)
	if args.fill_invert:
		list_snpinv = glob.glob(os.path.join('/storage1/shared/data/*',uID[0:2], uID[2:4], uID[4:6], '*_SNPinv.txt'))
		for fn in list_snpinv:
			clean_filled_prism_file(fn)
		
		
	#cleanup the resulting prism file	
	#add the part Johanna has written about removing lines from the merged file that only have uniprot/netsurfp features, i.e. we have no actual (gnomad/spliceAI,clinvar) data on them
	if args.cleanup:
	
		#adding the prism_type to column names, replacing old column names in the df and the metadata
		key_dic = {}
		new_columns = {}
		
		for key in all_merged_metadata['merged']:
			num = key.split("_")[1]
			name = all_merged_metadata['merged'][key].split('_')[1]
			for column in all_merged_metadata['columns']:
				if column.endswith(num):
					key_dic[column] = f'{name};{column}'
					new_columns[f'{name};{column}'] = all_merged_metadata['columns'][column]
		
		all_merged_df.rename(columns=key_dic, inplace=True)
		all_merged_metadata['columns'] = new_columns
		
		#list columns starting with uniprot or netsurfp
		non_relevant = [elem for elem in list(all_merged_df.columns) if elem.startswith('uniprot')] + [elem for elem in list(all_merged_df.columns) if elem.startswith('netsurfp')]

		#make a dataframe only containing the data cols, so without uniprot, netsurfp, general cols
		tmp_df = all_merged_df.copy()
		tmp_df = tmp_df[tmp_df.columns.difference(non_relevant+['variant','n_mut', 'aa_ref', 'resi', 'aa_var'])]
		
		#get the names of the rows that do not only have NAs (i.e. at least one non-NA col). dropna returns the retained df
		ret_idx = tmp_df.dropna(axis=0, how='all').index
		
		#subset orig df to only the rows to be retained
		good = all_merged_df.iloc[ret_idx,].copy().reset_index(drop=True)
		
		#drop columns that are only composed of NA's now. This can happen with uniprot feature cols that only had features for positions that were dropped because we have no other data for them
		drop_meta_columns = [key for key in all_merged_metadata['columns'].keys() if not key in good.columns]
		for key in drop_meta_columns:
			all_merged.metadata['columns'].pop(key)

		#only write this file if there are data rows left. If for example only a uniprot file and a netsurfp file got merged (this should not happen but it has ...), instead delete the merge file
		if good.empty:
			print(uID, ': the resulting dataframe after cleaning is empty. Removing the merge file:', merged_prism_file, file = sumOUT)
			os.remove(merged_prism_file)
			if args.fail:
				with open(fail_path, 'a') as fOUT:
					print(uID, file = fOUT)
		else:
			comment = ['cleaned of rows with only uniprot or netsurfp data',]
			write_prism(all_merged_metadata, good, merged_prism_file, comment = comment)
			#debug: write to diff filename for comparison
			#write_prism(all_merged_metadata, good, merged_prism_file+'.cleaned', comment = comment)
			#if args.verbose >= 1:
			#	print(str(datetime.now()))
			#at this point we assume writing and clean up has succeeded
			if args.typeslog:
				files_log = []
				if uniprot_file is not None:
					files_log.append('uniprot')
				if clinvar_file is not None:
					files_log.append('clinvar')	
				if gnomad_file is not None:
					files_log.append('gnomad')
				if spliceAI_file is not None:
					files_log.append('spliceai')
				print(uID, ','.join(files_log), len(files_log), sep = ';' , file = typeOUT)
			print(uID, ': success', file = sumOUT)	
		
	else:
		#if no cleanup was asked for we can still report which component files were merged
		if args.typeslog:
			files_log = []
			if uniprot_file is not None:
				files_log.append('uniprot')
			if clinvar_file is not None:
				files_log.append('clinvar')	
			if gnomad_file is not None:
				files_log.append('gnomad')
			if spliceAI_file is not None:
				files_log.append('spliceai')
			print(uID, ','.join(files_log), len(files_log), sep = ';' , file = typeOUT)
		print(uID, ': success', file = sumOUT)	
	#	if args.verbose >= 1:
	#		print(str(datetime.now()))


	#now, take care of the log file. If there were no problems, delete it
	if args.allow:
		logfile = '/'.join(merged_prism_file.split('/')[:-1])+'/'+uID+'.log'
		rm_log = 1
		with open(logfile,'r') as logIN:
			for line in logIN:
				if line.startswith('WARNING:'):
					rm_log = 0
					break
				elif re.match('New reference sequence will remove \d+ variants \(mismatch setting\: remove\;', line):
					rm_log = 0
					break
		
		if rm_log:
			os.remove(logfile)	

if args.typeslog:
	typeOUT.close()
sumOUT.close()

##################





