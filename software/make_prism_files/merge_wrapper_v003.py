
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
from Bio.SubsMat import MatrixInfo
from collections import namedtuple

import_path_base = '/storage1/hezscha/src/'
sys.path.insert(1, import_path_base + 'PRISM/prism/scripts/')
#from PrismData_HZ_allow_SNP_v4 import PrismParser, VariantData
from PrismData import PrismParser, VariantData

#parse arguments
################################################################################
#https://stackoverflow.com/questions/50021282/python-argparse-how-can-i-add-text-to-the-default-help-message
parser = argparse.ArgumentParser(epilog=textwrap.dedent('''\
         additional information:
         The options -uniprot, -ffs and -ff are mutually exlusive. Use only one of them. -seq only pairs with the special case of a single uniprot ID passed to -uniprot. Otherwise use -ffs if you have a list of uniprot IDs and associated sequences.
         '''))
parser.add_argument('-uniprot', dest="uniprot", help="The uniprot ID for which to merge files. The sequence of isoform 1 will be the authorative sequence for merging unless otherwise specified with -seq")
parser.add_argument('-ffs','--fromfileseq', dest="ffs", help="A file containing a header line specifying 'uniprot' and 'seq' columns. You can also run a bash loop over your list, calling this script once per line, and supply -uniprot and -seq as arguments.")
parser.add_argument('-ff','--fromfile', dest="ff", help="A file from which to read uniprot IDs (one per line).")
# ~ parser.add_argument('-seq', dest="seq", help="The sequence for which to merge prism files. If this arg is not used it will be isoform 1 of the given uniprotID. This only works when running on a single uniprot ID passed with -uniprot. If you have several uniprot IDs and the sequences of interest use -fromfile.")
# ~ parser.add_argument('-delim', dest="delim", default = ' ', help="The delimiter of the input file to fromfile. Assumed to be space if not given.")
parser.add_argument('-out_folder', dest='out_folder', help="Where the output files should be written. Default location is /storage1/shared/data/prism_merge/uniprot[0:2]/unipro[2:4]/uniprot[4:6]/")
parser.add_argument('-m', dest="mode", choices=['overwrite', 'leave'], default = 'leave', help="What do when the output file already exists. Leave (default) or overwrite")
parser.add_argument("-v","--verbose", action="count", default=0, help="Level of output, default zero is no output")
parser.add_argument('-allow', dest="allow", action='store_true', help="Allow merging with files that do not have the exact same sequence (but still the same uniprot ID). This feature is under development while I figure out how to select the closest sequence (by using the aligner in the parser so we don't have another layer of complexity).")
parser.add_argument('-check_uniprot', dest="check_uniprot", action='store_true', help="If there is no prism_uniprot file, check if this ID is on the fail list (if it's not we can try to make a prism_uniprot file to this ID).")
parser.add_argument('-cleanup', dest="cleanup", action='store_true', help="Reduce the resulting prism file to only rows that have values in the gnomad/spliceAI/clinvar columns since those are actual variants we have data for. Remove rows, i.e. variants for which we only have uniprot feature or netsurfp data.")
#parser.add_argument('-only', dest="only", help = "Only merge the specified file types. Comma separated list with NO white space, i.e. -only netsurfp,spliceai . You can request using full names(spliceAI, netsurfp, clinvar, uniprot, gnomad) or short names(sai, nsp, clin, uni, gnomad) or a mixture of both.")
#parser.add_argument('-req', dest="req", help = "Required file types for merging. If any of those is not available no merge file will be made. Comma separated list with NO white space, i.e. -req netsurfp,spliceai . You can request using full names(spliceAI, netsurfp, clinvar, uniprot, gnomad) or short names(sai, nsp, clin, uni, gnomad) or a mixture of both. As a special case -req all means all 5 types are required")
parser.add_argument('-fail', dest="fail", action='store_true', help="Consult a file listing uniprot IDs for which we tried to make this type of merge file before and failed. Used to limit the amount of seq lookups we do to uniprot since this is failing a lot now")
parser.add_argument('-min_id', dest="min_id",default=0.8,help="Minimum identity (1.00 = 100%) to authoritative sequence needed to qualify file for the merge. Default 0.8.")
parser.add_argument('-min_cov', dest="min_cov",default=0.1,help="Minimum coverage (1.00 = 100%) of authoritative sequence needed to qualify file for the merge. Default 0.1.")
parser.add_argument('-snp_af', dest="snp_af",default=0.01,help="Variants with allel frequencies greater than or equal to this will be regarded as SNPs, i.e. common variants. Default 0.01.")
parser.add_argument("--fill_invert", action='store_true', default = False, help="Use with --merge. If there are known SNP positions this will create the inverse variants.")
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

def output_file_exists(uID, outsuffix='.all.txt'):
	#define output file name and folder and check if exists
	if args.out_folder:
		if args.out_folder == '.':
			args.out_folder = os.getcwd()
		merged_prism_file = os.path.join(args.out_folder,'prism_merged_'+uID+outsuffix)
		#merged_prism_file = os.path.join(args.out_folder,file_name)
	else:	
		merged_prism_file = os.path.join(make_default_outfolder(uID),'prism_merged_'+uID+outsuffix)
		#merged_prism_file = os.path.join(make_default_outfolder(uID),file_name)
	
	if args.verbose >= 1:
		print('Output file:', merged_prism_file)
	
	if os.path.exists(merged_prism_file):
		if args.mode == 'leave':
			print(merged_prism_file, 'exists. Use -m overwrite to overwrite it.')
			return
		
		elif args.mode == 'overwrite':
		#else:
			return merged_prism_file
	else:
		return merged_prism_file		

def do_align(seq1,seq2, same_length = False, verbose = 0):#, return_scores = False, return_nr_align = False):
	
	#named tuple approach: 
	align_obj = namedtuple('align_obj',['n_gaps','subs','best_score','n_align','all_scores','aln_len', 'cov', 'id_s1', 'id_s2'])
	
	if same_length:
		if len(seq1) != len(seq2):
			if verbose:
				print('Seq 1 and 2 are not of equal length!')
			ret_tuple = align_obj(None, [], None, None,[], None, None, None, None)
			return(ret_tuple)
				
	n_gaps = 0
	subs = []
	align = pairwise2.align.globalds(seq1, seq2, MatrixInfo.blosum62, -3, -1)
	align_scores = np.array([a[2] for a in align])
	n_align = len(align)
	align = align[0] #only proceed with first alignment
	pos = 1
	for a, b in zip(align[0], align[1]):
		if a == "-" or b == "-": 
			n_gaps += 1
		elif a != b:
			subs.append((pos,a,b))
		pos += 1

	#add coverage and perc identity calc
	n_res = max(len(seq1), len(seq2))
	coverage = (n_res-n_gaps)/n_res
	identity_s1 = (len(seq1)-len(subs))/len(seq1)
	identity_s2 = (len(seq2)-len(subs))/len(seq2)
	
	ret_tuple = align_obj(n_gaps, subs, max(align_scores),n_align,align_scores, n_res, coverage, identity_s1, identity_s2)
	return(ret_tuple)

	#return(n_gaps,subs)
	#if it ever matters what the alignment score is and how many alignments were made I wrote some code below that can be expanded. For now I generally only care about the best alignment
	
	# ~ if return_scores:
		# ~ return()
	# ~ else:
		# ~ return(n_gaps,subs)

def get_uniprot_seq(uID):
	while True:
		try:
			r = requests.get('https://www.uniprot.org/uniprot/'+uID+'.fasta', headers={ "Content-Type" : "text/plain"})
			seq = ''.join(r.text.split('\n')[1:])
			#in rare cases entries have become obsolete and those will have an empty string as seq. I guess we are not making merge files for those
			if not seq:
				print('Returned sequence for isoform 1 of', uID, 'is empty. The entry may be obsolete. No merge file will be made since we have no authoritative sequence.')
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

#orig version where if there is only one file and allow is on it is automatically accepted. 
def find_prism_file_orig(uID, file_type, seq, v=args.verbose):
	#in contrast to the earlier version, find_prism_file_orig, this function only ever returns at max 1 file name
	
	if file_type not in prism_types:
		raise ValueError("Invalid prism file type. Chose one of: %s" % prism_types)
		
	path = os.path.join(base_prism_dir,file_type,uID[0:2], uID[2:4], uID[4:6])
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
	
	if len(listing) > 1:
		print('There is more than one' ,file_type, 'file for this uniprot ID:', uID)#, file = sys.stderr)
		if v > 1:
			print(listing)
		for prism_file in listing:
			meta, df = read_from_prism(prism_file)
			#debug
			#print(meta['protein']['sequence'])
			#debug
			if meta['protein']['sequence'] == seq:
				#if v > 1:
				print(prism_file, 'matches the authoritative sequence.')
				return prism_file
		
		else:
			if args.allow:
				max_score = 0
				best_match = None
				#todo! If args.allow, go through all files again and find the closest seq to the requested
				for prism_file in listing:
					meta, df = read_from_prism(prism_file)
					probe_seq = meta['protein']['sequence']
					target_seq = seq
				
					#this is the align function used in the Prism parser:
					align = pairwise2.align.globalds(target_seq.upper(), probe_seq.upper(), MatrixInfo.blosum62, -3, -1)
					align_scores = np.array([a[2] for a in align])
					if v > 1:
						print('Aligning', prism_file, 'to target seq.')
						print(len(align),'alignments found.')
						print('Scores are:', align_scores)
					
					if max(align_scores) > max_score:
						max_score = max(align_scores)
						best_match = prism_file
					
				#in the end, return the file that had the best alignment score, we'll try to merge that one
				if v > 1:
					print("None of the files match the authoritative sequence but allow has been used so one of the files will be merged anyway.")
					print('Best matching sequence is', best_match)
				return best_match
			
			else:
				print("None of the files match the authoritative sequence. Use -allow to enable merging of prism files with not the same sequence (but the same uniprot ID). This feature doesn't work yet.")
				return None
		
	elif (len(listing) == 1):
		#read in file to confirm seq is the desired/authoritative sequence
		meta, df = read_from_prism(listing[0])
		if meta['protein']['sequence'] == seq:
			return listing[0]
		#special case: if there is only one eligible prism file to that uniprot its seq is the closest to the req'd seq 
		elif args.allow:
			if v >= 1:
				print(listing[0], "does not match the authoritative sequence but allow has been used so it will be merged anyway.")
			return listing[0]
		else:
			if v > 1:
				print('The sequence in', listing[0], "does not match the authoritative sequence. Use -allow to enable merging of prism files with not the same sequence (but the same uniprot ID). This feature doesn't work yet.")
			return None	
			
	else:
		print('No ', file_type,' file found for uniprot ID:', uID)#, file = sys.stderr)
		return None

#new version: check seq identity since we align anyway and only qualify files that have at least x% identity with the target seq (right now 0.8) because at lower identities the parser will refuse. I can change that but it might not make sense to merge those. args.min_ID
def find_prism_file(uID, file_type, seq, v=args.verbose):
	#in contrast to the earlier version, find_prism_file_orig, this function only ever returns at max 1 file name
	
	if file_type not in prism_types:
		raise ValueError("Invalid prism file type. Chose one of: %s" % prism_types)
		
	path = os.path.join(base_prism_dir,file_type,uID[0:2], uID[2:4], uID[4:6])
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
		print('There is more than one' ,file_type, 'file for this uniprot ID:', uID)#, file = sys.stderr)
		if v > 1:
			print(listing)
		for prism_file in listing:
			meta, df = read_from_prism(prism_file)
			#debug
			#print(meta['protein']['sequence'])
			#debug
			if meta['protein']['sequence'] == seq:
				#if v > 1:
				print(prism_file, 'matches the authoritative sequence.')
				return prism_file
		
		else:
			if args.allow:
				max_score = 0
				best_match = None
				target_seq = seq
				#n_res_target = len(target_seq)
				#todo! If args.allow, go through all files again and find the closest seq to the requested
				for prism_file in listing:
					meta, df = read_from_prism(prism_file)
					probe_seq = meta['protein']['sequence']
					#n_res_data = len(probe_seq)
				
					#new version: switch to do_align
					align_obj = do_align(target_seq.upper(), probe_seq.upper())
					if v > 1:
						print('Aligning', prism_file, 'to target seq.')
						print(align_obj.n_align,'alignments found.')
						print('Scores are:', align_obj.all_scores)

					if align_obj.best_score > max_score:
						max_score = align_obj.best_score
						best_match = prism_file
						n_res = align_obj.aln_len
						coverage = align_obj.cov
						identity = align_obj.id_s1
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
	
	#disabled this because we also want an alignment if there is only one file so we can check coverage and identity	
	# ~ elif (len(listing) == 1):
		# ~ #read in file to confirm seq is the desired/authoritative sequence
		# ~ meta, df = read_from_prism(listing[0])
		# ~ if meta['protein']['sequence'] == seq:
			# ~ return listing[0]
		# ~ #special case: if there is only one eligible prism file to that uniprot its seq is the closest to the req'd seq 
		# ~ elif args.allow:
			# ~ if v >= 1:
				# ~ print(listing[0], "does not match the authoritative sequence but allow has been used so it will be merged anyway.")
			# ~ return listing[0]
		# ~ else:
			# ~ if v > 1:
				# ~ print('The sequence in', listing[0], "does not match the authoritative sequence. Use -allow to enable merging of prism files with not the same sequence (but the same uniprot ID). This feature doesn't work yet.")
			# ~ return None	
			
	else:
		print('No ', file_type,' file found for uniprot ID:', uID)#, file = sys.stderr)
		return None

def make_default_outfolder(uniprot_ID):
	#sometimes we cannot find the unirpot ID to a transcript. Put them in the none folder then
	if not uniprot_ID == 'None':
		if not os.path.exists(base_prism_dir+'prism_merge/'+ uniprot_ID[0:2]):
			os.makedirs(base_prism_dir+'prism_merge/'+ uniprot_ID[0:2])
		if not os.path.exists(base_prism_dir+'prism_merge/'+ uniprot_ID[0:2]+ '/' + uniprot_ID[2:4]):
			os.makedirs(base_prism_dir+'prism_merge/'+ uniprot_ID[0:2]+ '/' + uniprot_ID[2:4])
		if not os.path.exists(base_prism_dir+'prism_merge/'+ uniprot_ID[0:2]+ '/' + uniprot_ID[2:4] + '/' + uniprot_ID[4:6]):
			os.makedirs(base_prism_dir+'prism_merge/'+ uniprot_ID[0:2]+ '/' + uniprot_ID[2:4]+ '/' + uniprot_ID[4:6])
	
		return(base_prism_dir+'prism_merge/'+ uniprot_ID[0:2]+ '/' + uniprot_ID[2:4]+ '/' + uniprot_ID[4:6])
	else:
		if not os.path.exists(base_prism_dir+'prism_merge/None/'):
			os.makedirs(base_prism_dir+'prism_merge/None/')
		return(base_prism_dir+'prism_merge/None/')

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

#parse arguments
if (args.uniprot and args.ff) or (args.uniprot and args.ffs) or (args.ff and args.ffs):
	print('The options -uniprot, -ffs and -ff are mutually exlusive. Use only one of them.')
	sys.exit()

if args.verbose >= 1:
	print(str(datetime.now()))

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
	
	file_type = 'prism_clinvar'
	path = os.path.join(base_prism_dir,file_type,uID[0:2], uID[2:4], uID[4:6])
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
			path = os.path.join(base_prism_dir,file_type,uID[0:2], uID[2:4], uID[4:6])
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
				file_type = 'prism_gnomad'
				path = os.path.join(base_prism_dir,file_type,uID[0:2], uID[2:4], uID[4:6])
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

	else:
		if args.verbose > 1:
			print('No clinvar files, looking for gnomad file.')
		#else choose best gnomad file, which is just the one with the most entries
		file_type = 'prism_gnomad'
		path = os.path.join(base_prism_dir,file_type,uID[0:2], uID[2:4], uID[4:6])
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
			file_type = 'prism_gnomad'
			path = os.path.join(base_prism_dir,file_type,uID[0:2], uID[2:4], uID[4:6])
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

	if args.verbose > 1:
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
	
	#if we want to adapt SNP pos to the target seq do this
	if args.fill_invert:
		#if there is a gnomad file, find SNPS
		if gnomad_file is not None:
			
			#deal with gnomad file
			#load in file
			md, df = read_from_prism(gnomad_file)
			gnomad_seq = md['protein']['sequence']
			new_AF_tot = {}
			new_AF_ex = {}
			new_AF_g = {}
			add_rows = []
			
			#I do need to make another alignment to return the substitutions I want to try and fix. I think it would be more messy to try and return that while chosing the file
			#if the gnomad_seq is the target_seq it's pointless to align them. 
			if gnomad_seq == target_seq:
				subs = []
				#coverage = 1.0
				#identity = 1.0
				#pass #if they are identical there are no gaps or subs. We still get gnomad SNPs in case we need them for the uniprot and netsurfp files
				
			else: 
				align_obj = do_align(gnomad_seq, target_seq, verbose = args.verbose)
				#coverage = align_obj.cov
				#identity = align_obj.id_s2
				#n_gaps = align_obj.n_gaps
				subs = align_obj.subs
			
			#alignment quality of gnomad seq to target seq has already been checked when finding the gnomad file, so just proceed
			
			#these vars are candidates for being inverted since they are SNPs
			#they are also the candidates for being mapped onto uniprot and netsurfp seqs if there is sufficient agreement to the gnomad seq
			#mapped_SNPs[gnomad_file] = set()
			mapped_SNPs = {}
			#mapped_SNPs[gnomad_file] = {}
			#for SNP_var_name in [s.upper() for s in df.loc[df['SNP'] == True, 'variant']]:
			for SNP_var_name in [s.upper() for s in df.loc[df['AF_tot'] >= args.snp_af, 'variant']]:
				aa_ref = SNP_var_name[0]
				aa_alt = SNP_var_name[-1]
				aa_pos = int(SNP_var_name[1:-1])
				new_var = aa_alt+str(aa_pos)+aa_ref
				#calc inverse allel freqs
				new_AF_tot[new_var] = 1 - float(df.loc[df['variant'] == SNP_var_name, 'AF_tot'])
				new_AF_ex[new_var] = 1 - float(df.loc[df['variant'] == SNP_var_name, 'AF_ex'])
				new_AF_g[new_var] = 1 - float(df.loc[df['variant'] == SNP_var_name, 'AF_g'])
				#add to list of SNPs
				#mapped_SNPs[gnomad_file].add(SNP_var_name)
				if aa_pos in mapped_SNPs:
					mapped_SNPs[aa_pos].update([aa_ref, aa_alt])
				else:	
					mapped_SNPs[aa_pos] = set([aa_ref, aa_alt])
						
			#now we can go through the substiutions,i.e. diffs btw gnomad and target and find out which of them we can replace by inversion of SNPs
			for (pos,gnomad_aa,target_aa) in subs:
				if pos in mapped_SNPs and gnomad_aa in mapped_SNPs[pos] and target_aa in mapped_SNPs[pos]:
						new_var = target_aa+str(pos)+gnomad_aa
						old_var = gnomad_aa+str(pos)+target_aa
						#n_mut should be 1 since it's not a double/multi mutant. 
						new_row = [new_var] + ['NA']*(df.shape[1]-5) + [1,[gnomad_aa], [pos], [target_aa]]
						#print(new_row)
						#print('len(new_row):', len(new_row))
						#print('df shape before adding:', df.shape)
						df.loc[len(df)] = new_row
						#set SNP and AF cols. The proper row is now len -1 since we have already added the dummy row filled with NAs
						df.loc[len(df)-1, 'SNP'] = True
						df.loc[len(df)-1, 'AF_tot'] = new_AF_tot[df.loc[len(df)-1,'variant']]
						df.loc[len(df)-1, 'AF_ex'] = new_AF_ex[df.loc[len(df)-1,'variant']]
						df.loc[len(df)-1, 'AF_g'] = new_AF_g[df.loc[len(df)-1,'variant']]
						#drop the old line because otherwise the prism parser complains about two diff WTs at pos
						df.drop(df[df["variant"] == old_var].index, inplace=True)
						df.reset_index(drop=True, inplace=True)
						#fix the seq in the metadata
						md['protein']['sequence'] = md['protein']['sequence'][:pos-1] + target_aa + md['protein']['sequence'][pos:]
						#print('df shape after adding and removing orig line:', df.shape)
				
			#only do this if we have made changes
			if len(subs):
				#print(df.tail(n=10))
				#df.to_csv('/storage1/hezscha/genome_proteome_map/results/error_checking/gnomad_inv.csv')
				#write temp prism file and change pointer to that file
				new_gnomad = re.sub('.txt','_SNPinv.txt',gnomad_file)
				write_prism(md,df,new_gnomad)
				gnomad_file = new_gnomad
				#print('leaving for now')
				#sys.exit()
				
			#attempt to transfer SNPs to uniprot and netsurfp seqs
			for filename in [uniprot_file, netsurfp_file]:
				if filename is not None:
					md, df = read_from_prism(filename)
					probe_seq = md['protein']['sequence']
					
					#first, establish if SNP inv necessary, i.e. if there are mismatches to target seq
					if probe_seq == target_seq:
						continue #if they are identical we have nothing to do
					
					else:
						align_obj = do_align(probe_seq, target_seq, same_length = True, verbose = args.verbose)

						#following the same logic as above, if uniprot_file/netsurfp_file are not None their seq has passed alignment requirements to the target seq (wrt coverage and identity) so we do not need to confirm this again (since we have moved away from 'only same isoform')
						#SNP inv is only possible if seq and target seq are same isoform	
						#if align_obj.n_gaps is None or align_obj.n_gaps > 0 or len(align_obj.subs) > 10:
						#	continue
						
					
						#else, check if SNP transfer from gnomad is possible. We have not aligned these two previsouly
						align_gnomad = do_align(probe_seq, gnomad_seq, same_length = True, verbose = args.verbose)
						coverage = align_gnomad.cov
						#I'm using the identity of the gnomad seq here, normalized by the length of the gnomad seq, hence s2
						identity = align_gnomad.id_s2 
						if coverage < args.min_cov or identity < args.min_id:
							if v >= 1:
								print("Alignment between", filename, "and gnomad does not fullfill the requirements for minimum alignment coverage and identity so we cannot transfer SNPs from gnomad. It has coverage", round(coverage,2), "and identity", round(identity,2), ". Required are at least coverage", args.min_cov, 'and identity', args.min_id)
								continue
						
						#if align_gnomad.n_gaps is not None and align_gnomad.n_gaps == 0 and len(align_gnomad.subs) <= 10:
						else:
							#now, go through the subs btw seq and target seq and see if this pos is a SNP in gnomad. If yes we just invert the var name of that row in current file.
							for (pos,seq_aa,target_aa) in subs:
								#is this pos known as a SNP in gnomad and both the target_aa and seq_aa are known allowed aa's? Bascially the same check as above for gnomad SNP inv
								if pos in mapped_SNPs and seq_aa in mapped_SNPs[pos] and target_aa in mapped_SNPs[pos]:
									#in pos spec files vars are encoded only on their pos (assuming features apply to all 20 possible vars)
									new_var = target_aa+str(pos)+'='
									old_var = seq_aa+str(pos)+'='
									#print('old_var:', old_var)
									#print('new_var:', new_var)
									
									#swapp the variant name and aa_ref and aa_var columns
									df.loc[df['variant'] == old_var, 'variant'] = new_var
									df.loc[df['variant'] == new_var, 'aa_ref'] = [target_aa]
									#df.loc[df['variant'] == new_var, 'aa_var'] = ['=']
									
									#print(df.loc[df['variant'] == old_var])
									#print(df.loc[df['variant'] == new_var])

									#fix the seq in the metadata
									#print('target_aa:', target_aa)
									md['protein']['sequence'] = md['protein']['sequence'][:pos-1] + target_aa + md['protein']['sequence'][pos:]
									#print('df shape after adding and removing orig line:', df.shape)
				
							if len(subs):
								#print(md['protein']['sequence'][545])
								#print(df.tail(n=10))
								#print new file and change pointer
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
		print('Not enough files to merge, need at least 2.')
		continue
	#if only the netsurfp and uniprot files exist but we have asked for more types so they will be filled up and we have used -cleanup, the resulting df will always be empty since we remove all lines having only netsurfp and uniprot data. So just stop at this point is this is the case	
	elif len(files) == 2 and uniprot_file in files and netsurfp_file in files and uniprot_file.endswith('_filled.txt') and netsurfp_file.endswith('_filled.txt') and args.cleanup:
		print('Only prism_uniprot and prism_netsurfp files are available and cleanup was requested which will result in an empty dataframe so no merging is done.')
		if args.fail:
			with open(fail_path, 'a') as fOUT:
				print(uID, file = fOUT)
		continue
		
	else:	
		(all_merged_metadata, all_merged_df) = merge_prism(out=merged_prism_file,files=files,uID=uID,target_seq=target_seq)
		
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
			print(uID, ': the resulting dataframe after cleaning is empty. Removing the merge file:', merged_prism_file)
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
		
	#else:
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

##################





