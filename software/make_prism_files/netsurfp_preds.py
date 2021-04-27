
#go through a (human proteome) fasta file, check for each entry if there is a spliceAI file, if yes run netsurfp

import argparse
#import requests #url requests
#from requests import HTTPError
import sys
import re
import os
import glob
import subprocess

import_path_base = '/storage1/hezscha/src/'
sys.path.insert(1, import_path_base + 'PRISM/prism/scripts/')
from PrismData import PrismParser, VariantData

#parse arguments
################################################################################
parser = argparse.ArgumentParser()
parser.add_argument('-in', dest="infile", type = str, required=True, help="The list of uniprot IDs for which to make netsurfp predictions.")
parser.add_argument('-m', dest="mode", choices=['overwrite', 'leave'], default = 'leave', help="What do when the output file already exists. Leave (default) or overwrite")
parser.add_argument("-v","--verbose", action="count", default=0, help="Level of output, default zero is no output")
parser.add_argument('-spliceAI', dest="spliceAI", action='store_true', help="Only run if spliceAI file exists.")
parser.add_argument('-proteome', dest="proteome", action='store_true', help="Only run if uniprot ID in human proteome UP000005640_9606_ID.list (1 seq per gene, 20K proteins). I've been running on the chromosome proteome lists but those also include other proteins not in the human proteome 20K UP000005640_9606_ID.list.")
#parser.add_argument("--first_res_num", metavar="INT", help="In the output data, assign this number to the first given amino acid of the sequence")
args = parser.parse_args()
################################################################################

def call_netsurfp(uniprot_id,file_path, v=args.verbose):
	os.chdir('/storage1/hezscha/src/netsurfp')
	shell_command = ['/storage1/hezscha/src/netsurfp/build/netsurfp2', '--csv', '/storage1/shared/data/netsurfp_preds/csv/'+uniprot_id+'.csv', '--hhdb', '/storage1/shared/data/HH_databases/UniRef30_2020/UniRef30_2020_02', 'hhblits', '/storage1/hezscha/src/netsurfp/models/hhsuite.pb', file_path, '/storage1/shared/data/netsurfp_preds/']
	if v > 1:
		print(shell_command)
	subprocess.run(shell_command,
				   stdout=subprocess.PIPE, stderr=subprocess.PIPE,
				   check=True, text=True)

if args.proteome:
	proteome = set()
	with open('/storage1/shared/data/uniprot_datasets/human_proteome_UP000005640_9606.list','r') as IN:
		for line in IN:
			proteome.add(line.rstrip())

with open(args.infile, 'r') as IN:
	for line in IN:
		(seq_name, uID, seq) = line.rstrip().split(' ')
		
		#if output file exists, skip
		if os.path.exists('/storage1/shared/data/netsurfp_preds/csv/'+uID+'.csv'):
			print('netsurfp prediction for', uID, 'exists. Skipping to next protein.')
			continue

		#check existance of spliceAI file's directory and skip to next protein if doesn't exist
		if args.spliceAI:
			spliceAI_path = os.path.join('/storage1/shared/data/prism_spliceai', uID[0:2], uID[2:4], uID[4:6])
			#print(spliceAI_path)
			#if os.path.exists(spliceAI_path) and os.listdir(spliceAI_path):
			if os.path.exists(spliceAI_path) and glob.glob(os.path.join(spliceAI_path, 'prism_spliceai_XXX_'+uID+'-*')):
				pass
			else:
				print('No spliceAI prism file for', uID, ', skipping to next protein.')
				continue
		
		if args.proteome:
			if not uID in proteome:
				print(uID, 'is not part of the human proteome UP000005640_9606_ID, skipping to next protein.')
				continue
		
		#if we're not skipping print current ID to know where we are in going through the list
		print('Running netsurfp2 for', uID)

		#make tmp fasta
		tmp_fasta = os.path.join('/storage1/hezscha/genome_proteome_map/tmp_fasta',uID+'.fasta')
		tmp_out = open(tmp_fasta, 'w')
		print('>', uID, sep = '', file = tmp_out)
		for i in range(0,len(seq),60):
			print(seq[i:i+60], file = tmp_out)
		tmp_out.close()
		
		#run netsurfp
		call_netsurfp(uID, tmp_fasta)
		
		#delete fasta
		if os.path.isfile(tmp_fasta):
			os.remove(tmp_fasta)
		

#netsurfp2 --csv example/dhfr_human.csv --hhdb /storage1/shared/data/HH_databases/UniRef30_2020/UniRef30_2020_02 hhblits /storage1/hezscha/src/netsurfp/models/hhsuite.pb ../dhfr_human.fasta example/
