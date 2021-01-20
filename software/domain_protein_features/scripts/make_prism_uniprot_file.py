import argparse
import re
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
import urllib.parse
import urllib.request

# Third party imports
from Bio.PDB import PDBParser
from Bio.PDB.PDBList import PDBList 
from Bio.PDB.DSSP import DSSP
#from IPython.display import display
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd 
import seaborn as sns
from scipy import stats
#import sklearn.model_selection
import xmltodict
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

#local paths for HZ. Change this to your paths
#server	
path_base = '/storage1/hezscha/'
import_path_base = path_base + 'src/'

try:
	# The insertion index should be 1 because index 0 is this file
	#sys.path.insert(1, '/groups/sbinlab/tiemann/repos/PRISM/prism/scripts/')  # the type of path is string
	sys.path.insert(1, import_path_base + 'PRISM/prism/scripts/')
	# because the system path already have the absolute path to folder a
	# so it can recognize file_a.py while searching 
	from PrismData import PrismParser, VariantData#VariantParser, VariantData

	#sys.path.insert(1, '/groups/sbinlab/tiemann/repos/PRISM/PRISM/software/rosetta_ddG_pipeline')
	sys.path.insert(1, import_path_base + 'PRISM/software/rosetta_ddG_pipeline')
	from mp_prepare import mp_superpose_opm, mp_TMalign_opm, mp_span_from_pdb_octopus, mp_lipid_acc_resi, mp_span_from_pdb_dssp

	#sys.path.insert(1, '/groups/sbinlab/tiemann/repos/PRISM/MPRISM/helper')
	sys.path.insert(1, import_path_base + 'MPRISM/helper')
	from uniprot_search import extract_uniprot_fasta, extract_by_uniprot_fasta
	from restAPI_scripts import extract_from_semanticscholar, extract_from_biorxiv

	sys.path.insert(1, import_path_base + 'MPRISM/helper/plot/')
	import prism2heatmap

	sys.path.insert(1, import_path_base + 'PRISM/software/domain_protein_features/scripts/')
	from get_domains import extract_single_protein_pfam, extract_pfam_nested, extract_pfam_pdb_mapping, extract_pfam_all_release
	from get_db_features import get_disprot_disordered, get_mobidb_disordered
	from get_human_proteome import extract_uniprot_human_proteome_ids
	from get_uniprot_features import extract_uniprot_human_proteome_info, extract_uniprot_info
	from map_domains_features import extract_single_infos, map_pfam_pdb, get_domain_features_human_proteome

except (ModuleNotFoundError, ImportError) as e:
	logger.error("{} fileure".format(type(e)))
	print(e)
else:
	logger.info("Import succeeded")

#parse arguments
################################################################################
parser = argparse.ArgumentParser()
parser.add_argument('-uniprot', dest="uniprot", help="Comma separated list of uniprot IDs (no white space)")
parser.add_argument('-fromfile', dest="fromfile", help="A file from which to read uniprot IDs (one per line)")
parser.add_argument('-outdir', dest="outdir", help="Output directory for prism_uniprot files. Will be current working directory if not given.")
parser.add_argument('-m', dest="mode", choices=['overwrite', 'leave'], default = 'leave', help="What do when the output file already exists. Leave (default) or overwrite")
args = parser.parse_args()
################################################################################

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

def merge_prism(source, target, out):
	shell_command = ['python', '/lindorffgrp-isilon/tiemann/dev/repos/prism/scripts/PrismData.py',
					 source, target, '--merge', out]
	subprocess.run(shell_command,
				   stdout=subprocess.PIPE, stderr=subprocess.PIPE,
				   check=True, text=True)
	meta_data, dataframe = read_from_prism(out)
	return meta_data, dataframe

def download_pdb(pdb_id, output_dir):
	pdbl = PDBList()
	pdbl.retrieve_pdb_file(pdb_id, file_format='pdb', pdir=output_dir)
	output_name = os.path.join(output_dir, f'{pdb_id}.pdb')
	shutil.move(os.path.join(output_dir, f'pdb{pdb_id}.ent'), output_name)
	return(output_name)

#features are requested from the uniprot API using the terminology defined on this page: https://www.uniprot.org/help/uniprotkb_column_names (right side, 'Column names as displayed in URL')
#in the returned object, the names used for the features are as described here: https://web.expasy.org/docs/userman.html#FT_line
#as an example, we request feature(TRANSMEMBRANE) and get a string like this: TRANSMEM 35..58;  /note="Helical; Name=1"; TRANSMEM 72..95;  /note="Helical; Name=2"; TRANSMEM 107..129;  /note="Helical; Name=3"; ect
#mappings of the requested feature name and the returned feature name are in the two dicts below. Right now I'm only using return_to_request_feature

request_to_return_feature = {
	'feature(ACTIVE%20SITE)':'ACT_SITE', 
	'feature(BINDING%20SITE)':'BINDING', 
    'feature(DNA%20BINDING)':'DNA_BIND', 
    'feature(METAL%20BINDING)':'METAL',
    'feature(NP%20BIND)':'NP_BIND', 
    'feature(SITE)':'SITE',
    #topologies
    #INTRAMEM - Extent of a region located in a membrane without crossing it.
    'feature(TRANSMEMBRANE)':'TRANSMEM', 
    'feature(INTRAMEMBRANE)':'INTRAMEM', 
    'feature(TOPOLOGICAL%20DOMAIN)':'TOPO_DOM',
    #PTMs
    #CHAIN - Extent of a polypeptide chain in the mature protein.
    #PEPTIDE - Extent of a released active peptide.
    'feature(LIPIDATION)':'LIPID', 
    'feature(MODIFIED%20RESIDUE)':'MOD_RES',
    'feature(GLYCOSYLATION)':'CARBOHYD',
    'feature(INITIATOR%20METHIONINE)':'INIT_MET', 
    'feature(PEPTIDE)':'PEPTIDE',
    'feature(SIGNAL)': 'SIGNAL', 
    'feature(TRANSIT)':'TRANSIT',
    'feature(CROSS%20LINK)':'CROSSLNK', 
    'feature(DISULFIDE%20BOND)':'DISULFID',
    #I don't think we need these
    #'feature(PROPEPTIDE)', 'feature(CHAIN)', 
    #secondary structure
    'feature(BETA%20STRAND)':'STRAND', 
    'feature(HELIX)':'HELIX', 
    'feature(TURN)':'TURN',
    #Domains: Aren't we getting those from pfam?
    #COMPOSITIONAL%20BIAS - Extent of a compositionally biased region. Example: /note="Glu-rich (acidic)"
    #MOTIF - Short (up to 20 amino acids) sequence motif of biological interest.
    #REGION - Extent of a region of interest in the sequence. Examples: /note="Zymogen activation region", /note="Possesses antibiotic activity"
    'feature(COILED%20COIL)':'COILED', 
    'feature(COMPOSITIONAL%20BIAS)':'COMPBIAS',
    'feature(MOTIF)':'MOTIF', 
    'feature(REGION)':'REGION', 
    'feature(ZINC%20FINGER)':'ZN_FING'
}

return_to_request_feature = {
	'ACT_SITE':'feature(ACTIVE%20SITE)', 
	'BINDING':'feature(BINDING%20SITE)', 
    'DNA_BIND':'feature(DNA%20BINDING)', 
    'METAL':'feature(METAL%20BINDING)',
    'NP_BIND':'feature(NP%20BIND)', 
    'SITE':'feature(SITE)',
    #topologies
    #INTRAMEM - Extent of a region located in a membrane without crossing it.
    'TRANSMEM':'feature(TRANSMEMBRANE)', 
    'INTRAMEM':'feature(INTRAMEMBRANE)', 
    'TOPO_DOM':'feature(TOPOLOGICAL%20DOMAIN)',
    #PTMs
    #CHAIN - Extent of a polypeptide chain in the mature protein.
    #PEPTIDE - Extent of a released active peptide.
    'LIPID':'feature(LIPIDATION)', 
    'MOD_RES':'feature(MODIFIED%20RESIDUE)',
    'CARBOHYD':'feature(GLYCOSYLATION)',
    'INIT_MET':'feature(INITIATOR%20METHIONINE)', 
    'PEPTIDE':'feature(PEPTIDE)',
    'SIGNAL':'feature(SIGNAL)', 
    'TRANSIT':'feature(TRANSIT)',
    'CROSSLNK':'feature(CROSS%20LINK)', 
    'DISULFID':'feature(DISULFIDE%20BOND)',
    #I don't think we need these
    #'feature(PROPEPTIDE)', 'feature(CHAIN)', 
    #secondary structure
    'STRAND':'feature(BETA%20STRAND)', 
    'HELIX':'feature(HELIX)', 
    'TURN':'feature(TURN)',
    #Domains: Aren't we getting those from pfam?
    #COMPOSITIONAL%20BIAS - Extent of a compositionally biased region. Example: /note="Glu-rich (acidic)"
    #MOTIF - Short (up to 20 amino acids) sequence motif of biological interest.
    #REGION - Extent of a region of interest in the sequence. Examples: /note="Zymogen activation region", /note="Possesses antibiotic activity"
    'COILED':'feature(COILED%20COIL)', 
    'COMPBIAS':'feature(COMPOSITIONAL%20BIAS)',
    'MOTIF':'feature(MOTIF)', 
    'REGION':'feature(REGION)', 
    'ZN_FING':'feature(ZINC%20FINGER)'
}

#in some features, the 'FROM' and 'TO' endpoints designate the two residues which are linked instead of a range the feature applies to! This is the case for: DISULFID, CROSSLNK
direct_link_features = {
    'feature(CROSS%20LINK)':'CROSSLNK', 
    'feature(DISULFIDE%20BOND)':'DISULFID', 
}

def make_uniprot_prism_files(uniprot_id, prism_file, version=1):
	
	'''
	
	#OBS!
	It turns out features that can have notes don't always do. I.e. p53 P04637 has this Dna binding feature without a note: DNA_BIND 102..292
	
	Feature that have several entries i.e. DISULFID 265..?274;  /evidence="ECO:0000255|PROSITE-ProRule:PRU00460"; DISULFID 267..295;  /evidence="ECO:0000255|PROSITE-ProRule:PRU00460"; DISULFID 297..306;  /evidence="ECO:0000255|PROSITE-ProRule:PRU00460"; DISULFID 309..329;  /evidence="ECO:0000255|PROSITE-ProRule:PRU00460"; DISULFID 332..341;  /evidence="ECO:0000255|PROSITE-ProRule:PRU00460"; DISULFID 334..359;  /evidence="ECO:0000255|PROSITE-ProRule:PRU00460"; DISULFID 362..371;  /evidence="ECO:0000255|PROSITE-ProRule:PRU00460" 
	are usually distinct from each other and we should probably write something different for each entry. If feature entries have names we can use that, otherwise enumerate. Actually I asked Johanna and she said we should enumerate even if they have names, so name + number
	
	
	#two conceptually different types of features:
	#1. features where you need the other fields to usefully describe the feature:
	#function
	#ACT_SITE - what kind? /note="Proton acceptor" /note="Charge relay system; for serine protease NS3; activity"
	#BINDING - to what? /note="Heme (covalent)"
	#DNA_BIND - which motif? /note="Homeobox" /note="H-T-H motif"
	#METAL - /note="Iron (heme axial ligand); shared with dimeric partner" 
	#NP_BIND - which compound is bound? /note="ATP" /note="FAD"
	#SITE - /note="Cleavage; by thrombin"
	
	#topologies
	#TOPO_DOM - which one? /note="Cytoplasmic" /note="Extracellular"
	
	#PTMs
	#LIPID - which one? /note="N-myristoyl glycine; by host"
	#MOD_RES -  /note="N-acetylalanine" /note="Phosphoserine; by CK1"
	#CARBOHYD - /note="N-linked (GlcNAc...) asparagine" /note="O-linked (Ara...) hydroxyproline"
	#PEPTIDE - which active peptide is this? /note="Arg-vasopressin" /note="Met-enkephalin"
	#TRANSIT - one might want to know where to: /note="Chloroplast" /note="Mitochondrion" ect
	#CROSSLNK - by what?  /note="Isoglutamyl cysteine thioester (Cys-Gln)"
	#OBS! CROSSLNK is also a direct link feature, so the two positions named are linked instead of describing a range!
	
	#regions
	#COMPBIAS - which one? /note="Glu-rich (acidic)"
	#MOTIF - which one? /note="Nuclear localization signal"
	#REGION - which region? /note="Zymogen activation region"
	
	#-> potential issue: /note can go over several 'lines', so there might be ; in between (see METAL, ACT_SITE)
	#-> might want to do regex instead of split
	
	#2. features that describe themselves:
	#function
	
	#topologies
	#TRANSMEM, INTRAMEM
	
	#PTMs
	#INIT_MET
	#SIGNAL - it doesn't say more info
	#DISULFID
	
	#secondary structure
	#HELIX, STRAND, TURN
	
	#regions
	#COILED - there doesn't seem to be more info in the example
	#ZN_FING - though they seem to have categories: /note="GATA-type" /note="NR C4-type"
	'''
	
	#compile regex for later
	left_missing = re.compile('\s\?[.]{2}\d+')
	right_missing = re.compile('\d+[.]{2}\?$')
	
	left_open = re.compile('\?(\d+)[.]{2}(\d+)')
	right_open = re.compile('(\d+)[.]{2}\?(\d+)')
	both_open = re.compile('\?(\d+)[.]{2}\?(\d+)')
	
	read_note = re.compile('\s(\d+)(?:\.\.(\d+))?(?:;\s+/note="([^"]+?)"(?:[;]|$))*')
	only_note = re.compile(';\s+/note="([^"]+?)"(?:[;]|$)')
	
	#self_features are features that explain themselves as listed above. We do not extract further info for these but we do enumerate them like TRANSMEM1, TRANSMEM2, ect
	self_features = [
		'TRANSMEM', 'INTRAMEM',
		'INIT_MET', 'SIGNAL', 'DISULFID',
		'HELIX', 'STRAND', 'TURN',
		'COILED', 'ZN_FING'
	]
	
	dlf = ['DISULFID', 'CROSSLNK']
	
	logger.info('Extract info from uniprot')

	uniprot_info_df = extract_uniprot_info(uniprot_id)
	sleep(1)
	#for debugging
	#print(uniprot_info_df.columns)

	logger.info('Extract domain from pfam')
	pfam_dict = extract_single_protein_pfam( uniprot_id, verbose=False )
	sleep(1)

	logger.info('Generate dataframe with info')
	
	#the order in this list imposes the order of columns in the resulting dataframe and prism_file. 
	#NOTE: Except for 'variant' and 'pfam_name' which are not retrieved from the uniprot API you must use the terms specified in return_to_request_feature otherwise it can't find them!
	variant_list = ['variant', 
					'pfam_name',
					'ACT_SITE',
					'BINDING',
					'DNA_BIND',
					'METAL',
					'NP_BIND',
					'SITE',
					'TRANSMEM', 
					'INTRAMEM',
					'TOPO_DOM', 
					'LIPID',
					'MOD_RES',
					'CARBOHYD',
					'INIT_MET',
					'PEPTIDE',
					'SIGNAL',
					'TRANSIT',
					'CROSSLNK',
					'DISULFID',
					'STRAND',
					'HELIX',
					'TURN',
					'COILED',
					'COMPBIAS',
					'MOTIF', 
					'REGION', 
					'ZN_FING']
	
	#for later haha
	DBs = ['disprot', 'mobidb']
	
	#some entries may have become obsolete like A0A087X1A0. The return then contains only the entry and entry name
	if uniprot_info_df['sequence'].empty:
		print("A0A087X1A0 is obsolete, can't obtain data.")
		return()
	
	#print(uniprot_info_df['sequence'])
	#print(uniprot_info_df['sequence'][0])
	#print(len(uniprot_info_df['sequence'][0]))
	#print(range(0, len(uniprot_info_df['sequence'][0])))
	
	#As far as I can see the dataframe row names have to start at 0 for the write_prism method to work.
	#So for assignements generally index = position-1 or the other way around, position = index+1
	output_df = pd.DataFrame(data='None', index=range(0, len(uniprot_info_df['sequence'][0])), columns=variant_list+DBs)
	#print(output_df.shape)
	#list of features that had uncertain placements, indicated by a ? in front of either start, end or both
	uncertain = [] 
	missing_placement = []
	#for features that are explained by their notes we want to count up per note to avoid confusion. I.e. if there are two active site with /note="Proton acceptor" the first one is Proton_acceptor|1 and the second one Proton_acceptor|2. To do that, we need to remember which notes exist per feature and count up
	existing_notes = {}
	
	
	for index, res in enumerate(uniprot_info_df['sequence'][0]):
		#need +1 here since enumerate starts at 0
		output_df.loc[index,'variant'] = res+str(index+1)+'='
	
	#pfam
	#return can be empty if this uniprot ID was not in pfam
	if pfam_dict is not None:
		for pfam in pfam_dict:
			#all positions between pfam['start'] and pfam['end'] should inherit the pfam['id']
			for index in range(int(pfam['start']), int(pfam['end'])+1):
				output_df.loc[index-1,'pfam_name'] = pfam['id']
		
	#features from the uniprot API
	for ft in variant_list:
		if ft == 'variant' or ft == 'pfam_name':
			continue
		#self explaining features just get their name plus a number as entries	
		elif ft in self_features:
			#count for enumerating instances of the feature
			count = 1
			#we use the mapping between the return name and the request name of features to address the correct field in uniprot_info_df
			if len(uniprot_info_df[return_to_request_feature[ft]][0]) > 0:
				ls = uniprot_info_df[return_to_request_feature[ft]][0].split(';')
				for i in range(len(ls)):
					if ls[i].startswith(' '+ft) or ls[i].startswith(ft):
						ls[i] = ls[i].split()[1]
						#is the placement uncertain and if yes in which position? If either of the positions is replaced entirely by ? we do not write this feature into the df and instead make an entry in the missing_placement row of the header
						
						#uncertain:
						if '?' in ls[i]:
							#missing:
							if left_missing.match(ls[i]) or right_missing.match(ls[i]):
								missing_placement.append(ft+'_'+ls[i])
								continue #skip to the next entry of this feature, i.e. the next ls[i]
						
							#left open:
							elif left_open.match(ls[i]):
								match_obj = left_open.match(ls[i])
								f_start = int(match_obj.groups()[0])
								f_end = int(match_obj.groups()[1])
								uncertain.append(ft+'_'+ls[i])
								
								#direct link or range feature?
								if ft in dlf:
									output_df.loc[f_start-1,ft] = '?'+ft+'|'+str(count)
									output_df.loc[f_end-1,ft] = ft+'|'+str(count)
								else:
									output_df.loc[f_start-1,ft] = '?'+ft+'|'+str(count)
									for index in range(f_start+1, f_end+1):
										output_df.loc[index-1,ft] = ft+'|'+str(count)
								
								count += 1
								
							#right open:
							elif right_open.match(ls[i]):
								match_obj = right_open.match(ls[i])
								f_start = int(match_obj.groups()[0])
								f_end = int(match_obj.groups()[1])
								uncertain.append(ft+'_'+ls[i])
								
								if ft in dlf:
									output_df.loc[f_start-1,ft] = ft+'|'+str(count)
									output_df.loc[f_end-1,ft] = '?'+ft+'|'+str(count)
								else:
									for index in range(f_start, f_end):
										output_df.loc[index-1,ft] = ft+'|'+str(count)
									output_df.loc[f_end-1,ft] = '?'+ft+'|'+str(count)	
								
								count += 1
								
							#both open
							elif both_open.match(ls[i]):
								match_obj = right_open.match(ls[i])
								f_start = int(match_obj.groups()[0])
								f_end = int(match_obj.groups()[1])
								uncertain.append(ft+'_'+ls[i])
								
								if ft in dlf:
									output_df.loc[f_start-1,ft] = '?'+ft+'|'+str(count)
									output_df.loc[f_end-1,ft] = '?'+ft+'|'+str(count)
								else:
									output_df.loc[f_start-1,ft] = '?'+ft+'|'+str(count)
									for index in range(f_start+1, f_end):
										output_df.loc[index-1,ft] = ft+'|'+str(count)
									output_df.loc[f_end-1,ft] = '?'+ft+'|'+str(count)
									
								count += 1	

							else:
								print('Could not determine position in', ls[i])
								sys.exit()
								
						#else the pos is certain. 
						#do we have two or one positions? Disulfid links can also be between chains, so there can be only a single position named (like position 30 links somewhere in another chain)
						else:
							if ft in dlf:
								if '..' in ls[i]:
									f_start = int(ls[i].split('..')[0])
									f_end = int(ls[i].split('..')[1])
									output_df.loc[f_start-1,ft] = ft+'|'+str(count)
									output_df.loc[f_end-1,ft] = ft+'|'+str(count)
									
									count += 1
								
								else:
									#actually we should read the note in this case because it describes where the link goes
									f_pos = int(ls[i])
									output_df.loc[f_pos-1,ft] = ft+'|'+str(count)
									
									count += 1
						
							else:
								if '..' in ls[i]:
									f_start = int(ls[i].split('..')[0])
									f_end = int(ls[i].split('..')[1])
								
									for index in range(f_start, f_end+1):
										output_df.loc[index-1,ft] = ft+'|'+str(count)
										
									count += 1
										
								else:
									f_pos = int(ls[i].split()[-1])
									output_df.loc[f_pos-1,ft] = ft+'|'+str(count)
									
									count += 1	
						
						
		#else it's an explained feature. explained features get what is in their notes as entry
		#we also want to enumerate these. The issue is we need to know how many i.e. Proton_acceptor there have been already so we cannot just use the iterator i.
		else:
			#init a subdict for this feature in notes dict
			existing_notes[ft] = {}
			#some of these features lack a note and for those we will use the no_note_count
			no_note_count = 1
			
			if len(uniprot_info_df[return_to_request_feature[ft]][0]) > 0:
				string_entry = uniprot_info_df[return_to_request_feature[ft]][0]
				ls = string_entry.split(ft)[1:] #the first item of the split will always be empty because the string starts with the feature name
				
				for item in ls:
					#might need separate case for uncertain placements becaues I need to know which case it is missing or open and which side
					if left_missing.match(item) or right_missing.match(item):
						missing_placement.append(ft+'_'+item.split(';'))
						continue #skip to the next entry of this feature, i.e. the next ls[i]
					
					#left open:
					elif left_open.match(item):
						match_obj = left_open.match(item)
						f_start = int(match_obj.groups()[0])
						f_end = int(match_obj.groups()[1])
						
						#need to read out note if it exists
						if '/note=' in item:
							match_obj = only_note.search(item)
							if match_obj:
								note = match_obj.groups()[1]
								if note in existing_notes[ft]:
									count = existing_notes[ft][note]+1
									existing_notes[ft][note] += 1
								else:
									count = 1
									existing_notes[ft][note] = 1
							
							else:
								print('Could not extract note from', item)
								sys.exit()	
						else:
							note = ft
							count = no_note_count
							no_note_count += 1
						
						uncertain.append(ft+'_'+item.split(';'))
						
						#direct link or range feature?
						if ft in dlf:
							output_df.loc[f_start-1,ft] = '?'+note+'|'+str(count)
							output_df.loc[f_end-1,ft] = note+'|'+str(count)
						else:
							output_df.loc[f_start-1,ft] = '?'+note+'|'+str(count)
							for index in range(f_start+1, f_end+1):
								output_df.loc[index-1,ft] = note+'|'+str(count)
						
					#right open:
					elif right_open.match(item):
						match_obj = right_open.match(item)
						f_start = int(match_obj.groups()[0])
						f_end = int(match_obj.groups()[1])
						
						#need to read out note if it exists
						if '/note=' in item:
							match_obj = only_note.search(item)
							if match_obj:
								note = match_obj.groups()[1]
								if note in existing_notes[ft]:
									count = existing_notes[ft][note]+1
									existing_notes[ft][note] += 1
								else:
									count = 1
									existing_notes[ft][note] = 1
							
							else:
								print('Could not extract note from', item)
								sys.exit()	
						else:
							note = ft
							count = no_note_count
							no_note_count += 1
						
						uncertain.append(ft+'_'+ls[i])
						
						if ft in dlf:
							output_df.loc[f_start-1,ft] = note+'|'+str(count)
							output_df.loc[f_end-1,ft] = '?'+note+'|'+str(count)
						else:
							for index in range(f_start, f_end):
								output_df.loc[index-1,ft] = note+'|'+str(count)
							output_df.loc[f_end-1,ft] = '?'+note+'|'+str(count)	
					
					#both open
					elif both_open.match(item):
						match_obj = right_open.match(item)
						f_start = int(match_obj.groups()[0])
						f_end = int(match_obj.groups()[1])
						
						#need to read out note if it exists
						if '/note=' in item:
							match_obj = only_note.search(item)
							if match_obj:
								note = match_obj.groups()[1]
								if note in existing_notes[ft]:
									count = existing_notes[ft][note]+1
									existing_notes[ft][note] += 1
								else:
									count = 1
									existing_notes[ft][note] = 1
							
							else:
								print('Could not extract note from', item)
								sys.exit()	
						else:
							note = ft
							count = no_note_count
							no_note_count += 1
						
						uncertain.append(ft+'_'+item)
						
						if ft in dlf:
							output_df.loc[f_start-1,ft] = '?'+note+'|'+str(count)
							output_df.loc[f_end-1,ft] = '?'+note+'|'+str(count)
						else:
							output_df.loc[f_start-1,ft] = '?'+note+'|'+str(count)
							for index in range(f_start+1, f_end):
								output_df.loc[index-1,ft] = note+'|'+str(count)
							output_df.loc[f_end-1,ft] = '?'+note+'|'+str(count) 
					
					#else the placement is not uncertain. Perform the 'normal' extraction except also getting a count added to the note
					else:
						#a feature that has a note field and a location 
						match_obj = read_note.search(item)

						#though features in this group may have a note, they don't necessarily do. I.e. p53 P04637 has this DNA binding feature without a note: 
						#DNA_BIND 102..292 
						#and those Metal binding features with notes: 
						#METAL 176;  /note="Zinc"; METAL 179;  /note="Zinc"; METAL 238;  /note="Zinc"; METAL 242;  /note="Zinc" 
						
						if match_obj:
							f_start = int(match_obj.groups()[0])
							
							#check if there was a note found:
							#in any case, you should not allow the note to have white space otherwise it will make problems with writing the whitespace separated prism files
							if match_obj.groups()[2] is None:
								#if there was no note to clarify, we will just save the feature's name instead as we do with features that don't have notes/do not require an explanation
								note = ft
								count = no_note_count
								no_note_count += 1
							else:
								note = re.sub(' ', '_', match_obj.groups()[2])
								if note in existing_notes[ft]:
									count = existing_notes[ft][note]+1
									existing_notes[ft][note] += 1
								else:
									count = 1
									existing_notes[ft][note] = 1
							
							#there are either one or two positions named: If it's one, the features applies to that position. If it's two, the feature applies to the range between them.
							#if the second match group is empty there was only one position
							if match_obj.groups()[1] is None:
								#assign the note to the single named position
								output_df.loc[f_start-1,ft] = note+'|'+str(count)
							
							else:
								f_end = int(match_obj.groups()[1])
								#if the feature is both one that needs explaining and a direct link feature (right now that's only crosslink but let's keep it general), only assign the note to the two named positions
								if ft in dlf:
									output_df.loc[f_start-1,ft] = note+'|'+str(count)
									output_df.loc[f_end-1,ft] = note+'|'+str(count)
							
								#else it's a range, fill every position between the two named including them with the note
								else:
									for index in range(f_start, f_end+1):
										output_df.loc[index-1,ft] = note+'|'+str(count)
						
						#can't process info with the current setup, print it to see what's up
						else:
							print("Couldn't process:")
							print(item)
							print('of', ft)
							print('whole string:', string_entry)

	#database references
	#the request to uniprot API returns the identifier of this protein in another DB and then you can go look it up there
	#we need this for disprot and mobi
	#example is p35, P04637(uniprot), DP00086(disprot)
	#actually for disprot as Johanna also wrote you can just directly use the uniprot ID without getting the disprot one
	disp_list = get_disprot_disordered(uniprot_id)
	sleep(1)
	#print('Disprot:', disp_list)
	if not disp_list == 'NA':
		count = 1
		#returns a list of lists: (start, end, type),(start, end, type)
		for sublist in disp_list:
			f_start = int(sublist[0])
			f_end = int(sublist[1])
			for index in range(f_start, f_end+1):
				output_df.loc[index-1,'disprot'] = sublist[2] + '|' + str(count)
			count += 1
	
	#mobi
	#mobi DB actually has an entry for p53/P04637 which you can get with:
	#r_mobi = requests.get('https://mobidb.bio.unipd.it/ws/P04637/consensus', headers={ "Content-Type" : "application/json"})
	#however, the uniprot api doesn't return a crossref (it does for disprot though), so we always query. I have re-written the query to return NA if the entry does not exist or there are no predicted disordered regions
	mobi_ls = get_mobidb_disordered(uniprot_id)
	sleep(1)
	#print('MobiDB', mobi_ls)
	if not mobi_ls == 'NA':
		count = 1
		#returns a list of lists but only with positions: (start,end),(start,end),...
		for sublist in mobi_ls:
			f_start = int(sublist[0])
			f_end = int(sublist[1])
			for index in range(f_start, f_end+1):
				output_df.loc[index-1,'mobidb'] = 'disordered' + '|' + str(count)
			count += 1	
	
	#for checking:
	#print('resulting df:')
	#print(output_df.head())
	#print(output_df)
	#print(output_df.size)
	#print(output_df.shape)
	
	#drop columns and rows that have only None
	#first turn None into proper NAs 
	output_df.replace(to_replace=['None'], value=np.nan, inplace=True)
	output_df.dropna(axis='columns', how='all', inplace=True)
	#when dropping rows with only NA ignore the first col which is the name of the var
	output_df.dropna(axis='index', how='all', subset = list(output_df.columns.values[1:]), inplace=True)
	#for some uniprot IDs like A0A0G2JQZ4 pfam and uniprot disagree on the length of the protein. So I will also remove any rows where the variant name is missing (those are created if pfam is longer than uniprot)
	output_df.dropna(axis='index', how='all', subset = ['variant'], inplace=True)
	
	#remember to reset the index after dropping rows so that it starts with 0 again, otherwise writing the prism file will fail
	output_df = output_df.reset_index(drop=True)
	
	#in rare cases the dataframe can be empty now like for A0A075B6U7.
	if output_df.empty:
		print('No features are available for', uniprot_id)
		return()
	
	#print('df after omitting NA rows and cols:')
	#print(output_df)
	#print(output_df.size)
	#print(output_df.shape)
	
	logger.info('Generate metadata')
	#column descriptions for prism header
	columns_dic = {}
	# ~ for elem in variant_list[1:]:
	for elem in output_df.columns[1:]:
		columns_dic[elem]= elem
	metadata = {
		"version": version,
		"protein": {
			"name": uniprot_info_df['Entry name'][0],
			"organism": uniprot_info_df['organism'][0],
			"uniprot": uniprot_info_df['Entry'][0],
			"sequence": uniprot_info_df['sequence'][0],
		},
		"uniprot": {
		"reviewed": uniprot_info_df['reviewed'][0],
		"uncertain_placement": ','.join(uncertain),
		"missing_placement": ','.join(missing_placement)
		},
		"columns": columns_dic,
	}

	comment = [ f"version {version} - {datetime.date(datetime.now())} - uniprot extract script",] 
	logger.info('Writing prism file')
	write_prism(metadata, output_df, prism_file, comment=comment)
	logger.info('uniprot prism files written!')
	return('Success')

#trying out the function make_uniprot_prism_files
#output_dir = path_base + 'pos_spec_prism_files/results/'
output_dir = args.outdir if args.outdir else os.getcwd()

#more examples
#uniprotIDs = ['P07550','P04637', 'P35520']
#uniprotIDs = ['P04637', 'Q9NTF0']

if args.uniprot:
	uniprotIDs = args.uniprot.split(',')
elif args.fromfile:
	uniprotIDs = []
	with open(args.fromfile) as IN:
		for line in IN:
			uniprotIDs.append(line.rstrip())
else:
	sys.exit('Please use either -uniprot to pass one or several comma separated IDs or -fromfile to read in a file with uniprot IDs.')	

for uniprot_id in uniprotIDs:
	print(uniprot_id)
	prism_file = os.path.join(output_dir, f'prism_uniprot_XXX_{uniprot_id}.txt')
	
	#if file exists, either exit or mkae new one anyway, depending on mode
	if os.path.exists(prism_file):
		if args.mode == 'leave':
			print(prism_file, 'already exists. If you wish to overwrite it use -m overwrite.')
		#there are only two modes, the other one is overwrite
		else:	
			ret = make_uniprot_prism_files(uniprot_id, prism_file, version=1)
			#try to reimport the file we just made to check for compatibility with the prism parser
			if ret == 'Success':
				meta_data, dataframe = read_from_prism(prism_file)
				logger.info('Prism file passed check')
	else:
		ret = make_uniprot_prism_files(uniprot_id, prism_file, version=1)
		#try to reimport the file we just made to check for compatibility with the prism parser
		if ret == 'Success':
			meta_data, dataframe = read_from_prism(prism_file)
			logger.info('Prism file passed check')

