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
parser.add_argument('-outdir', dest="outdir", help="Output directory for prism_uniprot files. Will be current working directory if not given.")
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
	
	self_features = [
		'TRANSMEM', 'INTRAMEM',
		'INIT_MET', 'SIGNAL', 'DISULFID',
		'HELIX', 'STRAND', 'TURN',
		'COILED', 'ZN_FING'
	]
	
	dlf = ['DISULFID', 'CROSSLNK']
	
	logger.info('Extract info from uniprot')

	uniprot_info_df = extract_uniprot_info(uniprot_id)
	#for debugging
	#print(uniprot_info_df.columns)

	logger.info('Extract domain from pfam')
	pfam_dict = extract_single_protein_pfam( uniprot_id, verbose=False )

	logger.info('Generate dataframe with info')
	
	#the order in this list imposes the order of columns in the resulting dataframe and prism_file. 
	#NOTE: Except for 'variant' and 'pfam_name' which are not retrieved from the uniprot API you must use the terms specified in return_to_request_feature otherwise it can't find them!
	# ~ variant_list = ['variant', 'pfam_name', 'TRANSMEM', 'TOPO_DOM', 
					 # ~ 'MOTIF', 'REGION', 'DISULFID','BINDING']
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

	#As far as I can see the dataframe row names have to start at 0 for the write_prism method to work.
	#So for assignements generally index = position-1 or the other way around, position = index+1
	output_df = pd.DataFrame(data='None', index=range(0, len(uniprot_info_df['sequence'][0])), columns=variant_list+DBs)
	
	for index, res in enumerate(uniprot_info_df['sequence'][0]):
		#need +1 here since enumerate starts at 0
		output_df.loc[index,'variant'] = res+str(index+1)+'='
	
	#pfam
	for pfam in pfam_dict:
		#all positions between pfam['start'] and pfam['end'] should inherit the pfam['id']
		for index in range(int(pfam['start']), int(pfam['end'])+1):
			output_df.loc[index-1,'pfam_name'] = pfam['id']
		
	#features from the uniprot API
	for ft in variant_list:
		if ft == 'variant' or ft == 'pfam_name':
			continue
		#self explaining features just get their name as entries	
		elif ft in self_features:
			#we use the mapping between the return name and the request name of features to address the correct field in uniprot_info_df
			if len(uniprot_info_df[return_to_request_feature[ft]][0]) > 0:
				ls = uniprot_info_df[return_to_request_feature[ft]][0].split(';')
				#is it a region ft or direct link feature?
				#direct link features: the positions indicated are not a range but the two positions linked by the feature
				if ft in dlf:
					for i in range(len(ls)):
						if ls[i].startswith(' '+ft) or ls[i].startswith(ft):
						#range or single position? Disulfid links can also be between chains, so there can be only a single position named (like position 30 links to position 30 in the duplicate of the chain)
							if '..' in ls[i]:
								f_start = int(ls[i].split()[1].split('..')[0])
								f_end = int(ls[i].split()[1].split('..')[1])
								output_df.loc[f_start-1,ft] = ft
								output_df.loc[f_end-1,ft] = ft
							else:
								f_pos = int(ls[i].split()[-1])
								output_df.loc[f_pos-1,ft] = ft
					
				#else it's not a direct link feature and all positions between Start..End including them will inherit the feature
				else:
					for i in range(len(ls)):
						if ls[i].startswith(' '+ft) or ls[i].startswith(ft):
							#is range given or single position?
							if '..' in ls[i]:
								f_start = int(ls[i].split()[1].split('..')[0])
								f_end = int(ls[i].split()[1].split('..')[1])
								for index in range(f_start, f_end+1):
									output_df.loc[index-1,ft] = ft
							else:
								f_pos = int(ls[i].split()[-1])
								output_df.loc[f_pos-1,ft] = ft

		#explained features get what is in their notes as entry
		else:
			if len(uniprot_info_df[return_to_request_feature[ft]][0]) > 0:
				string_entry = uniprot_info_df[return_to_request_feature[ft]][0]
				ls = string_entry.split(ft)[1:] #the first item of the split will always be empty because the string starts with the feature name
				
				for item in ls:
					#a feature that has a note field and a location 
					match_obj = re.search('\s(\d+)(?:\.\.(\d+))?(?:;\s+/note="([^"]+?)"(?:[;]|$))*', item)

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
						else:
							note = re.sub(' ', '_', match_obj.groups()[2])
						
						#there are either one or two positions named: If it's one, the features applies to that position. If it's two, the feature applies to the range between them.
						#if the second match group is empty there was only one position
						if match_obj.groups()[1] is None:
							#assign the note to the single named position
							output_df.loc[f_start-1,ft] = note
						
						else:
							f_end = int(match_obj.groups()[1])
							#if the feature is both one that needs explaining and a direct link feature (right now that's only crosslink but let's keep it general), only assign the note to the two named positions
							if ft in dlf:
								output_df.loc[f_start-1,ft] = note
								output_df.loc[f_end-1,ft] = note
						
							#else it's a range, fill every position between the two named including them with the note
							else:
								for index in range(f_start, f_end+1):
									output_df.loc[index-1,ft] = note
					
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
	#print('Disprot:', disp_list)
	if not disp_list == 'NA':
		#returns a list of lists: (start, end, type),(start, end, type)
		for sublist in disp_list:
			f_start = int(sublist[0])
			f_end = int(sublist[1])
			for index in range(f_start, f_end+1):
				output_df.loc[index-1,'disprot'] = sublist[2]
	
	#mobi
	#mobi DB actually has an entry for p53/P04637 which you can get with:
	#r_mobi = requests.get('https://mobidb.bio.unipd.it/ws/P04637/consensus', headers={ "Content-Type" : "application/json"})
	#however, the uniprot api doesn't return a crossref (it does for disprot though), so we always query. I have re-written the query to return NA if the entry does not exist or there are no predicted disordered regions
	mobi_ls = get_mobidb_disordered(uniprot_id)
	#print('MobiDB', mobi_ls)
	if not mobi_ls == 'NA':
		#returns a list of lists but only with positions: (start,end),(start,end),...
		for sublist in mobi_ls:
			f_start = int(sublist[0])
			f_end = int(sublist[1])
			for index in range(f_start, f_end+1):
				output_df.loc[index-1,'mobidb'] = 'disordered'
	
	#for checking:
	# ~ print('resulting df:')
	# ~ print(output_df.head())
	# ~ print(output_df.size)
	# ~ print(output_df.shape)
	
	#drop columns and rows that have only None
	#first turn None into proper NAs 
	output_df.replace(to_replace=['None'], value=np.nan, inplace=True)
	output_df.dropna(axis='columns', how='all', inplace=True)
	#when dropping rows with only NA ignore the first col which is the name of the var
	output_df.dropna(axis='index', how='all', subset = list(output_df.columns.values[1:]), inplace=True)
	
	# ~ print('df after omitting NA rows and cols:')
	# ~ print(output_df.head())
	# ~ print(output_df.size)
	# ~ print(output_df.shape)
	
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
		"uniprot": {},
		"columns": columns_dic,
	}

	comment = [ f"version {version} - {datetime.date(datetime.now())} - uniprot extract script",] 
	logger.info('Writing prism file')
	write_prism(metadata, output_df, prism_file, comment=comment)
	logger.info('uniprot prism files written!')

#trying out the function make_uniprot_prism_files
#output_dir = path_base + 'pos_spec_prism_files/results/'
output_dir = args.outdir if args.outdir else os.getcwd()

#more examples
#uniprotIDs = ['P07550','P04637', 'P35520']
#uniprotIDs = ['P04637', 'Q9NTF0']
uniprotIDs = args.uniprot.split(',')

for uniprot_id in uniprotIDs:
	print(uniprot_id)
	prism_file = os.path.join(output_dir, f'prism_uniprot_XXX_{uniprot_id}.txt')
	make_uniprot_prism_files(uniprot_id, prism_file, version=1)
	#try to reimport the file we just made to check for compatibility with the prism parser
	meta_data, dataframe = read_from_prism(prism_file)
	#print('\n\n')



