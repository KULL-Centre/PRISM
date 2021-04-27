"""get_db_features obtains information from different databases

Author: Johanna K.S. Tiemann

Date of last major changes: 2020-08-24

"""

# Standard library imports
import json
import urllib.request
import requests
from requests import HTTPError
import sys

def get_disprot_disordered(uniprot_id, server = 'https://disprot.org/api/'):
	"""
	Extracts disorder content (%) from disprot (https://disprot.org/) 
	or present uniprot ID
	"""
	ext = uniprot_id
	try:
		r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
		r.raise_for_status()
	except HTTPError as e:
		status_code = r.status_code
		#404 means there is no such entry
		if status_code == 404:
			return('NA')
		print('While looking up', uniprot_id, 'on disprot the following error occured',status_code)
		sys.exit()
	
	else:
		contents = r.json()
		if 'status' in contents.keys() and contents['status'] == 'not-exist':
			return('NA')
		else:
			disordered_region = []
			for elem in contents['disprot_consensus']['full']:
				disordered_region.append([elem['start'],elem['end'],elem['type']])
			return(disordered_region)

def get_mobidb_lite_disordered(uniprot_id, server = 'https://mobidb.bio.unipd.it/api/'):
	"""
	Extracts disorder content (%) from mobidb: https://mobidb.bio.unipd.it/help/apidoc
	for present uniprot_ID
	"""
	ext = "download?acc="+uniprot_id+"&projection=prediction-disorder-mobidb_lite&format=json"
	try:
		r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
		r.raise_for_status()
	except HTTPError as e:
		status_code = r.status_code
		print('While looking up', uniprot_id, 'on mobiDB the following error occured:',status_code)
		#P0DTE7 produced a code 500, internal server error. I checked and other uniprot IDs still work fine, so I'm not sure if it is a temporary issue. For now, I think we should return NA
		if status_code == 500:
			return('NA')
		sys.exit()
	
	else:
		contents = r.json()[0]
		#interestingly, some uniprot IDs that are predicted to have no disordered regions like P00374 have an entry called 'prediction-disorder-mobidb_lite' but it doesn't have any regions, while others (also predicted to have no disordered regions) just have an empty return.
		if not contents:
			return('NA')
		
		elif 'regions' in contents['prediction-disorder-mobidb_lite']:
			disordered_region = []
			for elem in contents['prediction-disorder-mobidb_lite']['regions']:
				#disordered_region.append([elem[0],elem[1]])
				disordered_region.append(elem)
			return(disordered_region)
			
		#if this key doesn't exist there are no disordered regions predicted by mobiDB-lite	
		else:
			return('NA')

def get_mobidb_consensus_priority(uniprot_id, server = 'https://mobidb.bio.unipd.it/api/'):
	"""
	Extracts disordered regions for present uniprot_ID with 'priority' as a consensus source, compare https://mobidb.bio.unipd.it/about/mobidb "Feature identifiers"
	"""
	
	#go through the evidence levels and take the highest one available
	evidence_levels = ['curated', 'derived', 'homology', 'prediction']
	
	for el in evidence_levels:
		
		#we're only interested in the regions field
		ext = "download?acc="+uniprot_id+"&projection="+el+"-disorder-priority.regions"
		try:
			r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
			r.raise_for_status()
		except HTTPError as e:
			status_code = r.status_code
			print('While looking up', uniprot_id, 'on mobiDB the following error occured:',status_code)
			#P0DTE7 produced a code 500, internal server error. I checked and other uniprot IDs still work fine, so I'm not sure if it is a temporary issue. For now, I think we should return NA
			if status_code == 500:
				return('NA')
			sys.exit()
			

		else:
			#download returns a list of jsons which in our case has only one entry since we're only asking for a single uniprot ID 
			contents = r.json()[0]
			#if this json is empty, the requested evidence level doesn't exist. Try with the next one
			if not contents:
				continue
				
			else:
				disordered_region = []
				for elem in contents[el+"-disorder-priority"]['regions']:
					#disordered_region.append([elem[0],elem[1]])
					disordered_region.append(elem)
				return(disordered_region)
				
	#if we went through all four evidence levels and haven't returned anything, return NA now
	return('NA')


def get_disprot_disordered_orig(uniprot_id):
    """
    Extracts disorder content (%) from disprot (https://disprot.org/) 
    or present uniprot ID
    """
    disprot_url = f'https://disprot.org/api/{uniprot_id}'
    request_header = {"format": "json"} 
    request = urllib.request.Request(disprot_url, headers=request_header)
    contents = urllib.request.urlopen(request).read()
    disordered_region = []
    if contents:
        all_json = json.loads(contents.decode("utf-8"))
        for elem in all_json['disprot_consensus']['full']:
            disordered_region.append([elem['start'],elem['end'],elem['type']])
    return disordered_region


def get_mobidb_disordered_orig(uniprot_id):
    """
    Extracts disorder content (%) from mobidb (https://mobidb.bio.unipd.it/) 
    for present uniprot_ID
    """
    mobidb_url = f'https://mobidb.bio.unipd.it/ws/{uniprot_id}/consensus'
    request_header = {"format": "json"} 
    request = urllib.request.Request(mobidb_url, headers=request_header)
    contents = urllib.request.urlopen(request).read()
    all_json = json.loads(contents.decode("utf-8"))
    disordered_region = []
    if 'db' in all_json['mobidb_consensus']['disorder'].keys():
        for elem in all_json['mobidb_consensus']['disorder']['db'][0]['regions']:
            if elem[2]!='S':
                disordered_region.append(elem)
    return disordered_region
